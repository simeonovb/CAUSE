process RESOLVE_TRAITS {
    tag "${keyword}"
    conda "${projectDir}/envs/gwas.yml"

    // retry transient API failures, but not exit status 3 (no trait matches)
    errorStrategy { task.exitStatus == 3 ? 'terminate' : 'retry' }
    maxRetries 3

    input:
    val keyword

    output:
    path "*.traits.tsv"

    script:
    def prefix = keyword.replaceAll(/\W+/, '_').toLowerCase()
    def quoted = "'" + keyword.replace("'", "'\\''") + "'"
    """
    resolve_traits.py --keyword ${quoted} --out ${prefix}.traits.tsv
    """
}

process FETCH_GWAS {
    tag "${trait}"
    conda "${projectDir}/envs/gwas.yml"
    publishDir "${params.outdir}/gwas_catalog/per_trait", mode: 'copy', overwrite: true

    // the API is a remote service, so retry transient failures,
    // but not exit status 3, which means the trait name is unknown
    errorStrategy { task.exitStatus == 3 ? 'terminate' : 'retry' }
    maxRetries 3
    maxForks 4

    input:
    val trait
    val allow_empty

    output:
    path "*.tsv", optional: true

    script:
    def prefix = trait.replaceAll(/\W+/, '_').toLowerCase()
    def quoted = "'" + trait.replace("'", "'\\''") + "'"
    def empty_flag = allow_empty ? '--allow-empty' : ''
    """
    fetch_gwas.py --trait ${quoted} --out ${prefix}.tsv ${empty_flag}
    """
}

process MERGE_GWAS {
    conda "${projectDir}/envs/gwas.yml"
    publishDir "${params.outdir}/gwas_catalog", mode: 'copy', overwrite: true

    input:
    path tsvs, stageAs: "?/*"

    output:
    path "gwas_catalog_traits.tsv"

    script:
    """
    #!/usr/bin/env python3
    import glob

    import pandas as pd

    files = sorted(glob.glob("*/*.tsv"))
    df = pd.concat([pd.read_csv(f, sep="\\t", dtype=str) for f in files])

    # an association can come back under several traits: keep it once
    # and list every trait that returned it
    traits = df.groupby("association_id")["query_trait"].agg(lambda s: "; ".join(sorted(set(s))))
    df = df.drop_duplicates("association_id").drop(columns="query_trait").merge(traits, on="association_id")
    df.to_csv("gwas_catalog_traits.tsv", sep="\\t", index=False)
    """
}

process EXTRACT_SNPS {
    conda "${projectDir}/envs/gwas.yml"
    publishDir "${params.outdir}/ld_proxy", mode: 'copy', overwrite: true

    input:
    path table
    val pop
    val genome_build

    output:
    path "query_snps.tsv", emit: snps
    path "skipped_variants.tsv", emit: skipped

    script:
    """
    extract_snps.py --table ${table} --out query_snps.tsv --skipped skipped_variants.tsv \\
        --pop '${pop}' --genome-build ${genome_build}
    """
}

process LD_PROXY {
    tag "${snps.size()} SNPs from ${snps[0]}"
    conda "${projectDir}/envs/gwas.yml"
    secret 'LDLINK_TOKEN'

    // retry transient API failures, but not exit status 3 (bad token or settings)
    errorStrategy { task.exitStatus == 3 ? 'terminate' : 'retry' }
    maxRetries 3
    // LDlink only allows one request at a time per token
    maxForks 1

    input:
    val snps
    val pop
    val genome_build

    output:
    path "ld_proxy.tsv.gz", emit: proxies
    path "failed.tsv", emit: failed

    script:
    """
    ld_proxy.py --snps ${snps.join(',')} --pop '${pop}' --genome-build ${genome_build} \\
        --window ${params.ld_window} --out ld_proxy.tsv.gz --failed failed.tsv
    """
}

process MERGE_LD {
    conda "${projectDir}/envs/gwas.yml"
    publishDir "${params.outdir}/ld_proxy", mode: 'copy', overwrite: true

    input:
    path chunks, stageAs: "chunks/?/*"
    path failed, stageAs: "failures/?/*"
    path snps
    val r2_min

    output:
    path "ld_proxies.tsv"
    path "ld_failed.tsv"

    script:
    """
    merge_ld.py --chunks chunks --failures failures --snps ${snps} --r2-min ${r2_min} \\
        --out ld_proxies.tsv --failed-out ld_failed.tsv
    """
}

// ask a question in the terminal; used for any setting not given as a flag
def ask(String question) {
    def console = System.console()
    if (!console) {
        error "Some settings are missing and there is no interactive terminal to ask in. Pass them as flags, e.g. --traits \"heart,cardiac\" --match contains --ld_pop EUR --ld_build grch38 --ld_r2_min 0.8"
    }
    return console.readLine('%s', question)
}

// make sure an LDlink token is stored for this user, asking for one (once) if not
def ensureLdToken() {
    def store = nextflow.secret.SecretsLoader.instance.load()
    if (store.getSecret('LDLINK_TOKEN')) {
        return
    }
    def console = System.console()
    if (!console) {
        error "No LDlink token is stored. Get one at https://ldlink.nih.gov/?tab=apiaccess and run: nextflow secrets set LDLINK_TOKEN <token>"
    }
    def typed = console.readPassword('%s', 'The LD proxy step needs a personal LDlink token.\nGet one (free, sent by email) at https://ldlink.nih.gov/?tab=apiaccess\nPaste your token (hidden as you type; it is saved so you are only asked once): ')
    def token = typed ? new String(typed).trim() : ''
    if (!token) {
        error "No LDlink token given"
    }

    // try the token on one small lookup before saving it
    def reply
    try {
        def query = "var=rs3&pop=CEU&r2_d=r2&window=1000&genome_build=grch37&token=${URLEncoder.encode(token, 'UTF-8')}"
        reply = new URL("https://ldlink.nih.gov/LDlinkRest/ldproxy?${query}").getText(connectTimeout: 30000, readTimeout: 120000)
    }
    catch (IOException e) {
        error "Could not reach LDlink to check the token: ${e.message}"
    }
    if (reply.toLowerCase().contains('api token')) {
        error "LDlink did not accept that token, so it was not saved. Check it and run the pipeline again"
    }
    store.putSecret('LDLINK_TOKEN', token)
    store.close()  // writes the secrets store to disk
    log.info "LDlink token saved. To change it later: nextflow secrets set LDLINK_TOKEN <token>"
}

// 1000 Genomes populations LDlink can compute LD in
def ldPopulations() {
    return [
        ALL: 'All populations below combined',
        AFR: 'African (all 7 below)',
        YRI: '  Yoruba in Ibadan, Nigeria',
        LWK: '  Luhya in Webuye, Kenya',
        GWD: '  Gambian in Western Gambia',
        MSL: '  Mende in Sierra Leone',
        ESN: '  Esan in Nigeria',
        ASW: '  Americans of African ancestry in SW USA',
        ACB: '  African Caribbeans in Barbados',
        AMR: 'Ad Mixed American (all 4 below)',
        MXL: '  Mexican ancestry from Los Angeles, USA',
        PUR: '  Puerto Ricans from Puerto Rico',
        CLM: '  Colombians from Medellin, Colombia',
        PEL: '  Peruvians from Lima, Peru',
        EAS: 'East Asian (all 5 below)',
        CHB: '  Han Chinese in Beijing, China',
        JPT: '  Japanese in Tokyo, Japan',
        CHS: '  Southern Han Chinese',
        CDX: '  Chinese Dai in Xishuangbanna, China',
        KHV: '  Kinh in Ho Chi Minh City, Vietnam',
        EUR: 'European (all 5 below)',
        CEU: '  Utah residents with Northern and Western European ancestry',
        TSI: '  Toscani in Italy',
        FIN: '  Finnish in Finland',
        GBR: '  British in England and Scotland',
        IBS: '  Iberian population in Spain',
        SAS: 'South Asian (all 5 below)',
        GIH: '  Gujarati Indians from Houston, USA',
        PJL: '  Punjabi from Lahore, Pakistan',
        BEB: '  Bengali from Bangladesh',
        STU: '  Sri Lankan Tamil from the UK',
        ITU: '  Indian Telugu from the UK',
    ]
}

// reference panels LDlink offers, in menu order
def ldBuilds() {
    return [
        grch37: 'GRCh37 / hg19, 1000 Genomes phase 3',
        grch38: 'GRCh38 / hg38, 1000 Genomes phase 3 lifted over',
        grch38_high_coverage: 'GRCh38 / hg38, 1000 Genomes 30x high coverage',
    ]
}

def askPop() {
    def menu = ldPopulations().collect { code, label -> "  ${code}  ${label}" }.join('\n')
    return ask("Which populations should LD be calculated in?\n${menu}\nType ALL for every population, or one or more codes, comma-separated (e.g. EUR or CEU,YRI): ")
}

def askBuild() {
    def menu = ldBuilds().collect { code, label -> code }.withIndex()
        .collect { code, i -> "  ${i + 1}) ${code}  (${ldBuilds()[code]})" }.join('\n')
    return ask("Which genome assembly?\n${menu}\nChoose 1, 2 or 3 [3]: ")
}

// turn "eur, yri" into the "EUR+YRI" form LDlink expects
def parsePop(answer) {
    def codes = (answer instanceof List ? answer.join(',') : answer.toString())
        .toUpperCase().tokenize(',+ ').unique()
    def unknown = codes.findAll { !ldPopulations().containsKey(it) }
    if (!codes || unknown) {
        error "Unknown population code(s): ${unknown.join(', ') ?: '(none given)'}. Choose from: ${ldPopulations().keySet().join(', ')}"
    }
    return codes.contains('ALL') ? 'ALL' : codes.join('+')
}

def askR2() {
    return ask('Which R2 threshold should proxies meet? Proxies with R2 below it are left out.\nType a number from 0 to 1 (0 keeps everything LDlink returns): ')
}

def parseR2(answer) {
    def text = answer.toString().trim().replace(',', '.')
    if (!text.isNumber() || text.toBigDecimal() < 0 || text.toBigDecimal() > 1) {
        error "R2 threshold must be a number from 0 to 1, got '${text}'"
    }
    return text.toBigDecimal()
}

def parseBuild(answer) {
    def builds = ldBuilds().keySet() as List
    def text = answer.toString().trim().toLowerCase() ?: '3'
    def build = text in ['1', '2', '3'] ? builds[text.toInteger() - 1] : text
    if (!(build in builds)) {
        error "Unknown genome assembly '${text}'. Choose from: ${builds.join(', ')}"
    }
    return build
}

def askMatch() {
    def answer = ask('How should traits be matched?\n  1) exact trait name (e.g. "Atrial fibrillation")\n  2) any trait containing a word (e.g. "heart", "cardio")\nChoose 1 or 2 [1]: ')
    return [(''): 'exact', '1': 'exact', '2': 'contains'].get(answer?.trim()?.toLowerCase(), answer?.trim()?.toLowerCase())
}

workflow {
    // match mode: --match, else exact when --traits was given, else ask
    match = params.match ?: (params.traits ? 'exact' : askMatch())
    if (!(match in ['exact', 'contains'])) {
        error "Unknown match mode '${match}'. Use 'exact' or 'contains'"
    }

    // accept a list from a config/params file, a comma-separated --traits string, or ask
    question = match == 'exact'
        ? 'Which traits do you want SNPs for? (exact names, comma-separated): '
        : 'Which words should the trait names contain? (comma-separated): '
    answer = params.traits ?: ask(question)
    terms = (answer instanceof List ? answer : answer.toString().tokenize(','))
        .collect { it.toString().trim() }
        .findAll { it }
        .unique { it.toLowerCase() }
    if (!terms) {
        error "No traits given"
    }

    // LD settings are asked up front too, so the run needs no further input
    ld_pop = parsePop(params.ld_pop ?: askPop())
    ld_build = parseBuild(params.ld_build ?: askBuild())
    ld_r2_min = parseR2(params.ld_r2_min != null ? params.ld_r2_min : askR2())
    ensureLdToken()
    log.info "LD proxies: populations ${ld_pop}, assembly ${ld_build}, R2 >= ${ld_r2_min}"

    if (match == 'contains') {
        log.info "Fetching GWAS Catalog associations for traits containing: ${terms.join(', ')}"
        RESOLVE_TRAITS(channel.fromList(terms))

        // keep a record of which traits each word matched
        RESOLVE_TRAITS.out.collectFile(
            name: 'matched_traits.tsv',
            keepHeader: true,
            storeDir: "${params.outdir}/gwas_catalog",
        )
        traits_ch = RESOLVE_TRAITS.out
            .splitCsv(sep: '\t', header: true)
            .map { row -> row.efo_trait }
            .unique { it.toLowerCase() }
    }
    else {
        log.info "Fetching GWAS Catalog associations for: ${terms.join(', ')}"
        traits_ch = channel.fromList(terms)
    }

    // a matched trait may have no associations; only an exact name must have some
    FETCH_GWAS(traits_ch, match == 'contains')
    MERGE_GWAS(FETCH_GWAS.out.collect())

    // look up LD proxies for every SNP, in chunks so an interrupted run can resume
    EXTRACT_SNPS(MERGE_GWAS.out, ld_pop, ld_build)
    snps_ch = EXTRACT_SNPS.out.snps.splitCsv(sep: '\t', header: true).map { row -> row.query_snp }
    if (params.ld_max_snps) {
        snps_ch = snps_ch.take(params.ld_max_snps.toString().toInteger())
    }
    LD_PROXY(snps_ch.collate(params.ld_chunk_size.toString().toInteger()), ld_pop, ld_build)
    MERGE_LD(LD_PROXY.out.proxies.collect(), LD_PROXY.out.failed.collect(), EXTRACT_SNPS.out.snps, ld_r2_min)
}
