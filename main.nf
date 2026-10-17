// the settings as answered, filled in by the workflow and read again when the run
// finishes. A map is used because the handler at the end does not see plain
// variables set inside the workflow
run_settings = [keep_genome: false]

process RESOLVE_TRAITS {
    tag "${keyword}"
    conda "${projectDir}/envs/gwas.yml"
    container "cause-nf/gwas:${workflow.manifest.version}"

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
    container "cause-nf/gwas:${workflow.manifest.version}"
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
    container "cause-nf/gwas:${workflow.manifest.version}"
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
    container "cause-nf/gwas:${workflow.manifest.version}"
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
    container "cause-nf/gwas:${workflow.manifest.version}"
    secret 'LDLINK_TOKEN'

    // retry transient API failures, but not exit status 3 (bad token or settings)
    // or a task killed because the run was stopped (129, 130, 143): retrying that
    // leaves a lookup running after the pipeline has gone, holding the token
    errorStrategy { task.exitStatus in [3, 129, 130, 143] ? 'terminate' : 'retry' }
    maxRetries 3
    // LDlink's API terms: "API access is limited to sequential requests only",
    // so one request at a time per token. Do not raise this.
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
    container "cause-nf/gwas:${workflow.manifest.version}"
    publishDir "${params.outdir}/ld_proxy", mode: 'copy', overwrite: true

    input:
    path chunks, stageAs: "chunks/?/*"
    path failed, stageAs: "failures/?/*"
    path snps
    val r2_min

    output:
    path "ld_proxies.tsv", emit: proxies
    path "ld_failed.tsv", emit: failed
    path "ld_snps.bed", emit: bed

    script:
    """
    merge_ld.py --chunks chunks --failures failures --snps ${snps} --r2-min ${r2_min} \\
        --out ld_proxies.tsv --failed-out ld_failed.tsv --bed-out ld_snps.bed
    """
}

process OVERLAP_PEAKS {
    tag "${name}"
    conda "${projectDir}/envs/gwas.yml"
    container "cause-nf/gwas:${workflow.manifest.version}"
    publishDir "${params.outdir}/atac_overlap", mode: 'copy', overwrite: true

    input:
    // staged under a fixed name, since peak file names can hold spaces or brackets
    tuple val(name), path(peaks, stageAs: 'peaks.bed')
    path proxies

    output:
    path "*.peaks_with_snps.bed", emit: peaks
    path "*.snps_in_peaks.tsv", emit: snps

    script:
    """
    overlap_peaks.py --peaks ${peaks} --snps ${proxies} \\
        --out-peaks ${name}.peaks_with_snps.bed --out-snps ${name}.snps_in_peaks.tsv
    """
}

process GET_FLASHFRY {
    conda "${projectDir}/envs/sgrna.yml"
    container "cause-nf/sgrna:${workflow.manifest.version}"
    // kept next to the genome, so it is downloaded once
    storeDir "${params.genome_dir}"

    output:
    path "flashfry.jar"

    script:
    """
    curl -sSL --retry 3 -o flashfry.jar ${params.flashfry_url}
    """
}

process DOWNLOAD_GENOME {
    tag "${build}"
    conda "${projectDir}/envs/sgrna.yml"
    container "cause-nf/sgrna:${workflow.manifest.version}"
    storeDir "${params.genome_dir}"
    debug params.progress

    input:
    val build

    output:
    tuple val(build), path("${build}.fa")

    script:
    // the assembly chosen for the LD lookup decides the genome, so guides are
    // always designed on the same coordinates the peaks and SNPs are on
    def ucsc = [grch37: 'hg19', grch38: 'hg38', grch38_high_coverage: 'hg38'][build]
    """
    echo "Downloading the ${ucsc} genome for ${build} (about 1 GB)"
    curl -sSL --retry 3 -o genome.fa.gz ${params.genome_url_base}/${ucsc}/bigZips/${ucsc}.fa.gz &
    download=\$!
    while kill -0 \$download 2>/dev/null; do
        sleep 60
        echo "  ...\$(du -h genome.fa.gz 2>/dev/null | cut -f1) downloaded"
    done
    wait \$download

    echo "Unpacking it (about 3 GB on disk)"
    gunzip -c genome.fa.gz > ${build}.fa
    rm genome.fa.gz
    """
}

process FLASHFRY_INDEX {
    tag "${build} ${enzyme}"
    conda "${projectDir}/envs/sgrna.yml"
    container "cause-nf/sgrna:${workflow.manifest.version}"
    storeDir "${params.genome_dir}"
    debug params.progress

    input:
    tuple val(build), path(fasta)
    path jar
    val enzyme

    output:
    tuple path("${build}.${enzyme}.db"), path("${build}.${enzyme}.db.header")

    script:
    """
    mkdir -p flashfry_tmp
    echo "Building the FlashFry index for ${build}, ${enzyme}."
    echo "This takes about half an hour and 20 GB of scratch space. It writes 16384 bins,"
    echo "and reports every thousandth one, so the count below is the progress."
    # FlashFry logs its progress on stderr, which is merged in so it can be seen live
    java -Xmx${params.flashfry_memory} -jar ${jar} index --reference ${fasta} \\
        --database ${build}.${enzyme}.db --tmpLocation flashfry_tmp --enzyme ${enzyme} 2>&1
    rm -rf flashfry_tmp
    """
}

process DESIGN_SGRNAS {
    tag "${name}"
    conda "${projectDir}/envs/sgrna.yml"
    container "cause-nf/sgrna:${workflow.manifest.version}"
    publishDir "${params.outdir}/sgrna", mode: 'copy', overwrite: true

    input:
    tuple val(name), path(peaks)
    tuple val(build), path(fasta)
    tuple path(database), path(header)
    path jar
    val enzyme
    val per_peak
    path snps

    output:
    path "${name}.sgrnas.tsv", emit: guides
    path "${name}.sgrnas.bed", emit: bed
    path "${name}.sgrna_diagnostics.tsv", emit: diagnostics

    script:
    """
    design_sgrnas.py --peaks ${peaks} --fasta ${fasta} --database ${database} --jar ${jar} \\
        --enzyme ${enzyme} --per-peak ${per_peak} \\
        --max-run ${params.sgrna_max_run} --max-1mm ${params.sgrna_max_1mm} \\
        --max-2mm ${params.sgrna_max_2mm} --min-hsu ${params.sgrna_min_hsu} \\
        --rank ${params.sgrna_rank} \\
        --min-specificity ${params.sgrna_min_specificity} \\
        --snps ${snps} --out ${name}.sgrnas.tsv --out-bed ${name}.sgrnas.bed \\
        --out-diagnostics ${name}.sgrna_diagnostics.tsv
    """
}

// ask a question in the terminal; used for any setting not given as a flag
def ask(String question) {
    def console = System.console()
    if (!console) {
        error "Some settings are missing and there is no interactive terminal to ask in. Pass them as flags, e.g. --traits \"heart,cardiac\" --match contains --ld_pop EUR --ld_build grch38 --ld_r2_min 0.8 --atac_peaks peaks.bed"
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

def askPeaks(build) {
    return ask("Which ATAC-seq peak file(s) should the SNPs be overlapped with?\nBED format, on the same assembly as the LD lookup (${build}).\nType one or more paths, comma-separated: ")
}

// turn "PM.bed, SN.bed" into a list of [name, file], checking each file exists
def parsePeaks(answer) {
    def paths = (answer instanceof List ? answer : answer.toString().tokenize(','))
        .collect { it.toString().trim().replaceAll(/^["']|["']$/, '') }
        .findAll { it }
    if (!paths) {
        error "No ATAC-seq peak file given"
    }
    def peaks = paths.collect { path ->
        def bed = new File(path.replaceFirst(/^~/, System.getProperty('user.home')))
        if (!bed.isFile()) {
            error "ATAC-seq peak file not found: ${path}"
        }
        // outputs are named after the file, e.g. PM_peaks.bed gives PM_peaks
        def name = bed.name.replaceFirst(/(?i)(\.bed|\.narrowPeak|\.broadPeak)?(\.gz)?$/, '').replaceAll(/[^\w.-]+/, '_')
        [name, bed.toPath().toAbsolutePath()]
    }
    def names = peaks.collect { it[0] }
    if (names.toUnique().size() < names.size()) {
        error "ATAC-seq peak files must have different file names, got: ${names.join(', ')}"
    }
    return peaks
}

// the nucleases offered in the menu; the other FlashFry profiles can still be
// typed in full, or passed with --sgrna_enzyme
def sgrnaEnzymes() {
    return [
        'spcas9-ngg-20': 'SpCas9, NGG PAM, 20 nt spacer. Scored for efficiency and specificity',
        'sacas9-nngrrt-21': 'SaCas9, NNGRRT PAM, 21 nt spacer. Ranked by genome uniqueness only',
    ]
}

def allSgrnaEnzymes() {
    return ['spcas9-ngg-20', 'spcas9-ngg-19', 'spcas9-nag-20', 'spcas9-ngg-nag-20',
            'spcas9-ngg-nag-19', 'spg-ngn-20', 'sacas9-nngrrt-21', 'sacas9-nngrrt-22',
            'cpf1-tttn-20', 'ascas12a-tttv-23', 'lbcas12a-tttv-23', 'encas12a-tttn-23',
            'cas12a-2xnls-tttn-23']
}

def askEnzyme() {
    def codes = sgrnaEnzymes().keySet() as List
    def menu = codes.withIndex().collect { code, i -> "  ${i + 1}) ${code}  (${sgrnaEnzymes()[code]})" }.join('\n')
    return ask("Which nuclease will the guides be used with? The PAM decides where guides can sit.\n${menu}\nChoose 1 or 2 [1]: ")
}

def parseEnzyme(answer) {
    def codes = sgrnaEnzymes().keySet() as List
    def text = answer.toString().trim().toLowerCase() ?: '1'
    def enzyme = text ==~ /\d+/ && text.toInteger() <= codes.size() ? codes[text.toInteger() - 1] : text
    if (!(enzyme in allSgrnaEnzymes())) {
        error "Unknown nuclease '${text}'. Choose from: ${allSgrnaEnzymes().join(', ')}"
    }
    return enzyme
}

def askKeepGenome() {
    return ask('The genome and its FlashFry index take about 6 GB, and building the index takes around half an hour.\n  1) keep them when the run finishes, so the next run needs neither\n  2) delete them when the run finishes, to save the space\nChoose 1 or 2 [1]: ')
}

def parseKeepGenome(answer) {
    if (answer instanceof Boolean) {
        return answer
    }
    def text = answer.toString().trim().toLowerCase() ?: '1'
    if (text in ['1', 'keep', 'k', 'yes', 'y', 'true']) {
        return true
    }
    if (text in ['2', 'delete', 'd', 'no', 'n', 'false']) {
        return false
    }
    error "Answer 1 to keep the genome or 2 to delete it, got '${text}'"
}

def askPerPeak() {
    return ask('How many sgRNAs per peak? They are spread evenly across the peak.\nType a whole number [4]: ')
}

def parsePerPeak(answer) {
    def text = answer.toString().trim() ?: '4'
    if (!(text ==~ /\d+/) || text.toInteger() < 1) {
        error "The number of sgRNAs per peak must be a whole number of 1 or more, got '${text}'"
    }
    return text.toInteger()
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

    // LD settings and peak files are asked up front too, so the run needs no further input
    ld_pop = parsePop(params.ld_pop ?: askPop())
    ld_build = parseBuild(params.ld_build ?: askBuild())
    ld_r2_min = parseR2(params.ld_r2_min != null ? params.ld_r2_min : askR2())
    peaks = parsePeaks(params.atac_peaks ?: askPeaks(ld_build))
    enzyme = parseEnzyme(params.sgrna_enzyme ?: askEnzyme())
    per_peak = parsePerPeak(params.sgrna_per_peak != null ? params.sgrna_per_peak : askPerPeak())
    run_settings.keep_genome = parseKeepGenome(params.keep_genome != null ? params.keep_genome : askKeepGenome())
    ensureLdToken()
    log.info "LD proxies: populations ${ld_pop}, assembly ${ld_build}, R2 >= ${ld_r2_min}"
    log.info "ATAC-seq peaks: ${peaks.collect { it[1].fileName }.join(', ')}"
    log.info "sgRNAs: ${per_peak} per peak, ${enzyme}, designed with FlashFry on ${ld_build}"
    log.info "The genome and its index will be ${run_settings.keep_genome ? 'kept in ' + params.genome_dir : 'deleted when the run finishes'}"

    // kept for the run summary written at the end
    run_settings.putAll([match: match, traits: terms.join('; '), ld_population: ld_pop,
                         genome_build: ld_build, ld_r2_min: ld_r2_min,
                         atac_peaks: peaks.collect { it[1] }.join('; '),
                         nuclease: enzyme, sgrnas_per_peak: per_peak])

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

    // keep the peaks holding at least one query SNP or proxy, one task per peak file
    OVERLAP_PEAKS(channel.fromList(peaks), MERGE_LD.out.proxies)

    // design CRISPRi guides in those peaks. The genome is downloaded and indexed
    // once per assembly and nuclease, then removed again unless --keep_genome
    named_peaks = OVERLAP_PEAKS.out.peaks
        .filter { bed -> bed.size() > 0 }
        .map { bed -> tuple(bed.name.replace('.peaks_with_snps.bed', ''), bed) }

    // the genome is only worth fetching and indexing if some peak holds a SNP.
    // These come out as value channels, so every peak file can use the same ones
    needed = named_peaks.map { ld_build }.first()
    jar = GET_FLASHFRY()
    genome = DOWNLOAD_GENOME(needed)
    database = FLASHFRY_INDEX(genome, jar, enzyme)
    DESIGN_SGRNAS(named_peaks, genome, database, jar, enzyme, per_peak, MERGE_LD.out.bed)
}

// the genome and its index are large, so they go once the guides are designed
workflow.onComplete {
    // a record of what produced these results, since many settings change them
    def report = file("${params.outdir}/pipeline_info")
    report.mkdirs()
    def git = { command ->
        try {
            def process = ["git", "-C", "${projectDir}", *command].execute()
            process.waitFor()
            return process.exitValue() == 0 ? process.text.trim() : null
        }
        catch (Exception ignored) {
            return null
        }
    }
    def lines = ["# ${workflow.manifest.name ?: 'pipeline'} run summary", ""]
    [
        'pipeline_version': workflow.manifest.version,
        'nextflow_version': nextflow.version,
        'git_commit'      : git(['rev-parse', '--short', 'HEAD']),
        'git_branch'      : git(['rev-parse', '--abbrev-ref', 'HEAD']),
        'git_uncommitted' : git(['status', '--porcelain'])?.with { it ? 'yes' : 'no' },
        'run_name'        : workflow.runName,
        'session_id'      : workflow.sessionId,
        'started'         : workflow.start,
        'completed'       : workflow.complete,
        'duration'        : workflow.duration,
        'succeeded'       : workflow.success,
        'exit_status'     : workflow.exitStatus,
        'command_line'    : workflow.commandLine,
        'work_dir'        : workflow.workDir,
    ].each { key, value -> lines << "${key}: ${value ?: '-'}" }

    lines << "" << "## settings as answered" << ""
    run_settings.sort().each { key, value -> lines << "${key}: ${value}" }

    lines << "" << "## every parameter, including defaults" << ""
    params.sort().each { key, value -> lines << "${key}: ${value}" }

    file("${report}/run_${workflow.runName}.txt").text = lines.join('\n') + '\n'

    if (run_settings.keep_genome || !workflow.success) {
        return
    }
    def cache = file(params.genome_dir)
    def large = cache.exists() ? cache.listFiles().findAll { it.name ==~ /.*\.(fa|fai|db|header)$/ } : []
    large.each { it.delete() }
    if (large) {
        log.info "Deleted the genome and its FlashFry index (${large.size()} files) to save space. " +
                 "Answer 1 to the genome question, or pass --keep_genome true, to keep them next time"
    }
}
