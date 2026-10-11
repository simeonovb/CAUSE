# CAUSE — Chromatin-Accessible Uncovering of SNP Effects

A Nextflow pipeline that goes from a trait to sgRNAs for its regulatory variants:

1. **Collect SNPs** for the traits you choose from the [GWAS Catalog](https://www.ebi.ac.uk/gwas/).
2. **Expand them** with the SNPs commonly inherited alongside (LD proxies), using [LDlink](https://ldlink.nih.gov/?tab=ldproxy).
3. **Keep the open chromatin** — the ATAC-seq peaks that contain at least one of those SNPs.
4. **Design CRISPRi sgRNAs** spread across each of those peaks, with [FlashFry](https://github.com/mckennalab/FlashFry) (MIT-licensed).

## Quick start

```bash
nextflow run main.nf -resume
```

The pipeline asks everything it needs up front, then runs without further input. To skip the questions, pass the answers instead:

```bash
nextflow run main.nf -resume \
  --traits "heart failure" --match exact \
  --ld_pop EUR --ld_build grch37 --ld_r2_min 0.8 \
  --atac_peaks "ATAC_brain_peaks.bed,ATAC_skin_peaks.bed" \
  --sgrna_enzyme spcas9-ngg-20 --sgrna_per_peak 4 --keep_genome true
```

Always use `-resume`: finished work is reused, so an interrupted run carries on where it stopped.

## What you need

- [Nextflow](https://www.nextflow.io/docs/latest/install.html) (tested with 25.10) and Java 17 or newer.
- [Conda](https://docs.conda.io/en/latest/miniconda.html). Packages install themselves on the first run. This version runs on conda only; see [containers](#containers) if you need Docker or Singularity.
- An LDlink token, free from <https://ldlink.nih.gov/?tab=apiaccess>. You are asked for it once and it is saved. Tokens are personal, so everyone running the pipeline registers their own.
- Your ATAC-seq peaks as BED files (plain or gzipped; MACS2 narrowPeak works too).
- About 25 GB of free disk for the sgRNA step, which downloads and indexes a genome.

## The questions

| Question | Options |
| --- | --- |
| How traits are matched | `exact` trait names, or any trait name that `contains` your words |
| Traits or words | comma-separated, e.g. `Atrial fibrillation, Heart rate` or `heart, cardio` |
| LDlink token | first run only |
| Populations | `ALL`, or 1000 Genomes codes such as `EUR` or `CEU,YRI`; the full list is shown |
| Genome assembly | `grch37`, `grch38` or `grch38_high_coverage` |
| R2 threshold | 0 to 1; proxies below it are dropped (0 keeps everything LDlink returns) |
| ATAC-seq peak files | one or more paths, comma-separated, e.g. one per cell type |
| Nuclease | SpCas9 or SaCas9 |
| sgRNAs per peak | a whole number, spread evenly across each peak |
| Keep the genome? | keep it for next time, or delete it when the run finishes |

Two things the pipeline cannot check for you:

- **Your peaks must be on the same assembly you pick for the LD lookup.** Mixing hg19 peaks with GRCh38 coordinates gives wrong overlaps and no error.
- **Exact trait names must match the GWAS Catalog's spelling.** If one is not found, the error suggests close names.

To set the token without being asked:

```bash
nextflow secrets set LDLINK_TOKEN <your token>
```

## Outputs

Everything lands in `results/`. For the peak and sgRNA steps there is one set of files per peak file, named after it, so `PM_peaks.bed` gives `PM_peaks.*`.

**`gwas_catalog/`** — what the GWAS Catalog returned

| File | Contents |
| --- | --- |
| `gwas_catalog_traits.tsv` | every association, one row each. `discovery_ancestry` says who was studied |
| `matched_traits.tsv` | in `contains` mode, which traits each word matched |
| `per_trait/` | one table per trait (files from earlier runs are kept) |

**`ld_proxy/`** — the SNPs to target

| File | Contents |
| --- | --- |
| `ld_proxies.tsv` | each query SNP and its proxies, with LDlink's columns. A SNP is listed as its own proxy at R2 1 |
| `ld_snps.bed` | the same SNPs as BED, for a browser or `bedtools`. Score is the highest R2 × 1000 |
| `query_snps.tsv` | the SNPs that were looked up |
| `ld_failed.tsv` | SNPs LDlink could not look up, with its reason |
| `skipped_variants.tsv` | variants given by position instead of an rsID, so their assembly is unknown |

**`atac_overlap/`** — the peaks worth targeting

| File | Contents |
| --- | --- |
| `<name>.peaks_with_snps.bed` | peaks containing at least one SNP, copied unchanged from your file |
| `<name>.snps_in_peaks.tsv` | which SNPs sit in each peak, one row per peak and SNP |

**`pipeline_info/`** — what produced this

| File | Contents |
| --- | --- |
| `run_<name>.txt` | the pipeline version, the Nextflow version, the git commit it ran from and whether anything was uncommitted, the command line, the settings as you answered them, and every parameter including defaults |

Many settings change the guides you get — the LD population, the R2 threshold, the nuclease, each filter — so this file is what lets you tell two result sets apart later. One is written per run, named after the run.

**`sgrna/`** — the guides

| File | Contents |
| --- | --- |
| `<name>.sgrnas.tsv` | the chosen guides: peak, bin, spacer, PAM, scores, and any trait SNP they cover |
| `<name>.sgrnas.bed` | the same guides as BED, with efficiency in the score column |
| `<name>.sgrna_diagnostics.tsv` | per peak, how many candidates survived each filter — read this when a peak gets few guides |

### Two columns people mix up

- **`gwas_ancestry`** — who the SNP was *found* in, from its GWAS studies.
- **`ld_population`** — the 1000 Genomes population *you* chose for the LD lookup.

Proxies are trustworthy when these agree. A SNP from a Korean study, with proxies computed in Europeans, may get proxies that are not inherited together in the people actually studied.

## How guides are chosen

Only peaks holding a trait SNP get guides. Each peak is split into as many equal bins as you asked for guides, and the best guide in each bin is taken, so they spread across the element rather than clustering. The `selection` column says `even` (one per bin) or `extra` (a bin had none, so the next best guide was used).

A guide must pass all of these:

| Filter | Why |
| --- | --- |
| sits fully inside the peak | guides are scored with a margin, but only kept if they target the element |
| unique in the genome | a second exact match silences another place too |
| no run of 4+ identical bases | hard to synthesise, and awkward to map in sequencing |
| GC not extreme | FlashFry's own warning |
| no off-target site differing by 1 or 2 bases | this is the off-target risk that matters for CRISPRi |
| MIT (Hsu 2013) score ≥ 50 | the familiar convention, kept as a safety net |

Survivors are ranked by Doench 2014 on-target efficiency, with CFD specificity breaking ties. All three scores are reported so you can re-sort as you like.

### Why the MIT score is not enough here

If you have designed guides before, you probably used MIT > 50 and nothing else. That works for cutting, but not for silencing.

**The MIT score is an average.** It squeezes every off-target site into one number. Most sites differ by 3 or 4 bases and are harmless, and because there are hundreds of them they dominate the number — so one genuinely dangerous site disappears into the crowd. It is a hotel's 4.3 average rating hiding the review about the broken lock.

**Cutting tolerates that; silencing does not.** An extra cut gets repaired. A dCas9-KRAB parked on a near-identical site switches off a second gene for the whole experiment, and one such site is enough to confuse the result.

Measured on a peak in this project: one guide had a site elsewhere in the genome matching with a **single base different**, and its MIT score was **55.8** — a comfortable pass. Nine guides scoring above 50 had a 1- or 2-mismatch site. So the pipeline counts those sites directly instead of trusting the summary, and keeps MIT ≥ 50 as a backstop.

CFD (Doench 2016) is the better of FlashFry's two scores, since it accounts for which bases are mismatched and not only where. But it also sums every predicted off-target, including the many weak false positives, so in repeat-rich peaks it collapses towards zero and would discard every guide. It is reported, not used as a gate; `--sgrna_min_specificity` turns it into one if you want.

### SaCas9 has no efficiency score

The Doench and Hsu models were trained on SpCas9 with an NGG PAM, and FlashFry will not apply them to anything else. With SaCas9 the efficiency and specificity columns are empty, the MIT filter does not apply, and guides are ranked by genome uniqueness alone (`mismatches_0_1_2_3_4`). The run says so when it happens.

The menu offers SpCas9 and SaCas9. FlashFry's other profiles (`spg-ngn-20`, the Cas12a ones) still work if passed with `--sgrna_enzyme`, and behave like SaCas9 in this respect.

### Guides on top of the SNP

A guide covering a trait SNP is flagged in `overlapping_snps`, with `bp_to_nearest_snp` giving the distance to the closest one. Whether you want that is your call: if your cell line carries the other allele, the guide binds less well.

## Licence

This pipeline is released under the MIT licence; see `LICENSE`. The tools it calls carry their own terms, listed under [credits](#credits).

## Version

The pipeline's own version is in the `manifest` block of `nextflow.config`, and `nextflow run main.nf -version`-style metadata ends up in each run's summary file. Tag releases in git with the same number, so a result folder can always be traced back to the code that made it.

## The genome and its index

The sgRNA step needs the genome your peaks are on, and an index FlashFry builds from it. Together they are the slow, bulky part of the pipeline: the genome is 1 GB compressed and 3 GB unpacked, the index is a few GB more, and building it needs around 15 GB of scratch space while it works.

The pipeline therefore **asks whether to keep them** when the run finishes. Keeping them (in `genome_cache/`) means later runs skip straight past both steps; deleting them frees the space but means building again next time. The index is tied to one assembly *and* one nuclease, so switching either one means a rebuild.

Nothing is downloaded if no peak contains a SNP — the genome step waits until there is something to design guides for.

The LD lookups are the other slow part, and they cannot be sped up by running more at once: LDlink states that "API access is limited to sequential requests only", so the pipeline sends one request at a time and waits for each result. Do not raise that limit.

How long they take depends on the trait and the population: a trait with thousands of SNPs looked up in `ALL` takes far longer than a small trait in one population. Use `--ld_max_snps 20` to try a large trait quickly.

LDlink also asks that you "contact support if you plan to make large volumes of API requests". A trait with thousands of SNPs means one request per SNP, so get in touch with them before running at that scale, and before several people start running the pipeline at once.

## Watching a run

The display shows one line per step with a count, such as `[4 of 67]`, so you can see how far the LD chunks or peak files have got.

Two steps cannot be followed that way, because each is a single task that runs for a long time: the genome download and the FlashFry index build. Those stream their own progress to the terminal, so you can see the download grow and the index work through the chromosomes and then its 16,384 bins. `--progress false` silences them.

The LD lookups also name each SNP as they go, but that is not streamed, since the step counter already shows how many chunks are done. To follow any step in detail, take its directory from the display and read its log:

```bash
tail -f work/95/dd85bf*/.command.log
```

Streamed output and the live table do overwrite each other a little. If that bothers you, `-ansi-log false` replaces the table with one line per task.

Afterwards, `nextflow log` lists runs with durations, and `-with-report report.html` breaks down time and memory per step.

## Other options

| Flag | Default | Meaning |
| --- | --- | --- |
| `--ld_window` | `500000` | bp searched either side of each SNP |
| `--ld_max_snps` | all | only look up the first N SNPs |
| `--ld_chunk_size` | `25` | SNPs per LDlink task |
| `--outdir` | `results` | where outputs go |
| `--sgrna_max_run` | `3` | longest run of one base allowed; `3` rejects AAAA, TTTT, CCCC, GGGG |
| `--sgrna_max_1mm` | `0` | 1-mismatch off-target sites allowed |
| `--sgrna_max_2mm` | `0` | 2-mismatch off-target sites allowed |
| `--sgrna_min_hsu` | `50` | lowest MIT score; `0` turns it off |
| `--sgrna_min_specificity` | `0.0` | optional CFD specificity floor; `0` turns it off |
| `--keep_genome` | asks | `true` keeps the genome and index, `false` deletes them at the end |
| `--genome_dir` | `genome_cache` | where the genome and index live |
| `--flashfry_memory` | `8g` | memory for the index build |
| `--progress` | `true` | stream the genome download and index build to the terminal |

## If you get no results

An empty overlap or guide table is usually real, not a bug:

- **No peaks contain a SNP.** Check that the LD population suits the trait's `gwas_ancestry`, and that peaks and `ld_build` use the same assembly. A trait with few SNPs often simply misses your peaks.
- **A peak got fewer guides than asked.** Open `<name>.sgrna_diagnostics.tsv`: it shows how many candidates each filter removed.
## Credits

This pipeline is a wrapper. The work it depends on belongs to other people, and if you publish results from it, please cite them.

**[GWAS Catalog](https://www.ebi.ac.uk/gwas/)** — where every SNP in this pipeline comes from, including the trait names, p-values and study ancestries.

> Cerezo M, Sollis E, Ji Y, et al. The NHGRI-EBI GWAS Catalog: standards for reusability, sustainability and diversity. *Nucleic Acids Research* 2025;53(D1):D998–D1005. doi:[10.1093/nar/gkae1070](https://doi.org/10.1093/nar/gkae1070)

The Catalog asks that you also cite the original study behind an association, and note the date you downloaded the data. Its own guidance is at <https://www.ebi.ac.uk/gwas/docs/citation>.

**[LDlink](https://ldlink.nih.gov/)** — the LDproxy module computes every proxy SNP, using 1000 Genomes as its reference panel. Each run needs your own free token. Its API terms say access "is limited to sequential requests only", which this pipeline respects, and ask that you "contact support if you plan to make large volumes of API requests".

> Machiela MJ, Chanock SJ. LDlink: a web-based application for exploring population-specific haplotype structure and linking correlated alleles of possible functional variants. *Bioinformatics* 2015;31(21):3555–3557. doi:[10.1093/bioinformatics/btv402](https://doi.org/10.1093/bioinformatics/btv402)

**[FlashFry](https://github.com/mckennalab/FlashFry)** — finds and scores every candidate guide. MIT licensed, which is why this pipeline can use it freely.

> McKenna A, Shendure J. FlashFry: a fast and flexible tool for large-scale CRISPR target design. *BMC Biology* 2018;16:74. doi:[10.1186/s12915-018-0545-0](https://doi.org/10.1186/s12915-018-0545-0)

### Also used

- **Genome assemblies** are downloaded from the [UCSC Genome Browser](https://genome.ucsc.edu/), which distributes the hg19 and hg38 sequences.
- **Reference panel** for all LD calculations is the [1000 Genomes Project](https://www.internationalgenome.org/), through LDlink.
- **Guide scores** reported in the sgRNA tables come from Doench et al., *Nat Biotechnol* 2014;32(12):1262–1267 (on-target efficiency), Doench et al., *Nat Biotechnol* 2016;34(2):184–191 (the CFD off-target score), and Hsu et al., *Nat Biotechnol* 2013;31(9):827–832 (the MIT specificity score). FlashFry implements all three.
