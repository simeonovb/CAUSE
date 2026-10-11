# GWAS SNPs and LD proxies pipeline

A Nextflow pipeline that:

1. downloads the SNP associations for the traits you choose from the [GWAS Catalog](https://www.ebi.ac.uk/gwas/), and
2. looks up, for each of those SNPs, the other SNPs that are commonly inherited with it (LD proxies), using [LDlink LDproxy](https://ldlink.nih.gov/?tab=ldproxy).

## What you need

- [Nextflow](https://www.nextflow.io/docs/latest/install.html) (tested with 25.10) and Java 17 or newer
- [Conda](https://docs.conda.io/en/latest/miniconda.html). The Python packages are installed automatically on the first run, which takes a few minutes.
- A personal LDlink token. Register at <https://ldlink.nih.gov/?tab=apiaccess>; the token is sent by email.

## Running it

```bash
nextflow run main.nf -resume
```

The pipeline asks for everything it needs before it starts work:

| Question | Options |
| --- | --- |
| How traits are matched | `exact` trait names, or any trait whose name `contains` the words you give |
| Traits or words | comma-separated, e.g. `Atrial fibrillation, Heart rate` or `heart, cardio` |
| LDlink token | asked on your first run only, then saved on your machine |
| Populations | `ALL`, or one or more 1000 Genomes codes such as `EUR` or `CEU,YRI`; the full list is shown |
| Genome assembly | `grch37`, `grch38` or `grch38_high_coverage` |
| R2 threshold | a number from 0 to 1; proxies below it are left out (0 keeps everything LDlink returns) |

Exact trait names must match the GWAS Catalog's name for the trait. If a name is not found, the error lists similar names.

To skip the questions, pass the answers as flags:

```bash
nextflow run main.nf -resume --traits "heart,cardio" --match contains --ld_pop EUR --ld_build grch38 --ld_r2_min 0.8
```

The token cannot be passed as a flag. To set or change it without being asked:

```bash
nextflow secrets set LDLINK_TOKEN <your token>
```

### Other options

| Flag | Default | Meaning |
| --- | --- | --- |
| `--ld_window` | `500000` | base pairs searched either side of each SNP |
| `--ld_max_snps` | all | only look up the first N SNPs, for a quick test |
| `--ld_chunk_size` | `25` | SNPs per LDlink task |
| `--outdir` | `results` | where outputs are written |

## How long it takes

LDlink allows one request at a time per token, and each SNP takes roughly 7 to 9 seconds. A thousand SNPs is about two hours. Use `--ld_max_snps 20` to try the pipeline out first.

Always run with `-resume`: finished work is reused, so an interrupted run continues where it stopped. To pull fresh data from the GWAS Catalog, run once without `-resume`.

## Outputs

In `results/gwas_catalog/`:

- `gwas_catalog_traits.tsv`: every association found, one row each. `query_trait` lists the trait or traits that returned it.
- `matched_traits.tsv`: in `contains` mode, the traits each word matched.
- `per_trait/`: one table per trait. Files from earlier runs are not removed.

In `results/ld_proxy/`:

- `ld_proxies.tsv`: each query SNP with its traits, the population and assembly used, and its proxies with LDlink's columns (position, alleles, MAF, D', R2 and others).
- `query_snps.tsv`: the SNPs that were looked up, with their traits.
- `ld_failed.tsv`: SNPs LDlink could not look up, with its reason.
- `skipped_variants.tsv`: variants the GWAS Catalog lists by position rather than rsID. These are not looked up, because their genome assembly is not stated.
