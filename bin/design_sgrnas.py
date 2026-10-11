#!/usr/bin/env python3
"""Design CRISPRi sgRNAs in the ATAC-seq peaks that hold a trait SNP.

Candidate guides are found and scored with FlashFry, filtered, then the best
ones are picked so that they sit evenly across each peak.
"""
import argparse
import re
import subprocess
import sys

import pandas as pd
from pyfaidx import Fasta

# the nuclease profiles FlashFry 2.0 can index, with the spacer length and
# whether the PAM sits 3' of the spacer (SpCas9-like) or 5' of it (Cas12a-like).
# spcas9_scores says whether the Doench/Hsu models apply: they were trained on
# SpCas9 with an NGG PAM, and FlashFry drops them for any other nuclease.
PROFILES = {
    "spcas9-ngg-20": dict(spacer=20, pam="NGG", pam_3prime=True, spcas9_scores=True),
    "spcas9-ngg-19": dict(spacer=19, pam="NGG", pam_3prime=True, spcas9_scores=True),
    "spcas9-nag-20": dict(spacer=20, pam="NAG", pam_3prime=True, spcas9_scores=True),
    "spcas9-ngg-nag-20": dict(spacer=20, pam="NGG/NAG", pam_3prime=True, spcas9_scores=True),
    "spcas9-ngg-nag-19": dict(spacer=19, pam="NGG/NAG", pam_3prime=True, spcas9_scores=True),
    "spg-ngn-20": dict(spacer=20, pam="NGN", pam_3prime=True, spcas9_scores=True),
    "sacas9-nngrrt-21": dict(spacer=21, pam="NNGRRT", pam_3prime=True, spcas9_scores=False),
    "sacas9-nngrrt-22": dict(spacer=22, pam="NNGRRT", pam_3prime=True, spcas9_scores=False),
    "cpf1-tttn-20": dict(spacer=20, pam="TTTN", pam_3prime=False, spcas9_scores=False),
    "ascas12a-tttv-23": dict(spacer=23, pam="TTTV", pam_3prime=False, spcas9_scores=False),
    "lbcas12a-tttv-23": dict(spacer=23, pam="TTTV", pam_3prime=False, spcas9_scores=False),
    "encas12a-tttn-23": dict(spacer=23, pam="TTTN", pam_3prime=False, spcas9_scores=False),
    "cas12a-2xnls-tttn-23": dict(spacer=23, pam="TTTN", pam_3prime=False, spcas9_scores=False),
}

COLUMNS = ["peak_chr", "peak_start", "peak_end", "peak_name", "bin", "selection",
           "guide_chr", "guide_start", "guide_end", "strand", "spacer", "pam",
           "efficiency_doench2014", "specificity_cfd", "specificity_hsu2013",
           "mismatches_0_1_2_3_4", "overlapping_snps", "bp_to_nearest_snp", "enzyme"]

# every peak gets a row here, so a peak with few or no guides can be explained
FUNNEL = ["candidates", "inside_peak", "unique_in_genome", "no_long_run", "gc_ok",
          "efficiency_known", "few_near_offtargets", "mit_score_ok", "specificity_ok", "selected"]

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--peaks", required=True, help="BED of the peaks to target")
parser.add_argument("--fasta", required=True, help="genome FASTA the peaks are on")
parser.add_argument("--database", required=True, help="FlashFry off-target database")
parser.add_argument("--jar", required=True, help="FlashFry jar")
parser.add_argument("--enzyme", required=True, choices=sorted(PROFILES))
parser.add_argument("--per-peak", type=int, required=True, help="guides wanted per peak")
parser.add_argument("--max-run", type=int, default=3,
                    help="longest run of one base allowed in the spacer, e.g. 3 rejects AAAA")
parser.add_argument("--max-1mm", type=int, default=0,
                    help="off-target sites with 1 mismatch allowed")
parser.add_argument("--max-2mm", type=int, default=0,
                    help="off-target sites with 2 mismatches allowed")
parser.add_argument("--rank", choices=["efficiency", "uniqueness"], default="efficiency",
                    help="how to order the guides that pass the filters. 'uniqueness' is the only "
                         "option for nucleases without an efficiency model, and using it for "
                         "SpCas9 too makes the two sets comparable")
parser.add_argument("--min-hsu", type=float, default=50.0,
                    help="lowest MIT (Hsu 2013) specificity score to accept, 0 to not use it")
parser.add_argument("--min-specificity", type=float, default=0.0,
                    help="lowest CFD specificity to accept, 0 to not use it as a filter")
parser.add_argument("--snps", help="BED of the trait SNPs, to flag guides that cover one")
parser.add_argument("--out", required=True)
parser.add_argument("--out-bed", required=True)
parser.add_argument("--out-diagnostics", required=True)
args = parser.parse_args()

profile = PROFILES[args.enzyme]
# Doench2014 reads 4 bases before and 3 after the target, so peaks are cut out
# with a margin; guides are still only kept if they lie inside the peak itself
MARGIN = 30
names = ["peak_chr", "peak_start", "peak_end", "peak_name"]
peaks = pd.read_csv(args.peaks, sep="\t", header=None, comment="#",
                    usecols=[0, 1, 2, 3], names=names, dtype={0: str, 3: str}) \
    if sum(1 for _ in open(args.peaks)) else pd.DataFrame(columns=names)
counts = {peak: dict.fromkeys(FUNNEL, 0) for peak in peaks.index}


def write_diagnostics():
    table = pd.DataFrame.from_dict(counts, orient="index").reindex(columns=FUNNEL).fillna(0)
    table = peaks.join(table).assign(width=peaks["peak_end"] - peaks["peak_start"])
    table.to_csv(args.out_diagnostics, sep="\t", index=False)


def finish(reason):
    """Write the outputs as they stand and stop; used when nothing is left to do."""
    pd.DataFrame(columns=COLUMNS).to_csv(args.out, sep="\t", index=False)
    open(args.out_bed, "w").close()
    write_diagnostics()
    print(reason)
    sys.exit(0)


if peaks.empty:
    finish("No peaks to design guides for")

genome = Fasta(args.fasta, as_raw=True, sequence_always_upper=True)


def contig(chrom):
    """Match 'chr1' in the peaks to '1' in the genome, or the other way round."""
    for name in (chrom, chrom.removeprefix("chr"), f"chr{chrom}"):
        if name in genome:
            return name
    return None


# one FASTA record per peak, named by its row so the results can be matched back
with open("peaks.fa", "w") as fh:
    kept = 0
    for row in peaks.itertuples():
        name = contig(row.peak_chr)
        if name is None:
            print(f"{row.peak_chr} is not in the genome, skipping {row.peak_name}", file=sys.stderr)
            continue
        start = max(0, row.peak_start - MARGIN)
        end = min(len(genome[name]), row.peak_end + MARGIN)
        fh.write(f">{row.Index}::{start}\n{genome[name][start:end]}\n")
        kept += 1
if not kept:
    finish("No peaks could be read from the genome")


def flashfry(command, *arguments):
    subprocess.run(["java", "-Xmx4g", "-jar", args.jar, command, *arguments], check=True)


flashfry("discover", "--database", args.database, "--fasta", "peaks.fa", "--output", "candidates.txt")
metrics = ["dangerous", "minot"]
if profile["spcas9_scores"]:
    metrics = ["doench2014ontarget", "doench2016cfd", "hsu2013", *metrics, "rank"]
flashfry("score", "--input", "candidates.txt", "--output", "scored.txt",
         "--database", args.database, "--scoringMetrics", ",".join(metrics))

guides = pd.read_csv("scored.txt", sep="\t", dtype=str)
if guides.empty:
    finish(f"FlashFry found no {args.enzyme} target site in any peak")

# FlashFry reports positions inside the sequence it was given, so put them back
# on the chromosome. A guide keeps the peak's row number in its contig name.
guides["peak"] = guides["contig"].str.split("::").str[0].astype(int)
offset = guides["contig"].str.split("::").str[1].astype(int)
guides["guide_start"] = offset + guides["start"].astype(int)
guides["guide_end"] = offset + guides["stop"].astype(int)
guides["strand"] = guides["orientation"].map({"FWD": "+", "RVS": "-"}).fillna(".")
pam_length = guides["target"].str.len() - profile["spacer"]
guides["spacer"] = [t[p:] if not profile["pam_3prime"] else t[:-p or None]
                    for t, p in zip(guides["target"], pam_length)]
guides["pam"] = [t[:p] if not profile["pam_3prime"] else t[-p:] if p else ""
                 for t, p in zip(guides["target"], pam_length)]

for column in ("Doench2014OnTarget", "DoenchCFD_specificityscore", "Hsu2013"):
    guides[column] = pd.to_numeric(guides.get(column), errors="coerce")
if "0-1-2-3-4_mismatch" not in guides:
    guides["0-1-2-3-4_mismatch"] = ""
mismatches = guides["0-1-2-3-4_mismatch"].str.split(",", expand=True).apply(pd.to_numeric, errors="coerce")
mismatches.columns = [f"mm{i}" for i in range(mismatches.shape[1])]
guides = pd.concat([guides, mismatches], axis=1)
guides = guides.merge(peaks, left_on="peak", right_index=True)


def record(stage, kept_guides):
    """Note how many guides each peak still has at this point of the filtering."""
    per_peak = kept_guides.groupby("peak").size()
    for peak in counts:
        counts[peak][stage] = int(per_peak.get(peak, 0))
    return kept_guides


guides = record("candidates", guides)
# the margin was only there for scoring: keep guides that sit inside the peak
guides = record("inside_peak", guides[(guides["guide_start"] >= guides["peak_start"])
                                      & (guides["guide_end"] <= guides["peak_end"])])
# a guide that matches more than one place in the genome targets all of them
unique = guides.get("dangerous_in_genome", pd.Series("IN_GENOME=1", index=guides.index))
guides = record("unique_in_genome", guides[unique == "IN_GENOME=1"])
# runs of one base make oligos hard to synthesise and reads hard to map, and a
# run of four or more T can cut transcription of the guide short
run = re.compile("|".join(f"(?:{base}{{{args.max_run + 1},}})" for base in "ACGT"))
guides = record("no_long_run", guides[~guides["spacer"].str.contains(run)])
if "dangerous_GC" in guides:
    guides = record("gc_ok", guides[guides["dangerous_GC"] == "NONE"])
else:
    guides = record("gc_ok", guides)

if profile["spcas9_scores"]:
    guides = record("efficiency_known", guides[guides["Doench2014OnTarget"].notna()])
else:
    guides = record("efficiency_known", guides)

# the specificity that matters for CRISPRi: dCas9 binds off-target sites that
# differ by only one or two bases, so those are counted rather than scored
near = (guides["mm1"].fillna(0) <= args.max_1mm) & (guides["mm2"].fillna(0) <= args.max_2mm)
guides = record("few_near_offtargets", guides[near])
# the MIT score (Hsu 2013) is the long-standing convention, where 50 is the
# usual cut-off. It rarely removes anything the off-target counts above have not
# already removed, but it does catch guides in very repetitive surroundings
if args.min_hsu > 0 and profile["spcas9_scores"]:
    guides = record("mit_score_ok", guides[guides["Hsu2013"] >= args.min_hsu])
else:
    guides = record("mit_score_ok", guides)
# the CFD specificity score is an optional extra filter: it adds up every
# predicted off-target, which drops guides in repeat-rich peaks very steeply
if args.min_specificity > 0:
    guides = record("specificity_ok", guides[guides["DoenchCFD_specificityscore"] >= args.min_specificity])
else:
    guides = record("specificity_ok", guides)

if guides.empty:
    finish("No guide passed the filters in any peak")

# an efficiency model exists for SpCas9 only, so ranking by uniqueness is both
# the fallback for other nucleases and the way to compare one nuclease to another
rank_by_efficiency = args.rank == "efficiency" and profile["spcas9_scores"]
if rank_by_efficiency:
    # efficiency decides the order, with specificity breaking ties
    order = guides.sort_values(["Doench2014OnTarget", "DoenchCFD_specificityscore"], ascending=False)
else:
    # the best guide is the one with the fewest near-matches elsewhere in the
    # genome; position breaks ties so the choice does not depend on row order
    order = guides.sort_values([*mismatches.columns[1:], "guide_start"], ascending=True)

snps = {}
if args.snps and sum(1 for _ in open(args.snps)):
    table = pd.read_csv(args.snps, sep="\t", header=None, usecols=[0, 1, 2, 3],
                        names=["chrom", "start", "end", "name"], dtype={0: str, 3: str})
    for chrom, rows in table.groupby("chrom"):
        snps[chrom] = list(zip(rows["start"], rows["end"], rows["name"]))


def snp_annotation(chrom, start, end):
    """rsIDs the guide covers, and how far away the nearest SNP is."""
    found = snps.get(chrom) or snps.get(f"chr{chrom}") or snps.get(chrom.removeprefix("chr")) or []
    covered = [name for s, e, name in found if s < end and start < e]
    distances = [0 if s < end and start < e else min(abs(s - end), abs(start - e)) for s, e, *_ in found]
    return "; ".join(sorted(set(covered))), min(distances, default="")


# spread the guides over each peak: split it into as many equal bins as guides
# wanted and take the best guide in each, then fill any empty bin with the next best
selected = []
for peak, candidates in order.groupby("peak", sort=True):
    row = peaks.loc[peak]
    width = max(row.peak_end - row.peak_start, 1)
    middle = (candidates["guide_start"] + candidates["guide_end"]) // 2
    bins = ((middle - row.peak_start) * args.per_peak // width).clip(0, args.per_peak - 1)
    candidates = candidates.assign(bin=bins)
    picked = candidates.drop_duplicates("bin")
    extra = candidates.drop(picked.index).head(args.per_peak - len(picked))
    chosen = pd.concat([picked.assign(selection="even"), extra.assign(selection="extra")])
    selected.append(chosen.sort_values(["bin", "guide_start"]))

guides = record("selected", pd.concat(selected))
annotation = [snp_annotation(c, s, e) for c, s, e
              in zip(guides["peak_chr"], guides["guide_start"], guides["guide_end"])]
guides["overlapping_snps"] = [a[0] for a in annotation]
guides["bp_to_nearest_snp"] = [a[1] for a in annotation]
guides["guide_chr"] = guides["peak_chr"]
guides["enzyme"] = args.enzyme
guides["bin"] = guides["bin"] + 1  # count bins from 1 in the output
guides = guides.rename(columns={"Doench2014OnTarget": "efficiency_doench2014",
                                "DoenchCFD_specificityscore": "specificity_cfd",
                                "Hsu2013": "specificity_hsu2013",
                                "0-1-2-3-4_mismatch": "mismatches_0_1_2_3_4"})
guides = guides.sort_values(["peak_chr", "peak_start", "bin"])
guides[COLUMNS].to_csv(args.out, sep="\t", index=False, float_format="%.4f")

bed = pd.DataFrame({
    "chrom": guides["peak_chr"],
    "start": guides["guide_start"],
    "end": guides["guide_end"],
    "name": guides["peak_name"] + "_sg" + guides["bin"].astype(str),
    # the score column holds the efficiency, or 0 where there is no model for it
    "score": (guides["efficiency_doench2014"].fillna(0) * 1000).round().astype(int),
    "strand": guides["strand"],
})
bed.to_csv(args.out_bed, sep="\t", index=False, header=False)
write_diagnostics()

short = sum(1 for peak in counts if counts[peak]["selected"] < args.per_peak)
print(f"{len(guides)} guides for {guides['peak_name'].nunique()} of {len(peaks)} peaks; "
      f"{(guides['overlapping_snps'] != '').sum()} cover a trait SNP")
if short:
    print(f"{short} peak(s) got fewer than {args.per_peak} guides; "
          f"the counts in {args.out_diagnostics} show which filter removed them", file=sys.stderr)
if not rank_by_efficiency:
    reason = ("FlashFry has no on-target efficiency model for " + args.enzyme
              if not profile["spcas9_scores"] else "you asked for --rank uniqueness")
    print(f"Note: {reason}, so the guides were ranked by genome uniqueness", file=sys.stderr)
