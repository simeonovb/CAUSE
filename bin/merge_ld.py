#!/usr/bin/env python3
"""Combine LDproxy chunks, keeping proxies at or above an R2 threshold."""
import argparse
import glob

import pandas as pd

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--chunks", required=True, help="directory of LDproxy chunk tables")
parser.add_argument("--failures", required=True, help="directory of failed-SNP tables")
parser.add_argument("--snps", required=True, help="query SNP table with their traits")
parser.add_argument("--r2-min", type=float, required=True)
parser.add_argument("--out", required=True)
parser.add_argument("--failed-out", required=True)
parser.add_argument("--bed-out", required=True, help="BED of every query SNP and proxy kept")
args = parser.parse_args()


def read_all(directory, pattern):
    files = sorted(glob.glob(f"{directory}/**/{pattern}", recursive=True))
    return pd.concat([pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False) for f in files])


proxies = read_all(args.chunks, "*.tsv.gz")
proxies = proxies[proxies["R2"].astype(float) >= args.r2_min]

# show each query SNP's traits next to it
snps = pd.read_csv(args.snps, sep="\t", dtype=str, keep_default_na=False)
proxies = snps.merge(proxies, on="query_snp", how="right")
proxies.to_csv(args.out, sep="\t", index=False)

# the same SNPs as a BED file, one line per SNP: BED starts count from 0, so a SNP
# at position p is p-1 to p. The score is its highest R2 with a query SNP, times 1000.
coord = proxies["Coord"].str.extract(r"^(?P<chrom>[^:]+):(?P<pos>\d+)$").dropna()
bed = pd.DataFrame({
    "chrom": coord["chrom"],
    "start": coord["pos"].astype(int) - 1,
    "end": coord["pos"].astype(int),
    # a few variants have no rsID; name those by position
    "name": proxies["RS_Number"].where(proxies["RS_Number"].str.startswith("rs"), proxies["Coord"]),
    "score": (proxies["R2"].astype(float) * 1000).round().astype(int),
    "strand": ".",
}).dropna(subset=["chrom"])
bed = bed.sort_values("score", ascending=False).drop_duplicates(["chrom", "start", "name"])
chrom_order = bed["chrom"].str.replace("chr", "").map(lambda c: (0, int(c)) if c.isdigit() else (1, c))
bed = bed.assign(order=chrom_order).sort_values(["order", "start"]).drop(columns="order")
bed.to_csv(args.bed_out, sep="\t", index=False, header=False)

failed = read_all(args.failures, "*.tsv")
failed.to_csv(args.failed_out, sep="\t", index=False)

print(f"{len(proxies)} proxy rows for {proxies['query_snp'].nunique()} query SNPs, "
      f"{len(failed)} SNPs failed")
