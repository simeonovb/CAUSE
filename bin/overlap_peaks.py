#!/usr/bin/env python3
"""Find the ATAC-seq peaks that contain at least one SNP."""
import argparse
import gzip
import sys

import pandas as pd


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--peaks", required=True, help="ATAC-seq peaks in BED format (may be gzipped)")
    parser.add_argument("--snps", required=True, help="LD proxy table with a Coord column like chr1:12345")
    parser.add_argument("--out-peaks", required=True, help="BED of the peaks containing a SNP")
    parser.add_argument("--out-snps", required=True, help="table of the SNPs in each of those peaks")
    return parser.parse_args(argv)


def chrom_key(names):
    """Chromosome names without the chr prefix, so 'chr1' and '1' match."""
    return names.str.replace(r"^chr", "", case=False, regex=True)


def read_peaks(path):
    """Return the peaks as a table, plus the file's lines so they can be written back unchanged."""
    # BED has no header row, but may start with track, browser or comment lines
    with open(path, "rb") as fh:
        opener = gzip.open if fh.read(2) == b"\x1f\x8b" else open
    with opener(path, "rt") as fh:
        lines = [line.rstrip("\r\n") for line in fh
                 if line.strip() and not line.startswith(("#", "track", "browser"))]

    fields = [line.split("\t") if "\t" in line else line.split() for line in lines]
    try:
        peaks = pd.DataFrame({
            "peak_chr": [f[0] for f in fields],
            "peak_start": [int(f[1]) for f in fields],
            "peak_end": [int(f[2]) for f in fields],
            "peak_name": [f[3] if len(f) > 3 else f"{f[0]}:{f[1]}-{f[2]}" for f in fields],
        })
    except (IndexError, ValueError):
        sys.exit("The peak file is not in BED format: each line needs a chromosome, a start and an end")
    if peaks.empty:
        sys.exit("The peak file has no peaks in it")
    peaks["chrom"] = chrom_key(peaks["peak_chr"])
    return peaks, lines


def read_snps(path):
    """Return the proxy table, with chrom and pos pulled out of its Coord column."""
    snps = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    coord = snps["Coord"].str.extract(r"^(?P<chrom>[^:]+):(?P<pos>\d+)$")
    snps = snps[coord["pos"].notna()].copy()
    snps["chrom"] = chrom_key(coord["chrom"].dropna())
    snps["pos"] = coord["pos"].dropna().astype(int)
    return snps


def overlap(peaks, snps):
    """Pair each peak with the SNP positions that sit inside it.

    BED starts count from 0 and SNP positions from 1, so a SNP at position p
    is inside a peak when start < p <= end.
    """
    pairs = []
    for chrom, positions in snps.groupby("chrom")["pos"].unique().items():
        chrom_peaks = peaks[peaks["chrom"] == chrom]
        starts, ends = chrom_peaks["peak_start"].to_numpy(), chrom_peaks["peak_end"].to_numpy()
        for pos in positions:
            for peak in chrom_peaks.index[(starts < pos) & (pos <= ends)]:
                pairs.append((peak, chrom, pos))
    return pd.DataFrame(pairs, columns=["peak", "chrom", "pos"])


def snps_in_peaks(pairs, peaks, snps):
    """One row per peak and SNP (and per query SNP the SNP is a proxy of)."""
    peak_columns = ["peak_chr", "peak_start", "peak_end", "peak_name"]
    table = (pairs.merge(peaks[peak_columns], left_on="peak", right_index=True)
             .merge(snps, on=["chrom", "pos"])
             .sort_values(["peak", "pos", "query_snp"]))
    return table.drop(columns=["peak", "chrom", "pos"])


def main(args):
    peaks, lines = read_peaks(args.peaks)
    snps = read_snps(args.snps)
    pairs = overlap(peaks, snps)

    # the peaks, as they were given, in their original order
    with open(args.out_peaks, "w") as fh:
        for peak in sorted(pairs["peak"].unique()):
            fh.write(lines[peak] + "\n")

    snps_in_peaks(pairs, peaks, snps).to_csv(args.out_snps, sep="\t", index=False)

    print(f"{pairs['peak'].nunique()} of {len(peaks)} peaks contain at least one SNP; "
          f"{len(pairs.drop_duplicates(['chrom', 'pos']))} of {len(snps.drop_duplicates(['chrom', 'pos']))} "
          f"SNP positions are in a peak")
    if not set(peaks["chrom"]) & set(snps["chrom"]):
        print("Warning: the peaks and the SNPs have no chromosome names in common", file=sys.stderr)


if __name__ == "__main__":
    main(parse_args())
