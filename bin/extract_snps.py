#!/usr/bin/env python3
"""Pull the unique SNP rsIDs out of the merged GWAS Catalog table."""
import argparse
import ast
import csv
import re

import pandas as pd

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--table", required=True)
parser.add_argument("--out", required=True)
parser.add_argument("--skipped", required=True)
parser.add_argument("--pop", required=True, help="population(s) LD will be calculated in")
parser.add_argument("--genome-build", required=True)
args = parser.parse_args()

columns = ["snp_allele", "efo_traits", "discovery_ancestry"]
df = pd.read_csv(args.table, sep="\t", dtype=str, usecols=columns)


def ancestry_group(label):
    """'12730 European (Germany, U.K.)' -> 'European'."""
    return re.sub(r"\s*\([^)]*\)$", "", re.sub(r"^[\d,]+\s+", "", label.strip()))


# snp_allele holds a list like [{'rs_id': 'rs123', 'effect_allele': 'A'}],
# efo_traits one like [{'efo_id': 'MONDO_0004981', 'efo_trait': 'atrial fibrillation'}]
# and discovery_ancestry one like ['12730 European (Germany, U.K.)']
traits_by_variant = {}
ancestry_by_variant = {}
for snp_cell, trait_cell, ancestry_cell in df[columns].dropna(subset=["snp_allele"]).itertuples(index=False):
    traits = [t["efo_trait"] for t in ast.literal_eval(trait_cell)] if isinstance(trait_cell, str) else []
    ancestries = [ancestry_group(a) for a in ast.literal_eval(ancestry_cell)] if isinstance(ancestry_cell, str) else []
    for entry in ast.literal_eval(snp_cell):
        variant = entry["rs_id"].strip()
        traits_by_variant.setdefault(variant, set()).update(traits)
        ancestry_by_variant.setdefault(variant, set()).update(ancestries)

# LD proxy lookups need rsIDs; chr:pos identifiers have no stated genome build
rsids = sorted(v for v in traits_by_variant if re.fullmatch(r"rs\d+", v))
skipped = sorted(set(traits_by_variant) - set(rsids))


def joined_traits(variant):
    return "; ".join(sorted(traits_by_variant[variant], key=str.lower))


def joined_ancestry(variant):
    return "; ".join(sorted(ancestry_by_variant[variant], key=str.lower)) or "not reported"


with open(args.out, "w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
    # gwas_ancestry is who the SNP was found in; ld_population is who its proxies are looked up in
    writer.writerow(["query_snp", "traits", "gwas_ancestry", "ld_population", "genome_build"])
    for rsid in rsids:
        writer.writerow([rsid, joined_traits(rsid), joined_ancestry(rsid), args.pop, args.genome_build])

with open(args.skipped, "w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
    writer.writerow(["variant", "traits", "reason"])
    for v in skipped:
        writer.writerow([v, joined_traits(v), "not an rsID"])

print(f"{len(rsids)} rsIDs, {len(skipped)} variants skipped")
