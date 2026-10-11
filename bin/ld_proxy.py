#!/usr/bin/env python3
"""Query LDlink LDproxy for a list of SNPs, one request at a time."""
import argparse
import csv
import io
import os
import sys
import time

import pandas as pd
import requests

API_URL = "https://ldlink.nih.gov/LDlinkRest/ldproxy"
COLUMNS = ["RS_Number", "Coord", "Alleles", "MAF", "Distance", "Dprime", "R2",
           "Correlated_Alleles", "FORGEdb", "RegulomeDB", "Function"]
BAD_SETTINGS = 3  # exit status telling the pipeline not to retry
# errors about the request as a whole rather than about one SNP
SETTINGS_ERRORS = ("api token", "ancestral population", "genome build")

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--snps", required=True, help="comma-separated rsIDs")
parser.add_argument("--pop", required=True, help="population codes joined with +")
parser.add_argument("--genome-build", required=True)
parser.add_argument("--window", type=int, default=500000)
parser.add_argument("--out", required=True)
parser.add_argument("--failed", required=True)
args = parser.parse_args()

token = os.environ.get("LDLINK_TOKEN")
if not token:
    print("LDLINK_TOKEN is not set. Run: nextflow secrets set LDLINK_TOKEN <your token>", file=sys.stderr)
    sys.exit(BAD_SETTINGS)


def ld_proxy(snp):
    """Return (table, None) on success or (None, message) when LDlink rejects the SNP."""
    query = {"var": snp, "pop": args.pop, "r2_d": "r2", "window": args.window,
             "genome_build": args.genome_build, "token": token}
    for attempt in range(1, 6):
        r = requests.get(API_URL, params=query, timeout=300)
        # LDlink reports errors as JSON, with HTTP 200
        error = r.json().get("error", r.text) if r.text.lstrip().startswith("{") else None
        busy = r.status_code in (429, 500, 502, 503, 504) or (error and "concurrent" in error.lower())
        if busy:
            time.sleep(15 * attempt)
            continue
        r.raise_for_status()
        if error is None:
            return pd.read_csv(io.StringIO(r.text), sep="\t", dtype=str, keep_default_na=False), None
        if any(s in error.lower() for s in SETTINGS_ERRORS):
            print(error, file=sys.stderr)
            sys.exit(BAD_SETTINGS)
        return None, error.strip()
    raise SystemExit(f"LDlink kept refusing requests for {snp}")


tables, failures = [], []
snp_list = args.snps.split(",")
for number, snp in enumerate(snp_list, start=1):
    # printed as each SNP is done, so a long run visibly moves along
    print(f"  {number} of {len(snp_list)} SNPs in this chunk: {snp}", flush=True)
    table, error = ld_proxy(snp)
    if error:
        failures.append([snp, error])
        continue
    table.insert(0, "query_snp", snp)
    tables.append(table)

proxies = pd.concat(tables) if tables else pd.DataFrame(columns=["query_snp"] + COLUMNS)
proxies.to_csv(args.out, sep="\t", index=False)

with open(args.failed, "w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
    writer.writerow(["query_snp", "error"])
    writer.writerows(failures)

print(f"{len(tables)} SNPs queried, {len(failures)} failed")
