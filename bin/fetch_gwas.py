#!/usr/bin/env python3
"""Download all GWAS Catalog associations for one trait as a TSV."""
import argparse
import sys

import pandas as pd
import requests

API_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2/associations"
TRAITS_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2/efo-traits"
UNKNOWN_TRAIT = 3  # exit status telling the pipeline not to retry

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--trait", required=True)
parser.add_argument("--out", required=True)
parser.add_argument("--allow-empty", action="store_true",
                    help="exit quietly, writing no file, when the trait has no associations")
args = parser.parse_args()

records = []
url = API_URL
params = {"efo_trait": args.trait, "show_child_traits": "true", "size": 200}
while url:
    r = requests.get(url, params=params, timeout=60)
    r.raise_for_status()
    data = r.json()
    for a in data.get("_embedded", {}).get("associations", []):
        a["query_trait"] = args.trait  # remember which trait this came from
        records.append(a)
    url = data.get("_links", {}).get("next", {}).get("href")
    params = None  # the next link already carries the query string

if not records and args.allow_empty:
    print(f"{args.trait}: 0 associations, skipped")
    sys.exit(0)

if not records:
    # the API only matches exact ontology trait names, so suggest similar ones
    msg = f"No associations found for trait '{args.trait}'."
    r = requests.get(TRAITS_URL, params={"efo_trait": args.trait, "size": 20}, timeout=60)
    names = [t["efo_trait"] for t in r.json().get("_embedded", {}).get("efo_traits", [])] if r.ok else []
    if names:
        msg += " The GWAS Catalog needs the exact trait name. Did you mean: " + "; ".join(names) + "?"
    print(msg, file=sys.stderr)
    sys.exit(UNKNOWN_TRAIT)

df = pd.json_normalize(records)
df.to_csv(args.out, sep="\t", index=False)
print(f"{args.trait}: {len(df)} associations")
