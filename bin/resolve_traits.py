#!/usr/bin/env python3
"""List every GWAS Catalog trait whose name contains a keyword."""
import argparse
import csv
import sys

import requests

TRAITS_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2/efo-traits"
NO_MATCH = 3  # exit status telling the pipeline not to retry


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--keyword", required=True)
    parser.add_argument("--out", required=True)
    return parser.parse_args(argv)


def fetch_traits(keyword):
    """Every EFO trait whose name contains the keyword, following the API's paging."""
    traits = []
    url = TRAITS_URL
    params = {"efo_trait": keyword, "size": 200}
    while url:
        r = requests.get(url, params=params, timeout=60)
        r.raise_for_status()
        data = r.json()
        traits.extend(data.get("_embedded", {}).get("efo_traits", []))
        url = data.get("_links", {}).get("next", {}).get("href")
        params = None  # the next link already carries the query string
    return traits


def write_traits(path, keyword, traits):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["keyword", "efo_trait", "efo_id"])
        for t in traits:
            writer.writerow([keyword, t["efo_trait"], t["efo_id"]])


def main(args):
    traits = fetch_traits(args.keyword)
    if not traits:
        print(f"No GWAS Catalog trait name contains '{args.keyword}'.", file=sys.stderr)
        sys.exit(NO_MATCH)

    write_traits(args.out, args.keyword, traits)
    print(f"{args.keyword}: {len(traits)} matching traits")


if __name__ == "__main__":
    main(parse_args())
