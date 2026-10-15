#!/usr/bin/env python3
"""Download all GWAS Catalog associations for one trait as a TSV."""
import argparse
import sys
import time

import pandas as pd
import requests

API_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2/associations"
TRAITS_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2/efo-traits"
STUDIES_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2/studies"
UNKNOWN_TRAIT = 3  # exit status telling the pipeline not to retry

session = requests.Session()


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trait", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--allow-empty", action="store_true",
                        help="exit quietly, writing no file, when the trait has no associations")
    return parser.parse_args(argv)


def get(url, params=None):
    """GET that waits and tries again when the GWAS Catalog is slow or busy."""
    for attempt in range(1, 6):
        try:
            r = session.get(url, params=params, timeout=120)
            if r.status_code not in (429, 500, 502, 503, 504):
                return r
            problem = f"HTTP {r.status_code}"
        except (requests.ConnectionError, requests.Timeout) as e:
            problem = type(e).__name__
        print(f"{problem} from the GWAS Catalog, attempt {attempt} of 5", file=sys.stderr)
        time.sleep(10 * attempt)
    sys.exit(f"The GWAS Catalog kept failing for {url}")


def fetch_associations(trait):
    """Every association for the trait, following the API's paging."""
    records = []
    url = API_URL
    params = {"efo_trait": trait, "show_child_traits": "true", "size": 200}
    while url:
        r = get(url, params)
        r.raise_for_status()
        data = r.json()
        for a in data.get("_embedded", {}).get("associations", []):
            a["query_trait"] = trait  # remember which trait this came from
            records.append(a)
        url = data.get("_links", {}).get("next", {}).get("href")
        params = None  # the next link already carries the query string
    return records


def suggest_traits(trait):
    """Trait names close to one the API did not recognise."""
    r = get(TRAITS_URL, {"efo_trait": trait, "size": 20})
    return [t["efo_trait"] for t in r.json().get("_embedded", {}).get("efo_traits", [])] if r.ok else []


def fetch_ancestry(accessions):
    """Who was studied, per study accession, e.g.
    ['12730 European (Germany, U.K.)', '1606 African American or Afro-Caribbean (U.S.)'].
    An association does not carry this, so its study has to be asked."""
    ancestry = {}
    for accession in accessions:
        r = get(f"{STUDIES_URL}/{accession}")
        if r.status_code == 404:
            ancestry[accession] = []
            continue
        r.raise_for_status()
        ancestry[accession] = r.json().get("discovery_ancestry") or []
    return ancestry


def main(args):
    records = fetch_associations(args.trait)

    if not records and args.allow_empty:
        print(f"{args.trait}: 0 associations, skipped")
        sys.exit(0)

    if not records:
        # the API only matches exact ontology trait names, so suggest similar ones
        msg = f"No associations found for trait '{args.trait}'."
        names = suggest_traits(args.trait)
        if names:
            msg += " The GWAS Catalog needs the exact trait name. Did you mean: " + "; ".join(names) + "?"
        print(msg, file=sys.stderr)
        sys.exit(UNKNOWN_TRAIT)

    df = pd.json_normalize(records)
    ancestry = fetch_ancestry(df["accession_id"].dropna().unique())
    df["discovery_ancestry"] = df["accession_id"].map(ancestry)

    df.to_csv(args.out, sep="\t", index=False)
    print(f"{args.trait}: {len(df)} associations from {len(ancestry)} studies")


if __name__ == "__main__":
    main(parse_args())
