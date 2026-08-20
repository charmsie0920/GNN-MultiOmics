"""Fetch canonical SMILES for GDSC's screened drugs from PubChem, once, and cache to disk.

`03_graph_construction.py` needs a SMILES string per GDSC drug to build RDKit
Morgan fingerprints for the `drug` node type. Rather than hitting the live
PubChem API at graph-build time (fragile — no guarantee the API is reachable
whenever the graph script runs) or downloading all of PubChem (unnecessary --
GDSC only screens ~600 distinct compounds), this script resolves each drug
exactly once against PubChem's PUG REST API and writes the result to a local
CSV that `03_graph_construction.py` reads like any other local file.

The PubChem query is deduped on `DRUG_NAME`, not `DRUG_ID`: GDSC lists 621
`DRUG_ID`s but only 542 unique `DRUG_NAME`s (71 names are re-screened under a
second `DRUG_ID`, e.g. Erlotinib -> IDs 1 and 1168). Querying by name and then
fanning the resolved SMILES/CID out to every `DRUG_ID` sharing that name keeps
the API call count at the unique-compound count instead of the row count.
"""

from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

GDSC_COMPOUNDS_PATH = Path("data/raw/gdsc/screened_compounds_rel_8.5.csv")
OUTPUT_PATH = Path("data/raw/pubchem/gdsc_drug_smiles.csv")
OUTPUT_FIELDS = ["drug_id", "drug_name", "canonical_smiles", "pubchem_cid"]

PUBCHEM_PROPERTY_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/SMILES/JSON"
)
# PubChem's documented anonymous rate limit is ~5 requests/second.
REQUEST_INTERVAL_SECONDS = 0.21
REQUEST_TIMEOUT_SECONDS = 15
MAX_RETRIES = 3


def read_gdsc_drugs(path: Path) -> List[Dict[str, str]]:
    """Read GDSC's screened-compounds table as a list of row dicts."""

    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def split_synonyms(raw: str) -> List[str]:
    """Split GDSC's `SYNONYMS` column (comma-separated) into a clean list."""

    if not raw:
        return []
    return [part.strip() for part in raw.split(",") if part.strip()]


def query_pubchem_smiles(name: str) -> Optional[Tuple[str, str]]:
    """Query PubChem PUG REST for a compound name; return (smiles, cid) or None.

    Retries on transient errors (e.g. PubChem's 503 "server busy") but treats a
    404 (name not found) as a normal, non-retryable miss.
    """

    url = PUBCHEM_PROPERTY_URL.format(name=urllib.parse.quote(name, safe=""))
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            with urllib.request.urlopen(url, timeout=REQUEST_TIMEOUT_SECONDS) as response:
                payload = json.loads(response.read().decode("utf-8"))
            props = payload["PropertyTable"]["Properties"][0]
            return str(props["SMILES"]), str(props["CID"])
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return None
            if attempt == MAX_RETRIES:
                return None
            time.sleep(REQUEST_INTERVAL_SECONDS * attempt)
        except (urllib.error.URLError, TimeoutError, KeyError, IndexError, json.JSONDecodeError):
            if attempt == MAX_RETRIES:
                return None
            time.sleep(REQUEST_INTERVAL_SECONDS * attempt)
    return None


def resolve_drug_name(name: str, synonyms: Sequence[str]) -> Optional[Tuple[str, str]]:
    """Resolve one drug name to (smiles, cid), falling back through its synonyms."""

    result = query_pubchem_smiles(name)
    time.sleep(REQUEST_INTERVAL_SECONDS)
    if result is not None:
        return result
    for synonym in synonyms:
        result = query_pubchem_smiles(synonym)
        time.sleep(REQUEST_INTERVAL_SECONDS)
        if result is not None:
            return result
    return None


def load_existing_cache(path: Path) -> Dict[str, Dict[str, str]]:
    """Load a previously written output CSV, keyed by `drug_name`, for resumability."""

    if not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["drug_name"]: row for row in csv.DictReader(handle)}


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Resolve SMILES for every unique GDSC drug name and write the cached CSV."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gdsc-compounds", type=Path, default=GDSC_COMPOUNDS_PATH)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args(argv)

    if not args.gdsc_compounds.exists():
        raise FileNotFoundError(f"Missing GDSC compounds file: {args.gdsc_compounds}")

    rows = read_gdsc_drugs(args.gdsc_compounds)

    # name -> list of (drug_id, drug_name) rows sharing that name, in first-seen order.
    rows_by_name: Dict[str, List[Dict[str, str]]] = {}
    for row in rows:
        rows_by_name.setdefault(row["DRUG_NAME"], []).append(row)

    existing_by_name = load_existing_cache(args.output)

    resolved_by_name: Dict[str, Tuple[str, str]] = {}
    unresolved_names: List[str] = []
    queried_count = 0

    for name, name_rows in rows_by_name.items():
        cached = existing_by_name.get(name)
        if cached is not None and cached.get("canonical_smiles"):
            resolved_by_name[name] = (cached["canonical_smiles"], cached["pubchem_cid"])
            continue
        if cached is not None:
            # Previously attempted and unresolved; skip re-querying on resume.
            unresolved_names.append(name)
            continue

        synonyms = split_synonyms(name_rows[0].get("SYNONYMS", ""))
        result = resolve_drug_name(name, synonyms)
        queried_count += 1
        if result is None:
            unresolved_names.append(name)
        else:
            resolved_by_name[name] = result

    output_rows: List[Dict[str, str]] = []
    for name, name_rows in rows_by_name.items():
        smiles, cid = resolved_by_name.get(name, ("", ""))
        for row in name_rows:
            output_rows.append(
                {
                    "drug_id": row["DRUG_ID"],
                    "drug_name": name,
                    "canonical_smiles": smiles,
                    "pubchem_cid": cid,
                }
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS)
        writer.writeheader()
        writer.writerows(output_rows)

    output_size_kb = args.output.stat().st_size / 1024
    print(f"Unique drug names in GDSC: {len(rows_by_name)}")
    print(f"PubChem queries issued this run: {queried_count}")
    print(f"Names resolved: {len(resolved_by_name)}")
    print(f"Names unresolved: {len(unresolved_names)}")
    if unresolved_names:
        print("Unresolved drug names:")
        for name in unresolved_names:
            print(f"  - {name}")
    print(f"Output rows written: {len(output_rows)} -> {args.output}")
    print(f"Output file size: {output_size_kb:.1f} KB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
