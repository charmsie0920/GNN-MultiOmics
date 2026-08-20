"""Download the STRING v12.0 human protein alias file, once.

`03_graph_construction.py` needs to map GDSC's free-text gene-symbol drug
targets (e.g. `EGFR`) onto STRING's Ensembl protein IDs (`ENSP*`) so a
`(drug, targets, protein)` edge can be added between a resolved target and the
matching PPI node. No local alias file exists yet; this script downloads
`9606.protein.aliases.v12.0.txt.gz`, matching the v12.0 STRING release already
used for `data/raw/string/9606.protein.links.v12.0.txt.gz` (see
`ingest_and_align.py`).
"""

from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path
from typing import Optional, Sequence

STRING_ALIASES_URL = (
    "https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz"
)
OUTPUT_PATH = Path("data/raw/string/9606.protein.aliases.v12.0.txt.gz")
REQUEST_TIMEOUT_SECONDS = 60


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Download the STRING alias file to `data/raw/string/` if not already present."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", default=STRING_ALIASES_URL)
    parser.add_argument("--output", type=Path, default=OUTPUT_PATH)
    args = parser.parse_args(argv)

    args.output.parent.mkdir(parents=True, exist_ok=True)

    if args.output.exists():
        print(f"Already present, skipping download: {args.output}")
    else:
        with urllib.request.urlopen(args.url, timeout=REQUEST_TIMEOUT_SECONDS) as response:
            args.output.write_bytes(response.read())

    size_mb = args.output.stat().st_size / (1024 * 1024)
    print(f"Downloaded file: {args.output}")
    print(f"Downloaded file size: {size_mb:.1f} MB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
