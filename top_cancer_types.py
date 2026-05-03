from __future__ import annotations

from pathlib import Path
import csv
from collections import Counter

#Use your own dataset files in the LocalDataset directory
LOCAL_DATASET_DIR = Path(__file__).resolve().parent / "LocalDataset"

GDSC_PATH = LOCAL_DATASET_DIR / "GDSC_DATASET.csv"
CCLE_PATH = LOCAL_DATASET_DIR / "CCLE_GlobalChromatinProfiling_20181130.csv"

GDSC_CANCER_FIELDS = [
    "Cancer Type (matching TCGA label)",
    "GDSC Tissue descriptor 1",
    "TCGA_DESC",
]

CCLE_CELL_LINE_FIELD = "CellLineName"


def normalize_label(label: str) -> str:
    value = label.strip()
    if not value:
        return ""
    return value.replace(" ", "_").upper()


def pick_first(row: dict[str, str], fields: list[str]) -> str:
    for field in fields:
        value = (row.get(field) or "").strip()
        if value:
            return value
    return ""


def extract_ccle_tissue(cell_line: str) -> str:
    parts = [part for part in cell_line.strip().split("_") if part]
    if len(parts) < 2:
        return ""
    return "_".join(parts[1:])


def count_gdsc_types(dataset_path: Path) -> Counter[str]:
    counts: Counter[str] = Counter()
    with dataset_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            label = pick_first(row, GDSC_CANCER_FIELDS)
            label = normalize_label(label)
            if label:
                counts[label] += 1
    return counts


def count_ccle_types(dataset_path: Path) -> Counter[str]:
    counts: Counter[str] = Counter()
    with dataset_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            cell_line = (row.get(CCLE_CELL_LINE_FIELD) or "").strip()
            tissue = extract_ccle_tissue(cell_line)
            tissue = normalize_label(tissue)
            if tissue:
                counts[tissue] += 1
    return counts


def main() -> None:
    missing = [
        path.name
        for path in [GDSC_PATH, CCLE_PATH]
        if not path.exists()
    ]
    if missing:
        print("Missing dataset files:", ", ".join(missing))
        return

    gdsc_counts = count_gdsc_types(GDSC_PATH)
    ccle_counts = count_ccle_types(CCLE_PATH)

    combined = Counter()
    combined.update(gdsc_counts)
    combined.update(ccle_counts)

    print("Cancer type sources:")
    print(f"- GDSC fields: {GDSC_CANCER_FIELDS}")
    print(
        "- CCLE derived from CellLineName suffix (text after first underscore)"
    )

    print("\nTop 3 combined cancer types:")
    for label, count in combined.most_common(3):
        print(f"- {label}: {count}")


if __name__ == "__main__":
    main()
