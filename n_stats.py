from __future__ import annotations

from pathlib import Path
import csv

#Use your own dataset files in the LocalDataset directory
LOCAL_DATASET_DIR = Path(__file__).resolve().parent / "LocalDataset"

DATASET_FILES = {
    "GDSC": LOCAL_DATASET_DIR / "GDSC_DATASET.csv",
    "GDSC2": LOCAL_DATASET_DIR / "GDSC2-dataset.csv",
}

DRUG_FIELDS = ["DRUG_ID", "DRUG_NAME"]
CELL_FIELDS = ["COSMIC_ID", "CELL_LINE_NAME"]


def pick_first(row: dict[str, str], fields: list[str]) -> str:
    for field in fields:
        value = (row.get(field) or "").strip()
        if value:
            return value
    return ""


def collect_stats(dataset_path: Path) -> tuple[set[str], set[str], set[tuple[str, str]]]:
    unique_drugs: set[str] = set()
    unique_cells: set[str] = set()
    unique_pairs: set[tuple[str, str]] = set()

    with dataset_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            drug = pick_first(row, DRUG_FIELDS)
            cell = pick_first(row, CELL_FIELDS)

            if drug:
                unique_drugs.add(drug)
            if cell:
                unique_cells.add(cell)
            if drug and cell:
                unique_pairs.add((cell, drug))

    return unique_drugs, unique_cells, unique_pairs


def main() -> None:
    missing = [name for name, path in DATASET_FILES.items() if not path.exists()]
    if missing:
        print("Missing dataset files:", ", ".join(missing))
        return

    per_dataset: dict[str, tuple[set[str], set[str], set[tuple[str, str]]]] = {}

    for dataset_name, dataset_path in DATASET_FILES.items():
        per_dataset[dataset_name] = collect_stats(dataset_path)

    print("N statistics per dataset:")
    for dataset_name in DATASET_FILES:
        drugs, cells, pairs = per_dataset[dataset_name]
        print(
            f"- {dataset_name}: drugs={len(drugs)}, cells={len(cells)}, pairs={len(pairs)}"
        )

    combined_drugs: set[str] = set()
    combined_cells: set[str] = set()
    combined_pairs: set[tuple[str, str]] = set()

    for drugs, cells, pairs in per_dataset.values():
        combined_drugs |= drugs
        combined_cells |= cells
        combined_pairs |= pairs

    print(
        "\nCombined (GDSC + GDSC2): "
        f"drugs={len(combined_drugs)}, cells={len(combined_cells)}, pairs={len(combined_pairs)}"
    )


if __name__ == "__main__":
    main()
