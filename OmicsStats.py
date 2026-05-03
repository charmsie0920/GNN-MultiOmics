from __future__ import annotations

from pathlib import Path
import csv


LOCAL_DATASET_DIR = Path(__file__).resolve().parent / "LocalDataset"

DATASET_FILES = {
	"GDSC": LOCAL_DATASET_DIR / "GDSC_DATASET.csv",
	"GDSC2": LOCAL_DATASET_DIR / "GDSC2-dataset.csv",
}

MERGE_KEYS = ["COSMIC_ID", "DRUG_ID"]
JOIN_TYPES = ("inner", "left", "right", "outer")

# Omics types we look for in column names, case-insensitive.
OMICS_KEYWORDS = {
	"copy_number": ["cna", "copy number"],
	"gene_expression": ["gene expression", "expression"],
	"methylation": ["methylation", "methyl"],
	"mirna": ["mirna"],
	"rppa": ["rppa"],
	"metabolomics": ["metabolomics", "metabolite"],
	"chromatin": ["chromatin", "atac", "chip"],
}

MISSING_VALUES = {"", "n", "no", "0", "na", "nan"}


def detect_omics_columns(headers: list[str]) -> dict[str, set[str]]:
	matches: dict[str, set[str]] = {}
	for column in headers:
		column_lower = column.lower()
		for omics_type, keywords in OMICS_KEYWORDS.items():
			if any(keyword in column_lower for keyword in keywords):
				matches.setdefault(omics_type, set()).add(column)
	return matches



def has_signal(value: str) -> bool:
	return value.strip().lower() not in MISSING_VALUES



def normalize_value(value: str | None) -> str:
	return (value or "").strip()


def make_key(row: dict[str, str], key_fields: list[str]) -> tuple[str, ...] | None:
	key_parts: list[str] = []
	for field in key_fields:
		value = normalize_value(row.get(field))
		if not value:
			return None
		key_parts.append(value)
	return tuple(key_parts)


def build_index(
	dataset_path: Path,
	key_fields: list[str],
) -> tuple[list[str], dict[tuple[str, ...], list[dict[str, str]]], list[str]]:
	with dataset_path.open("r", encoding="utf-8", newline="") as handle:
		reader = csv.DictReader(handle)
		headers = reader.fieldnames or []
		missing_fields = [field for field in key_fields if field not in headers]
		if missing_fields:
			return headers, {}, missing_fields

		index: dict[tuple[str, ...], list[dict[str, str]]] = {}
		for row in reader:
			key = make_key(row, key_fields)
			if key is None:
				continue
			index.setdefault(key, []).append(row)

	return headers, index, []


def build_merged_headers(
	left_headers: list[str],
	right_headers: list[str],
	right_prefix: str,
) -> tuple[list[str], dict[str, str], dict[str, str]]:
	merged_headers: list[str] = []
	left_map: dict[str, str] = {}
	right_map: dict[str, str] = {}
	existing: set[str] = set()

	for header in left_headers:
		left_map[header] = header
		merged_headers.append(header)
		existing.add(header)

	for header in right_headers:
		if header in existing:
			output_name = f"{right_prefix}{header}"
		else:
			output_name = header
		right_map[header] = output_name
		merged_headers.append(output_name)
		existing.add(output_name)

	return merged_headers, left_map, right_map


def iter_joined_rows(
	left_index: dict[tuple[str, ...], list[dict[str, str]]],
	right_index: dict[tuple[str, ...], list[dict[str, str]]],
	join_type: str,
	left_headers: list[str],
	right_headers: list[str],
	left_map: dict[str, str],
	right_map: dict[str, str],
):
	left_keys = set(left_index)
	right_keys = set(right_index)

	if join_type == "inner":
		keys = left_keys & right_keys
	elif join_type == "left":
		keys = left_keys
	elif join_type == "right":
		keys = right_keys
	else:
		keys = left_keys | right_keys

	for key in keys:
		left_rows = left_index.get(key) or [None]
		right_rows = right_index.get(key) or [None]
		for left_row in left_rows:
			for right_row in right_rows:
				merged_row: dict[str, str] = {}
				for header in left_headers:
					merged_row[left_map[header]] = (
						left_row.get(header, "") if left_row else ""
					)
				for header in right_headers:
					merged_row[right_map[header]] = (
						right_row.get(header, "") if right_row else ""
					)
				yield merged_row


def detect_omics_usage_from_rows(
	headers: list[str],
	rows_iter,
) -> tuple[dict[str, set[str]], set[str], dict[str, int], int]:
	omics_columns = detect_omics_columns(headers)
	used_omics: set[str] = set()
	counts = {omics_type: 0 for omics_type in omics_columns}
	row_count = 0

	if not omics_columns:
		return omics_columns, used_omics, counts, row_count

	for row in rows_iter:
		row_count += 1
		for omics_type, columns in omics_columns.items():
			if any(has_signal(row.get(column, "")) for column in columns):
				counts[omics_type] += 1
				used_omics.add(omics_type)

	return omics_columns, used_omics, counts, row_count


def detect_omics_usage_in_file(
	dataset_path: Path,
) -> tuple[dict[str, set[str]], set[str], dict[str, int], int]:
	with dataset_path.open("r", encoding="utf-8", newline="") as handle:
		reader = csv.DictReader(handle)
		if reader.fieldnames is None:
			return {}, set(), {}, 0

		return detect_omics_usage_from_rows(reader.fieldnames, reader)


def main() -> None:
	missing = [name for name, path in DATASET_FILES.items() if not path.exists()]
	if missing:
		print("Missing dataset files:", ", ".join(missing))
		return

	results: dict[str, set[str]] = {}
	column_matches: dict[str, dict[str, set[str]]] = {}
	omics_counts: dict[str, dict[str, int]] = {}

	for dataset_name, dataset_path in DATASET_FILES.items():
		omics_columns, used_omics, counts, _ = detect_omics_usage_in_file(dataset_path)
		column_matches[dataset_name] = omics_columns
		results[dataset_name] = used_omics
		omics_counts[dataset_name] = counts

	union_omics = set().union(*results.values())
	intersection_omics = set.intersection(*results.values()) if results else set()

	print("Omics types detected by dataset:")
	for dataset_name in DATASET_FILES:
		omics_list = sorted(results[dataset_name])
		print(f"- {dataset_name}: {len(omics_list)} -> {omics_list}")

	print("\nOmics types by column matches (before value check):")
	for dataset_name in DATASET_FILES:
		matched = {
			omics: sorted(columns)
			for omics, columns in column_matches[dataset_name].items()
		}
		print(f"- {dataset_name}: {matched}")

	print("\nOmics mention counts (rows with signal):")
	for dataset_name in DATASET_FILES:
		counts = omics_counts[dataset_name]
		if counts:
			count_list = ", ".join(
				f"{omics}: {counts[omics]}" for omics in sorted(counts)
			)
		else:
			count_list = "none"
		print(f"- {dataset_name}: {count_list}")

	print("\nComparison summary:")
	print(f"- Union of omics types: {len(union_omics)} -> {sorted(union_omics)}")
	print(
		f"- Intersection of omics types: {len(intersection_omics)} -> {sorted(intersection_omics)}"
	)

	left_path = DATASET_FILES["GDSC"]
	right_path = DATASET_FILES["GDSC2"]
	left_headers, left_index, left_missing = build_index(left_path, MERGE_KEYS)
	right_headers, right_index, right_missing = build_index(right_path, MERGE_KEYS)

	if left_missing or right_missing:
		missing_message = ", ".join(sorted(set(left_missing + right_missing)))
		print("\nMerge skipped; missing columns:", missing_message)
		return

	merged_headers, left_map, right_map = build_merged_headers(
		left_headers,
		right_headers,
		right_prefix="GDSC2_",
	)

	print("\nMerged omics summary (COSMIC_ID + DRUG_ID):")
	for join_type in JOIN_TYPES:
		merged_rows = iter_joined_rows(
			left_index,
			right_index,
			join_type,
			left_headers,
			right_headers,
			left_map,
			right_map,
		)
		_, used_omics, counts, row_count = detect_omics_usage_from_rows(
			merged_headers,
			merged_rows,
		)
		omics_list = sorted(used_omics)
		print(f"- Join {join_type}: rows={row_count}")
		print(f"  Omics types: {len(omics_list)} -> {omics_list}")
		if counts:
			count_list = ", ".join(
				f"{omics}: {counts[omics]}" for omics in sorted(counts)
			)
		else:
			count_list = "none"
		print(f"  Omics mention counts: {count_list}")

	print("\nJoin guidance:")
	print("- inner: overlap only; best for strict comparability")
	print("- left: keep all GDSC rows; best when GDSC is primary")
	print("- right: keep all GDSC2 rows; best when GDSC2 is primary")
	print("- outer: keep all rows; best for coverage, includes unmatched")


if __name__ == "__main__":
	main()