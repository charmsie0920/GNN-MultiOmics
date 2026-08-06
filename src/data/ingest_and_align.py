"""Ingest and align the raw GDSC and CCLE/DepMap inputs.

This module builds a first-pass aligned bundle from the uploaded raw data:

- GDSC2 fitted dose-response rows are standardized against the model crosswalk.
- CCLE/DepMap RNA-seq, CNV, and mutation archives are normalized to shared IDs.
- STRING presence is validated and summarized for later graph construction.

The script is dependency-light and uses only the Python standard library so it
can run in the current workspace without pandas.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import zipfile
import xml.etree.ElementTree as ET
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple


RAW_ROOT = Path("data/raw")
OUTPUT_ROOT = Path("data/processed/aligned")

MODEL_LIST_PATH = RAW_ROOT / "ccle_depmap" / "model_list_20260724.csv"
CELL_LINES_PATH = RAW_ROOT / "ccle_depmap" / "Cell_lines_annotations_20181226.txt"
RNASEQ_ZIP_PATH = RAW_ROOT / "ccle_depmap" / "rnaseq_all_20260323.zip"
MUTATIONS_ZIP_PATH = RAW_ROOT / "ccle_depmap" / "mutations_all_20260724.zip"
CNV_ZIP_PATH = RAW_ROOT / "ccle_depmap" / "cnv_summary_20260316.zip"
GDSC_XLSX_PATH = RAW_ROOT / "gdsc" / "GDSC2_fitted_dose_response_27Oct23.xlsx"
STRING_PATH = RAW_ROOT / "string" / "9606.protein.links.v12.0.txt.gz"


@dataclass(frozen=True)
class CrosswalkRecord:
    """Unified model metadata used to standardize identifiers."""

    sanger_model_id: str
    depmap_id: str
    cell_line_name: str
    broad_id: str
    ccle_id: str
    cosmic_id: str
    rrid: str
    tissue: str
    cancer_type: str


def normalize_key(value: Optional[str]) -> str:
    """Return a canonical uppercase key for matching identifiers and names."""

    if value is None:
        return ""
    return re.sub(r"[^A-Z0-9]+", "", value.strip().upper())


def strip_value(value: Optional[str]) -> str:
    """Coerce empty values to an empty string and trim whitespace."""

    if value is None:
        return ""
    return str(value).strip()


def ensure_parent(path: Path) -> None:
    """Create the parent directory for a file path if needed."""

    path.parent.mkdir(parents=True, exist_ok=True)


def write_csv(path: Path, rows: Iterable[Dict[str, str]], fieldnames: Sequence[str]) -> int:
    """Write dictionaries to CSV and return the number of rows written."""

    ensure_parent(path)
    count = 0
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: strip_value(row.get(field)) for field in fieldnames})
            count += 1
    return count


def read_csv_rows(path: Path) -> Iterator[Dict[str, str]]:
    """Yield rows from a plain CSV file as dictionaries."""

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            yield {key: strip_value(value) for key, value in row.items()}


def read_tab_rows(path: Path) -> Iterator[Dict[str, str]]:
    """Yield rows from a tab-delimited text file as dictionaries."""

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            yield {key: strip_value(value) for key, value in row.items()}


def read_csv_from_zip(zip_path: Path) -> Iterator[Dict[str, str]]:
    """Stream a single CSV file from a ZIP archive as dictionaries."""

    with zipfile.ZipFile(zip_path) as archive:
        members = [name for name in archive.namelist() if not name.endswith("/")]
        if not members:
            return
        with archive.open(members[0]) as raw_handle:
            text_handle = (line.decode("utf-8", "ignore") for line in raw_handle)
            reader = csv.DictReader(text_handle)
            for row in reader:
                yield {key: strip_value(value) for key, value in row.items()}


def column_index(reference: str) -> int:
    """Convert an Excel cell reference such as 'AB12' to a zero-based column index."""

    match = re.match(r"([A-Z]+)", reference)
    if not match:
        return 0
    letters = match.group(1)
    index = 0
    for char in letters:
        index = index * 26 + (ord(char) - ord("A") + 1)
    return index - 1


def xlsx_first_sheet_path(archive: zipfile.ZipFile) -> str:
    """Return the worksheet path for the first sheet in an XLSX workbook."""

    workbook_ns = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    rel_ns = {"r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}
    workbook_root = ET.fromstring(archive.read("xl/workbook.xml"))
    first_sheet = workbook_root.find("a:sheets/a:sheet", workbook_ns)
    if first_sheet is None:
        raise ValueError("Workbook does not contain any sheets")
    rel_id = first_sheet.attrib.get(f"{{{rel_ns['r']}}}id")
    rel_root = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    for rel in rel_root:
        if rel.attrib.get("Id") == rel_id:
            target = rel.attrib["Target"].lstrip("/")
            if not target.startswith("xl/"):
                target = f"xl/{target}"
            return target
    raise ValueError("Could not resolve the first worksheet path")


def xlsx_shared_strings(archive: zipfile.ZipFile) -> List[str]:
    """Load shared strings from an XLSX archive when present."""

    if "xl/sharedStrings.xml" not in archive.namelist():
        return []
    ns = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
    values: List[str] = []
    for item in root.findall("a:si", ns):
        values.append("".join(node.text or "" for node in item.iterfind(".//a:t", ns)))
    return values


def xlsx_row_values(cell_elements: List[ET.Element], shared_strings: Sequence[str]) -> Dict[int, str]:
    """Convert a list of Excel cell elements into a column-indexed mapping."""

    ns = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    values: Dict[int, str] = {}
    for cell in cell_elements:
        ref = cell.attrib.get("r", "A1")
        col = column_index(ref)
        cell_type = cell.attrib.get("t")
        if cell_type == "inlineStr":
            text = "".join(node.text or "" for node in cell.findall(".//a:t", ns))
        else:
            raw_value = cell.findtext("a:v", default="", namespaces=ns)
            if cell_type == "s" and raw_value.isdigit() and int(raw_value) < len(shared_strings):
                text = shared_strings[int(raw_value)]
            else:
                text = raw_value
        values[col] = strip_value(text)
    return values


def iter_xlsx_rows(xlsx_path: Path) -> Iterator[Dict[str, str]]:
    """Yield rows from the first sheet of an XLSX workbook as dictionaries."""

    ns = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with zipfile.ZipFile(xlsx_path) as archive:
        sheet_path = xlsx_first_sheet_path(archive)
        shared_strings = xlsx_shared_strings(archive)
        headers: List[str] = []
        with archive.open(sheet_path) as sheet_handle:
            context = ET.iterparse(sheet_handle, events=("end",))
            for _, element in context:
                if element.tag.endswith("row"):
                    cells = list(element.findall("a:c", ns))
                    indexed_values = xlsx_row_values(cells, shared_strings)
                    if not headers:
                        max_index = max(indexed_values.keys(), default=-1)
                        headers = ["" for _ in range(max_index + 1)]
                        for index, value in indexed_values.items():
                            headers[index] = value
                    else:
                        row = {headers[index]: value for index, value in indexed_values.items() if index < len(headers) and headers[index]}
                        yield {key: strip_value(value) for key, value in row.items()}
                    element.clear()


def load_model_crosswalk(model_list_path: Path, cell_lines_path: Path) -> Tuple[List[CrosswalkRecord], Dict[str, CrosswalkRecord], Dict[str, CrosswalkRecord]]:
    """Build crosswalk records keyed by Sanger model ID and DepMap ID."""

    annotations_by_ccle: Dict[str, Dict[str, str]] = {}
    annotations_by_name: Dict[str, Dict[str, str]] = {}
    if cell_lines_path.exists():
        for row in read_tab_rows(cell_lines_path):
            ccle_key = strip_value(row.get("CCLE_ID"))
            depmap_key = strip_value(row.get("depMapID"))
            name_key = strip_value(row.get("Name"))
            if ccle_key:
                annotations_by_ccle[normalize_key(ccle_key)] = row
            if depmap_key:
                annotations_by_ccle[normalize_key(depmap_key)] = row
            if name_key:
                annotations_by_name[normalize_key(name_key)] = row

    records: List[CrosswalkRecord] = []
    by_sanger: Dict[str, CrosswalkRecord] = {}
    by_depmap: Dict[str, CrosswalkRecord] = {}

    for row in read_csv_rows(model_list_path):
        sanger_id = strip_value(row.get("model_id"))
        depmap_id = strip_value(row.get("BROAD_ID"))
        cell_line_name = strip_value(row.get("model_name"))
        ccle_id = strip_value(row.get("CCLE_ID"))
        cosmic_id = strip_value(row.get("COSMIC_ID"))
        rrid = strip_value(row.get("RRID"))
        tissue = strip_value(row.get("tissue"))
        cancer_type = strip_value(row.get("cancer_type"))

        annotation = annotations_by_ccle.get(normalize_key(ccle_id)) or annotations_by_name.get(normalize_key(cell_line_name))
        if annotation:
            if not depmap_id:
                depmap_id = strip_value(annotation.get("depMapID"))
            if not cell_line_name:
                cell_line_name = strip_value(annotation.get("Name"))
            if not ccle_id:
                ccle_id = strip_value(annotation.get("CCLE_ID"))

        record = CrosswalkRecord(
            sanger_model_id=sanger_id,
            depmap_id=depmap_id,
            cell_line_name=cell_line_name,
            broad_id=depmap_id,
            ccle_id=ccle_id,
            cosmic_id=cosmic_id,
            rrid=rrid,
            tissue=tissue,
            cancer_type=cancer_type,
        )
        records.append(record)
        if sanger_id:
            by_sanger[normalize_key(sanger_id)] = record
        if depmap_id:
            by_depmap[normalize_key(depmap_id)] = record

    return records, by_sanger, by_depmap


def load_gdsc_responses(xlsx_path: Path, by_sanger: Dict[str, CrosswalkRecord]) -> Tuple[List[Dict[str, str]], Dict[str, int]]:
    """Standardize GDSC2 response rows using the model crosswalk."""

    aligned_rows: List[Dict[str, str]] = []
    missing: Dict[str, int] = defaultdict(int)
    for row in iter_xlsx_rows(xlsx_path):
        sanger_id = strip_value(row.get("SANGER_MODEL_ID"))
        crosswalk = by_sanger.get(normalize_key(sanger_id))
        if crosswalk is None:
            missing[sanger_id or "__missing__"] += 1
            crosswalk = CrosswalkRecord(
                sanger_model_id=sanger_id,
                depmap_id="",
                cell_line_name=strip_value(row.get("CELL_LINE_NAME")),
                broad_id="",
                ccle_id="",
                cosmic_id="",
                rrid="",
                tissue="",
                cancer_type=strip_value(row.get("CANCER_TYPE")),
            )

        aligned_rows.append(
            {
                "record_type": "gdsc2_response",
                "standard_model_id": crosswalk.depmap_id or crosswalk.sanger_model_id,
                "depmap_id": crosswalk.depmap_id,
                "sanger_model_id": crosswalk.sanger_model_id,
                "cell_line_name": crosswalk.cell_line_name,
                "broad_id": crosswalk.broad_id,
                "ccle_id": crosswalk.ccle_id,
                "cosmic_id": crosswalk.cosmic_id,
                "rrid": crosswalk.rrid,
                "tissue": crosswalk.tissue,
                "cancer_type": crosswalk.cancer_type,
                "dataset": strip_value(row.get("DATASET")),
                "nlme_result_id": strip_value(row.get("NLME_RESULT_ID")),
                "nlme_curve_id": strip_value(row.get("NLME_CURVE_ID")),
                "gdsc_cell_line_name": strip_value(row.get("CELL_LINE_NAME")),
                "gdsc_sanger_model_id": sanger_id,
                "drug_id": strip_value(row.get("DRUG_ID")),
                "drug_name": strip_value(row.get("DRUG_NAME")),
                "putative_target": strip_value(row.get("PUTATIVE_TARGET")),
                "pathway_name": strip_value(row.get("PATHWAY_NAME")),
                "min_conc": strip_value(row.get("MIN_CONC")),
                "max_conc": strip_value(row.get("MAX_CONC")),
                "ln_ic50": strip_value(row.get("LN_IC50")),
                "auc": strip_value(row.get("AUC")),
                "rmse": strip_value(row.get("RMSE")),
                "z_score": strip_value(row.get("Z_SCORE")),
            }
        )
    return aligned_rows, dict(missing)


def stream_omics_table(
    zip_path: Path, source_name: str, record_type: str, counts_by_model: Dict[str, int]
) -> Iterator[Dict[str, str]]:
    """Stream a zipped omics CSV row by row, normalizing model identifiers.

    Yields rows instead of collecting them so callers (i.e. `write_csv`) can stream
    straight to disk. The RNA-seq archive alone unpacks to ~79M rows; materializing
    that as a list of dicts before writing needs on the order of 70GB of RAM, which
    reliably OOMs both local machines and Colab. `counts_by_model` is mutated as a
    side effect during iteration so the caller can still get per-model row counts
    without a second, memory-heavy pass.
    """

    for row in read_csv_from_zip(zip_path):
        standard_model_id = strip_value(row.get("model_id"))
        counts_by_model[normalize_key(standard_model_id)] += 1

        base_row = {
            "record_type": record_type,
            "source_name": source_name,
            "standard_model_id": standard_model_id,
            "depmap_id": standard_model_id,
            "model_name": strip_value(row.get("model_name")),
        }

        if record_type == "rnaseq":
            base_row.update(
                {
                    "dataset_id": strip_value(row.get("dataset_id")),
                    "gene_id": strip_value(row.get("gene_id")),
                    "gene_symbol": strip_value(row.get("gene_symbol")),
                    "ensembl_gene_id": strip_value(row.get("ensembl_gene_id")),
                    "htseq_read_count": strip_value(row.get("htseq_read_count")),
                    "rsem_expected_count": strip_value(row.get("rsem_expected_count")),
                    "rsem_fpkm": strip_value(row.get("rsem_fpkm")),
                    "rsem_tpm": strip_value(row.get("rsem_tpm")),
                    "htseq_fpkm": strip_value(row.get("htseq_fpkm")),
                    "data_source": strip_value(row.get("data_source")),
                    "duplicate": strip_value(row.get("duplicate")),
                }
            )
        elif record_type == "cnv":
            base_row.update(
                {
                    "symbol": strip_value(row.get("symbol")),
                    "gene_id": strip_value(row.get("gene_id")),
                    "total_copy_number": strip_value(row.get("total_copy_number")),
                    "cn_category": strip_value(row.get("cn_category")),
                    "data_type": strip_value(row.get("data_type")),
                    "source": strip_value(row.get("source")),
                }
            )
        elif record_type == "mutation":
            base_row.update(
                {
                    "gene_symbol": strip_value(row.get("gene_symbol")),
                    "ensembl_gene_id": strip_value(row.get("ensembl_gene_id")),
                    "transcript_id": strip_value(row.get("transcript_id")),
                    "protein_mutation": strip_value(row.get("protein_mutation")),
                    "rna_mutation": strip_value(row.get("rna_mutation")),
                    "cdna_mutation": strip_value(row.get("cdna_mutation")),
                    "chromosome": strip_value(row.get("chromosome")),
                    "position": strip_value(row.get("position")),
                    "reference": strip_value(row.get("reference")),
                    "alternative": strip_value(row.get("alternative")),
                    "cancer_driver": strip_value(row.get("cancer_driver")),
                    "cancer_predisposition_variant": strip_value(row.get("cancer_predisposition_variant")),
                    "effect": strip_value(row.get("effect")),
                    "vaf": strip_value(row.get("vaf")),
                    "coding": strip_value(row.get("coding")),
                    "source": strip_value(row.get("source")),
                    "gene_id": strip_value(row.get("gene_id")),
                    "found_in_matched_tumour": strip_value(row.get("found_in_matched_tumour")),
                }
            )
        yield base_row


def load_string_summary(path: Path) -> Dict[str, str]:
    """Collect a small summary for the STRING interaction archive."""

    summary: Dict[str, str] = {"path": str(path), "exists": str(path.exists()).lower()}
    if not path.exists():
        return summary

    with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as handle:
        header = handle.readline().strip()
        first_data = handle.readline().strip()
        summary["header"] = header
        summary["first_data_row"] = first_data
    return summary


def summarize_availability(models: Iterable[CrosswalkRecord], rnaseq_counts: Dict[str, int], cnv_counts: Dict[str, int], mutation_counts: Dict[str, int]) -> List[Dict[str, str]]:
    """Build a per-model availability table for the aligned omics sources."""

    rows: List[Dict[str, str]] = []
    for model in models:
        key = normalize_key(model.depmap_id or model.sanger_model_id)
        rows.append(
            {
                "sanger_model_id": model.sanger_model_id,
                "depmap_id": model.depmap_id,
                "cell_line_name": model.cell_line_name,
                "rnaseq_rows": str(rnaseq_counts.get(key, rnaseq_counts.get(model.depmap_id, 0))),
                "cnv_rows": str(cnv_counts.get(key, cnv_counts.get(model.depmap_id, 0))),
                "mutation_rows": str(mutation_counts.get(key, mutation_counts.get(model.depmap_id, 0))),
            }
        )
    return rows


def build_response_master(
    response_rows: List[Dict[str, str]],
    availability_rows: List[Dict[str, str]],
) -> List[Dict[str, str]]:
    """Add omics availability flags to the aligned GDSC response table."""

    availability_by_depmap = {strip_value(row.get("depmap_id")): row for row in availability_rows}
    master_rows: List[Dict[str, str]] = []
    for row in response_rows:
        depmap_id = strip_value(row.get("depmap_id"))
        availability = availability_by_depmap.get(depmap_id, {})
        enriched = dict(row)
        enriched.update(
            {
                "rnaseq_rows": strip_value(availability.get("rnaseq_rows", "0")),
                "cnv_rows": strip_value(availability.get("cnv_rows", "0")),
                "mutation_rows": strip_value(availability.get("mutation_rows", "0")),
            }
        )
        master_rows.append(enriched)
    return master_rows


def verify_inputs() -> None:
    """Ensure that the required source files are present before processing."""

    required = [MODEL_LIST_PATH, RNASEQ_ZIP_PATH, MUTATIONS_ZIP_PATH, CNV_ZIP_PATH, GDSC_XLSX_PATH, STRING_PATH]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required input files: " + ", ".join(missing))


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the first-pass ingest and alignment workflow."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", type=Path, default=RAW_ROOT, help="Root directory containing the organized raw subfolders.")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_ROOT, help="Directory for aligned outputs.")
    args = parser.parse_args(argv)

    verify_inputs()

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    crosswalk_records, by_sanger, _by_depmap = load_model_crosswalk(MODEL_LIST_PATH, CELL_LINES_PATH)
    response_rows, missing_responses = load_gdsc_responses(GDSC_XLSX_PATH, by_sanger)

    model_crosswalk_path = output_dir / "model_crosswalk.csv"
    response_master_path = output_dir / "gdsc2_response_master.csv"
    rnaseq_path = output_dir / "rnaseq_aligned.csv"
    cnv_path = output_dir / "cnv_aligned.csv"
    mutation_path = output_dir / "mutations_aligned.csv"
    availability_path = output_dir / "omics_availability_by_model.csv"
    string_summary_path = output_dir / "string_summary.json"
    summary_path = output_dir / "alignment_summary.json"

    write_csv(
        model_crosswalk_path,
        (
            {
                "sanger_model_id": record.sanger_model_id,
                "depmap_id": record.depmap_id,
                "cell_line_name": record.cell_line_name,
                "broad_id": record.broad_id,
                "ccle_id": record.ccle_id,
                "cosmic_id": record.cosmic_id,
                "rrid": record.rrid,
                "tissue": record.tissue,
                "cancer_type": record.cancer_type,
            }
            for record in crosswalk_records
        ),
        ["sanger_model_id", "depmap_id", "cell_line_name", "broad_id", "ccle_id", "cosmic_id", "rrid", "tissue", "cancer_type"],
    )

    # Stream each omics archive straight to disk row by row (rnaseq alone is ~79M
    # rows; collecting it into a list first needs ~70GB of RAM and reliably OOMs).
    # counts_by_model is filled in as a side effect of the streamed write.
    rnaseq_counts: Dict[str, int] = defaultdict(int)
    rnaseq_row_count = write_csv(
        rnaseq_path,
        stream_omics_table(RNASEQ_ZIP_PATH, "rnaseq_all_20260323", "rnaseq", rnaseq_counts),
        [
            "record_type",
            "source_name",
            "standard_model_id",
            "depmap_id",
            "model_name",
            "dataset_id",
            "gene_id",
            "gene_symbol",
            "ensembl_gene_id",
            "htseq_read_count",
            "rsem_expected_count",
            "rsem_fpkm",
            "rsem_tpm",
            "htseq_fpkm",
            "data_source",
            "duplicate",
        ],
    )
    cnv_counts: Dict[str, int] = defaultdict(int)
    cnv_row_count = write_csv(
        cnv_path,
        stream_omics_table(CNV_ZIP_PATH, "cnv_summary_20260316", "cnv", cnv_counts),
        [
            "record_type",
            "source_name",
            "standard_model_id",
            "depmap_id",
            "model_name",
            "symbol",
            "gene_id",
            "total_copy_number",
            "cn_category",
            "data_type",
            "source",
        ],
    )
    mutation_counts: Dict[str, int] = defaultdict(int)
    mutation_row_count = write_csv(
        mutation_path,
        stream_omics_table(MUTATIONS_ZIP_PATH, "mutations_all_20260724", "mutation", mutation_counts),
        [
            "record_type",
            "source_name",
            "standard_model_id",
            "depmap_id",
            "model_name",
            "gene_symbol",
            "ensembl_gene_id",
            "transcript_id",
            "protein_mutation",
            "rna_mutation",
            "cdna_mutation",
            "chromosome",
            "position",
            "reference",
            "alternative",
            "cancer_driver",
            "cancer_predisposition_variant",
            "effect",
            "vaf",
            "coding",
            "source",
            "gene_id",
            "found_in_matched_tumour",
        ],
    )

    availability_rows = summarize_availability(crosswalk_records, dict(rnaseq_counts), dict(cnv_counts), dict(mutation_counts))
    response_master_rows = build_response_master(response_rows, availability_rows)

    response_fields = [
        "record_type",
        "standard_model_id",
        "depmap_id",
        "sanger_model_id",
        "cell_line_name",
        "broad_id",
        "ccle_id",
        "cosmic_id",
        "rrid",
        "tissue",
        "cancer_type",
        "dataset",
        "nlme_result_id",
        "nlme_curve_id",
        "gdsc_cell_line_name",
        "gdsc_sanger_model_id",
        "drug_id",
        "drug_name",
        "putative_target",
        "pathway_name",
        "min_conc",
        "max_conc",
        "ln_ic50",
        "auc",
        "rmse",
        "z_score",
        "rnaseq_rows",
        "cnv_rows",
        "mutation_rows",
    ]
    write_csv(response_master_path, response_master_rows, response_fields)

    write_csv(
        availability_path,
        availability_rows,
        ["sanger_model_id", "depmap_id", "cell_line_name", "rnaseq_rows", "cnv_rows", "mutation_rows"],
    )

    with string_summary_path.open("w", encoding="utf-8") as handle:
        json.dump(load_string_summary(STRING_PATH), handle, indent=2, sort_keys=True)

    summary = {
        "inputs": {
            "model_list": str(MODEL_LIST_PATH),
            "cell_lines": str(CELL_LINES_PATH),
            "gdsc2_response": str(GDSC_XLSX_PATH),
            "rnaseq": str(RNASEQ_ZIP_PATH),
            "cnv": str(CNV_ZIP_PATH),
            "mutations": str(MUTATIONS_ZIP_PATH),
            "string": str(STRING_PATH),
        },
        "row_counts": {
            "model_crosswalk": len(crosswalk_records),
            "gdsc2_response_master": len(response_master_rows),
            "rnaseq_aligned": rnaseq_row_count,
            "cnv_aligned": cnv_row_count,
            "mutations_aligned": mutation_row_count,
            "omics_availability_by_model": len(availability_rows),
        },
        "missing_gdsc_model_ids": missing_responses,
        "outputs": {
            "model_crosswalk": str(model_crosswalk_path),
            "gdsc2_response_master": str(response_master_path),
            "rnaseq_aligned": str(rnaseq_path),
            "cnv_aligned": str(cnv_path),
            "mutations_aligned": str(mutation_path),
            "omics_availability_by_model": str(availability_path),
            "string_summary": str(string_summary_path),
        },
    }
    with summary_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)

    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())