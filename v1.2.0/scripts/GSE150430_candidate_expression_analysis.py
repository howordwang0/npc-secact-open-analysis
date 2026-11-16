#!/usr/bin/env python3
"""Summarize candidate and WNT receptor-context expression by annotated cell type.

The input is the author-processed GSE150430 gene-by-cell matrix. This script
uses the original cell-type labels, summarizes each patient-cell-type group,
and then reports medians across patients. It does not recluster cells, perform
count-based differential expression, or infer ligand-receptor interactions.

Example:
  python GSE150430_candidate_expression_analysis.py \
    --matrix GSE150430_npc_scRNA_hg19_processed_data.txt.gz \
    --outdir GSE150430_candidate_expression_results
"""
from __future__ import annotations
import argparse
import csv
import gzip
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path

MYELOID_LABELS = {
    "Macrophage": "Macrophage",
    "Monocyte": "Monocyte",
    "DC1 (CLEC9A+)": "DC1_CLEC9A+",
    "DC2 (CD1C+)": "DC2_CD1C+",
    "DC3 (CCR7+)": "DC3_CCR7+",
    "pDC": "pDC",
}
GENES = ["WNT7A", "VSTM2L", *[f"FZD{i}" for i in range(1, 11)],
         "LRP5", "LRP6", "ROR1", "ROR2"]
PATIENTS = {f"P{i:02d}" for i in range(2, 16)}
MIN_CELLS = 20


def open_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("rt", encoding="utf-8", newline="")


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True,
                        help="Author-processed GSE150430 gene-by-cell matrix (plain text or .gz).")
    parser.add_argument("--outdir", type=Path, required=True,
                        help="Directory for reproducible TSV/JSON outputs.")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    with open_text(args.matrix) as handle:
        barcodes = next(csv.reader(handle, delimiter="\t"))
        annotation = next(csv.reader(handle, delimiter="\t"))
    if annotation[0] != "Cell_type" or len(annotation) != len(barcodes) + 1:
        raise ValueError("Unexpected matrix header/Cell_type annotation structure")
    cell_types = annotation[1:]
    patient_ids = [barcode.split("_", 1)[0] for barcode in barcodes]
    keep = [patient in PATIENTS for patient in patient_ids]
    if sum(keep) == 0:
        raise ValueError("No P02-P15 patient cells found in the matrix")
    selected_patients = sorted(set(p for p, use in zip(patient_ids, keep) if use))
    if selected_patients != sorted(PATIENTS):
        raise ValueError(f"Expected patient IDs P02-P15, observed {selected_patients}")

    group_n = Counter((p, ct) for p, ct, use in zip(patient_ids, cell_types, keep) if use)
    accum = defaultdict(lambda: [0.0, 0])
    found = set()
    with open_text(args.matrix) as handle:
        next(handle)
        next(handle)
        for line_number, line in enumerate(handle, start=3):
            row = next(csv.reader([line], delimiter="\t"))
            if len(row) != len(barcodes) + 1:
                raise ValueError(f"Line {line_number} has {len(row)-1} values; expected {len(barcodes)}")
            gene = row[0].strip().upper()
            if gene not in GENES:
                continue
            if gene in found:
                raise ValueError(f"Duplicate target gene row: {gene}")
            found.add(gene)
            for i, value in enumerate(row[1:]):
                if not keep[i]:
                    continue
                x = float(value)
                if not math.isfinite(x) or x < 0:
                    raise ValueError(f"Invalid processed expression for {gene} on line {line_number}")
                a = accum[(patient_ids[i], cell_types[i], gene)]
                a[0] += x
                if x > 0:
                    a[1] += 1
    missing = sorted(set(GENES) - found)
    if missing:
        raise ValueError(f"Target genes absent from matrix: {missing}")

    patient_rows = []
    celltype_patients = defaultdict(lambda: defaultdict(list))
    for (patient, cell_type), n_cells in sorted(group_n.items()):
        if n_cells < MIN_CELLS:
            continue
        for gene in GENES:
            total, detected = accum[(patient, cell_type, gene)]
            mean = total / n_cells
            detection = detected / n_cells
            row = {
                "patient_id": patient,
                "cell_type": cell_type,
                "myeloid_subtype": MYELOID_LABELS.get(cell_type, ""),
                "gene": gene,
                "n_cells": n_cells,
                "n_detected": detected,
                "detection_fraction": detection,
                "mean_processed_expression": mean,
            }
            patient_rows.append(row)
            celltype_patients[(cell_type, gene)][patient] = row

    summary_rows = []
    for (cell_type, gene), patients in sorted(celltype_patients.items()):
        rows = list(patients.values())
        detections = [r["detection_fraction"] for r in rows]
        means = [r["mean_processed_expression"] for r in rows]
        sizes = [r["n_cells"] for r in rows]
        summary_rows.append({
            "cell_type": cell_type,
            "myeloid_subtype": MYELOID_LABELS.get(cell_type, ""),
            "gene": gene,
            "n_patients_with_at_least_20_cells": len(rows),
            "patients_with_detected_expression": sum(x > 0 for x in detections),
            "median_patient_detection_fraction": statistics.median(detections),
            "median_patient_mean_processed_expression": statistics.median(means),
            "minimum_cells_per_patient_group": min(sizes),
            "maximum_cells_per_patient_group": max(sizes),
            "minimum_cell_threshold": MIN_CELLS,
        })
    myeloid_rows = [r for r in summary_rows if r["cell_type"] in MYELOID_LABELS]

    write_tsv(args.outdir / "patient_level_celltype_expression.tsv", patient_rows,
              ["patient_id", "cell_type", "myeloid_subtype", "gene", "n_cells", "n_detected",
               "detection_fraction", "mean_processed_expression"])
    write_tsv(args.outdir / "celltype_expression_summary.tsv", summary_rows,
              ["cell_type", "myeloid_subtype", "gene", "n_patients_with_at_least_20_cells",
               "patients_with_detected_expression", "median_patient_detection_fraction",
               "median_patient_mean_processed_expression", "minimum_cells_per_patient_group",
               "maximum_cells_per_patient_group", "minimum_cell_threshold"])
    write_tsv(args.outdir / "myeloid_subtype_expression_summary.tsv", myeloid_rows,
              ["cell_type", "myeloid_subtype", "gene", "n_patients_with_at_least_20_cells",
               "patients_with_detected_expression", "median_patient_detection_fraction",
               "median_patient_mean_processed_expression", "minimum_cells_per_patient_group",
               "maximum_cells_per_patient_group", "minimum_cell_threshold"])

    metadata = {
        "dataset": "GSE150430",
        "input_matrix": args.matrix.name,
        "cell_count_in_matrix": len(barcodes),
        "original_cell_type_count": len(set(cell_types)),
        "patients_in_analysis": selected_patients,
        "patient_selection": "NPC patient IDs P02-P15; excludes normal N01 and P01 to match the existing 14-patient NPC compartment analysis",
        "minimum_cells_per_patient_cell_type": MIN_CELLS,
        "genes": GENES,
        "myeloid_annotations": MYELOID_LABELS,
        "expression_definition": "nonzero means detected; arithmetic mean of author-processed matrix values",
        "summary_definition": "median across patient-cell-type groups with at least 20 cells",
        "interpretation": "descriptive transcript context only; no raw-count differential expression, reclustering, protein abundance, or ligand-receptor interaction inference",
    }
    (args.outdir / "analysis_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps({"patients": selected_patients, "target_genes_found": sorted(found),
                      "patient_level_rows": len(patient_rows),
                      "celltype_summary_rows": len(summary_rows), "outdir": str(args.outdir.resolve())}, indent=2))


if __name__ == "__main__":
    main()
