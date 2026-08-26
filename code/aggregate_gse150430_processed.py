from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from audit_q1q2_enhancement_inputs import broad_lineage


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("matrix", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    with gzip.open(args.matrix, "rt", encoding="utf-8", newline="") as handle:
        barcodes = next(csv.reader(handle, delimiter="\t"))
        annotation = next(csv.reader(handle, delimiter="\t"))
    if annotation[0] != "Cell_type" or len(annotation) != len(barcodes) + 1:
        raise ValueError("Unexpected GSE150430 header/annotation structure")
    labels = annotation[1:]
    patients = np.asarray([barcode.split("_", 1)[0] for barcode in barcodes])
    lineages = np.asarray([broad_lineage(label) for label in labels])

    base_keys = sorted(set(zip(patients, lineages)))
    base_lookup = {key: index for index, key in enumerate(base_keys)}
    base_codes = np.asarray(
        [base_lookup[(patient, lineage)] for patient, lineage in zip(patients, lineages)],
        dtype=np.int32,
    )
    base_counts = np.bincount(base_codes, minlength=len(base_keys)).astype(np.int64)

    output_groups: list[dict] = []
    source_indices: list[np.ndarray] = []
    for index, ((patient, lineage), n_cells) in enumerate(zip(base_keys, base_counts)):
        output_groups.append(
            {
                "group_id": f"{patient}|{lineage}",
                "patient_id": patient,
                "sample_type": "normal" if patient == "N01" else "primary_NPC",
                "compartment": lineage,
                "aggregation_role": "author_broad_lineage",
                "n_cells": int(n_cells),
            }
        )
        source_indices.append(np.asarray([index], dtype=int))

    all_patients = sorted(set(patients))
    derived_definitions = {
        "NonMyeloid": lambda lineage: lineage != "Myeloid",
        "Lymphoid": lambda lineage: lineage in {"T_NK", "B", "PlasmaB"},
        "Epithelial": lambda lineage: lineage
        in {"Malignant", "Nonmalignant_epithelial"},
    }
    for patient in all_patients:
        for compartment, include in derived_definitions.items():
            indices = np.asarray(
                [
                    index
                    for index, (current_patient, lineage) in enumerate(base_keys)
                    if current_patient == patient and include(lineage)
                ],
                dtype=int,
            )
            if not len(indices):
                continue
            n_cells = int(base_counts[indices].sum())
            output_groups.append(
                {
                    "group_id": f"{patient}|{compartment}",
                    "patient_id": patient,
                    "sample_type": "normal" if patient == "N01" else "primary_NPC",
                    "compartment": compartment,
                    "aggregation_role": "derived_pooled_compartment",
                    "n_cells": n_cells,
                }
            )
            source_indices.append(indices)

    manifest = pd.DataFrame(output_groups)
    manifest["eligible_50_cells"] = manifest["n_cells"] >= 50
    manifest["eligible_5_cells"] = manifest["n_cells"] >= 5
    if manifest["group_id"].duplicated().any():
        raise ValueError("Duplicate aggregation group IDs")

    gene_names: list[str] = []
    expression_rows: list[np.ndarray] = []
    minimum = np.inf
    maximum = -np.inf
    nonfinite_rows = []
    width_mismatches = []
    with gzip.open(args.matrix, "rt", encoding="utf-8", newline="") as handle:
        handle.readline()
        handle.readline()
        for line_number, line in enumerate(handle, start=3):
            separator = line.find("\t")
            if separator <= 0:
                raise ValueError(f"No gene/value separator on line {line_number}")
            gene = line[:separator].upper()
            values = np.fromstring(line[separator + 1 :], sep="\t", dtype=np.float64)
            if len(values) != len(barcodes):
                if len(width_mismatches) < 20:
                    width_mismatches.append(
                        {
                            "line": line_number,
                            "gene": gene,
                            "observed": len(values),
                            "expected": len(barcodes),
                        }
                    )
                continue
            if not np.isfinite(values).all():
                nonfinite_rows.append(gene)
            minimum = min(minimum, float(np.nanmin(values)))
            maximum = max(maximum, float(np.nanmax(values)))
            base_sums = np.bincount(
                base_codes, weights=values, minlength=len(base_keys)
            )
            means = np.asarray(
                [
                    base_sums[indices].sum() / base_counts[indices].sum()
                    for indices in source_indices
                ],
                dtype=np.float64,
            )
            gene_names.append(gene)
            expression_rows.append(means)
            if len(gene_names) % 2000 == 0:
                print(f"AGGREGATED genes={len(gene_names)}", flush=True)

    if width_mismatches or nonfinite_rows:
        raise ValueError(
            f"Invalid matrix rows: width={width_mismatches[:3]}, "
            f"nonfinite={nonfinite_rows[:3]}"
        )
    duplicate_counter = Counter(gene_names)
    uppercase_duplicate_symbols = sorted(
        gene for gene, count in duplicate_counter.items() if count > 1
    )
    expression = pd.DataFrame(
        np.vstack(expression_rows),
        index=pd.Index(gene_names, name="gene"),
        columns=manifest["group_id"],
    )
    if uppercase_duplicate_symbols:
        expression = expression.groupby(level=0, sort=False).mean()
    if (expression.to_numpy() < 0).any():
        raise ValueError("Released processed expression contains negative values")

    expression.to_csv(
        args.output_dir / "GSE150430_patient_compartment_mean_processed.tsv.gz",
        sep="\t",
        compression="gzip",
    )
    manifest.to_csv(
        args.output_dir / "GSE150430_patient_compartment_manifest.tsv",
        sep="\t",
        index=False,
    )
    run = {
        "dataset": "GSE150430",
        "input": str(args.matrix),
        "n_cells": len(barcodes),
        "n_source_gene_rows": len(gene_names),
        "n_uppercase_unique_genes": len(expression),
        "n_uppercase_duplicate_symbols": len(uppercase_duplicate_symbols),
        "uppercase_duplicate_symbols": uppercase_duplicate_symbols,
        "n_output_groups": len(manifest),
        "released_value_minimum": minimum,
        "released_value_maximum": maximum,
        "aggregation": "arithmetic mean of released processed cell values",
        "duplicate_gene_rule": "arithmetic mean after uppercase normalization",
        "normalization_after_aggregation": "none",
    }
    (args.output_dir / "GSE150430_aggregation_run.json").write_text(
        json.dumps(run, indent=2), encoding="utf-8"
    )
    print(json.dumps(run, indent=2), flush=True)


if __name__ == "__main__":
    main()
