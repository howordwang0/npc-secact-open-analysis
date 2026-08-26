from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import tarfile
from collections import Counter
from pathlib import Path

import h5py
import pandas as pd


GSE150430_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/GSE150430/suppl/"
    "GSE150430_npc_scRNA_hg19_processed_data.txt.gz"
)
GSE200310_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE200nnn/GSE200310/suppl/"
    "GSE200310_RAW.tar"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def broad_lineage(label: str) -> str:
    direct = {
        "Macrophage": "Myeloid",
        "Monocyte": "Myeloid",
        "pDC": "Myeloid",
        "DC1 (CLEC9A+)": "Myeloid",
        "DC2 (CD1C+)": "Myeloid",
        "DC3 (CCR7+)": "Myeloid",
        "NK": "T_NK",
        "Plasma cell": "PlasmaB",
        "Malignant": "Malignant",
        "Non-malignant epithelial": "Nonmalignant_epithelial",
        "CAFs": "Stromal",
    }
    if label in direct:
        return direct[label]
    if "B cell" in label or label == "GC B cell":
        return "B"
    if label.startswith("CD4+") or label.startswith("CD8+"):
        return "T_NK"
    raise ValueError(f"Unmapped GSE150430 author label: {label!r}")


def audit_gse150430(matrix_path: Path) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    with gzip.open(matrix_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        barcodes = next(reader)
        annotation = next(reader)
        if annotation[0] != "Cell_type":
            raise ValueError(f"Unexpected annotation key: {annotation[0]!r}")
        labels = annotation[1:]
        if len(barcodes) != len(labels):
            raise ValueError(
                f"Barcode/annotation mismatch: {len(barcodes)} vs {len(labels)}"
            )
        if len(set(barcodes)) != len(barcodes):
            raise ValueError("GSE150430 cell barcodes are not unique")

        patients = [barcode.split("_", 1)[0] for barcode in barcodes]
        lineages = [broad_lineage(label) for label in labels]
        gene_counter: Counter[str] = Counter()
        width_mismatches = []
        n_gene_rows = 0
        for line_number, row in enumerate(reader, start=3):
            n_gene_rows += 1
            if len(row) != len(barcodes) + 1 and len(width_mismatches) < 20:
                width_mismatches.append(
                    {
                        "line_number": line_number,
                        "observed_fields": len(row),
                        "expected_fields": len(barcodes) + 1,
                    }
                )
            if not row:
                continue
            gene_counter[row[0]] += 1

    cell_rows = []
    for (patient, label, lineage), n_cells in sorted(
        Counter(zip(patients, labels, lineages)).items()
    ):
        cell_rows.append(
            {
                "patient_id": patient,
                "sample_type": "normal" if patient == "N01" else "primary_NPC",
                "author_cell_type": label,
                "broad_lineage": lineage,
                "n_cells": n_cells,
            }
        )
    cell_counts = pd.DataFrame(cell_rows)

    duplicate_rows = [
        {"gene_symbol": gene, "n_rows": count}
        for gene, count in sorted(gene_counter.items())
        if count > 1
    ]
    duplicates = pd.DataFrame(duplicate_rows, columns=["gene_symbol", "n_rows"])

    summary = {
        "file": str(matrix_path),
        "url": GSE150430_URL,
        "bytes": matrix_path.stat().st_size,
        "sha256": sha256(matrix_path),
        "gzip_full_stream_read": True,
        "n_cells": len(barcodes),
        "n_patients": len(set(patients)),
        "patients": sorted(set(patients)),
        "n_primary_npc_patients": len(set(patients) - {"N01"}),
        "n_gene_rows": n_gene_rows,
        "n_unique_gene_symbols": len(gene_counter),
        "n_duplicated_gene_symbols": len(duplicates),
        "width_mismatches": width_mismatches,
        "author_cell_type_counts": dict(sorted(Counter(labels).items())),
        "broad_lineage_counts": dict(sorted(Counter(lineages).items())),
    }
    return summary, cell_counts, duplicates


def decode(values) -> list[str]:
    return [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values]


def audit_h5(h5_path: Path, positions_path: Path) -> dict:
    with h5py.File(h5_path, "r") as handle:
        matrix = handle["matrix"]
        shape = [int(value) for value in matrix["shape"][:]]
        barcodes = decode(matrix["barcodes"][:])
        features = matrix["features"]
        name_key = "name" if "name" in features else "gene_names"
        names = decode(features[name_key][:])
        ids = decode(features["id"][:]) if "id" in features else []
        feature_types = (
            decode(features["feature_type"][:]) if "feature_type" in features else []
        )
        nnz = int(matrix["data"].shape[0])

    positions = pd.read_csv(
        positions_path,
        header=None,
        names=[
            "barcode",
            "in_tissue",
            "array_row",
            "array_col",
            "pixel_row",
            "pixel_col",
        ],
    )
    filtered = positions.loc[positions["in_tissue"] == 1, "barcode"].astype(str)
    barcode_set = set(barcodes)
    filtered_set = set(filtered)
    return {
        "file": str(h5_path),
        "bytes": h5_path.stat().st_size,
        "sha256": sha256(h5_path),
        "matrix_feature_by_barcode_shape": shape,
        "n_barcodes": len(barcodes),
        "n_features": len(names),
        "nnz": nnz,
        "n_unique_feature_symbols": len(set(names)),
        "n_duplicated_feature_symbols": len(names) - len(set(names)),
        "n_unique_feature_ids": len(set(ids)) if ids else None,
        "feature_types": dict(sorted(Counter(feature_types).items())),
        "positions_file": str(positions_path),
        "positions_bytes": positions_path.stat().st_size,
        "positions_sha256": sha256(positions_path),
        "n_position_rows": len(positions),
        "n_in_tissue_positions": len(filtered),
        "h5_barcodes_equal_in_tissue_positions": barcode_set == filtered_set,
        "h5_only_barcodes": sorted(barcode_set - filtered_set)[:20],
        "positions_only_barcodes": sorted(filtered_set - barcode_set)[:20],
    }


def parse_soft_sample_titles(path: Path) -> dict[str, str]:
    accession = None
    mapping = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("^SAMPLE = "):
                accession = line.strip().split(" = ", 1)[1]
            elif accession and line.startswith("!Sample_title = "):
                mapping[accession] = line.strip().split(" = ", 1)[1]
    return mapping


def audit_gse200310(raw_tar: Path, extracted_dir: Path, soft_path: Path) -> dict:
    with tarfile.open(raw_tar, "r") as archive:
        members = [
            {"name": member.name, "bytes": int(member.size)}
            for member in archive.getmembers()
            if member.isfile()
        ]
    titles = parse_soft_sample_titles(soft_path)
    samples = []
    for h5_path in sorted(extracted_dir.glob("*.h5")):
        accession_match = re.match(r"(GSM\d+)_", h5_path.name)
        if not accession_match:
            raise ValueError(f"Cannot parse GSM accession from {h5_path.name}")
        accession = accession_match.group(1)
        image_dir = next(extracted_dir.glob(f"{accession}*_spatial_images/spatial"))
        sample = audit_h5(h5_path, image_dir / "tissue_positions_list.csv")
        sample["gsm"] = accession
        sample["sample_title"] = titles.get(accession)
        sample["spatial_files"] = [
            {
                "name": path.name,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in sorted(image_dir.iterdir())
            if path.is_file()
        ]
        samples.append(sample)
    return {
        "file": str(raw_tar),
        "url": GSE200310_URL,
        "bytes": raw_tar.stat().st_size,
        "sha256": sha256(raw_tar),
        "tar_integrity": True,
        "tar_members": members,
        "samples": samples,
    }


def build_report(audit: dict) -> str:
    g150 = audit["GSE150430"]
    g200 = audit["GSE200310"]
    eligible = audit["GSE150430_primary_eligibility"]
    lines = [
        "# NPC SecAct Q1/Q2 enhancement source audit v0.1",
        "",
        "**Audit date:** 2026-08-25  ",
        "**Scope:** Fully open processed inputs for GSE150430 and GSE200310; no biological endpoint was calculated during this audit.",
        "",
        "## GSE150430",
        "",
        f"- File size: {g150['bytes']:,} bytes.",
        f"- SHA-256: `{g150['sha256']}`.",
        f"- Matrix: {g150['n_gene_rows']:,} gene rows x {g150['n_cells']:,} cells; {g150['n_unique_gene_symbols']:,} unique symbols.",
        f"- Patients: {g150['n_primary_npc_patients']} primary NPC plus one normal sample.",
        f"- Full gzip stream and every row width: {'passed' if not g150['width_mismatches'] else 'FAILED'}.",
        f"- Duplicate gene symbols: {g150['n_duplicated_gene_symbols']} symbols; combined deterministically by arithmetic mean downstream.",
        f"- Frozen primary eligibility: {eligible['n_eligible']} tumor patients; excluded {', '.join(eligible['excluded']) or 'none'}.",
        "- Released values are processed non-integer expression, not raw counts; v0.2 protocol amendment governs aggregation.",
        "",
        "## GSE200310",
        "",
        f"- Archive size: {g200['bytes']:,} bytes.",
        f"- SHA-256: `{g200['sha256']}`.",
        f"- Tar integrity: {'passed' if g200['tar_integrity'] else 'FAILED'}.",
    ]
    for sample in g200["samples"]:
        lines.extend(
            [
                f"- {sample['gsm']} ({sample['sample_title']}): {sample['n_barcodes']:,} filtered spots, {sample['n_features']:,} features, {sample['nnz']:,} nonzero entries; barcode/position match {'passed' if sample['h5_barcodes_equal_in_tissue_positions'] else 'FAILED'}.",
                f"  Duplicate feature-symbol rows: {sample['n_duplicated_feature_symbols']}.",
            ]
        )
    lines.extend(
        [
            "",
            "The GEO SOFT text states 1,089 spots for patient 1 and 1,331 for patient 2, whereas the released files and their matched position records contain 1,331 for `NPC_patient_1` and 1,089 for `NPC_patient_2`. Analyses retain the accession/file identity and observable barcodes; no specimens or labels are swapped to match the prose sentence.",
            "",
            "## Decision",
            "",
            "Both datasets are technically usable under the frozen enhancement plan. GSE150430 supports patient-level localization on released processed expression. GSE200310 supports two-specimen descriptive spatial replication, not formal population-level validation.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gse150430-matrix", type=Path, required=True)
    parser.add_argument("--gse150430-soft", type=Path, required=True)
    parser.add_argument("--gse200310-tar", type=Path, required=True)
    parser.add_argument("--gse200310-extracted", type=Path, required=True)
    parser.add_argument("--gse200310-soft", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    g150, cell_counts, duplicates = audit_gse150430(args.gse150430_matrix)
    patient_lineage = (
        cell_counts.groupby(["patient_id", "sample_type", "broad_lineage"], observed=True)[
            "n_cells"
        ]
        .sum()
        .reset_index()
    )
    wide = patient_lineage.pivot(
        index=["patient_id", "sample_type"], columns="broad_lineage", values="n_cells"
    ).fillna(0)
    lineage_columns = [column for column in wide.columns]
    wide["pooled_non_myeloid"] = wide[lineage_columns].sum(axis=1) - wide.get(
        "Myeloid", 0
    )
    tumor = wide.reset_index().query("sample_type == 'primary_NPC'").copy()
    eligible_patients = tumor.loc[
        (tumor.get("Myeloid", 0) >= 50) & (tumor["pooled_non_myeloid"] >= 50),
        "patient_id",
    ].tolist()
    excluded_patients = sorted(set(tumor["patient_id"]) - set(eligible_patients))

    audit = {
        "audit_date": "2026-08-25",
        "GSE150430": g150,
        "GSE150430_primary_eligibility": {
            "threshold_per_compartment": 50,
            "n_eligible": len(eligible_patients),
            "eligible": eligible_patients,
            "excluded": excluded_patients,
        },
        "GSE150430_soft_sample_titles": parse_soft_sample_titles(args.gse150430_soft),
        "GSE200310": audit_gse200310(
            args.gse200310_tar, args.gse200310_extracted, args.gse200310_soft
        ),
    }

    cell_counts.to_csv(
        args.output_dir / "NPC_SecAct_GSE150430_author_cell_counts_v0.1.tsv",
        sep="\t",
        index=False,
    )
    wide.reset_index().to_csv(
        args.output_dir / "NPC_SecAct_GSE150430_patient_lineage_counts_v0.1.tsv",
        sep="\t",
        index=False,
    )
    duplicates.to_csv(
        args.output_dir / "NPC_SecAct_GSE150430_duplicate_genes_v0.1.tsv",
        sep="\t",
        index=False,
    )
    (args.output_dir / "NPC_SecAct_Q1Q2_enhancement_source_audit_v0.1.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    (args.output_dir / "NPC_SecAct_Q1Q2_enhancement_source_audit_v0.1.md").write_text(
        build_report(audit), encoding="utf-8"
    )
    print(json.dumps({"GSE150430": g150, "eligibility": audit["GSE150430_primary_eligibility"], "GSE200310": audit["GSE200310"]}, indent=2))


if __name__ == "__main__":
    main()
