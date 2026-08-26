from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
META_PATH = ROOT / "work" / "data_audit" / "GSE102349_samples.csv"
EXPR_PATH = ROOT / "work" / "data_audit" / "GSE102349_NPC_mRNA_processed.txt.gz"
OUT_DIR = ROOT / "outputs"
WORK_DIR = ROOT / "work" / "data_audit"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_patient(text: object) -> str | None:
    match = re.search(r"NPC\d+PT\d+T\d+", str(text))
    return match.group(0) if match else None


def main() -> None:
    meta = pd.read_csv(META_PATH, dtype=str, keep_default_na=False)
    expression_header = pd.read_csv(EXPR_PATH, sep="\t", compression="gzip", nrows=0).columns.tolist()
    expression_columns = expression_header[1:]
    meta["patient_id"] = meta["description"].map(parse_patient)

    endpoint_time = pd.to_numeric(meta["time_to_event"].replace({"": np.nan, "N/A": np.nan}), errors="coerce")
    event_text = meta["event"].str.strip()
    event = event_text.eq("Disease progression").astype("Int64")
    event[event_text.eq("Last follow-up")] = 0
    event[event_text.isin(["", "N/A", "Not evaluated"])] = pd.NA

    mapped = pd.DataFrame({"expression_column": expression_columns})
    mapped["patient_id"] = mapped["expression_column"].map(parse_patient)
    mapped = mapped.merge(
        meta[
            [
                "geo_accession",
                "patient_id",
                "title",
                "event",
                "time_to_event",
                "clinical_stage",
                "cohort",
            ]
        ].assign(numeric_time=endpoint_time, event_code=event),
        on="patient_id",
        how="left",
        validate="one_to_one",
    )
    mapped["expression_column_present"] = mapped["expression_column"].isin(expression_columns)

    # Exact accession overlap is only a lower-bound check; patient-source overlap
    # cannot be excluded from GEO metadata and is therefore reported as unknown.
    prior_accessions: dict[str, set[str]] = {}
    for accession, path in {
        "GSE132112": ROOT / "outputs" / "NPC_SecAct_GSE132112_frozen_manifest_v0.3.csv",
        "GSE103611": ROOT / "outputs" / "NPC_SecAct_GSE103611_frozen_manifest_v0.3.csv",
    }.items():
        if path.exists():
            frame = pd.read_csv(path, dtype=str)
            values = set(frame.astype(str).stack().tolist())
            prior_accessions[accession] = set(mapped["geo_accession"].dropna()) & values
        else:
            prior_accessions[accession] = set()

    counts = {
        "metadata_rows": int(len(meta)),
        "expression_columns": int(len(expression_columns)),
        "expression_columns_with_patient_mapping": int(mapped["patient_id"].notna().sum()),
        "expression_columns_without_patient_mapping": int(mapped["patient_id"].isna().sum()),
        "unique_patient_ids": int(meta["patient_id"].nunique(dropna=True)),
        "duplicate_patient_ids": int(meta["patient_id"].duplicated(keep=False).sum()),
        "numeric_time": int(endpoint_time.notna().sum()),
        "progression_events": int((event == 1).sum()),
        "censored_last_followup": int((event == 0).sum()),
        "missing_endpoint": int(event.isna().sum()),
        "positive_time_values": int((endpoint_time > 0).sum()),
        "nonpositive_time_values": int((endpoint_time <= 0).sum()),
    }
    stage_counts = meta["clinical_stage"].replace("", "Missing").value_counts(dropna=False).to_dict()
    cohort_counts = meta["cohort"].replace("", "Missing").value_counts(dropna=False).to_dict()
    paired_counts = meta["cohort"].replace("", "Missing").value_counts(dropna=False).to_dict()

    audit = {
        "dataset": "GSE102349",
        "audit_date": "2026-08-26",
        "metadata_path": str(META_PATH),
        "expression_path": str(EXPR_PATH),
        "sha256": {"metadata": sha256(META_PATH), "expression": sha256(EXPR_PATH)},
        "counts": counts,
        "clinical_stage_counts": {str(k): int(v) for k, v in stage_counts.items()},
        "cohort_counts": {str(k): int(v) for k, v in cohort_counts.items()},
        "paired_field_counts": {str(k): int(v) for k, v in paired_counts.items()},
        "exact_geo_accession_overlap_with_prior_manifests": {
            key: sorted(value) for key, value in prior_accessions.items()
        },
        "source_patient_overlap": "unknown from GEO metadata; no exact patient-level provenance field is released",
        "endpoint_unit_status": "requires source-publication reconciliation; analysis may use released time_to_event units with this limitation",
        "protocol_gates": {
            "one_to_one_expression_mapping": bool(
                len(expression_columns) == len(meta)
                and mapped["patient_id"].notna().all()
                and mapped["patient_id"].is_unique
            ),
            "exact_endpoint_counts": counts["numeric_time"] == 88
            and counts["progression_events"] == 16
            and counts["censored_last_followup"] == 72
            and counts["missing_endpoint"] == 25,
            "no_duplicate_patients": counts["duplicate_patient_ids"] == 0,
            "minimum_event_gate": counts["progression_events"] >= 10,
        },
    }
    OUT_DIR.mkdir(exist_ok=True)
    WORK_DIR.mkdir(exist_ok=True)
    mapped.to_csv(WORK_DIR / "GSE102349_patient_mapping_audit_v0.1.csv", index=False)
    (WORK_DIR / "GSE102349_protocol_audit_v0.1.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (OUT_DIR / "NPC_SecAct_GSE102349_protocol_audit_report_v0.1.md").write_text(
        "# GSE102349 protocol audit report v0.1\n\n"
        "Audit performed before any molecular-to-endpoint association.\n\n"
        f"- Metadata rows: {counts['metadata_rows']}\n"
        f"- Expression columns: {counts['expression_columns']}\n"
        f"- Unique patient IDs: {counts['unique_patient_ids']}\n"
        f"- Numeric time values: {counts['numeric_time']}\n"
        f"- Progression events: {counts['progression_events']}\n"
        f"- Censored last-follow-up records: {counts['censored_last_followup']}\n"
        f"- Missing endpoint records: {counts['missing_endpoint']}\n"
        f"- Clinical stage distribution: {json.dumps(audit['clinical_stage_counts'], ensure_ascii=False)}\n"
        f"- Cohort field distribution: {json.dumps(audit['cohort_counts'], ensure_ascii=False)}\n"
        f"- Paired-field distribution: {json.dumps(audit['paired_field_counts'], ensure_ascii=False)}\n\n"
        "## Gate status\n\n"
        f"- One-to-one expression/patient mapping: {audit['protocol_gates']['one_to_one_expression_mapping']}\n"
        f"- Exact endpoint-count reproduction: {audit['protocol_gates']['exact_endpoint_counts']}\n"
        f"- No duplicate patient IDs: {audit['protocol_gates']['no_duplicate_patients']}\n"
        f"- At least 10 progression events: {audit['protocol_gates']['minimum_event_gate']}\n\n"
        "The GEO metadata do not prove that underlying patients or institutions are independent of every prior cohort; this remains an explicit limitation.\n",
        encoding="utf-8",
    )
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
