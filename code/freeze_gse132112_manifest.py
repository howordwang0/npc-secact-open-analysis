from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd


def normalize(value: str) -> str:
    value = str(value).strip().lower()
    value = re.sub(r"\.(cel|txt|gz)$", "", value)
    return re.sub(r"[^a-z0-9]+", "", value)


def aliases(row: pd.Series) -> set[str]:
    result = set()
    for field in ("geo_accession", "title", "source_name_ch1"):
        value = str(row.get(field, "")).strip()
        if value and value.lower() != "nan":
            result.add(normalize(value))
            match = re.search(r"LA[-_ ]?NPC[-_ ]?(\d+)", value, flags=re.IGNORECASE)
            if match:
                result.add(f"lanpc{int(match.group(1))}")
    return result


def match_expression_columns(clinical: pd.DataFrame, columns: list[str]) -> list[int]:
    alias_to_rows: dict[str, list[int]] = {}
    for index, row in clinical.iterrows():
        for alias in aliases(row):
            alias_to_rows.setdefault(alias, []).append(index)

    matched_rows: list[int] = []
    unmatched: list[str] = []
    ambiguous: dict[str, list[int]] = {}
    for column in columns:
        column_aliases = {normalize(column)}
        match = re.search(r"LA[-_ ]?NPC[-_ ]?(\d+)", column, flags=re.IGNORECASE)
        if match:
            column_aliases.add(f"lanpc{int(match.group(1))}")
        candidates = sorted({index for alias in column_aliases for index in alias_to_rows.get(alias, [])})
        if len(candidates) == 1:
            matched_rows.append(candidates[0])
        elif not candidates:
            unmatched.append(column)
        else:
            ambiguous[column] = candidates
    if unmatched or ambiguous:
        raise ValueError(f"Expression/clinical ID match failed; unmatched={unmatched}, ambiguous={ambiguous}")
    if len(set(matched_rows)) != len(matched_rows):
        raise ValueError("More than one expression column matched the same clinical row")
    return matched_rows


def main() -> None:
    parser = argparse.ArgumentParser(description="Freeze the phenotype manifest for GSE132112.")
    parser.add_argument("expression_tsv_gz", type=Path)
    parser.add_argument("geo_samples_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()

    columns = pd.read_csv(args.expression_tsv_gz, sep="\t", nrows=0).columns.tolist()[1:]
    clinical = pd.read_csv(args.geo_samples_csv, dtype=str)
    row_indices = match_expression_columns(clinical, columns)
    matched = clinical.loc[row_indices].copy().reset_index(drop=True)
    matched.insert(0, "expression_sample_id", columns)

    response_map = {
        "complete": "CR",
        "complete response": "CR",
        "partial response": "PR",
        "stable disease": "SD",
    }
    matched["response_original"] = matched["response"].str.strip().str.lower().map(response_map)
    if matched["response_original"].isna().any():
        unknown = matched.loc[matched["response_original"].isna(), "response"].unique().tolist()
        raise ValueError(f"Unknown response labels: {unknown}")
    matched["primary_response_cr_pr_vs_sd"] = matched["response_original"].isin(["CR", "PR"]).astype(int)
    matched["sensitivity_cr_vs_non_cr"] = (matched["response_original"] == "CR").astype(int)
    manifest = matched[
        [
            "expression_sample_id",
            "geo_accession",
            "title",
            "response_original",
            "primary_response_cr_pr_vs_sd",
            "sensitivity_cr_vs_non_cr",
            "platform_id",
            "source_name_ch1",
        ]
    ]
    if manifest["geo_accession"].nunique() != len(manifest):
        raise ValueError("GEO accessions are not unique at the patient/sample level")

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(args.output_csv, index=False, encoding="utf-8-sig")
    digest = hashlib.sha256(args.output_csv.read_bytes()).hexdigest()
    summary = {
        "expression_file": str(args.expression_tsv_gz),
        "geo_samples_file": str(args.geo_samples_csv),
        "manifest_file": str(args.output_csv),
        "patients": len(manifest),
        "unique_gsm": int(manifest["geo_accession"].nunique()),
        "response_counts": manifest["response_original"].value_counts().sort_index().to_dict(),
        "primary_endpoint": "CR/PR versus SD",
        "sensitivity_endpoint": "CR versus PR/SD",
        "sha256": digest,
    }
    Path(f"{args.output_csv}.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
