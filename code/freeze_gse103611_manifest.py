from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description="Freeze the GSE103611 metastasis manifest.")
    parser.add_argument("expression_tsv_gz", type=Path)
    parser.add_argument("geo_samples_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()

    columns = pd.read_csv(args.expression_tsv_gz, sep="\t", nrows=0).columns.tolist()[1:]
    samples = pd.read_csv(args.geo_samples_csv, dtype=str).set_index("geo_accession")
    if set(columns) != set(samples.index):
        raise ValueError("Expression columns do not match GEO accessions exactly")
    manifest = samples.loc[columns].reset_index()
    label_map = {"yes": 1, "no": 0}
    manifest["metastasis_binary"] = manifest["metastasis_status"].str.lower().map(label_map)
    if manifest["metastasis_binary"].isna().any():
        raise ValueError("Unknown metastasis status")
    manifest.insert(0, "expression_sample_id", columns)
    manifest = manifest[["expression_sample_id", "geo_accession", "title", "metastasis_status", "metastasis_binary", "platform_id"]]
    if manifest["geo_accession"].nunique() != 48:
        raise ValueError("Expected 48 unique patients")

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(args.output_csv, index=False, encoding="utf-8-sig")
    digest = hashlib.sha256(args.output_csv.read_bytes()).hexdigest()
    summary = {
        "patients": len(manifest),
        "metastasis_counts": manifest["metastasis_status"].value_counts().to_dict(),
        "sample_order_blocks": {"metastasis": "1-24", "non_metastasis": "25-48"},
        "sha256": digest,
    }
    Path(f"{args.output_csv}.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
