from __future__ import annotations

import argparse
import json
import sqlite3
from pathlib import Path

import pandas as pd
from secactpy import load_signature


def main() -> None:
    parser = argparse.ArgumentParser(description="Map HTA 2.0 core PSR features to unique gene symbols.")
    parser.add_argument("design_sqlite", type=Path)
    parser.add_argument("transcript_to_symbol_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()

    transcript_map = pd.read_csv(args.transcript_to_symbol_csv, dtype=str)
    transcript_map["is_multiple"] = pd.to_numeric(transcript_map["is_multiple"], errors="coerce").fillna(1).astype(int)
    transcript_map = transcript_map.dropna(subset=["probe_id", "symbol"])
    symbol_counts = transcript_map.groupby("probe_id")["symbol"].nunique()
    multiple_flags = transcript_map.groupby("probe_id")["is_multiple"].max()
    eligible_transcripts = symbol_counts.index[(symbol_counts == 1) & (multiple_flags == 0)]
    transcript_map = transcript_map[transcript_map["probe_id"].isin(eligible_transcripts)].drop_duplicates("probe_id")
    transcript_to_symbol = transcript_map.set_index("probe_id")["symbol"]

    connection = sqlite3.connect(args.design_sqlite)
    design = pd.read_sql_query(
        """
        SELECT DISTINCT f.man_fsetid AS probe_id, m.meta_fsetid AS transcript_cluster_id
        FROM featureSet AS f
        JOIN core_mps AS m USING (fsetid)
        WHERE m.meta_fsetid LIKE 'TC%'
        """,
        connection,
    )
    connection.close()
    transcript_counts = design.groupby("probe_id")["transcript_cluster_id"].nunique()
    unique_features = transcript_counts.index[transcript_counts == 1]
    design_unique = design[design["probe_id"].isin(unique_features)].drop_duplicates("probe_id")
    design_unique["symbol"] = design_unique["transcript_cluster_id"].map(transcript_to_symbol)
    mapped = design_unique.dropna(subset=["symbol"]).copy()
    mapped["is_multiple"] = 0
    mapped = mapped[["probe_id", "transcript_cluster_id", "symbol", "is_multiple"]].sort_values("probe_id")

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    mapped.to_csv(args.output_csv, index=False)
    signature_genes = set(load_signature("secact").index.astype(str))
    mapped_symbols = set(mapped["symbol"])
    summary = {
        "design_sqlite": str(args.design_sqlite),
        "transcript_annotation": str(args.transcript_to_symbol_csv),
        "core_feature_transcript_pairs": int(len(design)),
        "core_feature_ids": int(design["probe_id"].nunique()),
        "ambiguous_core_feature_ids_excluded": int((transcript_counts > 1).sum()),
        "unambiguous_core_features_with_gene_symbol": int(len(mapped)),
        "unique_transcript_clusters_with_gene_symbol": int(mapped["transcript_cluster_id"].nunique()),
        "unique_gene_symbols": int(mapped["symbol"].nunique()),
        "secact_signature_symbols_mappable": int(len(mapped_symbols & signature_genes)),
        "secact_signature_coverage_percent": round(100 * len(mapped_symbols & signature_genes) / len(signature_genes), 3),
        "output_csv": str(args.output_csv),
        "rule": "core_mps only; exclude PSRs assigned to more than one TC; exclude TCs with multiple gene symbols",
    }
    Path(f"{args.output_csv}.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
