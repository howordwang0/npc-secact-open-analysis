from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


SCENARIOS = {
    "hard_qc": lambda obs: np.ones(len(obs), dtype=bool),
    "doublet_excluded": lambda obs: ~obs["fixed_top_3p3pct_doublet"].to_numpy(bool),
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("qc_dir", type=Path)
    parser.add_argument("cell_metadata", type=Path)
    parser.add_argument("annotation_tsv", type=Path)
    parser.add_argument("work_output_dir", type=Path)
    parser.add_argument("cell_count_output", type=Path)
    args = parser.parse_args()
    args.work_output_dir.mkdir(parents=True, exist_ok=True)
    args.cell_count_output.parent.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(args.cell_metadata, sep="\t", index_col=0)
    annotation = pd.read_csv(args.annotation_tsv, sep="\t")
    annotation = annotation.loc[annotation["include_in_lineage_secact"]].copy()
    cluster_to_type = dict(
        zip(annotation["cluster"].astype(int), annotation["broad_cell_type"])
    )

    paths = sorted(args.qc_dir.glob("*.hard_qc.h5ad"))
    if len(paths) != 11:
        raise ValueError(f"Expected 11 hard-QC h5ad files, found {len(paths)}")

    first = ad.read_h5ad(paths[0], backed="r")
    gene_symbols = first.var["gene_symbol"].astype(str).to_numpy()
    first.file.close()
    valid_gene = pd.notna(gene_symbols) & (gene_symbols != "") & (gene_symbols != "nan")
    unique_symbols = np.array(sorted(set(gene_symbols[valid_gene])))
    symbol_to_code = {symbol: idx for idx, symbol in enumerate(unique_symbols)}
    feature_codes = np.array(
        [symbol_to_code.get(symbol, -1) for symbol in gene_symbols], dtype=np.int64
    )

    matrices = {scenario: {} for scenario in SCENARIOS}
    cell_rows = []
    sample_rows = []

    for path in paths:
        adata = ad.read_h5ad(path)
        if not np.array_equal(gene_symbols, adata.var["gene_symbol"].astype(str).to_numpy()):
            raise ValueError(f"Gene symbol order mismatch in {path.name}")
        sample_meta = metadata.reindex(adata.obs_names)
        if sample_meta["leiden_0p5"].isna().any():
            raise ValueError(f"Missing cluster metadata for cells in {path.name}")
        cell_types = sample_meta["leiden_0p5"].astype(int).map(cluster_to_type)
        matrix = adata.X.tocsr() if sparse.issparse(adata.X) else sparse.csr_matrix(adata.X)
        sample_id = str(adata.obs["sample_id"].iloc[0])
        patient_id = str(adata.obs["patient_id"].iloc[0])
        timepoint = str(adata.obs["timepoint"].iloc[0])
        sample_rows.append(
            {"sample_id": sample_id, "patient_id": patient_id, "timepoint": timepoint}
        )

        for scenario, scenario_filter in SCENARIOS.items():
            retained = scenario_filter(adata.obs)
            for cell_type in sorted(annotation["broad_cell_type"].unique()):
                selected = retained & (cell_types.to_numpy() == cell_type)
                n_cells = int(selected.sum())
                group_id = f"{sample_id}|{cell_type}"
                if n_cells:
                    feature_counts = np.asarray(matrix[selected].sum(axis=0)).ravel()
                    symbol_counts = np.bincount(
                        feature_codes[valid_gene],
                        weights=feature_counts[valid_gene],
                        minlength=len(unique_symbols),
                    )
                else:
                    symbol_counts = np.zeros(len(unique_symbols), dtype=float)
                matrices[scenario][group_id] = symbol_counts
                cell_rows.append(
                    {
                        "scenario": scenario,
                        "patient_id": patient_id,
                        "timepoint": timepoint,
                        "sample_id": sample_id,
                        "broad_cell_type": cell_type,
                        "group_id": group_id,
                        "n_cells": n_cells,
                        "total_umis": int(symbol_counts.sum()),
                        "eligible_50_cells": n_cells >= 50,
                        "eligible_20_cells": n_cells >= 20,
                    }
                )
        print(f"AGGREGATED {sample_id}", flush=True)

    for scenario, group_counts in matrices.items():
        counts = pd.DataFrame(group_counts, index=unique_symbols)
        library_sizes = counts.sum(axis=0)
        if (library_sizes == 0).any():
            raise ValueError(f"Zero-count pseudobulk group in {scenario}")
        log2cpm = np.log2(counts.divide(library_sizes, axis=1) * 1_000_000 + 1)
        counts.to_csv(
            args.work_output_dir / f"{scenario}_gene_counts.tsv.gz",
            sep="\t",
            compression="gzip",
        )
        log2cpm.to_csv(
            args.work_output_dir / f"{scenario}_log2cpm.tsv.gz",
            sep="\t",
            compression="gzip",
        )

    cell_counts = pd.DataFrame(cell_rows)
    cell_counts.to_csv(args.cell_count_output, sep="\t", index=False)
    sample_manifest = pd.DataFrame(sample_rows).drop_duplicates()
    sample_manifest.to_csv(
        args.work_output_dir / "sample_manifest.tsv", sep="\t", index=False
    )
    run = {
        "n_samples": int(sample_manifest.shape[0]),
        "n_patients": int(sample_manifest["patient_id"].nunique()),
        "paired_patients": ["P2", "P3", "P4", "P5", "P6"],
        "n_unique_gene_symbols": int(len(unique_symbols)),
        "scenarios": list(SCENARIOS),
        "normalization": "log2(CPM+1), denominator includes all nonempty gene symbols",
        "ambiguous_clusters_excluded": [10, 14],
    }
    (args.work_output_dir / "pseudobulk_run.json").write_text(
        json.dumps(run, indent=2), encoding="utf-8"
    )
    print("PSEUDOBULK COMPLETE", flush=True)


if __name__ == "__main__":
    main()
