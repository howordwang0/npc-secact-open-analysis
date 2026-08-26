from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("qc_dir", type=Path)
    parser.add_argument("cell_metadata", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(args.cell_metadata, sep="\t", index_col=0)
    clusters = sorted(metadata["leiden_0p5"].astype(int).unique())
    cluster_to_row = {cluster: row for row, cluster in enumerate(clusters)}

    paths = sorted(args.qc_dir.glob("*.hard_qc.h5ad"))
    if len(paths) != 11:
        raise ValueError(f"Expected 11 hard-QC h5ad files, found {len(paths)}")

    counts_sum = None
    detected_cells = None
    cluster_cells = np.zeros(len(clusters), dtype=np.int64)
    var_names = None

    for path in paths:
        adata = ad.read_h5ad(path)
        if var_names is None:
            var_names = adata.var_names.to_numpy().copy()
            counts_sum = np.zeros((len(clusters), adata.n_vars), dtype=np.float64)
            detected_cells = np.zeros((len(clusters), adata.n_vars), dtype=np.int64)
        elif not np.array_equal(var_names, adata.var_names.to_numpy()):
            raise ValueError(f"Feature order mismatch in {path.name}")

        sample_meta = metadata.reindex(adata.obs_names)
        if sample_meta["leiden_0p5"].isna().any():
            raise ValueError(f"Missing cluster metadata for cells in {path.name}")
        sample_clusters = sample_meta["leiden_0p5"].astype(int).to_numpy()
        matrix = adata.X.tocsr() if sparse.issparse(adata.X) else sparse.csr_matrix(adata.X)

        for cluster in np.unique(sample_clusters):
            row = cluster_to_row[int(cluster)]
            mask = sample_clusters == cluster
            cluster_matrix = matrix[mask]
            counts_sum[row] += np.asarray(cluster_matrix.sum(axis=0)).ravel()
            detected_cells[row] += np.asarray((cluster_matrix > 0).sum(axis=0)).ravel()
            cluster_cells[row] += int(mask.sum())

        print(f"AGGREGATED {path.stem}: cells={adata.n_obs}", flush=True)

    assert counts_sum is not None
    assert detected_cells is not None
    assert var_names is not None
    if int(cluster_cells.sum()) != len(metadata):
        raise ValueError("Aggregated cell total does not match cell metadata")

    total_counts = counts_sum.sum(axis=1)
    cpm = counts_sum / total_counts[:, None] * 1_000_000
    total_gene_counts = counts_sum.sum(axis=0)
    total_umis = total_counts.sum()
    rows = []
    top_rows = []

    for row, cluster in enumerate(clusters):
        other_counts = total_gene_counts - counts_sum[row]
        other_umis = total_umis - total_counts[row]
        other_cpm = other_counts / other_umis * 1_000_000
        log2fc = np.log2((cpm[row] + 0.5) / (other_cpm + 0.5))
        pct_cluster = detected_cells[row] / cluster_cells[row]
        other_detected = detected_cells.sum(axis=0) - detected_cells[row]
        other_cells = cluster_cells.sum() - cluster_cells[row]
        pct_other = other_detected / other_cells
        pct_diff = pct_cluster - pct_other

        eligible = (pct_cluster >= 0.10) & (pct_diff >= 0.05)
        order = np.lexsort((-pct_diff, -log2fc))
        order = [idx for idx in order if eligible[idx]][:30]
        for rank, idx in enumerate(order, start=1):
            top_rows.append(
                {
                    "cluster": cluster,
                    "rank": rank,
                    "gene": var_names[idx],
                    "log2fc_cluster_vs_rest": log2fc[idx],
                    "pct_cluster": pct_cluster[idx],
                    "pct_rest": pct_other[idx],
                    "cpm_cluster": cpm[row, idx],
                    "cpm_rest": other_cpm[idx],
                }
            )
        rows.append(
            {
                "cluster": cluster,
                "n_cells": cluster_cells[row],
                "total_umis": int(total_counts[row]),
                "n_top_markers_reported": len(order),
            }
        )

    pd.DataFrame(rows).to_csv(
        args.output_dir / "cluster_raw_count_summary.tsv", sep="\t", index=False
    )
    pd.DataFrame(top_rows).to_csv(
        args.output_dir / "cluster_top_markers_descriptive.tsv", sep="\t", index=False
    )
    print("MARKER AUDIT COMPLETE", flush=True)


if __name__ == "__main__":
    main()
