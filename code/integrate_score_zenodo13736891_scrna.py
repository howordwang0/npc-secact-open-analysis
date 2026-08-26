from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce


SEED = 20260823

LINEAGE_MARKERS = {
    "T_NK": ["CD3D", "CD3E", "CD3G", "TRAC", "CD2", "NKG7", "CCL5", "IL32"],
    "B": ["MS4A1", "CD79A", "CD79B", "BANK1", "CD74", "CD22", "VPREB3"],
    "Myeloid": ["LYZ", "AIF1", "TYROBP", "FCER1G", "CTSS", "CST3", "S100A8", "S100A9", "C1QA", "C1QB"],
    "PlasmaB": ["MZB1", "JCHAIN", "XBP1", "DERL3", "FKBP11", "IGHG1", "IGKC"],
    "Epi_Cancer": ["EPCAM", "KRT19", "KRT18", "KRT8", "KRT5", "KRT14", "FGFBP1", "CRABP2"],
    "Fibroblasts": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRB", "COL6A1", "COL6A2", "POSTN"],
    "Endothelial": ["VWF", "PECAM1", "KDR", "EMCN", "CLDN5", "RAMP2", "RNASE1", "ACKR1"],
    "Mast": ["TPSB2", "TPSAB1", "CPA3", "MS4A2", "KIT", "HPGDS", "HDC"],
    "RBC": ["HBB", "HBA1", "HBA2", "ALAS2", "AHSP"],
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("qc_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    paths = sorted(args.qc_dir.glob("*.hard_qc.h5ad"))
    if len(paths) != 11:
        raise ValueError(f"Expected 11 hard-QC h5ad files, found {len(paths)}")

    objects = [ad.read_h5ad(path) for path in paths]
    adata = ad.concat(objects, join="inner", merge="same", index_unique=None)
    if not adata.obs_names.is_unique:
        raise ValueError("Cell names are not unique after concatenation")
    if adata.n_vars != 60623:
        raise ValueError(f"Expected 60,623 features, found {adata.n_vars}")
    print(f"CONCAT cells={adata.n_obs} features={adata.n_vars}", flush=True)

    sc.pp.normalize_total(adata, target_sum=10000)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=2000,
        flavor="seurat",
        batch_key="sample_id",
        inplace=True,
    )
    hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(hvg, max_value=10)
    sc.tl.pca(hvg, n_comps=50, svd_solver="arpack", random_state=SEED)
    sce.pp.harmony_integrate(
        hvg,
        key="sample_id",
        basis="X_pca",
        adjusted_basis="X_pca_harmony",
        max_iter_harmony=20,
        random_state=SEED,
    )
    sc.pp.neighbors(
        hvg,
        n_neighbors=15,
        n_pcs=20,
        use_rep="X_pca_harmony",
        random_state=SEED,
    )
    sc.tl.leiden(
        hvg,
        resolution=0.5,
        random_state=SEED,
        key_added="leiden_0p5",
        flavor="leidenalg",
        n_iterations=-1,
    )
    sc.tl.umap(hvg, random_state=SEED)
    adata.obs["leiden_0p5"] = hvg.obs["leiden_0p5"].astype(str)
    adata.obsm["X_pca_harmony"] = hvg.obsm["X_pca_harmony"][:, :20]
    adata.obsm["X_umap"] = hvg.obsm["X_umap"]
    print(f"CLUSTERS n={adata.obs['leiden_0p5'].nunique()}", flush=True)

    all_markers = sorted({gene for genes in LINEAGE_MARKERS.values() for gene in genes})
    present_markers = [gene for gene in all_markers if gene in adata.var_names]
    missing_markers = sorted(set(all_markers) - set(present_markers))
    marker_values = pd.DataFrame(
        adata[:, present_markers].X.toarray(),
        index=adata.obs_names,
        columns=present_markers,
    )
    marker_values["cluster"] = adata.obs["leiden_0p5"].to_numpy()
    cluster_mean = marker_values.groupby("cluster", observed=True).mean()
    cluster_pct = marker_values.assign(
        **{gene: marker_values[gene] > 0 for gene in present_markers}
    ).groupby("cluster", observed=True)[present_markers].mean()

    gene_sd = cluster_mean.std(axis=0, ddof=0).replace(0, np.nan)
    cluster_z = (cluster_mean - cluster_mean.mean(axis=0)) / gene_sd
    lineage_scores = pd.DataFrame(index=cluster_mean.index)
    for lineage, genes in LINEAGE_MARKERS.items():
        usable = [gene for gene in genes if gene in cluster_z.columns]
        lineage_scores[lineage] = cluster_z[usable].mean(axis=1)
    ordered_scores = np.sort(lineage_scores.to_numpy(), axis=1)
    lineage_scores["best_lineage"] = lineage_scores.idxmax(axis=1)
    lineage_scores["score_margin"] = ordered_scores[:, -1] - ordered_scores[:, -2]
    lineage_scores["n_cells"] = adata.obs.groupby("leiden_0p5", observed=True).size()

    composition = (
        adata.obs.groupby(
            ["leiden_0p5", "patient_id", "timepoint", "sample_id"],
            observed=True,
        )
        .size()
        .rename("n_cells")
        .reset_index()
    )
    cell_metadata = adata.obs[
        [
            "patient_id",
            "timepoint",
            "sample_id",
            "n_genes_by_counts",
            "total_counts",
            "pct_counts_mt",
            "doublet_score",
            "predicted_doublet",
            "fixed_top_3p3pct_doublet",
            "leiden_0p5",
        ]
    ].copy()
    cell_metadata["UMAP_1"] = adata.obsm["X_umap"][:, 0]
    cell_metadata["UMAP_2"] = adata.obsm["X_umap"][:, 1]

    cluster_mean.to_csv(args.output_dir / "cluster_marker_logmean.tsv", sep="\t")
    cluster_pct.to_csv(args.output_dir / "cluster_marker_fraction.tsv", sep="\t")
    lineage_scores.to_csv(args.output_dir / "cluster_lineage_scores.tsv", sep="\t")
    composition.to_csv(args.output_dir / "cluster_sample_composition.tsv", sep="\t", index=False)
    cell_metadata.to_csv(
        args.output_dir / "cell_metadata_preannotation.tsv.gz",
        sep="\t",
        compression="gzip",
    )
    hvg.write_h5ad(args.output_dir / "integrated_hvg.h5ad", compression="gzip")

    run = {
        "seed": SEED,
        "scanpy_version": version("scanpy"),
        "anndata_version": version("anndata"),
        "harmonypy_version": version("harmonypy"),
        "igraph_version": version("igraph"),
        "leidenalg_version": version("leidenalg"),
        "n_cells": adata.n_obs,
        "n_features": adata.n_vars,
        "n_hvg": hvg.n_vars,
        "n_clusters": int(adata.obs["leiden_0p5"].nunique()),
        "missing_markers": missing_markers,
        "lineage_markers": LINEAGE_MARKERS,
    }
    (args.output_dir / "integration_run.json").write_text(
        json.dumps(run, indent=2),
        encoding="utf-8",
    )
    print("INTEGRATION COMPLETE", flush=True)


if __name__ == "__main__":
    main()
