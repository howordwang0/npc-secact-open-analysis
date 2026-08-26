from __future__ import annotations

import argparse
import csv
import json
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import scanpy as sc


SEED = 20260823


def write_summary(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("expected_counts", type=Path)
    parser.add_argument("extracted_dir", type=Path)
    parser.add_argument("h5ad_dir", type=Path)
    parser.add_argument("summary", type=Path)
    parser.add_argument("run_metadata", type=Path)
    args = parser.parse_args()

    with args.expected_counts.open(encoding="utf-8-sig", newline="") as handle:
        expected = {
            row["archive_sample_id"]: int(row["total_cells"])
            for row in csv.DictReader(handle)
        }
    with args.manifest.open(encoding="utf-8-sig", newline="") as handle:
        manifest = list(csv.DictReader(handle))

    args.h5ad_dir.mkdir(parents=True, exist_ok=True)
    summaries: list[dict[str, object]] = []
    for item in manifest:
        sample_id = item["sample_id"]
        print(f"QC START {sample_id}", flush=True)
        sample_dir = args.extracted_dir / sample_id
        adata = sc.read_10x_mtx(
            sample_dir,
            var_names="gene_symbols",
            make_unique=False,
            gex_only=True,
            cache=False,
        )
        adata.var["gene_symbol"] = adata.var_names.astype(str)
        adata.var_names_make_unique()
        adata.obs_names = [f"{sample_id}_{barcode}" for barcode in adata.obs_names]
        adata.obs["sample_id"] = sample_id
        adata.obs["patient_id"] = item["patient_id"]
        adata.obs["timepoint"] = item["timepoint"]
        adata.var["mt"] = adata.var["gene_symbol"].str.upper().str.startswith("MT-")
        sc.pp.calculate_qc_metrics(
            adata,
            qc_vars=["mt"],
            percent_top=None,
            log1p=False,
            inplace=True,
        )

        hard_mask = (
            adata.obs["n_genes_by_counts"].between(201, 5999, inclusive="both")
            & adata.obs["total_counts"].between(1001, 24999, inclusive="both")
            & (adata.obs["pct_counts_mt"] <= 15.0)
        )
        raw_cells = adata.n_obs
        adata = adata[hard_mask].copy()
        hard_qc_cells = adata.n_obs

        sc.pp.scrublet(
            adata,
            sim_doublet_ratio=2.0,
            expected_doublet_rate=0.033,
            stdev_doublet_rate=0.02,
            synthetic_doublet_umi_subsampling=1.0,
            knn_dist_metric="euclidean",
            normalize_variance=True,
            log_transform=False,
            mean_center=True,
            n_prin_comps=30,
            use_approx_neighbors=None,
            threshold=None,
            verbose=True,
            random_state=SEED,
        )
        auto_predicted_doublets = int(adata.obs["predicted_doublet"].sum())
        fixed_doublet_count = int(round(hard_qc_cells * 0.033))
        score_order = adata.obs["doublet_score"].sort_values(
            ascending=False,
            kind="mergesort",
        ).index
        adata.obs["fixed_top_3p3pct_doublet"] = False
        adata.obs.loc[score_order[:fixed_doublet_count], "fixed_top_3p3pct_doublet"] = True
        fixed_score_cutoff = float(
            adata.obs.loc[score_order[fixed_doublet_count - 1], "doublet_score"]
        ) if fixed_doublet_count else None
        main_retained = hard_qc_cells
        sensitivity_retained = hard_qc_cells - fixed_doublet_count
        author_retained = expected[sample_id]
        adata.write_h5ad(args.h5ad_dir / f"{sample_id}.hard_qc.h5ad", compression="gzip")

        summaries.append(
            {
                "patient_id": item["patient_id"],
                "timepoint": item["timepoint"],
                "sample_id": sample_id,
                "raw_barcodes": raw_cells,
                "hard_qc_cells": hard_qc_cells,
                "scrublet_auto_threshold": adata.uns["scrublet"]["threshold"],
                "scrublet_auto_predicted_doublets": auto_predicted_doublets,
                "fixed_top_3p3pct_doublets": fixed_doublet_count,
                "fixed_top_3p3pct_score_cutoff": fixed_score_cutoff,
                "main_retained_cells": main_retained,
                "sensitivity_retained_cells": sensitivity_retained,
                "author_retained_cells": author_retained,
                "main_minus_author": main_retained - author_retained,
                "sensitivity_minus_author": sensitivity_retained - author_retained,
            }
        )
        write_summary(summaries, args.summary)
        print(
            f"QC VERIFIED {sample_id}: raw={raw_cells}, hard={hard_qc_cells}, "
            f"auto_doublets={auto_predicted_doublets}, fixed_doublets={fixed_doublet_count}, "
            f"main={main_retained}, sensitivity={sensitivity_retained}, author={author_retained}",
            flush=True,
        )

    metadata = {
        "seed": SEED,
        "scanpy_version": version("scanpy"),
        "anndata_version": version("anndata"),
        "qc": {
            "n_genes_by_counts": [201, 5999],
            "total_counts": [1001, 24999],
            "pct_counts_mt_max": 15.0,
        },
        "scrublet": {
            "sim_doublet_ratio": 2.0,
            "expected_doublet_rate": 0.033,
            "stdev_doublet_rate": 0.02,
            "synthetic_doublet_umi_subsampling": 1.0,
            "knn_dist_metric": "euclidean",
            "normalize_variance": True,
            "log_transform": False,
            "mean_center": True,
            "n_prin_comps": 30,
            "use_approx_neighbors": None,
            "threshold": None,
        },
        "doublet_analysis": {
            "main": "hard_qc_without_doublet_exclusion",
            "sensitivity": "exclude_top_3.3_percent_scrublet_scores_per_sample",
            "fixed_fraction": 0.033,
            "fixed_count_rule": "round(n_hard_qc * 0.033)",
        },
    }
    args.run_metadata.parent.mkdir(parents=True, exist_ok=True)
    args.run_metadata.write_text(json.dumps(metadata, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
