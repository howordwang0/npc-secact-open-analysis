from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats
from secactpy import load_signature, secact_activity_inference

from analyze_gse206245_visium_spatial import (
    enable_process_gsl_permutation_cache,
    group_level_matrix,
    standardize_rows,
)
from run_cytosig_benchmark import infer as cytosig_infer
from run_direct_expression_benchmark import direct_set_scores, read_definitions


SEED = 20260823
N_BOOTSTRAP = 2000
SAMPLES = [
    {
        "gsm": "GSM6030901",
        "sample": "NPC_patient_1",
        "h5": "GSM6030901_patient_1_filtered_feature_bc_matrix.h5",
        "spatial_dir": "GSM6030901_patient1_spatial_images/spatial",
        "expected_spots": 1331,
    },
    {
        "gsm": "GSM6030902",
        "sample": "NPC_patient_2",
        "h5": "GSM6030902_patient_2_filtered_feature_bc_matrix.h5",
        "spatial_dir": "GSM6030902_patient2_spatial_images/spatial",
        "expected_spots": 1089,
    },
]


def decode(values) -> list[str]:
    return [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values]


def read_h5_slice(root: Path, metadata: dict) -> dict:
    h5_path = root / metadata["h5"]
    with h5py.File(h5_path, "r") as handle:
        matrix = handle["matrix"]
        shape = tuple(int(value) for value in matrix["shape"][:])
        counts = sparse.csc_matrix(
            (
                matrix["data"][:],
                matrix["indices"][:],
                matrix["indptr"][:],
            ),
            shape=shape,
        ).tocsr()
        barcodes = np.asarray(decode(matrix["barcodes"][:]))
        features = matrix["features"]
        key = "name" if "name" in features else "gene_names"
        symbols = pd.Series(decode(features[key][:])).str.upper()

    unique_symbols = pd.Index(pd.unique(symbols), dtype=str)
    codes = pd.Categorical(symbols, categories=unique_symbols).codes
    collapse = sparse.csr_matrix(
        (
            np.ones(len(codes), dtype=np.float64),
            (codes, np.arange(len(codes))),
        ),
        shape=(len(unique_symbols), len(codes)),
    )
    counts = (collapse @ counts).tocsr()
    if counts.shape[1] != metadata["expected_spots"]:
        raise ValueError(f"Unexpected spot count for {metadata['sample']}: {counts.shape}")

    position_path = root / metadata["spatial_dir"] / "tissue_positions_list.csv"
    positions = pd.read_csv(
        position_path,
        header=None,
        names=[
            "barcode",
            "in_tissue",
            "array_row",
            "array_col",
            "pixel_row",
            "pixel_col",
        ],
    ).set_index("barcode")
    positions = positions.loc[barcodes]
    if not positions["in_tissue"].eq(1).all():
        raise ValueError(f"Filtered H5 includes out-of-tissue spots for {metadata['sample']}")
    totals = np.asarray(counts.sum(axis=0)).ravel().astype(float)
    detected = np.asarray(counts.getnnz(axis=0)).ravel().astype(int)
    if (totals <= 0).any():
        raise ValueError(f"Zero-library spot in {metadata['sample']}")
    return {
        "genes": unique_symbols,
        "counts": counts,
        "barcodes": barcodes,
        "positions": positions,
        "total_counts": totals,
        "detected_features": detected,
    }


def block_bootstrap(
    score: np.ndarray,
    context: np.ndarray,
    rows: np.ndarray,
    cols: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    def bin_four(values: np.ndarray) -> np.ndarray:
        low = float(values.min())
        high = float(values.max())
        if high == low:
            return np.zeros(len(values), dtype=int)
        return np.minimum(3, np.floor((values - low) / (high - low) * 4).astype(int))

    block_ids = bin_four(rows) * 4 + bin_four(cols)
    blocks = [np.flatnonzero(block_ids == block) for block in sorted(set(block_ids))]
    result = np.empty(N_BOOTSTRAP, dtype=float)
    for iteration in range(N_BOOTSTRAP):
        selected = rng.integers(0, len(blocks), size=len(blocks))
        indices = np.concatenate([blocks[index] for index in selected])
        result[iteration] = stats.spearmanr(score[indices], context[indices]).statistic
    return result


def toroidal_shift_null(
    score: np.ndarray,
    context: np.ndarray,
    rows: np.ndarray,
    cols: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    row_min, row_max = int(rows.min()), int(rows.max())
    col_min, col_max = int(cols.min()), int(cols.max())
    height = row_max - row_min + 1
    width = col_max - col_min + 1
    coordinate_to_index = {
        (int(row), int(col)): index for index, (row, col) in enumerate(zip(rows, cols))
    }
    threshold = int(np.ceil(0.70 * len(score)))
    null_values = []
    overlap_values = []
    for row_shift in range(height):
        for col_shift in range(width):
            if row_shift == 0 and col_shift == 0:
                continue
            source_indices = []
            target_indices = []
            for source, (row, col) in enumerate(zip(rows, cols)):
                target = (
                    row_min + (int(row) - row_min + row_shift) % height,
                    col_min + (int(col) - col_min + col_shift) % width,
                )
                if target in coordinate_to_index:
                    source_indices.append(source)
                    target_indices.append(coordinate_to_index[target])
            if len(source_indices) < threshold:
                continue
            rho = stats.spearmanr(
                score[np.asarray(source_indices)], context[np.asarray(target_indices)]
            ).statistic
            null_values.append(float(rho))
            overlap_values.append(len(source_indices))
    return np.asarray(null_values), np.asarray(overlap_values)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("extracted_dir", type=Path)
    parser.add_argument("dynamic_sets", type=Path)
    parser.add_argument("group_members", type=Path)
    parser.add_argument("myeloid_markers", type=Path)
    parser.add_argument("cytosig_summary", type=Path)
    parser.add_argument("work_output_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--n-rand", type=int, default=1000)
    parser.add_argument(
        "--samples",
        nargs="+",
        choices=[metadata["sample"] for metadata in SAMPLES],
    )
    args = parser.parse_args()
    args.work_output_dir.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    enable_process_gsl_permutation_cache()

    dynamic = pd.read_csv(args.dynamic_sets, sep="\t")
    set_members = {
        label: set(subset["signature_group"].astype(str))
        for label, subset in dynamic.groupby("set_label")
    }
    direct_sets, direct_members = read_definitions(args.dynamic_sets, args.group_members)
    marker_table = pd.read_csv(args.myeloid_markers, sep="\t")
    myeloid_markers = marker_table["gene"].astype(str).str.upper().tolist()
    cytosig_summary = pd.read_csv(args.cytosig_summary, sep="\t")
    support_flag = cytosig_summary["passes_frozen_spatial_support"].map(
        lambda value: str(value).lower() == "true"
    )
    supported_cytokines = cytosig_summary.loc[support_flag, "cytokine"].astype(str).tolist()
    if len(supported_cytokines) != 29:
        raise ValueError(f"Expected 29 frozen supported cytokines, found {len(supported_cytokines)}")

    secact_signature = load_signature("secact")
    secact_signature.index = secact_signature.index.astype(str).str.upper()
    cytosig_signature = load_signature("cytosig")
    cytosig_signature.index = cytosig_signature.index.astype(str).str.upper()
    direct_genes = set().union(
        *(set(direct_members[group]) for groups in direct_sets.values() for group in groups)
    )
    needed_genes = (
        set(secact_signature.index)
        | set(cytosig_signature.index)
        | set(myeloid_markers)
        | direct_genes
    )

    score_frames = []
    summary_rows = []
    cytokine_rows = []
    null_rows = []
    coverage_rows = []
    rng = np.random.default_rng(SEED)

    samples_to_run = SAMPLES
    if args.samples:
        requested = set(args.samples)
        samples_to_run = [
            metadata for metadata in SAMPLES if metadata["sample"] in requested
        ]
    for metadata in samples_to_run:
        sample = metadata["sample"]
        print(f"LOAD START {sample}", flush=True)
        released = read_h5_slice(args.extracted_dir, metadata)
        genes = released["genes"]
        lookup = pd.Series(np.arange(len(genes)), index=genes)
        observed = sorted(set(genes) & needed_genes)
        row_indices = lookup.loc[observed].to_numpy(int)
        counts = released["counts"][row_indices].toarray().astype(np.float64)
        log2cpm = np.log2(
            counts / released["total_counts"][None, :] * 1_000_000 + 1
        )
        expression = pd.DataFrame(log2cpm, index=observed, columns=released["barcodes"])

        secact_input = expression.loc[expression.index.intersection(secact_signature.index)]
        secact_checkpoint = args.work_output_dir / f"{metadata['gsm']}_secact_group_activity.tsv.gz"
        if secact_checkpoint.exists():
            group_activity = pd.read_csv(secact_checkpoint, sep="\t", index_col=0)
            if group_activity.columns.tolist() != expression.columns.tolist():
                raise ValueError(f"Invalid SecAct checkpoint columns for {sample}")
            print(f"SECACT CHECKPOINT {sample}", flush=True)
        else:
            print(f"SECACT START {sample}", flush=True)
            secact_result = secact_activity_inference(
                secact_input,
                is_differential=False,
                sig_matrix="secact",
                n_rand=args.n_rand,
                seed=SEED,
                rng_method="gsl",
                sort_genes=True,
                verbose=True,
            )
            group_activity = group_level_matrix(secact_result["zscore"])
            group_activity.to_csv(
                secact_checkpoint, sep="\t", compression="gzip"
            )
        standardized_activity, omitted_groups = standardize_rows(group_activity)
        score = pd.DataFrame(index=expression.columns)
        for label, members in set_members.items():
            retained = sorted(set(standardized_activity.index) & members)
            if len(retained) != len(members):
                raise ValueError(f"Invariant/missing frozen groups in {sample}, {label}")
            score[label] = standardized_activity.loc[retained].mean(axis=0)

        marker_expression = expression.loc[
            expression.index.intersection(myeloid_markers)
        ]
        standardized_markers, omitted_markers = standardize_rows(marker_expression)
        if len(standardized_markers) < 16:
            raise ValueError(f"Myeloid-context coverage failed for {sample}")
        score["myeloid_context"] = standardized_markers.mean(axis=0)

        direct, direct_coverage = direct_set_scores(
            expression, direct_sets, direct_members, excluded_genes=set(myeloid_markers)
        )
        direct = direct.rename(
            columns={
                "dynamic_up_set": "direct_dynamic_up_set",
                "dynamic_down_set": "direct_dynamic_down_set",
            }
        )
        score = score.join(direct)

        cytosig_checkpoint = args.work_output_dir / f"{metadata['gsm']}_cytosig_activity.tsv.gz"
        if cytosig_checkpoint.exists():
            cytosig = pd.read_csv(cytosig_checkpoint, sep="\t", index_col=0)
            if cytosig.columns.tolist() != expression.columns.tolist():
                raise ValueError(f"Invalid CytoSig checkpoint columns for {sample}")
            print(f"CYTOSIG CHECKPOINT {sample}", flush=True)
        else:
            print(f"CYTOSIG START {sample}", flush=True)
            cytosig = cytosig_infer(expression, cytosig_signature)
            cytosig.to_csv(cytosig_checkpoint, sep="\t", compression="gzip")
        missing_cytokines = set(supported_cytokines) - set(cytosig.index)
        if missing_cytokines:
            raise ValueError(f"Frozen CytoSig cytokines missing: {missing_cytokines}")
        standardized_cytosig, omitted_cytokines = standardize_rows(cytosig)
        score["cytosig_supported_set"] = standardized_cytosig.loc[
            supported_cytokines
        ].mean(axis=0)

        positions = released["positions"]
        rows = positions["array_row"].to_numpy(int)
        cols = positions["array_col"].to_numpy(int)
        context = score["myeloid_context"].to_numpy(float)
        programs = [
            "dynamic_up_set",
            "dynamic_down_set",
            "direct_dynamic_up_set",
            "direct_dynamic_down_set",
            "cytosig_supported_set",
        ]
        for program in programs:
            values = score[program].to_numpy(float)
            observed_rho = float(stats.spearmanr(values, context).statistic)
            bootstrap = block_bootstrap(values, context, rows, cols, rng)
            shift_null, overlaps = toroidal_shift_null(values, context, rows, cols)
            empirical_p = (
                float((1 + (np.abs(shift_null) >= abs(observed_rho)).sum()) / (1 + len(shift_null)))
                if len(shift_null)
                else np.nan
            )
            summary_rows.append(
                {
                    **metadata,
                    "program": program,
                    "n_spots": len(score),
                    "spearman_rho": observed_rho,
                    "spatial_block_bootstrap_ci95_lower": float(np.quantile(bootstrap, 0.025)),
                    "spatial_block_bootstrap_ci95_upper": float(np.quantile(bootstrap, 0.975)),
                    "n_accepted_toroidal_shifts": len(shift_null),
                    "minimum_shift_overlap_spots": int(overlaps.min()) if len(overlaps) else 0,
                    "toroidal_shift_p_two_sided": empirical_p,
                }
            )
            for iteration, value in enumerate(shift_null):
                null_rows.append(
                    {
                        "sample": sample,
                        "program": program,
                        "shift_index": iteration,
                        "spearman_rho_null": value,
                        "overlap_spots": int(overlaps[iteration]),
                    }
                )

        for cytokine in cytosig.index:
            cytokine_rows.append(
                {
                    **metadata,
                    "cytokine": cytokine,
                    "spearman_rho": float(
                        stats.spearmanr(cytosig.loc[cytokine], context).statistic
                    ),
                    "in_frozen_29_set": cytokine in supported_cytokines,
                }
            )

        score.insert(0, "barcode", score.index)
        score.insert(0, "sample", sample)
        score["array_row"] = rows
        score["array_col"] = cols
        score["total_umi"] = released["total_counts"]
        score["detected_genes"] = released["detected_features"]
        score_frames.append(score.reset_index(drop=True))
        for _, row in direct_coverage.iterrows():
            coverage_rows.append(
                {
                    "sample": sample,
                    "method": "direct_expression",
                    **row.to_dict(),
                }
            )
        coverage_rows.extend(
            [
                {
                    "sample": sample,
                    "method": "SecAct",
                    "set_label": "signature_genes",
                    "total_groups": len(secact_signature),
                    "represented_groups": len(secact_input),
                    "represented_genes": len(secact_input),
                },
                {
                    "sample": sample,
                    "method": "myeloid_context",
                    "set_label": "markers",
                    "total_groups": 20,
                    "represented_groups": len(standardized_markers),
                    "represented_genes": len(standardized_markers),
                },
                {
                    "sample": sample,
                    "method": "CytoSig",
                    "set_label": "signature_genes",
                    "total_groups": len(cytosig_signature),
                    "represented_groups": len(expression.index.intersection(cytosig_signature.index)),
                    "represented_genes": len(expression.index.intersection(cytosig_signature.index)),
                },
            ]
        )
        print(f"COMPLETE {sample}", flush=True)

    scores = pd.concat(score_frames, ignore_index=True)
    summary = pd.DataFrame(summary_rows)
    cytokines = pd.DataFrame(cytokine_rows)
    nulls = pd.DataFrame(null_rows)
    coverage = pd.DataFrame(coverage_rows)
    scores.to_csv(
        args.work_output_dir / "GSE200310_spot_scores.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    summary.to_csv(
        args.output_dir / "NPC_SecAct_GSE200310_external_spatial_summary_v0.1.tsv",
        sep="\t",
        index=False,
    )
    cytokines.to_csv(
        args.output_dir / "NPC_SecAct_GSE200310_CytoSig_correlations_v0.1.tsv",
        sep="\t",
        index=False,
    )
    nulls.to_csv(
        args.work_output_dir / "GSE200310_toroidal_shift_nulls.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    coverage.to_csv(
        args.output_dir / "NPC_SecAct_GSE200310_score_coverage_v0.1.tsv",
        sep="\t",
        index=False,
    )

    up = summary.loc[summary["program"] == "dynamic_up_set"].copy()
    direct_up = summary.loc[summary["program"] == "direct_dynamic_up_set"].copy()
    concordant = bool(len(up) == 2 and (up["spearman_rho"] > 0).all())
    report = f"""# GSE200310 external spatial replication report v0.1

**Analysis date:** 2026-08-25  
**Role:** Locked descriptive replication in two fully open primary NPC Visium specimens.

## Frozen treatment-up SecAct program

| Specimen | Spots | Spearman rho | Spatial-block 95% CI | Toroidal-shift P |
|---|---:|---:|---:|---:|
{chr(10).join(f"| {row.sample} | {int(row.n_spots):,} | {row.spearman_rho:.3f} | {row.spatial_block_bootstrap_ci95_lower:.3f} to {row.spatial_block_bootstrap_ci95_upper:.3f} | {row.toroidal_shift_p_two_sided:.4f} |" for row in up.itertuples())}

Both specimens positive under the frozen descriptive concordance rule: {'yes' if concordant else 'no'}.

## Direct-expression comparator

The corresponding direct treatment-up correlations were {', '.join(f"{row.sample}: {row.spearman_rho:.3f}" for row in direct_up.itertuples())}. These are retained alongside SecAct and determine whether any spatial association is method-specific.

## Interpretation boundary

With two specimens, this analysis is descriptive external replication rather than formal population validation. It does not establish secretion, cellular source, treatment effect, clinical prediction, or SecAct superiority. Treatment-down, direct-expression, frozen 29-cytokine CytoSig, all 43 individual cytokines, and spatial-null results are retained in the accompanying tables.
"""
    (args.output_dir / "NPC_SecAct_GSE200310_external_spatial_report_v0.1.md").write_text(
        report, encoding="utf-8"
    )
    run = {
        "dataset": "GSE200310",
        "secactpy_version": version("secactpy"),
        "seed": SEED,
        "n_rand": args.n_rand,
        "n_spatial_block_bootstrap": N_BOOTSTRAP,
        "n_samples": len(samples_to_run),
        "frozen_cytosig_set_size": len(supported_cytokines),
        "secact_up_positive_both": concordant,
    }
    (args.work_output_dir / "GSE200310_external_spatial_run.json").write_text(
        json.dumps(run, indent=2), encoding="utf-8"
    )
    print(summary.to_string(index=False), flush=True)
    print(json.dumps(run, indent=2), flush=True)


if __name__ == "__main__":
    main()
