from __future__ import annotations

import argparse
import itertools
import json
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from secactpy import load_signature, secact_activity_inference
from secactpy.inference import group_signatures


SEED = 20260823
PAIRED_PATIENTS = ["P2", "P3", "P4", "P5", "P6"]
ANALYSES = [
    ("hard_qc", "B", 50, "main"),
    ("hard_qc", "Myeloid", 50, "main"),
    ("hard_qc", "T_NK", 50, "main"),
    ("doublet_excluded", "B", 50, "doublet_sensitivity"),
    ("doublet_excluded", "Myeloid", 50, "doublet_sensitivity"),
    ("doublet_excluded", "T_NK", 50, "doublet_sensitivity"),
    ("hard_qc", "PlasmaB", 20, "cell_count_sensitivity"),
    ("doublet_excluded", "PlasmaB", 20, "combined_sensitivity"),
]


def group_level_matrix(expanded_activity: pd.DataFrame) -> pd.DataFrame:
    grouped_signature = group_signatures(load_signature("secact"), cor_threshold=0.9)
    activity_rows = []
    for group_name in grouped_signature.columns:
        representative = group_name.split("|")[0]
        activity_rows.append(expanded_activity.loc[representative].rename(group_name))
    return pd.DataFrame(activity_rows)


def exact_sign_flip_p(differences: np.ndarray) -> float:
    signs = np.asarray(list(itertools.product((-1.0, 1.0), repeat=len(differences))))
    permuted = signs @ differences / len(differences)
    observed = abs(differences.mean())
    return float((np.abs(permuted) >= observed - 1e-12).mean())


def summarize_set(
    scores: pd.DataFrame,
    set_label: str,
    expected_direction: int,
    scenario: str,
    cell_type: str,
    threshold: int,
    analysis_role: str,
) -> dict:
    wide = scores.pivot(index="patient_id", columns="timepoint", values=set_label)
    differences = (wide["OT"] - wide["BT"]).to_numpy(float)
    mean_difference = float(differences.mean())
    sd = float(differences.std(ddof=1))
    cohen_dz = mean_difference / sd if sd else np.nan
    correction = 1 - 3 / (4 * (len(differences) - 1) - 1)
    loo_means = np.array(
        [np.delete(differences, idx).mean() for idx in range(len(differences))]
    )
    observed_direction = int(np.sign(mean_difference))
    return {
        "scenario": scenario,
        "analysis_role": analysis_role,
        "broad_cell_type": cell_type,
        "cell_threshold": threshold,
        "set_label": set_label,
        "expected_direction": expected_direction,
        "n_paired_patients": len(differences),
        "mean_OT_minus_BT": mean_difference,
        "median_OT_minus_BT": float(np.median(differences)),
        "cohen_dz": cohen_dz,
        "hedges_gz": correction * cohen_dz,
        "patients_expected_direction": int(
            (np.sign(differences) == expected_direction).sum()
        ),
        "mean_matches_expected_direction": observed_direction == expected_direction,
        "leave_one_out_matches_expected_direction": bool(
            np.all(np.sign(loo_means) == expected_direction)
        ),
        "exact_sign_flip_p_two_sided": exact_sign_flip_p(differences),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("pseudobulk_dir", type=Path)
    parser.add_argument("cell_counts_tsv", type=Path)
    parser.add_argument("dynamic_sets_tsv", type=Path)
    parser.add_argument("work_output_dir", type=Path)
    parser.add_argument("score_output_tsv", type=Path)
    parser.add_argument("summary_output_tsv", type=Path)
    parser.add_argument("--n-rand", type=int, default=1000)
    parser.add_argument("--rng-method", choices=["srand", "gsl", "numpy"], default="gsl")
    args = parser.parse_args()
    args.work_output_dir.mkdir(parents=True, exist_ok=True)
    args.score_output_tsv.parent.mkdir(parents=True, exist_ok=True)

    cell_counts = pd.read_csv(args.cell_counts_tsv, sep="\t")
    dynamic_sets = pd.read_csv(args.dynamic_sets_tsv, sep="\t")
    expected = {"dynamic_up_set": 1, "dynamic_down_set": -1}
    set_members = {
        label: set(group["signature_group"].astype(str))
        for label, group in dynamic_sets.groupby("set_label")
    }
    if set(set_members) != set(expected):
        raise ValueError(f"Unexpected dynamic set labels: {sorted(set_members)}")

    signature_genes = set(load_signature("secact").index.astype(str))
    expression_cache = {}
    score_rows = []
    summary_rows = []
    run_rows = []

    for scenario, cell_type, threshold, analysis_role in ANALYSES:
        if scenario not in expression_cache:
            expression_cache[scenario] = pd.read_csv(
                args.pseudobulk_dir / f"{scenario}_log2cpm.tsv.gz",
                sep="\t",
                index_col=0,
            )
        all_expression = expression_cache[scenario]
        current_counts = cell_counts.loc[
            (cell_counts["scenario"] == scenario)
            & (cell_counts["broad_cell_type"] == cell_type)
            & (cell_counts["patient_id"].isin(PAIRED_PATIENTS))
        ].copy()
        count_wide = current_counts.pivot(
            index="patient_id", columns="timepoint", values="n_cells"
        )
        eligible_patients = [
            patient
            for patient in PAIRED_PATIENTS
            if patient in count_wide.index
            and count_wide.loc[patient, ["BT", "OT"]].min() >= threshold
        ]
        if len(eligible_patients) < 4:
            raise ValueError(
                f"Frozen analysis unexpectedly has fewer than four pairs: "
                f"{scenario}, {cell_type}, threshold={threshold}"
            )

        selected_rows = current_counts.loc[
            current_counts["patient_id"].isin(eligible_patients)
        ].sort_values(["patient_id", "timepoint"])
        selected_columns = selected_rows["group_id"].tolist()
        expression = all_expression.loc[
            all_expression.index.astype(str).isin(signature_genes), selected_columns
        ].copy()
        expression.columns = selected_rows["sample_id"].tolist()
        if expression.isna().any().any() or (expression.to_numpy() < 0).any():
            raise ValueError(f"Invalid expression values for {scenario}, {cell_type}")

        print(
            f"SECACT START scenario={scenario} cell_type={cell_type} "
            f"threshold={threshold} pairs={len(eligible_patients)}",
            flush=True,
        )
        result = secact_activity_inference(
            expression,
            is_differential=False,
            sig_matrix="secact",
            n_rand=args.n_rand,
            seed=SEED,
            rng_method=args.rng_method,
            sort_genes=True,
            verbose=True,
        )
        group_activity = group_level_matrix(result["zscore"])
        missing_groups = set().union(*set_members.values()) - set(group_activity.index)
        if missing_groups:
            raise ValueError(f"Frozen groups missing from SecAct output: {missing_groups}")

        row_sd = group_activity.std(axis=1, ddof=1)
        if (row_sd == 0).any():
            raise ValueError("Zero-variance SecAct signature group encountered")
        standardized = group_activity.sub(group_activity.mean(axis=1), axis=0).div(
            row_sd, axis=0
        )
        sample_scores = pd.DataFrame(index=standardized.columns)
        for set_label, members in set_members.items():
            sample_scores[set_label] = standardized.loc[sorted(members)].mean(axis=0)
        sample_lookup = selected_rows.set_index("sample_id")
        sample_scores["patient_id"] = sample_lookup.loc[
            sample_scores.index, "patient_id"
        ].to_numpy()
        sample_scores["timepoint"] = sample_lookup.loc[
            sample_scores.index, "timepoint"
        ].to_numpy()
        sample_scores["sample_id"] = sample_scores.index
        sample_scores["scenario"] = scenario
        sample_scores["analysis_role"] = analysis_role
        sample_scores["broad_cell_type"] = cell_type
        sample_scores["cell_threshold"] = threshold
        score_rows.append(sample_scores.reset_index(drop=True))

        for set_label, direction in expected.items():
            summary_rows.append(
                summarize_set(
                    sample_scores,
                    set_label,
                    direction,
                    scenario,
                    cell_type,
                    threshold,
                    analysis_role,
                )
            )

        analysis_dir = args.work_output_dir / f"{scenario}__{cell_type}__ge{threshold}"
        analysis_dir.mkdir(parents=True, exist_ok=True)
        expression.to_csv(
            analysis_dir / "secact_input_log2cpm.tsv.gz", sep="\t", compression="gzip"
        )
        group_activity.to_csv(
            analysis_dir / "secact_group_activity.tsv.gz", sep="\t", compression="gzip"
        )
        result["zscore"].to_csv(
            analysis_dir / "secact_zscore.tsv.gz", sep="\t", compression="gzip"
        )
        run_rows.append(
            {
                "scenario": scenario,
                "analysis_role": analysis_role,
                "broad_cell_type": cell_type,
                "cell_threshold": threshold,
                "n_pairs": len(eligible_patients),
                "patients": eligible_patients,
                "signature_gene_coverage": int(expression.shape[0]),
                "signature_groups": int(group_activity.shape[0]),
            }
        )
        print(f"SECACT COMPLETE scenario={scenario} cell_type={cell_type}", flush=True)

    scores = pd.concat(score_rows, ignore_index=True)
    summaries = pd.DataFrame(summary_rows)
    scores.to_csv(args.score_output_tsv, sep="\t", index=False)
    summaries.to_csv(args.summary_output_tsv, sep="\t", index=False)
    run = {
        "dataset": "Zenodo 13736891 scRNA-seq",
        "secactpy_version": version("secactpy"),
        "seed": SEED,
        "n_rand": args.n_rand,
        "rng_method": args.rng_method,
        "activity_standardization": "within cell type and signature group across eligible paired samples; sample SD",
        "set_aggregation": {label: len(members) for label, members in set_members.items()},
        "analyses": run_rows,
    }
    (args.work_output_dir / "scrna_dynamic_secact_run.json").write_text(
        json.dumps(run, indent=2), encoding="utf-8"
    )
    print("SCRNA DYNAMIC SECACT COMPLETE", flush=True)


if __name__ == "__main__":
    main()
