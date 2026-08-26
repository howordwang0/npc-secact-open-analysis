from __future__ import annotations

import argparse
import itertools
import json
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from secactpy import load_signature, secact_activity_inference
from secactpy.inference import group_signatures


SEED = 20260823


def group_level_matrix(expanded_activity: pd.DataFrame) -> pd.DataFrame:
    grouped = group_signatures(load_signature("secact"), cor_threshold=0.9)
    return pd.DataFrame(
        [
            expanded_activity.loc[group.split("|")[0]].rename(group)
            for group in grouped.columns
        ]
    )


def standardize_rows(values: pd.DataFrame) -> pd.DataFrame:
    sd = values.std(axis=1, ddof=1)
    retained = sd.index[sd.gt(0) & np.isfinite(sd)]
    return values.loc[retained].sub(values.loc[retained].mean(axis=1), axis=0).div(
        sd.loc[retained], axis=0
    )


def read_direct_definitions(
    dynamic_sets_path: Path, group_members_path: Path
) -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    dynamic = pd.read_csv(dynamic_sets_path, sep="\t")
    members = pd.read_csv(group_members_path, sep="\t")
    sets = {
        label: subset["signature_group"].astype(str).tolist()
        for label, subset in dynamic.groupby("set_label", sort=False)
    }
    group_members = {
        group: sorted(set(subset["protein"].astype(str).str.upper()))
        for group, subset in members.groupby("signature_group", sort=False)
    }
    return sets, group_members


def direct_scores(
    expression: pd.DataFrame,
    sets: dict[str, list[str]],
    group_members: dict[str, list[str]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    standardized = standardize_rows(expression)
    score_columns = {}
    coverage = []
    for label, groups in sets.items():
        rows = []
        genes_observed = set()
        for group in groups:
            genes = [gene for gene in group_members[group] if gene in standardized.index]
            if genes:
                rows.append(standardized.loc[genes].mean(axis=0).rename(group))
                genes_observed.update(genes)
        if not rows:
            raise ValueError(f"No direct-expression groups represented for {label}")
        score_columns[f"direct_{label}"] = pd.DataFrame(rows).mean(axis=0)
        coverage.append(
            {
                "set_label": label,
                "total_groups": len(groups),
                "represented_groups": len(rows),
                "represented_genes": len(genes_observed),
            }
        )
    return pd.DataFrame(score_columns), pd.DataFrame(coverage)


def exact_signed_rank_interval(
    differences: np.ndarray, alpha: float = 0.05
) -> tuple[float, float, float]:
    values = np.asarray(differences, dtype=float)
    walsh = np.sort(
        np.asarray(
            [(values[i] + values[j]) / 2 for i in range(len(values)) for j in range(i, len(values))]
        )
    )
    estimate = float(np.median(walsh))
    rank_sums = np.asarray(
        [
            sum((index + 1) for index, sign in enumerate(signs) if sign > 0)
            for signs in itertools.product((-1, 1), repeat=len(values))
        ],
        dtype=int,
    )
    critical = [
        value
        for value in range(int(rank_sums.max()) + 1)
        if np.mean(rank_sums <= value) <= alpha / 2
    ]
    if not critical:
        return estimate, float("-inf"), float("inf")
    k = max(critical)
    lower = float(walsh[k])
    upper = float(walsh[len(walsh) - k - 1])
    return estimate, lower, upper


def summarize_difference(
    scores: pd.DataFrame,
    score_column: str,
    set_label: str,
    analysis: str,
    expected_direction: int | None,
) -> tuple[dict, pd.DataFrame]:
    current = scores.loc[
        scores["compartment"].isin(["Myeloid", "NonMyeloid"]),
        ["patient_id", "compartment", score_column],
    ]
    wide = current.pivot(index="patient_id", columns="compartment", values=score_column)
    if wide.isna().any().any():
        raise ValueError(f"Incomplete primary pairs for {score_column}")
    differences = (wide["Myeloid"] - wide["NonMyeloid"]).to_numpy(float)
    estimate, lower, upper = exact_signed_rank_interval(differences)
    loo_medians = np.asarray(
        [np.median(np.delete(differences, index)) for index in range(len(differences))]
    )
    positive = int((differences > 0).sum())
    wilcoxon = stats.wilcoxon(
        differences, alternative="two-sided", zero_method="wilcox", method="exact"
    )
    sign_test = stats.binomtest(positive, len(differences), 0.5, alternative="two-sided")
    mean = float(differences.mean())
    sd = float(differences.std(ddof=1))
    correction = 1 - 3 / (4 * (len(differences) - 1) - 1)
    summary = {
        "analysis": analysis,
        "set_label": set_label,
        "score_column": score_column,
        "contrast": "Myeloid_minus_pooled_nonmyeloid",
        "n_patients": len(differences),
        "mean_difference": mean,
        "median_difference": float(np.median(differences)),
        "hodges_lehmann_pseudomedian": estimate,
        "exact_signed_rank_ci95_lower": lower,
        "exact_signed_rank_ci95_upper": upper,
        "cohen_dz": mean / sd if sd else np.nan,
        "hedges_gz": correction * mean / sd if sd else np.nan,
        "positive_patients": positive,
        "positive_fraction": positive / len(differences),
        "minimum_leave_one_out_median": float(loo_medians.min()),
        "all_leave_one_out_medians_positive": bool((loo_medians > 0).all()),
        "exact_wilcoxon_p_two_sided": float(wilcoxon.pvalue),
        "exact_sign_p_two_sided": float(sign_test.pvalue),
        "expected_direction": expected_direction,
    }
    summary["passes_frozen_primary_support"] = bool(
        analysis == "SecAct"
        and set_label == "dynamic_up_set"
        and len(differences) >= 8
        and np.median(differences) > 0
        and positive >= int(np.ceil(0.75 * len(differences)))
        and (loo_medians > 0).all()
        and wilcoxon.pvalue < 0.05
    )
    patient = wide.reset_index().rename(
        columns={
            "Myeloid": "myeloid_score",
            "NonMyeloid": "pooled_nonmyeloid_score",
        }
    )
    patient["difference_myeloid_minus_nonmyeloid"] = differences
    patient.insert(1, "analysis", analysis)
    patient.insert(2, "set_label", set_label)
    return summary, patient


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("expression", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("dynamic_sets", type=Path)
    parser.add_argument("group_members", type=Path)
    parser.add_argument("work_output_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--n-rand", type=int, default=1000)
    parser.add_argument("--rng-method", choices=["srand", "gsl", "numpy"], default="gsl")
    args = parser.parse_args()
    args.work_output_dir.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    expression = pd.read_csv(args.expression, sep="\t", index_col=0)
    expression.index = expression.index.astype(str).str.upper()
    manifest = pd.read_csv(args.manifest, sep="\t")
    manifest = manifest.loc[
        manifest["sample_type"].eq("primary_NPC")
        & manifest["eligible_50_cells"]
        & (
            manifest["aggregation_role"].eq("author_broad_lineage")
            | manifest["compartment"].eq("NonMyeloid")
        )
    ].copy()
    pair_counts = manifest.loc[
        manifest["compartment"].isin(["Myeloid", "NonMyeloid"])
    ].groupby("patient_id")["compartment"].nunique()
    eligible_patients = sorted(pair_counts.index[pair_counts.eq(2)])
    if len(eligible_patients) != 14:
        raise ValueError(f"Expected 14 frozen eligible patients, found {eligible_patients}")
    manifest = manifest.loc[manifest["patient_id"].isin(eligible_patients)].copy()
    columns = manifest["group_id"].tolist()
    expression = expression.loc[:, columns]
    if expression.isna().any().any() or (expression.to_numpy() < 0).any():
        raise ValueError("Invalid GSE150430 processed-expression input")

    dynamic = pd.read_csv(args.dynamic_sets, sep="\t")
    set_members = {
        label: set(subset["signature_group"].astype(str))
        for label, subset in dynamic.groupby("set_label")
    }
    if {label: len(members) for label, members in set_members.items()} != {
        "dynamic_down_set": 93,
        "dynamic_up_set": 70,
    }:
        raise ValueError("Frozen dynamic set membership changed")

    signature_genes = set(load_signature("secact").index.astype(str).str.upper())
    secact_input = expression.loc[expression.index.intersection(signature_genes)]
    print(
        f"SECACT START genes={secact_input.shape[0]} groups={secact_input.shape[1]}",
        flush=True,
    )
    result = secact_activity_inference(
        secact_input,
        is_differential=False,
        sig_matrix="secact",
        n_rand=args.n_rand,
        seed=SEED,
        rng_method=args.rng_method,
        sort_genes=True,
        verbose=True,
    )
    group_activity = group_level_matrix(result["zscore"])
    standardized = standardize_rows(group_activity)
    missing = set().union(*set_members.values()) - set(standardized.index)
    if missing:
        raise ValueError(f"Frozen SecAct groups missing or invariant: {sorted(missing)}")
    scores = pd.DataFrame(index=columns)
    for label, members in set_members.items():
        scores[label] = standardized.loc[sorted(members)].mean(axis=0)

    direct_sets, direct_members = read_direct_definitions(
        args.dynamic_sets, args.group_members
    )
    direct, direct_coverage = direct_scores(expression, direct_sets, direct_members)
    scores = scores.join(direct)
    scores["group_id"] = scores.index
    scores = scores.merge(manifest, on="group_id", validate="one_to_one")

    summaries = []
    patient_tables = []
    for analysis, prefix in (("SecAct", ""), ("Direct_expression", "direct_")):
        for label, expected in (("dynamic_up_set", 1), ("dynamic_down_set", None)):
            summary, patients = summarize_difference(
                scores, f"{prefix}{label}", label, analysis, expected
            )
            summaries.append(summary)
            patient_tables.append(patients)
    summary_table = pd.DataFrame(summaries)
    patient_table = pd.concat(patient_tables, ignore_index=True)

    scores.to_csv(
        args.output_dir / "NPC_SecAct_GSE150430_compartment_scores_v0.1.tsv",
        sep="\t",
        index=False,
    )
    patient_table.to_csv(
        args.output_dir / "NPC_SecAct_GSE150430_primary_patient_differences_v0.1.tsv",
        sep="\t",
        index=False,
    )
    summary_table.to_csv(
        args.output_dir / "NPC_SecAct_GSE150430_localization_summary_v0.1.tsv",
        sep="\t",
        index=False,
    )
    direct_coverage.to_csv(
        args.output_dir / "NPC_SecAct_GSE150430_direct_expression_coverage_v0.1.tsv",
        sep="\t",
        index=False,
    )
    group_activity.to_csv(
        args.work_output_dir / "GSE150430_secact_group_activity.tsv.gz",
        sep="\t",
        compression="gzip",
    )
    result["zscore"].to_csv(
        args.work_output_dir / "GSE150430_secact_zscore.tsv.gz",
        sep="\t",
        compression="gzip",
    )
    run = {
        "dataset": "GSE150430",
        "secactpy_version": version("secactpy"),
        "seed": SEED,
        "n_rand": args.n_rand,
        "rng_method": args.rng_method,
        "n_primary_patients": len(eligible_patients),
        "eligible_patients": eligible_patients,
        "n_scoring_domain_compartments": len(manifest),
        "secact_signature_gene_coverage": int(secact_input.shape[0]),
        "released_expression_scale": "author processed nonnegative values; patient-compartment arithmetic mean",
    }
    (args.work_output_dir / "GSE150430_localization_run.json").write_text(
        json.dumps(run, indent=2), encoding="utf-8"
    )

    primary = summary_table.loc[
        (summary_table["analysis"] == "SecAct")
        & (summary_table["set_label"] == "dynamic_up_set")
    ].iloc[0]
    report = f"""# GSE150430 independent single-cell localization report v0.1

**Analysis date:** 2026-08-25  
**Role:** Prespecified independent baseline localization of the frozen paired-bulk SecAct programs.

## Primary result

Fourteen primary NPC patients met the frozen 50-cell requirement. The frozen treatment-up SecAct score was compared within patient between myeloid cells and all pooled non-myeloid cells using patient-compartment mean released processed expression.

- Median myeloid-minus-non-myeloid difference: {primary['median_difference']:.4f}.
- Hodges-Lehmann pseudomedian (exact signed-rank 95% interval): {primary['hodges_lehmann_pseudomedian']:.4f} ({primary['exact_signed_rank_ci95_lower']:.4f}, {primary['exact_signed_rank_ci95_upper']:.4f}).
- Positive patients: {int(primary['positive_patients'])}/14.
- Two-sided exact Wilcoxon P: {primary['exact_wilcoxon_p_two_sided']:.5g}.
- All leave-one-patient-out medians positive: {'yes' if primary['all_leave_one_out_medians_positive'] else 'no'}.
- Frozen independent-localization support rule: {'passed' if primary['passes_frozen_primary_support'] else 'not passed'}.

## Interpretation boundary

This independent cohort localizes a transcriptome-inferred program on author-processed single-cell expression. It does not test a treatment change, prove secretion or cellular source, or validate a clinical predictor. The treatment-down and direct-expression comparator results are retained in the accompanying tables regardless of direction.
"""
    (args.output_dir / "NPC_SecAct_GSE150430_localization_report_v0.1.md").write_text(
        report, encoding="utf-8"
    )
    print(summary_table.to_string(index=False), flush=True)
    print(json.dumps(run, indent=2), flush=True)


if __name__ == "__main__":
    main()
