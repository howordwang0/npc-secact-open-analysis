from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from secactpy import load_signature, secact_activity_inference
from secactpy.inference import group_signatures
from sklearn.metrics import roc_auc_score
from statsmodels.stats.multitest import multipletests


SEED = 20260823


def hedges_g(case: np.ndarray, control: np.ndarray) -> float:
    n1, n0 = len(case), len(control)
    variance = ((n1 - 1) * case.var(ddof=1) + (n0 - 1) * control.var(ddof=1)) / (n1 + n0 - 2)
    if variance <= 0 or not np.isfinite(variance):
        return np.nan
    correction = 1 - 3 / (4 * (n1 + n0) - 9)
    return float(correction * (case.mean() - control.mean()) / np.sqrt(variance))


def stratified_bootstrap_g(
    case_matrix: np.ndarray,
    control_matrix: np.ndarray,
    iterations: int,
    rng: np.random.Generator,
) -> np.ndarray:
    n_features = case_matrix.shape[0]
    result = np.empty((iterations, n_features), dtype=np.float32)
    for index in range(iterations):
        sampled_case = case_matrix[:, rng.integers(0, case_matrix.shape[1], case_matrix.shape[1])]
        sampled_control = control_matrix[:, rng.integers(0, control_matrix.shape[1], control_matrix.shape[1])]
        n1, n0 = sampled_case.shape[1], sampled_control.shape[1]
        pooled = ((n1 - 1) * sampled_case.var(axis=1, ddof=1) + (n0 - 1) * sampled_control.var(axis=1, ddof=1)) / (n1 + n0 - 2)
        correction = 1 - 3 / (4 * (n1 + n0) - 9)
        with np.errstate(divide="ignore", invalid="ignore"):
            result[index] = correction * (sampled_case.mean(axis=1) - sampled_control.mean(axis=1)) / np.sqrt(pooled)
    return result


def compare_groups(
    feature_matrix: pd.DataFrame,
    labels: pd.Series,
    positive_name: str,
    negative_name: str,
    bootstrap_iterations: int,
    seed_offset: int,
) -> pd.DataFrame:
    labels = labels.loc[feature_matrix.columns]
    positive = feature_matrix.loc[:, labels == 1].to_numpy(dtype=float)
    negative = feature_matrix.loc[:, labels == 0].to_numpy(dtype=float)
    rows = []
    for index, feature in enumerate(feature_matrix.index):
        case = positive[index]
        control = negative[index]
        welch = stats.ttest_ind(case, control, equal_var=False, nan_policy="omit")
        mann_whitney = stats.mannwhitneyu(case, control, alternative="two-sided")
        rows.append(
            {
                "feature": feature,
                "n_positive": int(np.isfinite(case).sum()),
                "n_negative": int(np.isfinite(control).sum()),
                "positive_group": positive_name,
                "negative_group": negative_name,
                "mean_positive": float(np.nanmean(case)),
                "mean_negative": float(np.nanmean(control)),
                "mean_difference": float(np.nanmean(case) - np.nanmean(control)),
                "hedges_g": hedges_g(case[np.isfinite(case)], control[np.isfinite(control)]),
                "welch_t": float(welch.statistic),
                "welch_p": float(welch.pvalue),
                "mann_whitney_p": float(mann_whitney.pvalue),
                "apparent_auc": float(roc_auc_score(labels, feature_matrix.loc[feature].to_numpy(dtype=float))),
            }
        )
    result = pd.DataFrame(rows).set_index("feature")
    result["welch_fdr_bh"] = multipletests(result["welch_p"].fillna(1.0), method="fdr_bh")[1]

    boot = stratified_bootstrap_g(
        positive,
        negative,
        bootstrap_iterations,
        np.random.default_rng(SEED + seed_offset),
    )
    result["hedges_g_boot_ci_low"] = np.nanquantile(boot, 0.025, axis=0)
    result["hedges_g_boot_ci_high"] = np.nanquantile(boot, 0.975, axis=0)
    observed_sign = np.sign(result["hedges_g"].to_numpy())
    result["bootstrap_same_direction_fraction"] = np.nanmean(np.sign(boot) == observed_sign, axis=0)
    return result.sort_values(["welch_fdr_bh", "welch_p", "hedges_g"], ascending=[True, True, False])


def group_level_matrix(expanded_activity: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    grouped_signature = group_signatures(load_signature("secact"), cor_threshold=0.9)
    rows = []
    member_rows = []
    for group_name in grouped_signature.columns:
        members = group_name.split("|")
        representative = members[0]
        rows.append(expanded_activity.loc[representative].rename(group_name))
        for member in members:
            member_rows.append({"signature_group": group_name, "protein": member, "representative": representative})
    return pd.DataFrame(rows), pd.DataFrame(member_rows)


def ligand_expression_comparator(
    expression: pd.DataFrame,
    group_members: pd.DataFrame,
    sample_order: list[str],
) -> pd.DataFrame:
    expression = expression.loc[:, sample_order]
    row_mean = expression.mean(axis=1)
    row_sd = expression.std(axis=1, ddof=1).replace(0, np.nan)
    z_expression = expression.subtract(row_mean, axis=0).divide(row_sd, axis=0)
    rows = []
    for group_name, block in group_members.groupby("signature_group", sort=False):
        observed = [protein for protein in block["protein"] if protein in z_expression.index]
        if observed:
            rows.append(z_expression.loc[observed].mean(axis=0).rename(group_name))
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the phenotype-linked GSE132112 SecAct discovery analysis.")
    parser.add_argument("expression_tsv_gz", type=Path)
    parser.add_argument("manifest_csv", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--n-rand", type=int, default=1000)
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--rng-method", choices=["srand", "gsl", "numpy"], default="gsl")
    args = parser.parse_args()

    expression = pd.read_csv(args.expression_tsv_gz, sep="\t", index_col=0)
    expression.index = expression.index.astype(str)
    manifest = pd.read_csv(args.manifest_csv)
    if manifest["expression_sample_id"].duplicated().any():
        raise ValueError("Manifest sample IDs are not unique")
    if set(manifest["expression_sample_id"]) != set(expression.columns):
        missing_manifest = sorted(set(expression.columns) - set(manifest["expression_sample_id"]))
        missing_expression = sorted(set(manifest["expression_sample_id"]) - set(expression.columns))
        raise ValueError(f"Manifest mismatch: missing_manifest={missing_manifest}; missing_expression={missing_expression}")
    sample_order = manifest["expression_sample_id"].tolist()
    expression = expression.loc[:, sample_order]
    if expression.isna().any().any():
        expression = expression.apply(lambda row: row.fillna(row.median()), axis=1)

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
    activity = result["zscore"]
    group_activity, group_members = group_level_matrix(activity)
    primary_labels = manifest.set_index("expression_sample_id")["primary_response_cr_pr_vs_sd"].astype(int)
    sensitivity_labels = manifest.set_index("expression_sample_id")["sensitivity_cr_vs_non_cr"].astype(int)

    primary = compare_groups(group_activity, primary_labels, "CR/PR", "SD", args.bootstrap, 0)
    sensitivity = compare_groups(group_activity, sensitivity_labels, "CR", "PR/SD", args.bootstrap, 1)
    primary["candidate_longlist"] = (
        (primary["welch_fdr_bh"] < 0.05)
        & (primary["hedges_g"].abs() >= 0.5)
        & (primary["bootstrap_same_direction_fraction"] >= 0.90)
    )

    ligand = ligand_expression_comparator(expression, group_members, sample_order)
    ligand_primary = compare_groups(ligand, primary_labels, "CR/PR", "SD", args.bootstrap, 2)
    shared_groups = primary.index.intersection(ligand_primary.index)
    primary.loc[shared_groups, "activity_ligand_expression_spearman"] = [
        stats.spearmanr(group_activity.loc[group_name], ligand.loc[group_name], nan_policy="omit").statistic
        for group_name in shared_groups
    ]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name in ("beta", "se", "zscore", "pvalue"):
        result[name].to_csv(args.output_dir / f"GSE132112_secact_{name}.tsv.gz", sep="\t", compression="gzip")
    group_activity.to_csv(args.output_dir / "GSE132112_secact_group_activity.tsv.gz", sep="\t", compression="gzip")
    group_members.to_csv(args.output_dir / "GSE132112_secact_group_members.tsv", sep="\t", index=False)
    primary.to_csv(args.output_dir / "GSE132112_primary_CRPR_vs_SD.tsv", sep="\t")
    sensitivity.to_csv(args.output_dir / "GSE132112_sensitivity_CR_vs_nonCR.tsv", sep="\t")
    ligand_primary.to_csv(args.output_dir / "GSE132112_ligand_expression_primary.tsv", sep="\t")

    summary = {
        "dataset": "GSE132112",
        "patients": int(len(manifest)),
        "response_counts": manifest["response_original"].value_counts().sort_index().to_dict(),
        "gene_expression_rows": int(expression.shape[0]),
        "secact_output_proteins": int(activity.shape[0]),
        "independent_signature_groups_tested": int(group_activity.shape[0]),
        "primary_endpoint": "CR/PR versus SD",
        "sensitivity_endpoint": "CR versus PR/SD",
        "n_rand": args.n_rand,
        "rng_method": args.rng_method,
        "seed": SEED,
        "bootstrap_iterations": args.bootstrap,
        "multiple_testing_family": "657 grouped SecAct signatures; Benjamini-Hochberg",
        "candidate_longlist_rule": "primary FDR < 0.05, absolute Hedges g >= 0.5, bootstrap same-direction fraction >= 0.90",
        "candidate_groups": int(primary["candidate_longlist"].sum()),
        "apparent_auc_warning": "Descriptive discovery-set statistic only; not an internally or externally validated performance estimate",
        "covariate_adjustment": "None; GEO metadata expose response but no complete baseline clinical covariates",
    }
    (args.output_dir / "GSE132112_discovery_run.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
