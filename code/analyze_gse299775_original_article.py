from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import re
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from secactpy import secact_activity_inference
from sklearn.decomposition import PCA
from sklearn.metrics import roc_auc_score


SEED = 20260823
EXPECTED_UP_GROUPS = 70


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def hedges_g(case: np.ndarray, control: np.ndarray) -> float:
    case = np.asarray(case, dtype=float)
    control = np.asarray(control, dtype=float)
    n1, n0 = len(case), len(control)
    pooled = ((n1 - 1) * case.var(ddof=1) + (n0 - 1) * control.var(ddof=1)) / (n1 + n0 - 2)
    if pooled <= 0 or not np.isfinite(pooled):
        return np.nan
    correction = 1 - 3 / (4 * (n1 + n0) - 9)
    return float(correction * (case.mean() - control.mean()) / np.sqrt(pooled))


def exact_mean_difference_p(values: np.ndarray, n_case: int) -> tuple[float, int]:
    values = np.asarray(values, dtype=float)
    observed = values[:n_case].mean() - values[n_case:].mean()
    total_sum = values.sum()
    extreme = 0
    total = 0
    for indexes in itertools.combinations(range(len(values)), n_case):
        case_sum = values[list(indexes)].sum()
        difference = case_sum / n_case - (total_sum - case_sum) / (len(values) - n_case)
        extreme += abs(difference) >= abs(observed) - 1e-12
        total += 1
    return extreme / total, total


def bootstrap_effects(
    secact_case: np.ndarray,
    secact_control: np.ndarray,
    direct_case: np.ndarray,
    direct_control: np.ndarray,
    iterations: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(SEED + 299775)
    rows = []
    for _ in range(iterations):
        case_index = rng.integers(0, len(secact_case), len(secact_case))
        control_index = rng.integers(0, len(secact_control), len(secact_control))
        secact_g = hedges_g(secact_case[case_index], secact_control[control_index])
        direct_g = hedges_g(direct_case[case_index], direct_control[control_index])
        rows.append((secact_g, direct_g, abs(secact_g) - abs(direct_g)))
    return pd.DataFrame(rows, columns=["secact_hedges_g", "direct_hedges_g", "absolute_effect_difference"])


def summarize_score(name: str, score: pd.Series, labels: pd.Series) -> dict[str, float | int | str]:
    case = score.loc[labels == 1].to_numpy(float)
    control = score.loc[labels == 0].to_numpy(float)
    ordered = np.concatenate([case, control])
    exact_p, permutations = exact_mean_difference_p(ordered, len(case))
    welch = stats.ttest_ind(case, control, equal_var=False)
    mann = stats.mannwhitneyu(case, control, alternative="two-sided")
    return {
        "program": name,
        "positive_group": "MET",
        "negative_group": "nonMET",
        "n_positive": len(case),
        "n_negative": len(control),
        "mean_positive": float(case.mean()),
        "mean_negative": float(control.mean()),
        "mean_difference": float(case.mean() - control.mean()),
        "hedges_g": hedges_g(case, control),
        "welch_t": float(welch.statistic),
        "welch_p_two_sided": float(welch.pvalue),
        "mann_whitney_u": float(mann.statistic),
        "mann_whitney_p_two_sided": float(mann.pvalue),
        "exact_permutation_p_two_sided": float(exact_p),
        "exact_permutations": permutations,
        "apparent_auc_MET_as_positive": float(roc_auc_score(labels, score)),
        "apparent_auc_distance_from_0_5": float(abs(roc_auc_score(labels, score) - 0.5)),
    }


def roi_id(title: str) -> str:
    return str(title).split(",", 1)[0].strip()


def patient_id(title: str) -> str:
    match = re.search(r"\b(NPC|NP)-(\d+)_", str(title))
    if not match:
        raise ValueError(f"Cannot parse patient ID: {title}")
    return f"{match.group(1)}-{match.group(2)}"


def plate_fields(identifier: str) -> tuple[str, int]:
    match = re.search(r"-H-([A-Z])(\d{2})$", identifier)
    if not match:
        raise ValueError(f"Cannot parse plate position: {identifier}")
    return match.group(1), int(match.group(2))


def build_manifest(samples: pd.DataFrame, expression_columns: list[str]) -> pd.DataFrame:
    manifest = samples[["geo_accession", "title", "cell_type"]].copy()
    manifest["sample_id"] = manifest["title"].map(roi_id)
    manifest["patient_id"] = manifest["title"].map(patient_id)
    manifest[["plate_row", "well_column"]] = pd.DataFrame(
        manifest["sample_id"].map(plate_fields).tolist(), index=manifest.index
    )
    manifest["acquisition_order"] = np.arange(1, len(manifest) + 1)
    manifest["is_npc"] = manifest["cell_type"].isin(["MET", "nonMET"])
    manifest["met_binary"] = manifest["cell_type"].map({"MET": 1, "nonMET": 0})
    if manifest["sample_id"].tolist() != expression_columns:
        raise ValueError("Sample metadata order does not match the locked expression matrix")
    npc = manifest.loc[manifest["is_npc"]]
    if len(npc) != 24 or npc["patient_id"].nunique() != 24:
        raise ValueError("Expected 24 independent NPC patients")
    if npc["cell_type"].value_counts().to_dict() != {"MET": 14, "nonMET": 10}:
        raise ValueError("Unexpected MET/nonMET counts")
    return manifest


def batch_audit(expression_log: pd.DataFrame, manifest: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    npc = manifest.loc[manifest["is_npc"]].set_index("sample_id")
    labels = npc["met_binary"].astype(int)
    rows: list[dict] = []
    for variable in ["well_column", "acquisition_order"]:
        case = npc.loc[labels == 1, variable].to_numpy(float)
        control = npc.loc[labels == 0, variable].to_numpy(float)
        welch = stats.ttest_ind(case, control, equal_var=False)
        mann = stats.mannwhitneyu(case, control, alternative="two-sided")
        rows.append({
            "audit_variable": variable,
            "MET_mean": case.mean(),
            "nonMET_mean": control.mean(),
            "hedges_g": hedges_g(case, control),
            "welch_p": welch.pvalue,
            "mann_whitney_p": mann.pvalue,
        })
    variable_genes = expression_log.loc[:, npc.index].var(axis=1).sort_values(ascending=False).head(2000).index
    matrix = expression_log.loc[variable_genes, npc.index].T
    matrix = matrix.sub(matrix.mean(axis=0), axis=1).div(matrix.std(axis=0, ddof=1).replace(0, np.nan), axis=1).fillna(0)
    pca = PCA(n_components=3, random_state=SEED).fit_transform(matrix)
    for index in range(3):
        case = pca[labels.to_numpy() == 1, index]
        control = pca[labels.to_numpy() == 0, index]
        welch = stats.ttest_ind(case, control, equal_var=False)
        mann = stats.mannwhitneyu(case, control, alternative="two-sided")
        rows.append({
            "audit_variable": f"expression_PC{index + 1}",
            "MET_mean": case.mean(),
            "nonMET_mean": control.mean(),
            "hedges_g": hedges_g(case, control),
            "welch_p": welch.pvalue,
            "mann_whitney_p": mann.pvalue,
        })
    row_table = pd.crosstab(npc["plate_row"], npc["cell_type"])
    row_p = stats.chi2_contingency(row_table, correction=False).pvalue
    summary = {
        "plate_row_by_phenotype": row_table.to_dict(),
        "plate_row_chi_square_p_descriptive": float(row_p),
        "note": "GEO reports upstream batch-effect removal and Q3 normalization; no outcome-selected adjustment is applied.",
    }
    return pd.DataFrame(rows), summary


def grouped_activity(activity: pd.DataFrame, membership: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for group, block in membership.groupby("signature_group", sort=False):
        representative = block["representative"].iloc[0]
        if representative in activity.index:
            rows.append(activity.loc[representative].rename(group))
    return pd.DataFrame(rows)


def composite_scores(
    group_values: pd.DataFrame,
    up_groups: list[str],
    npc_samples: list[str],
) -> tuple[pd.Series, pd.DataFrame, pd.DataFrame]:
    represented = [group for group in up_groups if group in group_values.index]
    if len(represented) < 63:
        raise ValueError(f"Only {len(represented)}/70 frozen groups represented")
    npc = group_values.loc[represented, npc_samples]
    means = npc.mean(axis=1)
    standard_deviations = npc.std(axis=1, ddof=1).replace(0, np.nan)
    retained = standard_deviations.dropna().index.tolist()
    standardized = group_values.loc[retained].sub(means.loc[retained], axis=0).div(standard_deviations.loc[retained], axis=0)
    score = standardized.mean(axis=0)
    coverage = pd.DataFrame({
        "signature_group": up_groups,
        "represented": [group in group_values.index for group in up_groups],
        "nonzero_npc_variance": [group in retained for group in up_groups],
    })
    scaling = pd.DataFrame({"npc_mean": means, "npc_sd": standard_deviations})
    return score, coverage, scaling


def direct_group_values(
    expression_log: pd.DataFrame,
    membership: pd.DataFrame,
    npc_samples: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    npc = expression_log.loc[:, npc_samples]
    means = npc.mean(axis=1)
    standard_deviations = npc.std(axis=1, ddof=1).replace(0, np.nan)
    gene_z = expression_log.sub(means, axis=0).div(standard_deviations, axis=0)
    rows = []
    coverage_rows = []
    for group, block in membership.groupby("signature_group", sort=False):
        members = block["protein"].astype(str).tolist()
        observed = sorted(set(members) & set(gene_z.index))
        usable = [gene for gene in observed if np.isfinite(standard_deviations.get(gene, np.nan))]
        coverage_rows.append({
            "signature_group": group,
            "members_total": len(members),
            "members_observed": len(observed),
            "members_nonzero_npc_variance": len(usable),
            "observed_genes": "|".join(observed),
        })
        if usable:
            rows.append(gene_z.loc[usable].mean(axis=0).rename(group))
    return pd.DataFrame(rows), pd.DataFrame(coverage_rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the frozen GSE299775 Original Article extension.")
    parser.add_argument("expression", type=Path)
    parser.add_argument("samples", type=Path)
    parser.add_argument("dynamic_sets", type=Path)
    parser.add_argument("membership", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--n-rand", type=int, default=1000)
    parser.add_argument("--bootstrap", type=int, default=2000)
    args = parser.parse_args()

    expression = pd.read_csv(args.expression, sep="\t", index_col=0)
    expression.index = expression.index.astype(str)
    if expression.index.duplicated().any():
        expression = expression.groupby(level=0, sort=True).median()
    if expression.isna().any().any() or (expression < 0).any().any():
        raise ValueError("The locked matrix must be complete and non-negative")
    expression_log = np.log2(expression.astype(float) + 1)
    samples = pd.read_csv(args.samples)
    manifest = build_manifest(samples, expression.columns.tolist())
    dynamic = pd.read_csv(args.dynamic_sets, sep="\t")
    membership = pd.read_csv(args.membership, sep="\t")
    up_groups = dynamic.loc[dynamic["set_label"].eq("dynamic_up_set"), "signature_group"].astype(str).tolist()
    if len(up_groups) != EXPECTED_UP_GROUPS or len(set(up_groups)) != EXPECTED_UP_GROUPS:
        raise ValueError("The locked treatment-up set is not exactly 70 unique groups")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    audit, audit_summary = batch_audit(expression_log, manifest)
    audit.to_csv(args.output_dir / "GSE299775_batch_order_audit.tsv", sep="\t", index=False)
    (args.output_dir / "GSE299775_batch_order_audit.json").write_text(
        json.dumps(audit_summary, indent=2), encoding="utf-8"
    )

    inferred = secact_activity_inference(
        expression_log,
        is_differential=False,
        sig_matrix="secact",
        n_rand=args.n_rand,
        seed=SEED,
        rng_method="gsl",
        sort_genes=True,
        verbose=True,
    )
    for name in ["beta", "se", "zscore", "pvalue"]:
        inferred[name].to_csv(
            args.output_dir / f"GSE299775_secact_{name}.tsv.gz", sep="\t", compression="gzip"
        )
    all_group_activity = grouped_activity(inferred["zscore"], membership)
    all_group_activity.to_csv(
        args.output_dir / "GSE299775_secact_group_activity.tsv.gz", sep="\t", compression="gzip"
    )

    npc_manifest = manifest.loc[manifest["is_npc"]].set_index("sample_id")
    npc_samples = npc_manifest.index.tolist()
    labels = npc_manifest["met_binary"].astype(int)
    secact_score, secact_coverage, secact_scaling = composite_scores(all_group_activity, up_groups, npc_samples)
    direct_groups, direct_gene_coverage = direct_group_values(expression_log, membership, npc_samples)
    direct_score, direct_coverage, direct_scaling = composite_scores(direct_groups, up_groups, npc_samples)

    scores = manifest.set_index("sample_id")[["geo_accession", "patient_id", "cell_type", "is_npc", "plate_row", "well_column", "acquisition_order"]].copy()
    scores["secact_treatment_up_composite"] = secact_score
    scores["direct_treatment_up_composite"] = direct_score
    scores.to_csv(args.output_dir / "GSE299775_patient_level_program_scores.tsv", sep="\t")
    secact_coverage.assign(modality="SecAct").to_csv(
        args.output_dir / "GSE299775_secact_program_coverage.tsv", sep="\t", index=False
    )
    direct_coverage.merge(direct_gene_coverage, on="signature_group", how="left").assign(modality="direct_expression").to_csv(
        args.output_dir / "GSE299775_direct_program_coverage.tsv", sep="\t", index=False
    )
    secact_scaling.to_csv(args.output_dir / "GSE299775_secact_npc_scaling.tsv", sep="\t")
    direct_scaling.to_csv(args.output_dir / "GSE299775_direct_npc_scaling.tsv", sep="\t")

    secact_npc = scores.loc[npc_samples, "secact_treatment_up_composite"]
    direct_npc = scores.loc[npc_samples, "direct_treatment_up_composite"]
    primary_rows = [
        summarize_score("SecAct_treatment_up_composite", secact_npc, labels),
        summarize_score("direct_expression_treatment_up_composite", direct_npc, labels),
    ]
    bootstrap = bootstrap_effects(
        secact_npc.loc[labels == 1].to_numpy(float),
        secact_npc.loc[labels == 0].to_numpy(float),
        direct_npc.loc[labels == 1].to_numpy(float),
        direct_npc.loc[labels == 0].to_numpy(float),
        args.bootstrap,
    )
    primary = pd.DataFrame(primary_rows).set_index("program")
    primary.loc["SecAct_treatment_up_composite", "hedges_g_boot_ci_low"] = bootstrap["secact_hedges_g"].quantile(0.025)
    primary.loc["SecAct_treatment_up_composite", "hedges_g_boot_ci_high"] = bootstrap["secact_hedges_g"].quantile(0.975)
    primary.loc["direct_expression_treatment_up_composite", "hedges_g_boot_ci_low"] = bootstrap["direct_hedges_g"].quantile(0.025)
    primary.loc["direct_expression_treatment_up_composite", "hedges_g_boot_ci_high"] = bootstrap["direct_hedges_g"].quantile(0.975)
    primary.to_csv(args.output_dir / "GSE299775_primary_and_benchmark_statistics.tsv", sep="\t")
    bootstrap.to_csv(args.output_dir / "GSE299775_paired_bootstrap_effects.tsv.gz", sep="\t", index=False, compression="gzip")

    normal = scores.loc[~scores["is_npc"]].copy()
    normal_donor = normal.groupby("patient_id")[["secact_treatment_up_composite", "direct_treatment_up_composite"]].mean()
    normal_donor.to_csv(args.output_dir / "GSE299775_normal_donor_descriptive_scores.tsv", sep="\t")

    secact_row = primary.loc["SecAct_treatment_up_composite"]
    incremental_ci = bootstrap["absolute_effect_difference"].quantile([0.025, 0.975])
    supported = bool(
        secact_row["exact_permutation_p_two_sided"] < 0.05
        and secact_row["hedges_g_boot_ci_low"] * secact_row["hedges_g_boot_ci_high"] > 0
    )
    incremental = bool(incremental_ci.iloc[0] > 0)
    run = {
        "dataset": "GSE299775",
        "analysis": "frozen 70-group treatment-up composite in PanCK-positive GeoMx WTA",
        "protocol": "outputs/NPC_SecAct_Original_Article_GSE299775_frozen_protocol_v0.1.md",
        "patients_primary": 24,
        "MET": 14,
        "nonMET": 10,
        "normal_donors_descriptive": 3,
        "normal_rois_descriptive": 4,
        "n_rand": args.n_rand,
        "rng_method": "gsl",
        "seed": SEED,
        "bootstrap_iterations": args.bootstrap,
        "secactpy_version": version("secactpy"),
        "input_sha256": {str(path): sha256(path) for path in [args.expression, args.samples, args.dynamic_sets, args.membership]},
        "frozen_groups": len(up_groups),
        "secact_groups_retained": int(secact_coverage["nonzero_npc_variance"].sum()),
        "direct_groups_retained": int(direct_coverage["nonzero_npc_variance"].sum()),
        "primary_association_supported": supported,
        "absolute_effect_difference": float(abs(primary.iloc[0]["hedges_g"]) - abs(primary.iloc[1]["hedges_g"])),
        "absolute_effect_difference_boot_ci": [float(incremental_ci.iloc[0]), float(incremental_ci.iloc[1])],
        "secact_incremental_separation_supported": incremental,
        "interpretation_boundary": "Association and benchmarking only; no prospective prediction, cutoff, mechanism, or clinical utility claim.",
    }
    (args.output_dir / "GSE299775_original_article_run.json").write_text(json.dumps(run, indent=2), encoding="utf-8")
    print(json.dumps(run, indent=2))


if __name__ == "__main__":
    main()
