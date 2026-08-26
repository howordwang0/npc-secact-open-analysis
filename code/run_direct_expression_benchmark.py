from __future__ import annotations

import argparse
import gzip
import itertools
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.io import mmread


SEED = 20260823
UP = "#00796B"
DOWN = "#D55E00"
DIRECT = "#6A51A3"
SECACT = "#238B45"


def exact_sign_flip_p(differences: np.ndarray) -> float:
    signs = np.asarray(list(itertools.product((-1.0, 1.0), repeat=len(differences))))
    permuted = signs @ differences / len(differences)
    observed = abs(float(differences.mean()))
    return float((np.abs(permuted) >= observed - 1e-12).mean())


def paired_summary(differences: np.ndarray, expected_direction: int) -> dict:
    differences = np.asarray(differences, dtype=float)
    n = len(differences)
    mean = float(differences.mean())
    sd = float(differences.std(ddof=1))
    cohen_dz = mean / sd if sd else np.nan
    correction = 1 - 3 / (4 * (n - 1) - 1)
    loo = np.asarray([np.delete(differences, index).mean() for index in range(n)])
    return {
        "n_patients": n,
        "mean_difference": mean,
        "median_difference": float(np.median(differences)),
        "cohen_dz": cohen_dz,
        "hedges_gz": correction * cohen_dz,
        "expected_direction": expected_direction,
        "patients_expected_direction": int((np.sign(differences) == expected_direction).sum()),
        "leave_one_out_expected_direction": bool(np.all(np.sign(loo) == expected_direction)),
        "exact_sign_flip_p_two_sided": exact_sign_flip_p(differences),
    }


def one_sample_spatial_summary(values: np.ndarray) -> dict:
    values = np.asarray(values, dtype=float)
    loo_medians = np.asarray([np.median(np.delete(values, index)) for index in range(len(values))])
    positive = int((values > 0).sum())
    wilcoxon = stats.wilcoxon(values, alternative="two-sided", method="exact")
    sign = stats.binomtest(positive, len(values), 0.5, alternative="two-sided")
    return {
        "n_patients": len(values),
        "median_rho": float(np.median(values)),
        "mean_rho": float(values.mean()),
        "minimum_rho": float(values.min()),
        "maximum_rho": float(values.max()),
        "positive_patients": positive,
        "minimum_leave_one_out_median_rho": float(loo_medians.min()),
        "exact_wilcoxon_p_two_sided": float(wilcoxon.pvalue),
        "exact_sign_p_two_sided": float(sign.pvalue),
    }


def read_definitions(dynamic_sets_path: Path, group_members_path: Path) -> tuple[dict[str, list[str]], dict[str, list[str]]]:
    dynamic = pd.read_csv(dynamic_sets_path, sep="\t")
    members = pd.read_csv(group_members_path, sep="\t")
    sets = {
        set_label: subset["signature_group"].tolist()
        for set_label, subset in dynamic.groupby("set_label", sort=False)
    }
    group_members = {
        group: sorted(set(subset["protein"].astype(str).str.upper()))
        for group, subset in members.groupby("signature_group", sort=False)
    }
    return sets, group_members


def collapse_expression(expression: pd.DataFrame, reducer: str) -> pd.DataFrame:
    expression = expression.copy()
    expression.index = expression.index.astype(str).str.upper()
    if reducer == "median":
        return expression.groupby(level=0, sort=False).median()
    if reducer == "sum":
        return expression.groupby(level=0, sort=False).sum()
    raise ValueError(reducer)


def direct_set_scores(
    expression: pd.DataFrame,
    sets: dict[str, list[str]],
    group_members: dict[str, list[str]],
    excluded_genes: set[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    excluded = {gene.upper() for gene in (excluded_genes or set())}
    expression = collapse_expression(expression, "median")
    sd = expression.std(axis=1, ddof=1)
    valid = sd[(sd > 0) & np.isfinite(sd)].index
    standardized = expression.loc[valid].sub(expression.loc[valid].mean(axis=1), axis=0).div(sd.loc[valid], axis=0)

    scores: dict[str, pd.Series] = {}
    coverage_rows = []
    for set_label, groups in sets.items():
        group_scores = []
        for group in groups:
            eligible = [
                gene for gene in group_members[group]
                if gene not in excluded and gene in standardized.index
            ]
            if eligible:
                group_scores.append(standardized.loc[eligible].mean(axis=0).rename(group))
        if not group_scores:
            raise RuntimeError(f"No represented groups for {set_label}")
        group_frame = pd.DataFrame(group_scores)
        scores[set_label] = group_frame.mean(axis=0)
        coverage_rows.append(
            {
                "set_label": set_label,
                "total_groups": len(groups),
                "represented_groups": len(group_scores),
                "represented_fraction": len(group_scores) / len(groups),
            }
        )
    return pd.DataFrame(scores), pd.DataFrame(coverage_rows)


def load_bulk_expression(path: Path) -> pd.DataFrame:
    bulk = pd.read_csv(path, sep="\t", low_memory=False)
    sample_columns = [column for column in bulk if re.fullmatch(r"Patient_\d+_(?:BT|OT)\d+", column)]
    bulk = bulk.loc[bulk["gene_symbol"].notna() & (bulk["gene_symbol"].astype(str).str.len() > 0)]
    expression = bulk.set_index("gene_symbol")[sample_columns].apply(pd.to_numeric, errors="coerce")
    if expression.isna().any().any() or (expression < 0).any().any():
        raise ValueError("Bulk expression contains missing or negative values")
    expression = collapse_expression(expression, "median")
    return np.log2(expression + 1)


def run_bulk(
    expression_path: Path,
    sets: dict[str, list[str]],
    group_members: dict[str, list[str]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    expression = load_bulk_expression(expression_path)
    scores, coverage = direct_set_scores(expression, sets, group_members)
    patients = [12, 13, 19, 26, 28, 29, 33, 35, 39]
    rows = []
    summaries = []
    for set_label, expected in (("dynamic_up_set", 1), ("dynamic_down_set", -1)):
        differences = []
        for patient in patients:
            bt = float(scores.loc[f"Patient_{patient}_BT{patient}", set_label])
            ot = float(scores.loc[f"Patient_{patient}_OT{patient}", set_label])
            differences.append(ot - bt)
            rows.extend(
                [
                    {"patient": f"P{patient}", "timepoint": "BT", "set_label": set_label, "score": bt},
                    {"patient": f"P{patient}", "timepoint": "OT", "set_label": set_label, "score": ot},
                ]
            )
        summary = paired_summary(np.asarray(differences), expected)
        summary.update({"modality": "paired_bulk", "scenario": "main", "cell_type": "all", "set_label": set_label})
        summaries.append(summary)
    coverage.insert(0, "modality", "paired_bulk")
    coverage.insert(1, "scenario", "main")
    return pd.DataFrame(rows), pd.DataFrame(summaries), coverage


def run_scrna(
    pseudobulk_dir: Path,
    sets: dict[str, list[str]],
    group_members: dict[str, list[str]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    analyses = [
        ("hard_qc", "B", "main"),
        ("hard_qc", "Myeloid", "main"),
        ("hard_qc", "T_NK", "main"),
        ("doublet_excluded", "B", "doublet_sensitivity"),
        ("doublet_excluded", "Myeloid", "doublet_sensitivity"),
        ("doublet_excluded", "T_NK", "doublet_sensitivity"),
        ("hard_qc", "PlasmaB", "cell_count_sensitivity"),
        ("doublet_excluded", "PlasmaB", "combined_sensitivity"),
    ]
    score_rows = []
    summaries = []
    coverage_rows = []
    for matrix_scenario, cell_type, analysis_role in analyses:
        expression = pd.read_csv(
            pseudobulk_dir / f"{matrix_scenario}_log2cpm.tsv.gz",
            sep="\t",
            index_col=0,
        )
        columns = [
            f"Patient_{patient}_{timepoint}{patient}|{cell_type}"
            for patient in range(2, 7)
            for timepoint in ("BT", "OT")
        ]
        expression = expression.loc[:, columns]
        scores, coverage = direct_set_scores(expression, sets, group_members)
        for row in coverage.to_dict("records"):
            row.update({"modality": "scrna", "scenario": analysis_role, "cell_type": cell_type})
            coverage_rows.append(row)
        for set_label, expected in (("dynamic_up_set", 1), ("dynamic_down_set", -1)):
            differences = []
            for patient in range(2, 7):
                bt_column = f"Patient_{patient}_BT{patient}|{cell_type}"
                ot_column = f"Patient_{patient}_OT{patient}|{cell_type}"
                bt = float(scores.loc[bt_column, set_label])
                ot = float(scores.loc[ot_column, set_label])
                differences.append(ot - bt)
                score_rows.extend(
                    [
                        {"scenario": analysis_role, "cell_type": cell_type, "patient": f"P{patient}", "timepoint": "BT", "set_label": set_label, "score": bt},
                        {"scenario": analysis_role, "cell_type": cell_type, "patient": f"P{patient}", "timepoint": "OT", "set_label": set_label, "score": ot},
                    ]
                )
            summary = paired_summary(np.asarray(differences), expected)
            summary.update({"modality": "scrna", "scenario": analysis_role, "cell_type": cell_type, "set_label": set_label})
            summaries.append(summary)
    return pd.DataFrame(score_rows), pd.DataFrame(summaries), pd.DataFrame(coverage_rows)


def read_visium_expression(
    data_dir: Path,
    gsm: str,
    sample: str,
    selected_genes: set[str],
) -> pd.DataFrame:
    prefix = f"{gsm}_NPC_{sample}"
    with gzip.open(data_dir / f"{prefix}_features.tsv.gz", "rt", encoding="utf-8") as handle:
        features = pd.read_csv(handle, sep="\t", header=None)
    with gzip.open(data_dir / f"{prefix}_barcodes.tsv.gz", "rt", encoding="utf-8") as handle:
        barcodes = pd.read_csv(handle, sep="\t", header=None)[0].astype(str).tolist()
    with gzip.open(data_dir / f"{prefix}_matrix.mtx.gz", "rb") as handle:
        matrix = mmread(handle).tocsr()
    symbols = features.iloc[:, 1].astype(str).str.upper()
    total = np.asarray(matrix.sum(axis=0)).ravel().astype(float)
    if (total <= 0).any():
        raise ValueError(f"Nonpositive library size in {sample}")
    rows = {}
    for gene in sorted(selected_genes):
        indexes = np.flatnonzero(symbols.to_numpy() == gene)
        if len(indexes):
            counts = np.asarray(matrix[indexes, :].sum(axis=0)).ravel().astype(float)
            rows[gene] = np.log2(counts / total * 1_000_000 + 1)
    return pd.DataFrame(rows, index=barcodes).T


def run_visium(
    data_dir: Path,
    run_json_path: Path,
    spot_scores_path: Path,
    patient_secact_path: Path,
    marker_path: Path,
    sets: dict[str, list[str]],
    group_members: dict[str, list[str]],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    run = json.loads(run_json_path.read_text(encoding="utf-8"))
    spots = pd.read_csv(spot_scores_path, sep="\t")
    markers = set(pd.read_csv(marker_path, sep="\t")["gene"].astype(str).str.upper())
    selected_genes = {
        gene
        for groups in sets.values()
        for group in groups
        for gene in group_members[group]
        if gene not in markers
    }
    score_rows = []
    patient_rows = []
    coverage_rows = []
    for sample_info in run["samples"]:
        sample = sample_info["sample"]
        gsm = sample_info["gsm"]
        current = spots.loc[spots["sample"] == sample].set_index("barcode")
        expression = read_visium_expression(data_dir, gsm, sample, selected_genes)
        missing_barcodes = current.index.difference(expression.columns)
        if len(missing_barcodes):
            raise ValueError(f"Missing {len(missing_barcodes)} score barcodes in {sample}")
        for scenario, mask in (
            ("source_filtered", np.ones(len(current), dtype=bool)),
            ("fixed_spot_qc", current["passes_fixed_spot_qc"].to_numpy(bool)),
        ):
            barcodes = current.index[mask]
            scores, coverage = direct_set_scores(
                expression.loc[:, barcodes],
                sets,
                group_members,
                excluded_genes=markers,
            )
            for row in coverage.to_dict("records"):
                row.update({"modality": "visium", "scenario": scenario, "cell_type": sample})
                coverage_rows.append(row)
            aligned = current.loc[barcodes]
            for set_label in ("dynamic_up_set", "dynamic_down_set"):
                rho = float(stats.spearmanr(scores[set_label], aligned["myeloid_context"]).statistic)
                patient_rows.append(
                    {
                        "sample": sample,
                        "response": sample_info["response"],
                        "scenario": scenario,
                        "set_label": set_label,
                        "direct_expression_rho": rho,
                        "n_spots": len(barcodes),
                    }
                )
            sample_scores = scores.rename(
                columns={
                    "dynamic_up_set": "direct_dynamic_up_set",
                    "dynamic_down_set": "direct_dynamic_down_set",
                }
            )
            sample_scores.insert(0, "barcode", sample_scores.index)
            sample_scores.insert(0, "scenario", scenario)
            sample_scores.insert(0, "sample", sample)
            score_rows.append(sample_scores.reset_index(drop=True))

    patient = pd.DataFrame(patient_rows)
    summaries = []
    for scenario in ("source_filtered", "fixed_spot_qc"):
        for set_label in ("dynamic_up_set", "dynamic_down_set"):
            values = patient.loc[
                (patient["scenario"] == scenario) & (patient["set_label"] == set_label),
                "direct_expression_rho",
            ].to_numpy(float)
            summary = one_sample_spatial_summary(values)
            summary.update({"modality": "visium", "scenario": scenario, "cell_type": "all", "set_label": set_label})
            summaries.append(summary)

    secact = pd.read_csv(patient_secact_path, sep="\t")
    comparisons = []
    for scenario in ("source_filtered", "fixed_spot_qc"):
        secact_current = secact.loc[
            (secact["scenario"] == scenario) & (secact["set_label"] == "dynamic_up_set"),
            ["sample", "spearman_rho"],
        ]
        direct_current = patient.loc[
            (patient["scenario"] == scenario) & (patient["set_label"] == "dynamic_up_set"),
            ["sample", "direct_expression_rho"],
        ]
        merged = secact_current.merge(direct_current, on="sample", validate="one_to_one")
        differences = merged["spearman_rho"].to_numpy(float) - merged["direct_expression_rho"].to_numpy(float)
        test = stats.wilcoxon(differences, alternative="two-sided", method="exact")
        comparisons.append(
            {
                "scenario": scenario,
                "n_patients": len(merged),
                "median_secact_rho": float(merged["spearman_rho"].median()),
                "median_direct_expression_rho": float(merged["direct_expression_rho"].median()),
                "median_paired_difference_secact_minus_direct": float(np.median(differences)),
                "mean_paired_difference_secact_minus_direct": float(np.mean(differences)),
                "secact_higher_patients": int((differences > 0).sum()),
                "exact_paired_wilcoxon_p_two_sided": float(test.pvalue),
                "passes_incremental_rule": bool(np.median(differences) > 0 and test.pvalue < 0.05),
            }
        )
    return (
        pd.concat(score_rows, ignore_index=True),
        patient,
        pd.DataFrame(summaries),
        pd.DataFrame(coverage_rows),
        {"comparisons": comparisons},
    )


def plot_spatial_comparison(
    direct_patient: pd.DataFrame,
    secact_patient_path: Path,
    comparison: pd.DataFrame,
    output_prefix: Path,
) -> None:
    secact = pd.read_csv(secact_patient_path, sep="\t")
    secact = secact.loc[
        (secact["scenario"] == "source_filtered") & (secact["set_label"] == "dynamic_up_set"),
        ["sample", "spearman_rho", "response"],
    ]
    direct = direct_patient.loc[
        (direct_patient["scenario"] == "source_filtered") & (direct_patient["set_label"] == "dynamic_up_set"),
        ["sample", "direct_expression_rho"],
    ]
    merged = secact.merge(direct, on="sample", validate="one_to_one")
    order = ["ST5", "ST6", "ST7", "ST8", "ST9", "ST10", "ST11", "ST12", "ST16", "ST17", "ST18", "ST19"]
    merged["order"] = merged["sample"].map({sample: index for index, sample in enumerate(order)})
    merged = merged.sort_values("order")
    comp = comparison.loc[comparison["scenario"] == "source_filtered"].iloc[0]

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    ax = axes[0]
    for index, row in enumerate(merged.itertuples(index=False)):
        ax.plot([index, index], [row.direct_expression_rho, row.spearman_rho], color="#BDBDBD", linewidth=1)
    ax.scatter(np.arange(len(merged)), merged["direct_expression_rho"], color=DIRECT, s=42, label="Direct ligand expression", zorder=3)
    ax.scatter(np.arange(len(merged)), merged["spearman_rho"], color=SECACT, s=42, label="SecAct", zorder=3)
    ax.axhline(0, color="#333333", linewidth=0.8)
    ax.set_xticks(np.arange(len(merged)), merged["sample"], rotation=45, ha="right")
    ax.set_ylabel("Within-patient Spearman rho")
    ax.set_title("A  Spatial association by patient", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    ax.scatter(merged["direct_expression_rho"], merged["spearman_rho"], color="#3F3F3F", s=46)
    low = min(merged["direct_expression_rho"].min(), merged["spearman_rho"].min()) - 0.05
    high = max(merged["direct_expression_rho"].max(), merged["spearman_rho"].max()) + 0.05
    ax.plot([low, high], [low, high], linestyle="--", color="#777777", linewidth=1)
    for row in merged.itertuples(index=False):
        ax.annotate(row.sample, (row.direct_expression_rho, row.spearman_rho), xytext=(3, 3), textcoords="offset points", fontsize=7)
    ax.set_xlim(low, high)
    ax.set_ylim(low, high)
    ax.set_xlabel("Direct ligand-expression rho")
    ax.set_ylabel("SecAct rho")
    ax.set_title("B  Frozen paired benchmark", loc="left", fontweight="bold")
    ax.text(
        0.04,
        0.96,
        f"median difference={comp['median_paired_difference_secact_minus_direct']:.3f}\n"
        f"exact paired P={comp['exact_paired_wilcoxon_p_two_sided']:.4f}",
        transform=ax.transAxes,
        va="top",
        fontsize=8.5,
    )
    ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("SecAct versus direct ligand expression in GSE206245", fontsize=14, fontweight="bold")
    fig.tight_layout()
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{output_prefix}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{output_prefix}.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(
    output_path: Path,
    bulk_summary: pd.DataFrame,
    scrna_summary: pd.DataFrame,
    spatial_summary: pd.DataFrame,
    comparison: pd.DataFrame,
    coverage: pd.DataFrame,
) -> None:
    bulk_up = bulk_summary.loc[bulk_summary["set_label"] == "dynamic_up_set"].iloc[0]
    bulk_down = bulk_summary.loc[bulk_summary["set_label"] == "dynamic_down_set"].iloc[0]
    scrna_myeloid_up = scrna_summary.loc[
        (scrna_summary["scenario"] == "main")
        & (scrna_summary["cell_type"] == "Myeloid")
        & (scrna_summary["set_label"] == "dynamic_up_set")
    ].iloc[0]
    visium_up = spatial_summary.loc[
        (spatial_summary["scenario"] == "source_filtered")
        & (spatial_summary["set_label"] == "dynamic_up_set")
    ].iloc[0]
    visium_up_qc = spatial_summary.loc[
        (spatial_summary["scenario"] == "fixed_spot_qc")
        & (spatial_summary["set_label"] == "dynamic_up_set")
    ].iloc[0]
    comp = comparison.loc[comparison["scenario"] == "source_filtered"].iloc[0]
    comp_qc = comparison.loc[comparison["scenario"] == "fixed_spot_qc"].iloc[0]
    bulk_coverage = coverage.loc[(coverage["modality"] == "paired_bulk")]
    visium_coverage = coverage.loc[(coverage["modality"] == "visium") & (coverage["scenario"] == "source_filtered")]

    incremental_text = (
        "met the frozen incremental rule"
        if bool(comp["passes_incremental_rule"])
        else "did not meet the frozen incremental rule"
    )
    report = f"""# Direct ligand-expression benchmark report v0.1

## Result status

The spatial benchmark **{incremental_text}**. In the 12-patient source-filtered
Visium analysis, the median myeloid-context correlation was
{comp['median_secact_rho']:.3f} for SecAct and
{comp['median_direct_expression_rho']:.3f} for the direct ligand-expression
score. The median paired difference (SecAct minus direct expression) was
{comp['median_paired_difference_secact_minus_direct']:.3f}; SecAct was higher
in {int(comp['secact_higher_patients'])}/12 patients, and the exact paired
Wilcoxon P value was {comp['exact_paired_wilcoxon_p_two_sided']:.4f}.

## Coverage

- Paired bulk direct-expression coverage: treatment-up
  {int(bulk_coverage.loc[bulk_coverage['set_label'] == 'dynamic_up_set', 'represented_groups'].iloc[0])}/70 groups;
  treatment-down {int(bulk_coverage.loc[bulk_coverage['set_label'] == 'dynamic_down_set', 'represented_groups'].iloc[0])}/93 groups.
- Source-filtered Visium direct-expression coverage was assessed separately in
  every section after excluding any shared myeloid-context gene. The minimum
  represented counts were
  {int(visium_coverage.loc[visium_coverage['set_label'] == 'dynamic_up_set', 'represented_groups'].min())}/70 treatment-up groups and
  {int(visium_coverage.loc[visium_coverage['set_label'] == 'dynamic_down_set', 'represented_groups'].min())}/93 treatment-down groups.

## Paired bulk descriptive benchmark

| Frozen set | Mean OT-BT | Hedges gz | Patients in expected direction | Exact sign-flip P |
|---|---:|---:|---:|---:|
| Treatment-up direct expression | {bulk_up['mean_difference']:.3f} | {bulk_up['hedges_gz']:.3f} | {int(bulk_up['patients_expected_direction'])}/9 | {bulk_up['exact_sign_flip_p_two_sided']:.4f} |
| Treatment-down direct expression | {bulk_down['mean_difference']:.3f} | {bulk_down['hedges_gz']:.3f} | {int(bulk_down['patients_expected_direction'])}/9 | {bulk_down['exact_sign_flip_p_two_sided']:.4f} |

Because the frozen sets were selected from SecAct behavior in this same cohort,
these direct-expression contrasts are descriptive and cannot establish
incremental value.

## Single-cell myeloid benchmark

For the main hard-QC myeloid treatment-up comparison, the direct-expression
mean OT-BT difference was {scrna_myeloid_up['mean_difference']:.3f}, Hedges gz
was {scrna_myeloid_up['hedges_gz']:.3f},
{int(scrna_myeloid_up['patients_expected_direction'])}/5 patients changed in
the expected direction, and the exact sign-flip P value was
{scrna_myeloid_up['exact_sign_flip_p_two_sided']:.4f}. This must be interpreted
beside, not as a replacement for, the frozen SecAct result (4/5 patients,
Hedges gz 0.267, exact P 0.50).

## Visium direct-expression association

| Scenario | Median direct-expression rho | Positive patients | Wilcoxon P | Median SecAct-minus-direct | Paired P | Incremental rule |
|---|---:|---:|---:|---:|---:|---|
| Source-filtered | {visium_up['median_rho']:.3f} | {int(visium_up['positive_patients'])}/12 | {visium_up['exact_wilcoxon_p_two_sided']:.4f} | {comp['median_paired_difference_secact_minus_direct']:.3f} | {comp['exact_paired_wilcoxon_p_two_sided']:.4f} | {'yes' if comp['passes_incremental_rule'] else 'no'} |
| Fixed spot QC | {visium_up_qc['median_rho']:.3f} | {int(visium_up_qc['positive_patients'])}/12 | {visium_up_qc['exact_wilcoxon_p_two_sided']:.4f} | {comp_qc['median_paired_difference_secact_minus_direct']:.3f} | {comp_qc['exact_paired_wilcoxon_p_two_sided']:.4f} | {'yes' if comp_qc['passes_incremental_rule'] else 'no'} |

## Interpretation

This benchmark compares one frozen, group-equal-weight direct-expression
construction with SecAct. It does not compare proteins, secretion, receptor
engagement, or causality. Failure of the incremental rule requires the
manuscript to state that SecAct's added spatial association over direct ligand
expression was not demonstrated. Passing it would support only incremental
transcriptomic association in this cohort, not a clinical or mechanistic claim.
"""
    output_path.write_text(report, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("dynamic_sets", type=Path)
    parser.add_argument("group_members", type=Path)
    parser.add_argument("bulk_expression", type=Path)
    parser.add_argument("pseudobulk_dir", type=Path)
    parser.add_argument("visium_data_dir", type=Path)
    parser.add_argument("visium_run_json", type=Path)
    parser.add_argument("visium_spot_scores", type=Path)
    parser.add_argument("visium_patient_secact", type=Path)
    parser.add_argument("myeloid_markers", type=Path)
    parser.add_argument("work_output_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("figure_dir", type=Path)
    args = parser.parse_args()
    args.work_output_dir.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    sets, group_members = read_definitions(args.dynamic_sets, args.group_members)
    bulk_scores, bulk_summary, bulk_coverage = run_bulk(args.bulk_expression, sets, group_members)
    scrna_scores, scrna_summary, scrna_coverage = run_scrna(args.pseudobulk_dir, sets, group_members)
    visium_scores, visium_patient, visium_summary, visium_coverage, comparison_data = run_visium(
        args.visium_data_dir,
        args.visium_run_json,
        args.visium_spot_scores,
        args.visium_patient_secact,
        args.myeloid_markers,
        sets,
        group_members,
    )
    comparison = pd.DataFrame(comparison_data["comparisons"])
    coverage = pd.concat([bulk_coverage, scrna_coverage, visium_coverage], ignore_index=True)
    summary = pd.concat([bulk_summary, scrna_summary, visium_summary], ignore_index=True, sort=False)

    bulk_scores.to_csv(args.work_output_dir / "bulk_direct_expression_scores.tsv", sep="\t", index=False)
    scrna_scores.to_csv(args.work_output_dir / "scrna_direct_expression_scores.tsv", sep="\t", index=False)
    visium_scores.to_csv(args.work_output_dir / "visium_direct_expression_spot_scores.tsv.gz", sep="\t", index=False, compression="gzip")
    visium_patient.to_csv(args.output_dir / "NPC_SecAct_direct_expression_visium_patient_correlations_v0.1.tsv", sep="\t", index=False)
    summary.to_csv(args.output_dir / "NPC_SecAct_direct_expression_benchmark_summary_v0.1.tsv", sep="\t", index=False)
    coverage.to_csv(args.output_dir / "NPC_SecAct_direct_expression_benchmark_coverage_v0.1.tsv", sep="\t", index=False)
    comparison.to_csv(args.output_dir / "NPC_SecAct_direct_expression_incremental_comparison_v0.1.tsv", sep="\t", index=False)

    write_report(
        args.output_dir / "NPC_SecAct_direct_expression_benchmark_report_v0.1.md",
        bulk_summary,
        scrna_summary,
        visium_summary,
        comparison,
        coverage,
    )
    plot_spatial_comparison(
        visium_patient,
        args.visium_patient_secact,
        comparison,
        args.figure_dir / "NPC_SecAct_direct_expression_spatial_benchmark_v0.1",
    )
    print("DIRECT EXPRESSION BENCHMARK COMPLETE", flush=True)


if __name__ == "__main__":
    main()
