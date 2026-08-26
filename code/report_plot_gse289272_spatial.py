from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def save_figure(fig: plt.Figure, prefix: Path) -> None:
    fig.savefig(f"{prefix}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{prefix}.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("patient_scores_tsv", type=Path)
    parser.add_argument("summary_tsv", type=Path)
    parser.add_argument("roi_qc_tsv", type=Path)
    parser.add_argument("run_json", type=Path)
    parser.add_argument("figure_dir", type=Path)
    parser.add_argument("report_output", type=Path)
    args = parser.parse_args()
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    args.report_output.parent.mkdir(parents=True, exist_ok=True)

    patients = pd.read_csv(args.patient_scores_tsv, sep="\t")
    summary = pd.read_csv(args.summary_tsv, sep="\t")
    roi_qc = pd.read_csv(args.roi_qc_tsv, sep="\t")
    run = json.loads(args.run_json.read_text(encoding="utf-8"))
    main_scores = patients.loc[
        (patients["aggregation"] == "geomean")
        & (patients["roi_scenario"] == "source_released")
    ].copy()
    main = summary.loc[
        (summary["aggregation"] == "geomean")
        & (summary["roi_scenario"] == "source_released")
    ].iloc[0]

    sns.set_theme(style="whitegrid", context="notebook")
    fig, ax = plt.subplots(figsize=(7.5, 6.2))
    sns.scatterplot(
        data=main_scores,
        x="independent_myeloid_score",
        y="dynamic_up_set",
        hue="Tissue",
        style="Tissue",
        palette={"NPC": "#B33A3A", "CN": "#267A78"},
        s=85,
        ax=ax,
    )
    npc = main_scores.loc[main_scores["Tissue"] == "NPC"]
    sns.regplot(
        data=npc,
        x="independent_myeloid_score",
        y="dynamic_up_set",
        scatter=False,
        ci=None,
        color="#B33A3A",
        line_kws={"linewidth": 1.4, "linestyle": "--"},
        ax=ax,
    )
    for row in main_scores.itertuples(index=False):
        ax.annotate(
            row.PatientID,
            (row.independent_myeloid_score, row.dynamic_up_set),
            xytext=(4, 3),
            textcoords="offset points",
            fontsize=8,
        )
    ax.set(
        xlabel="Independent myeloid marker score",
        ylabel="Partial bulk-up direct-expression score",
        title=(
            "GSE289272 patient-level TLS context\n"
            f"NPC Spearman rho = {main['npc_patient_spearman_rho']:.3f}, "
            f"exact P = {main['npc_patient_exact_permutation_p_two_sided']:.3f}"
        ),
    )
    fig.tight_layout()
    save_figure(
        fig,
        args.figure_dir / "NPC_SecAct_GSE289272_patient_myeloid_context_v0.1",
    )

    pass_by_tissue = (
        roi_qc.groupby("Tissue")
        .agg(
            rois=("roi_id", "size"),
            passing=("passes_geomxtools_default_qc", "sum"),
            patients=("PatientID", "nunique"),
        )
        .reset_index()
    )
    qc_text = "; ".join(
        f"{row.Tissue}: {int(row.passing)}/{int(row.rois)} ROIs"
        for row in pass_by_tissue.itertuples(index=False)
    )
    report = f"""# GSE289272 GeoMx spatial contextualization report v0.1

## Result status

The analysis **did not meet** the frozen spatial-support rule. Across nine NPC
patients, the partial bulk-up direct-expression score was inversely associated
with the independent myeloid marker score (Spearman rho
{main['npc_patient_spearman_rho']:.3f}; exact two-sided permutation P
{main['npc_patient_exact_permutation_p_two_sided']:.4f}). The smallest
leave-one-patient-out rho was {main['npc_leave_one_out_min_rho']:.3f}, so the
required positive direction was not achieved.

## Data and processing

- 55 TLS ROIs from 12 patients were reconstructed from DCC files: 40 NPC ROIs
  from nine patients and 15 inflammatory-control ROIs from three patients.
- DCC files omitted zero-count RTS identifiers; all observed identifiers were
  valid PKC entries and missing identifiers were restored as zero. The maximum
  omitted count was {run['max_omitted_zero_probe_ids_per_dcc']} probes per ROI.
- Only {run['low_ratio_endogenous_probes']} endogenous probes failed the frozen
  GeoMxTools low-ratio rule.
- GeoMx default-QC sensitivity retained {run['n_rois_passing_geomxtools_default_qc']}/55
  ROIs ({qc_text}). All nine NPC patients remained represented.
- Q3 normalization and geometric-mean probe aggregation were used in the main
  analysis; median probe aggregation was prespecified as sensitivity.

## Coverage limitation

The CTA panel contains 1,779 endogenous targets but covers only 745/7,919
SecAct signature genes (9.408%). Full SecAct was therefore not run. The
fallback directly measured 35/70 frozen bulk-up groups and 25/93 frozen
bulk-down groups. These are partial direct-expression modules, not inferred
secreted-protein activities.

## Sensitivity analyses

| Probe aggregation | ROI rule | NPC patient rho | Exact permutation P | Minimum LOO rho | Spatial support |
|---|---|---:|---:|---:|---|
"""
    for row in summary.itertuples(index=False):
        report += (
            f"| {row.aggregation} | {row.roi_scenario} | "
            f"{row.npc_patient_spearman_rho:.3f} | "
            f"{row.npc_patient_exact_permutation_p_two_sided:.4f} | "
            f"{row.npc_leave_one_out_min_rho:.3f} | "
            f"{'yes' if row.passes_frozen_spatial_support else 'no'} |\n"
        )
    report += """

## Interpretation

This TLS-targeted GeoMx dataset does not corroborate a positive spatial
association between the scRNA-localized myeloid program and the partially
observed bulk-up genes. The negative patient-level association is itself not
evidence of antagonism because the panel is targeted, only half of the frozen
groups are observed, and every ROI was selected as TLS. ROI-level correlations
remain descriptive and were not used as independent replicates.

The result must be retained as a prespecified negative spatial test. It does
not invalidate the scRNA directional finding, but it prevents claiming spatial
confirmation from GSE289272. A broader-platform spatial dataset is required
for the remaining contextualization attempt.
"""
    args.report_output.write_text(report, encoding="utf-8")
    print(f"REPORT COMPLETE: {args.report_output}", flush=True)


if __name__ == "__main__":
    main()
