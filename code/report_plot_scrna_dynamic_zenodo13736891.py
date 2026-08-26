from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


PATIENT_COLORS = {
    "P2": "#1B6CA8",
    "P3": "#2A9D8F",
    "P4": "#E9C46A",
    "P5": "#F4A261",
    "P6": "#C44536",
}
CELL_TYPE_ORDER = ["T_NK", "B", "Myeloid", "PlasmaB"]


def save_figure(fig: plt.Figure, prefix: Path) -> None:
    fig.savefig(f"{prefix}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{prefix}.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("scores_tsv", type=Path)
    parser.add_argument("summary_tsv", type=Path)
    parser.add_argument("cell_counts_tsv", type=Path)
    parser.add_argument("figure_dir", type=Path)
    parser.add_argument("decision_output", type=Path)
    parser.add_argument("report_output", type=Path)
    args = parser.parse_args()
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    args.decision_output.parent.mkdir(parents=True, exist_ok=True)

    scores = pd.read_csv(args.scores_tsv, sep="\t")
    summary = pd.read_csv(args.summary_tsv, sep="\t")
    counts = pd.read_csv(args.cell_counts_tsv, sep="\t")
    expected = {"dynamic_up_set": 1, "dynamic_down_set": -1}

    decisions = []
    main_rows = summary.loc[summary["analysis_role"] == "main"]
    for _, row in main_rows.iterrows():
        sensitivity = summary.loc[
            (summary["scenario"] == "doublet_excluded")
            & (summary["broad_cell_type"] == row["broad_cell_type"])
            & (summary["set_label"] == row["set_label"])
        ].iloc[0]
        direction = expected[row["set_label"]]
        direction_consistent = (
            np.sign(row["mean_OT_minus_BT"])
            == np.sign(sensitivity["mean_OT_minus_BT"])
        )
        passed = bool(
            row["mean_matches_expected_direction"]
            and row["patients_expected_direction"] >= 4
            and row["leave_one_out_matches_expected_direction"]
            and sensitivity["mean_matches_expected_direction"]
            and direction_consistent
        )
        decisions.append(
            {
                "broad_cell_type": row["broad_cell_type"],
                "set_label": row["set_label"],
                "expected_direction": direction,
                "main_mean_OT_minus_BT": row["mean_OT_minus_BT"],
                "main_hedges_gz": row["hedges_gz"],
                "main_patients_expected_direction": int(
                    row["patients_expected_direction"]
                ),
                "main_leave_one_out_expected": bool(
                    row["leave_one_out_matches_expected_direction"]
                ),
                "doublet_sensitivity_mean_OT_minus_BT": sensitivity[
                    "mean_OT_minus_BT"
                ],
                "doublet_sensitivity_expected": bool(
                    sensitivity["mean_matches_expected_direction"]
                ),
                "passes_frozen_directional_support": passed,
            }
        )
    decisions = pd.DataFrame(decisions)
    decisions.to_csv(args.decision_output, sep="\t", index=False)

    sns.set_theme(style="whitegrid", context="notebook")
    plot_scores = scores.loc[
        ((scores["scenario"] == "hard_qc") & (scores["cell_threshold"] == 50))
        | (
            (scores["scenario"] == "hard_qc")
            & (scores["broad_cell_type"] == "PlasmaB")
            & (scores["cell_threshold"] == 20)
        )
    ].copy()
    fig, axes = plt.subplots(2, 4, figsize=(15.5, 7.5), sharex=True)
    for row_idx, set_label in enumerate(["dynamic_up_set", "dynamic_down_set"]):
        for col_idx, cell_type in enumerate(CELL_TYPE_ORDER):
            ax = axes[row_idx, col_idx]
            block = plot_scores.loc[plot_scores["broad_cell_type"] == cell_type]
            for patient, patient_block in block.groupby("patient_id"):
                patient_block = patient_block.set_index("timepoint").loc[["BT", "OT"]]
                ax.plot(
                    [0, 1],
                    patient_block[set_label],
                    marker="o",
                    linewidth=1.5,
                    markersize=4.5,
                    color=PATIENT_COLORS[patient],
                    alpha=0.9,
                    label=patient,
                )
            ax.axhline(0, color="#777777", linewidth=0.8, linestyle=":")
            ax.set_xticks([0, 1], ["BT", "OT"])
            threshold = 20 if cell_type == "PlasmaB" else 50
            if row_idx == 0:
                ax.set_title(f"{cell_type} (>= {threshold} cells)")
            if col_idx == 0:
                label = "Bulk-up set score" if row_idx == 0 else "Bulk-down set score"
                ax.set_ylabel(label)
            else:
                ax.set_ylabel("")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.suptitle("Patient-level SecAct trajectories by broad cell type", y=0.99)
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.955),
        ncol=5,
        frameon=False,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.89))
    save_figure(
        fig,
        args.figure_dir / "NPC_SecAct_Zenodo13736891_scrna_dynamic_trajectories_v0.1",
    )

    paired_counts = counts.loc[
        (counts["scenario"] == "hard_qc") & (counts["patient_id"] != "P1")
    ].pivot_table(
        index=["broad_cell_type", "patient_id"],
        columns="timepoint",
        values="n_cells",
    )
    paired_counts["pair_min"] = paired_counts[["BT", "OT"]].min(axis=1)
    heat = paired_counts["pair_min"].unstack("patient_id")
    heat = heat.reindex(
        [
            "T_NK",
            "B",
            "Myeloid",
            "PlasmaB",
            "Epi_Cancer",
            "Fibroblasts",
            "Endothelial",
            "Mast",
        ]
    )
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    sns.heatmap(
        np.log10(heat + 1),
        annot=heat.astype(int),
        fmt="d",
        cmap="YlGnBu",
        cbar_kws={"label": "log10(minimum BT/OT cell count + 1)"},
        linewidths=0.5,
        linecolor="white",
        ax=ax,
    )
    ax.set(xlabel="Patient", ylabel="Broad cell type")
    ax.set_title("Paired pseudobulk eligibility (numbers show minimum BT/OT cells)")
    fig.tight_layout()
    save_figure(
        fig,
        args.figure_dir / "NPC_SecAct_Zenodo13736891_scrna_pair_cell_counts_v0.1",
    )

    supported = decisions.loc[decisions["passes_frozen_directional_support"]]
    support_text = ", ".join(
        f"{row.broad_cell_type}/{row.set_label}"
        for row in supported.itertuples(index=False)
    ) or "none"
    report = f"""# Zenodo 13736891 scRNA-seq SecAct source-localization report v0.1

## Result status

The frozen directional-support rule was met by: **{support_text}**. This is a
directional source-localization result, not conventional statistical
significance and not evidence that the combined treatment acts through one
secreted protein.

## Data and QC

- 112,545 raw barcodes from 11 open 10x datasets were audited.
- 90,315 cells passed the frozen hard QC rules (80.25%).
- The fixed doublet sensitivity retained 87,336 cells, close to the source
  supplement total of 87,224.
- Seventeen Leiden 0.5 clusters were obtained. Cluster 10 (cycling mixed
  lymphocytes) and cluster 14 (T/B doublets) were excluded from lineage SecAct.
- SecAct coverage was 7,761/7,919 signature genes (98.01%).

## Patient-level eligibility

- T/NK, B, and Myeloid had 5/5 paired patients with at least 50 cells at both
  timepoints and entered the main analysis.
- PlasmaB had only 3/5 pairs at 50 cells but 5/5 at 20 cells, so it was retained
  only as a prespecified cell-count sensitivity analysis.
- Epithelial/cancer, fibroblast, endothelial, and mast groups did not reach
  four paired patients even at the 20-cell threshold and were not tested.

## Frozen-set results

| Cell type | Frozen set | Mean OT-BT | Hedges gz | Patients in expected direction | LOO expected | Doublet sensitivity expected | Frozen support |
|---|---|---:|---:|---:|---|---|---|
"""
    for row in decisions.itertuples(index=False):
        report += (
            f"| {row.broad_cell_type} | {row.set_label} | "
            f"{row.main_mean_OT_minus_BT:.3f} | {row.main_hedges_gz:.3f} | "
            f"{row.main_patients_expected_direction}/5 | "
            f"{'yes' if row.main_leave_one_out_expected else 'no'} | "
            f"{'yes' if row.doublet_sensitivity_expected else 'no'} | "
            f"{'yes' if row.passes_frozen_directional_support else 'no'} |\n"
        )
    report += """

## Interpretation

The bulk-up set showed modest concordant activity in Myeloid pseudobulk: 4/5
patients changed in the expected direction, all leave-one-patient-out means
remained positive, and the doublet-excluded analysis remained positive. The
effect was small (Hedges gz about 0.27) and the two-sided exact sign-flip P was
0.50. Therefore the defensible claim is limited to directional localization.

B and T/NK pseudobulk changed opposite to the bulk-up set, and no broad cell
type met the frozen support rule for the bulk-down set. This pattern suggests
that the bulk treatment signal is not a uniform immune-cell response and may
also reflect composition changes or cell types that lacked adequate paired
coverage. PlasmaB did not rescue either set at the 20-cell threshold.

## Go/No-Go decision

Proceed to spatial contextualization because one prespecified cell-type/set
combination met the directional rule. Do not describe the result as a
statistically significant mechanism, a response predictor, or proof of a
myeloid-derived secreted-protein program. Spatial ROIs and spots must be
aggregated at the patient level before inference.
"""
    args.report_output.write_text(report, encoding="utf-8")
    print(f"REPORT COMPLETE: {args.report_output}", flush=True)


if __name__ == "__main__":
    main()
