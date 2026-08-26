from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


DIRECT = "#6A51A3"
SECACT = "#238B45"
CYTOSIG = "#326B7D"


def natural_key(value: str) -> tuple[str, int]:
    match = re.match(r"([A-Za-z]+)(\d+)$", value)
    return (match.group(1), int(match.group(2))) if match else (value, 0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    out_dir = args.output_dir

    direct = pd.read_csv(out_dir / "NPC_SecAct_direct_expression_visium_patient_correlations_v0.1.tsv", sep="\t")
    secact = pd.read_csv(out_dir / "NPC_SecAct_GSE206245_patient_correlations_v0.2.tsv", sep="\t")
    comparison = pd.read_csv(out_dir / "NPC_SecAct_direct_expression_incremental_comparison_v0.1.tsv", sep="\t")
    bulk = pd.read_csv(out_dir / "NPC_SecAct_CytoSig_bulk_results_v0.1.tsv", sep="\t")
    scrna = pd.read_csv(out_dir / "NPC_SecAct_CytoSig_scrna_results_v0.1.tsv", sep="\t")
    cyt_patient = pd.read_csv(out_dir / "NPC_SecAct_CytoSig_visium_patient_correlations_v0.1.tsv", sep="\t")
    cyt_summary = pd.read_csv(out_dir / "NPC_SecAct_CytoSig_visium_summary_v0.1.tsv", sep="\t")

    direct = direct.query("scenario == 'source_filtered' and set_label == 'dynamic_up_set'")
    secact = secact.query("scenario == 'source_filtered' and set_label == 'dynamic_up_set'")
    paired = direct[["sample", "direct_expression_rho"]].merge(
        secact[["sample", "spearman_rho"]], on="sample", validate="one_to_one"
    )
    paired["order"] = paired["sample"].map(natural_key)
    paired = paired.sort_values("order")
    comp = comparison.query("scenario == 'source_filtered'").iloc[0]

    bulk = bulk.query("scenario == 'main' and cell_type == 'all'").set_index("cytokine")
    myeloid = scrna.query("scenario == 'main' and cell_type == 'Myeloid'").set_index("cytokine")
    common = sorted(set(bulk.index) & set(myeloid.index))
    bulk_effect = bulk.loc[common, "mean_difference_OT_minus_BT"]
    myeloid_effect = myeloid.loc[common, "mean_difference_OT_minus_BT"]
    concordance = stats.spearmanr(bulk_effect, myeloid_effect)

    source_summary = cyt_summary.query("scenario == 'source_filtered'").sort_values(
        "median_patient_rho", ascending=False
    )
    sample_order = sorted(cyt_patient["sample"].unique(), key=natural_key)
    matrix = (
        cyt_patient.query("scenario == 'source_filtered'")
        .pivot(index="cytokine", columns="sample", values="spearman_rho")
        .loc[source_summary["cytokine"], sample_order]
    )
    supported = source_summary.set_index("cytokine")["passes_frozen_spatial_support"].to_dict()
    ylabels = [f"{name}*" if supported[name] else name for name in matrix.index]

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.labelsize": 9.5,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8.5,
        "axes.linewidth": 0.8,
    })
    fig = plt.figure(figsize=(13.2, 13.2), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=[0.82, 1.35], width_ratios=[1, 1.45])
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])

    x = np.arange(len(paired))
    for index, row in enumerate(paired.itertuples(index=False)):
        ax_a.plot([index, index], [row.direct_expression_rho, row.spearman_rho], color="#BDBDBD", lw=1)
    ax_a.scatter(x, paired["direct_expression_rho"], color=DIRECT, s=42, label="Direct ligand expression", zorder=3)
    ax_a.scatter(x, paired["spearman_rho"], color=SECACT, s=42, label="SecAct", zorder=3)
    ax_a.axhline(0, color="#333333", lw=0.8)
    ax_a.set_xticks(x, paired["sample"], rotation=45, ha="right")
    ax_a.set_ylabel("Within-patient Spearman rho")
    ax_a.set_title("A  Spatial association by patient", loc="left", fontweight="bold")
    ax_a.legend(frameon=False, loc="upper right")

    ax_b.scatter(paired["direct_expression_rho"], paired["spearman_rho"], color="#3F3F3F", s=44)
    low = min(paired["direct_expression_rho"].min(), paired["spearman_rho"].min()) - 0.05
    high = max(paired["direct_expression_rho"].max(), paired["spearman_rho"].max()) + 0.05
    ax_b.plot([low, high], [low, high], ls="--", color="#777777", lw=1)
    label_offsets = {
        "ST7": (4, -11),
        "ST10": (-2, -15),
        "ST11": (4, -12),
        "ST12": (4, 5),
        "ST17": (-7, 10),
        "ST18": (-18, -11),
        "ST5": (4, -12),
    }
    for row in paired.itertuples(index=False):
        offset = label_offsets.get(row.sample, (4, 4))
        ax_b.annotate(
            row.sample,
            (row.direct_expression_rho, row.spearman_rho),
            xytext=offset,
            textcoords="offset points",
            fontsize=7.5,
        )
    ax_b.set_xlim(low, high)
    ax_b.set_ylim(low, high)
    ax_b.set_xlabel("Direct ligand-expression rho")
    ax_b.set_ylabel("SecAct rho")
    ax_b.set_title("B  Frozen paired spatial benchmark", loc="left", fontweight="bold")
    ax_b.text(
        0.03,
        0.97,
        f"median difference={comp.median_paired_difference_secact_minus_direct:.3f}\n"
        f"exact paired P={comp.exact_paired_wilcoxon_p_two_sided:.5f}",
        transform=ax_b.transAxes,
        ha="left",
        va="top",
    )

    ax_c.scatter(bulk_effect, myeloid_effect, s=32, color=CYTOSIG, alpha=0.9)
    ax_c.axhline(0, color="#999999", lw=0.8)
    ax_c.axvline(0, color="#999999", lw=0.8)
    ax_c.set_xlabel("Bulk mean OT-BT")
    ax_c.set_ylabel("Myeloid mean OT-BT")
    ax_c.set_title("C  CytoSig paired-effect concordance", loc="left", fontweight="bold")
    ax_c.text(
        0.03,
        0.97,
        f"Spearman rho={concordance.statistic:.3f}\nP={concordance.pvalue:.3f}",
        transform=ax_c.transAxes,
        ha="left",
        va="top",
    )

    image = ax_d.imshow(matrix.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-0.8, vmax=0.8)
    ax_d.set_xticks(range(matrix.shape[1]), matrix.columns, rotation=45, ha="right")
    ax_d.set_yticks(range(matrix.shape[0]), ylabels, fontsize=7)
    ax_d.set_title("D  CytoSig-myeloid spatial correlations", loc="left", fontweight="bold")
    ax_d.text(1.0, -0.09, "* met frozen positive spatial-support rule", transform=ax_d.transAxes, ha="right", va="top", fontsize=7.5)
    colorbar = fig.colorbar(image, ax=ax_d, label="Within-patient Spearman rho", shrink=0.72, pad=0.02)
    colorbar.ax.tick_params(labelsize=8)

    for ax in (ax_a, ax_b, ax_c):
        ax.spines[["top", "right"]].set_visible(False)

    figure_dir = out_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    prefix = figure_dir / "NPC_SecAct_Translational_Oncology_figure5_combined_benchmarks_v0.1"
    fig.savefig(f"{prefix}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(f"{prefix}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"WROTE {prefix}.png")
    print(f"WROTE {prefix}.pdf")


if __name__ == "__main__":
    main()
