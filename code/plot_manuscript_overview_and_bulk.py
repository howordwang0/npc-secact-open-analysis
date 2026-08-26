from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd


UP = "#00796B"
DOWN = "#D55E00"
NEUTRAL = "#B8B8B8"
STOP = "#B2182B"
DIRECTIONAL = "#E69F00"
SUPPORT = "#2E7D32"
INK = "#202020"


def save_figure(fig: plt.Figure, prefix: Path) -> None:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{prefix}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(f"{prefix}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def add_box(
    ax: plt.Axes,
    xy: tuple[float, float],
    width: float,
    height: float,
    title: str,
    body: str,
    edge: str,
    face: str = "#FFFFFF",
    title_size: float = 11,
    body_size: float = 9,
) -> None:
    x, y = xy
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.008,rounding_size=0.012",
        linewidth=1.8,
        edgecolor=edge,
        facecolor=face,
    )
    ax.add_patch(patch)
    ax.text(
        x + 0.025 * width,
        y + height * 0.78,
        title,
        ha="left",
        va="center",
        fontsize=title_size,
        fontweight="bold",
        color=INK,
    )
    ax.text(
        x + 0.025 * width,
        y + height * 0.43,
        body,
        ha="left",
        va="center",
        fontsize=body_size,
        color=INK,
        linespacing=1.35,
    )


def add_arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    color: str = "#666666",
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=14,
            linewidth=1.5,
            color=color,
            shrinkA=4,
            shrinkB=4,
        )
    )


def plot_overview(output_prefix: Path) -> None:
    fig, ax = plt.subplots(figsize=(15, 7.2))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(
        0.02,
        0.95,
        "Patient-level multimodal evidence flow",
        fontsize=18,
        fontweight="bold",
        color=INK,
        ha="left",
        va="top",
    )
    ax.text(
        0.02,
        0.895,
        "Cells, ROIs and spots were never counted as independent patients",
        fontsize=10,
        color="#555555",
        ha="left",
        va="top",
    )

    add_box(
        ax,
        (0.02, 0.58),
        0.205,
        0.23,
        "Clinical response screen",
        "GSE132112 | n=95\n0/657 groups passed\nminimum FDR 0.9985",
        STOP,
        "#FFF7F7",
    )
    add_box(
        ax,
        (0.02, 0.25),
        0.205,
        0.23,
        "Metastasis screen",
        "GSE103611 | n=48\n0/657 groups passed\nminimum FDR 0.6097",
        STOP,
        "#FFF7F7",
    )
    ax.text(
        0.1225,
        0.16,
        "Clinical biomarker route stopped",
        fontsize=9,
        color=STOP,
        fontweight="bold",
        ha="center",
    )

    add_box(
        ax,
        (0.29, 0.42),
        0.205,
        0.29,
        "Paired bulk dynamics",
        "Zenodo 13736891 | n=9\n70 treatment-up groups\n93 treatment-down groups\nBroad remodeling",
        UP,
        "#F4FBFA",
    )
    add_box(
        ax,
        (0.555, 0.42),
        0.19,
        0.29,
        "Paired single-cell",
        "5 paired patients\nMyeloid / treatment-up\n4/5 concordant; gz=0.267\nExact P=0.50",
        DIRECTIONAL,
        "#FFFBF0",
    )
    ax.text(
        0.65,
        0.36,
        "Directional localization",
        fontsize=9,
        color=DIRECTIONAL,
        fontweight="bold",
        ha="center",
    )

    add_box(
        ax,
        (0.805, 0.58),
        0.175,
        0.23,
        "GeoMx TLS",
        "9 NPC patients\n9.4% predictor coverage\nPartial expression rho=-0.667\nNo spatial support",
        STOP,
        "#FFF7F7",
    )
    add_box(
        ax,
        (0.805, 0.25),
        0.175,
        0.23,
        "Whole-transcriptome Visium",
        "12 patients | 23,320 spots\nMedian rho=0.343; 11/12 positive\nWilcoxon P=0.00098\nFormal spatial support",
        SUPPORT,
        "#F5FBF5",
        title_size=9.2,
        body_size=8.2,
    )

    add_arrow(ax, (0.495, 0.565), (0.555, 0.565))
    add_arrow(ax, (0.745, 0.565), (0.805, 0.695))
    add_arrow(ax, (0.745, 0.565), (0.805, 0.365))

    ax.text(0.3925, 0.75, "Treatment-associated programs", ha="center", fontsize=9, color=UP)
    ax.text(0.8925, 0.84, "Targeted panel: negative", ha="center", fontsize=9, color=STOP)
    ax.text(0.8925, 0.19, "Context supported, source not proven", ha="center", fontsize=9, color=SUPPORT)

    ax.text(
        0.29,
        0.12,
        "Defensible synthesis:",
        fontsize=10,
        fontweight="bold",
        color=INK,
        ha="left",
    )
    ax.text(
        0.29,
        0.075,
        "broad treatment-associated secreted-activity remodeling with modest myeloid localization\n"
        "and a reproducible myeloid-rich spatial context",
        fontsize=9.5,
        color=INK,
        ha="left",
        va="center",
        linespacing=1.3,
    )
    ax.text(
        0.29,
        0.018,
        "Not established: a single-protein mechanism, myeloid secretion, treatment causality,\n"
        "or a validated response predictor",
        fontsize=8.5,
        color=STOP,
        ha="left",
        va="center",
        linespacing=1.25,
    )

    save_figure(fig, output_prefix)


def row_zscore(frame: pd.DataFrame) -> pd.DataFrame:
    mean = frame.mean(axis=1)
    sd = frame.std(axis=1, ddof=1).replace(0, np.nan)
    return frame.sub(mean, axis=0).div(sd, axis=0)


def plot_bulk(
    results_path: Path,
    activity_path: Path,
    output_prefix: Path,
) -> None:
    results = pd.read_csv(results_path, sep="\t")
    activity = pd.read_csv(activity_path, sep="\t", index_col=0)
    activity.index.name = "signature_group"
    results = results.set_index("signature_group").loc[activity.index].reset_index()

    effect = results["hedges_gz"].to_numpy(float)
    pvalue = results["exact_sign_flip_p"].to_numpy(float)
    dynamic = results["exploratory_dynamic_group"].astype(bool).to_numpy()
    colors = np.full(len(results), NEUTRAL, dtype=object)
    colors[dynamic & (effect > 0)] = UP
    colors[dynamic & (effect < 0)] = DOWN

    activity_z = row_zscore(activity)
    up_groups = results.loc[dynamic & (effect > 0), "signature_group"].tolist()
    down_groups = results.loc[dynamic & (effect < 0), "signature_group"].tolist()

    up_score = activity_z.loc[up_groups].mean(axis=0)
    down_score = activity_z.loc[down_groups].mean(axis=0)

    patients = [12, 13, 19, 26, 28, 29, 33, 35, 39]
    bt_columns = [f"Patient_{patient}_BT{patient}" for patient in patients]
    ot_columns = [f"Patient_{patient}_OT{patient}" for patient in patients]
    ordered_columns = [column for pair in zip(bt_columns, ot_columns) for column in pair]

    selected_results = results.set_index("signature_group")
    heat_order = (
        selected_results.loc[up_groups].sort_values("hedges_gz", ascending=False).index.tolist()
        + selected_results.loc[down_groups].sort_values("hedges_gz", ascending=True).index.tolist()
    )
    heat = activity_z.loc[heat_order, ordered_columns]

    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )
    fig = plt.figure(figsize=(14, 10))
    grid = fig.add_gridspec(2, 3, height_ratios=[1, 1.45], hspace=0.35, wspace=0.32)
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[0, 2])
    ax_d = fig.add_subplot(grid[1, :])

    ax_a.scatter(effect[~dynamic], -np.log10(pvalue[~dynamic]), s=16, c=NEUTRAL, alpha=0.65, linewidths=0)
    ax_a.scatter(effect[dynamic & (effect > 0)], -np.log10(pvalue[dynamic & (effect > 0)]), s=24, c=UP, alpha=0.9, linewidths=0)
    ax_a.scatter(effect[dynamic & (effect < 0)], -np.log10(pvalue[dynamic & (effect < 0)]), s=24, c=DOWN, alpha=0.9, linewidths=0)
    ax_a.axvline(0, color="#555555", linewidth=0.8)
    ax_a.set_xlabel("Hedges gz (OT - BT)")
    ax_a.set_ylabel("-log10 exact sign-flip P")
    ax_a.set_title("A  All 657 independent groups", loc="left", fontweight="bold")
    ax_a.scatter([], [], s=28, c=UP, label="Treatment-up (n=70)")
    ax_a.scatter([], [], s=28, c=DOWN, label="Treatment-down (n=93)")
    ax_a.scatter([], [], s=22, c=NEUTRAL, label="Did not pass (n=494)")
    ax_a.legend(frameon=False, fontsize=8, loc="upper left")
    ax_a.spines[["top", "right"]].set_visible(False)

    for patient, bt_column, ot_column in zip(patients, bt_columns, ot_columns):
        ax_b.plot([0, 1], [up_score[bt_column], up_score[ot_column]], color="#8CC8BE", linewidth=1.2, alpha=0.85)
        ax_b.scatter([0, 1], [up_score[bt_column], up_score[ot_column]], color=UP, s=22, zorder=3)
    ax_b.set_xticks([0, 1], ["BT", "OT"])
    ax_b.set_xlim(-0.25, 1.25)
    ax_b.set_ylabel("Mean group-wise z score")
    ax_b.set_title("B  Frozen treatment-up program", loc="left", fontweight="bold")
    ax_b.spines[["top", "right"]].set_visible(False)

    for patient, bt_column, ot_column in zip(patients, bt_columns, ot_columns):
        ax_c.plot([0, 1], [down_score[bt_column], down_score[ot_column]], color="#F0A477", linewidth=1.2, alpha=0.85)
        ax_c.scatter([0, 1], [down_score[bt_column], down_score[ot_column]], color=DOWN, s=22, zorder=3)
    ax_c.set_xticks([0, 1], ["BT", "OT"])
    ax_c.set_xlim(-0.25, 1.25)
    ax_c.set_ylabel("Mean group-wise z score")
    ax_c.set_title("C  Frozen treatment-down program", loc="left", fontweight="bold")
    ax_c.spines[["top", "right"]].set_visible(False)

    norm = TwoSlopeNorm(vmin=-2.5, vcenter=0, vmax=2.5)
    image = ax_d.imshow(heat.to_numpy(float), aspect="auto", interpolation="nearest", cmap="RdBu_r", norm=norm)
    ax_d.axhline(len(up_groups) - 0.5, color="black", linewidth=1.2)
    labels = []
    for patient in patients:
        labels.extend([f"P{patient}\nBT", f"P{patient}\nOT"])
    ax_d.set_xticks(np.arange(len(labels)), labels)
    ax_d.set_yticks([len(up_groups) / 2, len(up_groups) + len(down_groups) / 2], ["70 treatment-up groups", "93 treatment-down groups"])
    ax_d.tick_params(axis="y", length=0)
    ax_d.set_title("D  Patient-paired activity profiles for all 163 frozen dynamic groups", loc="left", fontweight="bold")
    colorbar = fig.colorbar(image, ax=ax_d, fraction=0.018, pad=0.015)
    colorbar.set_label("Within-group z score")

    fig.suptitle(
        "Broad SecAct remodeling after one cycle of combination treatment",
        fontsize=16,
        fontweight="bold",
        y=0.99,
    )
    fig.text(
        0.5,
        0.01,
        "Inference unit: nine paired patients. Dynamic status required exact sign-flip FDR <0.10, |Hedges gz| >=1.0 and bootstrap same-direction fraction >=0.90.",
        ha="center",
        fontsize=8.5,
        color="#444444",
    )
    save_figure(fig, output_prefix)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("paired_results", type=Path)
    parser.add_argument("group_activity", type=Path)
    parser.add_argument("figure_dir", type=Path)
    args = parser.parse_args()

    plot_overview(args.figure_dir / "NPC_SecAct_manuscript_figure1_study_flow_v0.1")
    plot_bulk(
        args.paired_results,
        args.group_activity,
        args.figure_dir / "NPC_SecAct_manuscript_figure2_paired_bulk_v0.1",
    )
    print("MANUSCRIPT FIGURES COMPLETE", flush=True)


if __name__ == "__main__":
    main()
