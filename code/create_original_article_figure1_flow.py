from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


INK = "#202428"
RED = "#A8343A"
GREEN = "#26734D"
BLUE = "#246B9E"
AMBER = "#B06B16"
GRAY = "#68737C"


def box(
    axis: plt.Axes,
    x: float,
    y: float,
    width: float,
    height: float,
    title: str,
    body: str,
    edge: str,
    face: str,
    title_size: float = 9.2,
    body_size: float = 7.8,
) -> None:
    axis.add_patch(
        FancyBboxPatch(
            (x, y),
            width,
            height,
            boxstyle="round,pad=0.006,rounding_size=0.008",
            linewidth=1.5,
            edgecolor=edge,
            facecolor=face,
        )
    )
    axis.text(
        x + width * 0.04,
        y + height * 0.78,
        title,
        fontsize=title_size,
        fontweight="bold",
        color=INK,
        va="center",
    )
    axis.text(
        x + width * 0.04,
        y + height * 0.39,
        body,
        fontsize=body_size,
        color=INK,
        va="center",
        linespacing=1.28,
    )


def arrow(axis: plt.Axes, start: tuple[float, float], end: tuple[float, float]) -> None:
    axis.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=12,
            linewidth=1.2,
            color="#7A838A",
            shrinkA=3,
            shrinkB=3,
        )
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("figure_dir", type=Path)
    args = parser.parse_args()

    figure, axis = plt.subplots(figsize=(15.5, 8.8))
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")

    axis.text(0.02, 0.965, "Patient-level multimodal evidence flow", fontsize=18, fontweight="bold", color=INK, va="top")
    axis.text(0.02, 0.922, "Fully open data; cells, ROIs and spots were not treated as independent patients", fontsize=9.5, color="#566068", va="top")

    axis.text(0.025, 0.85, "Clinical routes", fontsize=10, fontweight="bold", color=RED)
    box(axis, 0.02, 0.60, 0.17, 0.20, "Response screen", "GSE132112 | n=95\n0/657 groups\nminimum FDR 0.9985", RED, "#FFF5F5")
    box(axis, 0.02, 0.33, 0.17, 0.20, "Metastasis screen", "GSE103611 | n=48\n0/657 groups\nminimum FDR 0.6097", RED, "#FFF5F5")
    box(axis, 0.02, 0.065, 0.17, 0.17, "Public time-to-event", "GSE102349 | n=88 eligible\n16 events; 25 missing\nNo association FDR<0.05\nExploratory only", AMBER, "#FFF9EC", title_size=8.7, body_size=7.2)
    axis.text(0.105, 0.275, "Biomarker route stopped", fontsize=8.3, color=RED, fontweight="bold", ha="center")

    box(axis, 0.225, 0.49, 0.18, 0.25, "Frozen treatment programs", "Paired bulk | n=9\n70 treatment-up groups\n93 treatment-down groups\nBroad bidirectional remodeling", BLUE, "#F3F8FC")
    arrow(axis, (0.19, 0.64), (0.225, 0.64))

    axis.text(0.445, 0.85, "Independent localization", fontsize=10, fontweight="bold", color=GREEN)
    box(axis, 0.44, 0.60, 0.17, 0.20, "Paired single-cell", "n=5 paired\nmyeloid up: 4/5\ngz=0.267; P=0.50\nDirectional", AMBER, "#FFF9EC")
    box(axis, 0.44, 0.33, 0.17, 0.20, "Independent single-cell", "GSE150430 | n=14\n14/14 positive\nmedian difference 0.735\nP=0.00012", GREEN, "#F3FAF6")
    arrow(axis, (0.405, 0.615), (0.44, 0.69))
    arrow(axis, (0.405, 0.615), (0.44, 0.43))

    axis.text(0.65, 0.85, "Spatial contextualization", fontsize=10, fontweight="bold", color=GREEN)
    box(axis, 0.645, 0.60, 0.155, 0.20, "GeoMx TLS", "n=9 NPC\n9.4% predictor coverage\npartial rho=-0.667\nNo support", RED, "#FFF5F5", title_size=8.8, body_size=7.5)
    box(axis, 0.645, 0.33, 0.155, 0.20, "Primary Visium", "GSE206245 | n=12\nmedian rho 0.343\n11/12; P=0.00098\nFormal support", GREEN, "#F3FAF6", title_size=8.8, body_size=7.5)
    box(axis, 0.825, 0.60, 0.155, 0.20, "Clinical spatial endpoint", "GSE299775 | n=24\nSecAct g=-0.738; P=0.0789\nDirect g=-1.023; P=0.0188\nNo SecAct increment", AMBER, "#FFF9EC", title_size=8.3, body_size=7.1)
    box(axis, 0.825, 0.33, 0.155, 0.20, "External Visium", "GSE200310 | n=2\nrho 0.237 and 0.362\nBoth positive\nDescriptive replication", GREEN, "#F3FAF6", title_size=8.8, body_size=7.5)
    arrow(axis, (0.61, 0.69), (0.645, 0.69))
    arrow(axis, (0.61, 0.43), (0.645, 0.43))
    arrow(axis, (0.80, 0.69), (0.825, 0.69))
    arrow(axis, (0.80, 0.43), (0.825, 0.43))

    axis.text(0.65, 0.265, "Prespecified robustness and comparators", fontsize=10, fontweight="bold", color=BLUE)
    box(axis, 0.645, 0.065, 0.105, 0.15, "Composition", "adjusted rho 0.243\n12/12; P=0.00049", BLUE, "#F3F8FC", title_size=8.2, body_size=7.0)
    box(axis, 0.76, 0.065, 0.105, 0.15, "Spatial null", "999 shifts\ncohort P=0.001", BLUE, "#F3F8FC", title_size=8.2, body_size=7.0)
    box(axis, 0.875, 0.065, 0.105, 0.15, "Matched programs", "1,000 controls\n98.9th percentile\nP=0.012", BLUE, "#F3F8FC", title_size=8.2, body_size=7.0)
    arrow(axis, (0.722, 0.33), (0.70, 0.215))
    arrow(axis, (0.722, 0.33), (0.812, 0.215))
    arrow(axis, (0.722, 0.33), (0.925, 0.215))

    axis.text(0.225, 0.245, "Final claim boundary", fontsize=10, fontweight="bold", color=INK)
    axis.text(
        0.225,
        0.185,
        "Supported: frozen treatment-up activity localizes to myeloid-rich\n"
        "cellular and spatial contexts and remains associated under measured\n"
        "composition, spatial and matched-program controls.",
        fontsize=8.5,
        color=GREEN,
        va="center",
        linespacing=1.35,
    )
    axis.text(
        0.225,
        0.095,
        "Not established: secretion, cellular source, treatment causality,\n"
        "a single mediator, clinical prediction or SecAct superiority.\n"
        "Direct expression was stronger in all 14 Visium specimens and\n"
        "associated in GSE299775, without a conclusive effect difference.",
        fontsize=8.2,
        color=RED,
        va="center",
        linespacing=1.35,
    )

    output = args.figure_dir / "NPC_SecAct_Original_Article_figure1_study_flow_v0.1"
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(f"{output}.png", dpi=300, bbox_inches="tight", facecolor="white")
    figure.savefig(f"{output}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(figure)
    print(f"WROTE {output}.png")
    print(f"WROTE {output}.pdf")


if __name__ == "__main__":
    main()
