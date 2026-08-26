from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "outputs"


def main() -> None:
    result = pd.read_csv(OUTPUT / "NPC_SecAct_GSE102349_survival_results_v0.1.tsv", sep="\t")
    result = result.loc[result["association"].isin(["up_secact", "down_secact", "up_direct", "down_direct"])].copy()
    labels = {
        "up_secact": "Treatment-up SecAct",
        "down_secact": "Treatment-down SecAct",
        "up_direct": "Treatment-up direct expression",
        "down_direct": "Treatment-down direct expression",
    }
    colors = {
        "up_secact": "#1B6CA8",
        "down_secact": "#6E7781",
        "up_direct": "#C6533D",
        "down_direct": "#8B5A2B",
    }
    y = np.arange(len(result))[::-1]
    fig, ax = plt.subplots(figsize=(7.6, 4.2), constrained_layout=True)
    for yi, (_, row) in zip(y, result.iterrows()):
        ax.plot([row["score_ci_low"], row["score_ci_high"]], [yi, yi], color=colors[row["association"]], linewidth=2.0)
        ax.scatter(row["score_hr"], yi, color=colors[row["association"]], edgecolor="white", linewidth=0.7, s=52, zorder=3)
        ax.text(
            1.02,
            yi,
            f"HR {row['score_hr']:.2f}  [{row['score_ci_low']:.2f}, {row['score_ci_high']:.2f}]\nP={row['score_p']:.3f}; FDR={row['fdr_bh']:.3f}",
            transform=ax.get_yaxis_transform(),
            va="center",
            fontsize=8,
        )
    ax.axvline(1.0, color="#3F474D", linewidth=0.9, linestyle="--")
    ax.set_xscale("log")
    ax.set_xticks([0.3, 0.5, 1.0, 2.0, 3.0])
    ax.set_xticklabels(["0.3", "0.5", "1", "2", "3"])
    ax.set_yticks(y, [labels[name] for name in result["association"]])
    ax.set_xlabel("Hazard ratio per 1 SD higher score (log scale)")
    ax.set_title("GSE102349 public time-to-event association (n=88; 16 events)", fontsize=11, fontweight="bold")
    ax.grid(axis="x", color="#D9DEE2", linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.text(
        0.0,
        -0.27,
        "Continuous prespecified scores; no cutoff or prediction model. No association passed BH-FDR < 0.05.",
        transform=ax.transAxes,
        fontsize=8,
        color="#4D5961",
    )
    base = OUTPUT / "figures" / "NPC_SecAct_Q1Q2_figure7_GSE102349_survival_v0.1"
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{base}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{base}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"WROTE {base}.png")
    print(f"WROTE {base}.pdf")


if __name__ == "__main__":
    main()
