from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


COLORS = {"MET": "#D97706", "nonMET": "#197278"}
METHOD_COLORS = {"SecAct": "#3A6EA5", "Direct expression": "#B24C63"}


def grouped_panel(ax: plt.Axes, frame: pd.DataFrame, column: str, title: str, annotation: str) -> None:
    rng = np.random.default_rng(20260823)
    positions = {"nonMET": 0, "MET": 1}
    for group in ["nonMET", "MET"]:
        values = frame.loc[frame["cell_type"].eq(group), column].to_numpy(float)
        x = positions[group] + rng.uniform(-0.09, 0.09, len(values))
        ax.scatter(x, values, s=48, color=COLORS[group], edgecolor="white", linewidth=0.6, zorder=3)
        q1, median, q3 = np.quantile(values, [0.25, 0.5, 0.75])
        ax.plot([positions[group] - 0.18, positions[group] + 0.18], [median, median], color="#222222", lw=2.0, zorder=4)
        ax.plot([positions[group], positions[group]], [q1, q3], color="#222222", lw=4.5, solid_capstyle="butt", zorder=2)
    ax.axhline(0, color="#A0A0A0", lw=0.8, ls="--")
    ax.set_xticks([0, 1], ["nonMET\n(n=10)", "MET\n(n=14)"])
    ax.set_ylabel("Frozen 70-group composite score")
    ax.set_title(title, loc="left", fontweight="bold")
    ax.text(0.03, 0.97, annotation, transform=ax.transAxes, ha="left", va="top", fontsize=8.5)
    ax.spines[["top", "right"]].set_visible(False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("analysis_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    scores = pd.read_csv(args.analysis_dir / "GSE299775_patient_level_program_scores.tsv", sep="\t")
    scores = scores.loc[scores["is_npc"]].copy()
    primary = pd.read_csv(
        args.analysis_dir / "GSE299775_primary_and_benchmark_statistics.tsv", sep="\t", index_col="program"
    )
    bootstrap = pd.read_csv(args.analysis_dir / "GSE299775_paired_bootstrap_effects.tsv.gz", sep="\t")
    secact = primary.loc["SecAct_treatment_up_composite"]
    direct = primary.loc["direct_expression_treatment_up_composite"]

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 10.5,
        "axes.labelsize": 9.5,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "axes.linewidth": 0.8,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "pdf.fonttype": 42,
    })
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.5), constrained_layout=True)

    grouped_panel(
        axes[0, 0], scores, "secact_treatment_up_composite",
        "A  SecAct composite",
        f"Hedges g={secact.hedges_g:.3f}\nexact P={secact.exact_permutation_p_two_sided:.4f}",
    )
    grouped_panel(
        axes[0, 1], scores, "direct_treatment_up_composite",
        "B  Direct-expression composite",
        f"Hedges g={direct.hedges_g:.3f}\nexact P={direct.exact_permutation_p_two_sided:.4f}",
    )

    ax = axes[1, 0]
    effects = [secact.hedges_g, direct.hedges_g]
    lows = [secact.hedges_g_boot_ci_low, direct.hedges_g_boot_ci_low]
    highs = [secact.hedges_g_boot_ci_high, direct.hedges_g_boot_ci_high]
    labels = ["SecAct", "Direct expression"]
    for y, label, effect, low, high in zip([1, 0], labels, effects, lows, highs):
        ax.plot([low, high], [y, y], color=METHOD_COLORS[label], lw=2.2)
        ax.scatter(effect, y, color=METHOD_COLORS[label], s=60, edgecolor="white", linewidth=0.6, zorder=3)
    ax.axvline(0, color="#555555", lw=0.9, ls="--")
    ax.set_yticks([1, 0], labels)
    ax.set_xlabel("Hedges g (MET minus nonMET), bootstrap 95% CI")
    delta_low, delta_high = bootstrap["absolute_effect_difference"].quantile([0.025, 0.975])
    delta = abs(secact.hedges_g) - abs(direct.hedges_g)
    ax.text(0.03, 0.05, f"|g| difference={delta:.3f}\n95% CI {delta_low:.3f} to {delta_high:.3f}", transform=ax.transAxes, ha="left", va="bottom", fontsize=8.5)
    ax.set_title("C  Prespecified effect-size benchmark", loc="left", fontweight="bold")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)

    ax = axes[1, 1]
    for group in ["nonMET", "MET"]:
        current = scores.loc[scores["cell_type"].eq(group)]
        ax.scatter(
            current["direct_treatment_up_composite"], current["secact_treatment_up_composite"],
            s=52, color=COLORS[group], edgecolor="white", linewidth=0.6, label=f"{group} (n={len(current)})",
        )
    rho, p_value = stats.spearmanr(
        scores["direct_treatment_up_composite"], scores["secact_treatment_up_composite"]
    )
    ax.axhline(0, color="#A0A0A0", lw=0.7, ls="--")
    ax.axvline(0, color="#A0A0A0", lw=0.7, ls="--")
    ax.set_xlabel("Direct-expression composite")
    ax.set_ylabel("SecAct composite")
    ax.set_title("D  Patient-level score concordance", loc="left", fontweight="bold")
    ax.text(0.03, 0.97, f"Spearman rho={rho:.3f}\nP={p_value:.5f}", transform=ax.transAxes, ha="left", va="top", fontsize=8.5)
    ax.legend(frameon=False, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle(
        "Frozen treatment-up program in GSE299775 PanCK-positive NPC compartments",
        fontsize=12, fontweight="bold",
    )
    prefix = args.output_dir / "NPC_SecAct_GSE299775_Original_Article_figure_v0.1"
    fig.savefig(f"{prefix}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(f"{prefix}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"WROTE {prefix}.png")
    print(f"WROTE {prefix}.pdf")


if __name__ == "__main__":
    main()
