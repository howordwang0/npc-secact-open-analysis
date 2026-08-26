from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


COLORS = {"CR": "#2A6F97", "PR": "#5FA8D3", "SD": "#C65D3B"}


def short_label(value: str, limit: int = 35) -> str:
    return value if len(value) <= limit else f"{value[: limit - 3]}..."


def save_figure(fig: plt.Figure, prefix: Path) -> None:
    fig.savefig(f"{prefix}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{prefix}.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Create publication-oriented GSE132112 discovery figures.")
    parser.add_argument("primary_stats", type=Path)
    parser.add_argument("group_activity", type=Path)
    parser.add_argument("manifest_csv", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()

    stats = pd.read_csv(args.primary_stats, sep="\t", index_col=0)
    activity = pd.read_csv(args.group_activity, sep="\t", index_col=0)
    manifest = pd.read_csv(args.manifest_csv).set_index("expression_sample_id")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="notebook")

    plot_stats = stats.copy()
    plot_stats["minus_log10_fdr"] = -np.log10(plot_stats["welch_fdr_bh"].clip(lower=np.finfo(float).tiny))
    fig, ax = plt.subplots(figsize=(8.2, 6.2))
    candidate = plot_stats["candidate_longlist"].astype(bool)
    ax.scatter(
        plot_stats.loc[~candidate, "hedges_g"],
        plot_stats.loc[~candidate, "minus_log10_fdr"],
        s=18,
        color="#9AA0A6",
        alpha=0.65,
        linewidth=0,
    )
    ax.scatter(
        plot_stats.loc[candidate, "hedges_g"],
        plot_stats.loc[candidate, "minus_log10_fdr"],
        s=30,
        color="#B33A3A",
        alpha=0.9,
        linewidth=0,
    )
    ax.axhline(-np.log10(0.05), color="#333333", linestyle="--", linewidth=1)
    ax.axvline(-0.5, color="#777777", linestyle=":", linewidth=1)
    ax.axvline(0.5, color="#777777", linestyle=":", linewidth=1)
    labels = plot_stats.loc[candidate].sort_values("welch_fdr_bh").head(10)
    for feature, row in labels.iterrows():
        ax.annotate(short_label(feature, 28), (row["hedges_g"], row["minus_log10_fdr"]), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.set(xlabel="Hedges g (CR/PR minus SD)", ylabel="-log10(BH FDR)", title="GSE132112 SecAct discovery")
    save_figure(fig, args.output_dir / "GSE132112_secact_volcano")

    candidate_features = stats.index[stats["candidate_longlist"].astype(bool)].tolist()
    if candidate_features:
        selected = candidate_features[: min(20, len(candidate_features))]
        selection_label = "pre-specified candidate longlist"
    else:
        selected = stats.sort_values("welch_p").head(20).index.tolist()
        selection_label = "top-ranked exploratory groups; no group passed the candidate rule"
    ordered_samples = manifest.sort_values("response_original").index.tolist()
    heat = activity.loc[selected, ordered_samples]
    heat = heat.subtract(heat.mean(axis=1), axis=0).divide(heat.std(axis=1, ddof=1).replace(0, np.nan), axis=0)
    fig, ax = plt.subplots(figsize=(12, max(4.5, 0.35 * len(selected))))
    sns.heatmap(
        heat,
        cmap="vlag",
        center=0,
        xticklabels=False,
        yticklabels=[short_label(item, 40) for item in heat.index],
        cbar_kws={"label": "Row-standardized activity"},
        ax=ax,
    )
    ax.set(xlabel=f"Patients ordered by response ({selection_label})", ylabel="SecAct signature group")
    save_figure(fig, args.output_dir / "GSE132112_secact_heatmap")

    top_features = (candidate_features if candidate_features else stats.index.tolist())[:6]
    long = (
        activity.loc[top_features]
        .T.join(manifest[["response_original"]])
        .reset_index(names="sample_id")
        .melt(id_vars=["sample_id", "response_original"], var_name="signature_group", value_name="activity_zscore")
    )
    long["signature_group"] = long["signature_group"].map(lambda value: short_label(value, 30))
    fig, axes = plt.subplots(2, 3, figsize=(13, 8), squeeze=False)
    for ax, (feature, block) in zip(axes.ravel(), long.groupby("signature_group", sort=False)):
        sns.boxplot(
            data=block,
            x="response_original",
            y="activity_zscore",
            order=["CR", "PR", "SD"],
            palette=COLORS,
            hue="response_original",
            legend=False,
            width=0.6,
            fliersize=0,
            ax=ax,
        )
        sns.stripplot(
            data=block,
            x="response_original",
            y="activity_zscore",
            order=["CR", "PR", "SD"],
            color="#222222",
            size=2.5,
            alpha=0.6,
            jitter=0.18,
            ax=ax,
        )
        ax.set(title=feature, xlabel="", ylabel="SecAct z-score")
    for ax in axes.ravel()[len(top_features) :]:
        ax.axis("off")
    fig.tight_layout()
    save_figure(fig, args.output_dir / "GSE132112_top_group_distributions")

    summary = {
        "candidate_groups": len(candidate_features),
        "heatmap_groups": selected,
        "distribution_groups": top_features,
        "no_candidate_fallback_used": not bool(candidate_features),
        "fallback_warning": "Top-ranked groups are exploratory and are not called candidates" if not candidate_features else None,
    }
    (args.output_dir / "GSE132112_figure_manifest.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
