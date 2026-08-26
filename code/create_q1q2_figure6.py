from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


COLORS = {
    "secact": "#1B6CA8",
    "direct": "#C6533D",
    "down": "#6E7781",
    "cytosig": "#2D8B57",
    "adjusted": "#8A5A9E",
    "null": "#B9C2C9",
}


def style_axis(axis: plt.Axes, panel: str) -> None:
    axis.text(
        -0.13,
        1.07,
        panel,
        transform=axis.transAxes,
        fontsize=12,
        fontweight="bold",
        va="top",
    )
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.tick_params(labelsize=8, width=0.8)
    axis.grid(axis="y", color="#D9DEE2", linewidth=0.6, alpha=0.7)
    axis.set_axisbelow(True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("work_dir", type=Path)
    args = parser.parse_args()

    localization = pd.read_csv(
        args.output_dir / "NPC_SecAct_GSE150430_primary_patient_differences_v0.1.tsv",
        sep="\t",
    )
    external = pd.read_csv(
        args.output_dir / "NPC_SecAct_GSE200310_external_spatial_summary_v0.1.tsv",
        sep="\t",
    )
    adjusted = pd.read_csv(
        args.output_dir
        / "NPC_SecAct_GSE206245_composition_spatial_robustness_v0.1.tsv",
        sep="\t",
    )
    random_programs = pd.read_csv(
        args.output_dir / "NPC_SecAct_GSE206245_matched_random_programs_v0.1.tsv",
        sep="\t",
    )
    random_summary = json.loads(
        (
            args.output_dir
            / "NPC_SecAct_GSE206245_matched_random_summary_v0.1.json"
        ).read_text(encoding="utf-8")
    )

    figure, axes = plt.subplots(2, 2, figsize=(11.2, 8.3), constrained_layout=True)

    axis = axes[0, 0]
    up = localization.loc[localization["set_label"].eq("dynamic_up_set")].copy()
    secact = up.loc[up["analysis"].eq("SecAct")].set_index("patient_id")
    direct = up.loc[up["analysis"].eq("Direct_expression")].set_index("patient_id")
    patients = sorted(set(secact.index) & set(direct.index))
    for patient in patients:
        values = [
            secact.loc[patient, "difference_myeloid_minus_nonmyeloid"],
            direct.loc[patient, "difference_myeloid_minus_nonmyeloid"],
        ]
        axis.plot([0, 1], values, color="#AAB2B8", linewidth=0.8, alpha=0.7, zorder=1)
    axis.scatter(
        np.zeros(len(patients)),
        secact.loc[patients, "difference_myeloid_minus_nonmyeloid"],
        s=27,
        color=COLORS["secact"],
        edgecolor="white",
        linewidth=0.5,
        zorder=2,
    )
    axis.scatter(
        np.ones(len(patients)),
        direct.loc[patients, "difference_myeloid_minus_nonmyeloid"],
        s=27,
        color=COLORS["direct"],
        edgecolor="white",
        linewidth=0.5,
        zorder=2,
    )
    axis.axhline(0, color="#3F474D", linewidth=0.8)
    axis.set_xticks([0, 1], ["SecAct", "Direct expression"])
    axis.set_ylabel("Myeloid minus pooled non-myeloid score")
    axis.set_title("Independent single-cell localization (n=14)", fontsize=10)
    axis.text(
        0.48,
        0.97,
        "SecAct: 14/14 positive\nmedian difference 0.735; P=0.00012",
        transform=axis.transAxes,
        fontsize=8,
        va="top",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 1.5},
    )
    style_axis(axis, "A")

    axis = axes[0, 1]
    program_order = [
        "dynamic_up_set",
        "direct_dynamic_up_set",
        "dynamic_down_set",
        "cytosig_supported_set",
    ]
    program_labels = ["SecAct up", "Direct up", "SecAct down", "CytoSig-29"]
    program_colors = [
        COLORS["secact"],
        COLORS["direct"],
        COLORS["down"],
        COLORS["cytosig"],
    ]
    for specimen_index, (sample, marker) in enumerate(
        [("NPC_patient_1", "o"), ("NPC_patient_2", "s")]
    ):
        specimen = external.loc[external["sample"].eq(sample)].set_index("program")
        values = specimen.loc[program_order, "spearman_rho"].to_numpy(float)
        x = np.arange(len(program_order)) + (-0.08 if specimen_index == 0 else 0.08)
        axis.plot(x, values, color="#9FA7AD", linewidth=0.7, alpha=0.55, zorder=1)
        for pos, value, color in zip(x, values, program_colors):
            axis.scatter(
                pos,
                value,
                marker=marker,
                s=38,
                color=color,
                edgecolor="white",
                linewidth=0.5,
                zorder=2,
            )
    axis.axhline(0, color="#3F474D", linewidth=0.8)
    axis.set_xticks(np.arange(len(program_order)), program_labels, rotation=20, ha="right")
    axis.set_ylabel("Spearman rho with myeloid context")
    axis.set_title("External spatial replication (n=2)", fontsize=10)
    axis.text(
        0.02,
        0.96,
        "circle: patient 1   square: patient 2",
        transform=axis.transAxes,
        fontsize=8,
        va="top",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 1.5},
    )
    style_axis(axis, "B")

    axis = axes[1, 0]
    for _, row in adjusted.iterrows():
        axis.plot(
            [0, 1],
            [row["unadjusted_rho"], row["composition_adjusted_rho"]],
            color="#AAB2B8",
            linewidth=0.8,
            alpha=0.75,
            zorder=1,
        )
    axis.scatter(
        np.zeros(len(adjusted)),
        adjusted["unadjusted_rho"],
        s=28,
        color=COLORS["secact"],
        edgecolor="white",
        linewidth=0.5,
        zorder=2,
    )
    axis.scatter(
        np.ones(len(adjusted)),
        adjusted["composition_adjusted_rho"],
        s=28,
        color=COLORS["adjusted"],
        edgecolor="white",
        linewidth=0.5,
        zorder=2,
    )
    axis.axhline(0, color="#3F474D", linewidth=0.8)
    axis.set_xticks([0, 1], ["Unadjusted", "Composition-proxy adjusted"])
    axis.set_ylabel("Patient-level Spearman rho")
    axis.set_title("GSE206245 composition robustness (n=12)", fontsize=10)
    axis.text(
        0.02,
        0.96,
        "adjusted median 0.243\n12/12 positive; P=0.00049",
        transform=axis.transAxes,
        fontsize=8,
        va="top",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 1.5},
    )
    style_axis(axis, "C")

    axis = axes[1, 1]
    null_values = random_programs["median_patient_rho"].to_numpy(float)
    observed = random_summary["observed_program"]["median_patient_rho"]
    axis.hist(
        null_values,
        bins=34,
        color=COLORS["null"],
        edgecolor="white",
        linewidth=0.5,
    )
    axis.axvline(observed, color=COLORS["direct"], linewidth=2.0)
    axis.set_xlabel("Median patient rho of matched random program")
    axis.set_ylabel("Number of programs")
    axis.set_title("Activity-matched random-program null", fontsize=10)
    axis.text(
        0.03,
        0.96,
        "observed=0.343\n98.9th percentile; P=0.012",
        transform=axis.transAxes,
        fontsize=8,
        va="top",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 1.5},
    )
    style_axis(axis, "D")

    figure.suptitle(
        "Independent localization and prespecified robustness analyses",
        fontsize=13,
        fontweight="bold",
    )
    output_base = args.output_dir / "figures" / "NPC_SecAct_Q1Q2_figure6_enhancement_v0.1"
    output_base.parent.mkdir(parents=True, exist_ok=True)
    output_png = Path(f"{output_base}.png")
    output_pdf = Path(f"{output_base}.pdf")
    figure.savefig(output_png, dpi=300, bbox_inches="tight")
    figure.savefig(output_pdf, bbox_inches="tight")
    plt.close(figure)
    print(f"WROTE {output_png}")
    print(f"WROTE {output_pdf}")


if __name__ == "__main__":
    main()
