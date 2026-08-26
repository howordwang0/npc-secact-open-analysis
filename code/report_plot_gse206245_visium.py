from __future__ import annotations

import argparse
import gzip
import io
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image


RESPONSE_COLORS = {
    "CR": "#2166AC",
    "PR": "#67A9CF",
    "SD": "#FDAE61",
    "PD": "#B2182B",
}


def save_figure(fig: plt.Figure, prefix: Path) -> None:
    fig.savefig(f"{prefix}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{prefix}.pdf", bbox_inches="tight")
    plt.close(fig)


def read_lowres(data_dir: Path, gsm: str, sample: str) -> tuple[Image.Image, float]:
    prefix = f"{gsm}_NPC_{sample}"
    image_path = data_dir / f"{prefix}_tissue_lowres_image.png.gz"
    scale_path = data_dir / f"{prefix}_scalefactors_json.json.gz"
    with gzip.open(image_path, "rb") as handle:
        image = Image.open(io.BytesIO(handle.read())).convert("RGB")
    with gzip.open(scale_path, "rt", encoding="utf-8") as handle:
        scale = float(json.load(handle)["tissue_lowres_scalef"])
    return image, scale


def robust_limits(values: pd.Series) -> tuple[float, float]:
    lower, upper = np.nanquantile(values.to_numpy(float), [0.02, 0.98])
    if lower == upper:
        lower -= 1
        upper += 1
    return float(lower), float(upper)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("spot_scores_tsv_gz", type=Path)
    parser.add_argument("patient_correlations_tsv", type=Path)
    parser.add_argument("summary_tsv", type=Path)
    parser.add_argument("coverage_tsv", type=Path)
    parser.add_argument("response_secondary_tsv", type=Path)
    parser.add_argument("run_json", type=Path)
    parser.add_argument("data_dir", type=Path)
    parser.add_argument("figure_dir", type=Path)
    parser.add_argument("report_output", type=Path)
    args = parser.parse_args()
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    args.report_output.parent.mkdir(parents=True, exist_ok=True)

    spots = pd.read_csv(args.spot_scores_tsv_gz, sep="\t")
    patients = pd.read_csv(args.patient_correlations_tsv, sep="\t")
    summary = pd.read_csv(args.summary_tsv, sep="\t")
    coverage = pd.read_csv(args.coverage_tsv, sep="\t")
    response = pd.read_csv(args.response_secondary_tsv, sep="\t")
    run = json.loads(args.run_json.read_text(encoding="utf-8"))

    main_patients = patients.loc[
        (patients["scenario"] == "source_filtered")
        & (patients["set_label"] == "dynamic_up_set")
    ].copy()
    qc_patients = patients.loc[
        (patients["scenario"] == "fixed_spot_qc")
        & (patients["set_label"] == "dynamic_up_set")
    ].copy()
    main = summary.loc[
        (summary["scenario"] == "source_filtered")
        & (summary["set_label"] == "dynamic_up_set")
    ].iloc[0]
    sensitivity = summary.loc[
        (summary["scenario"] == "fixed_spot_qc")
        & (summary["set_label"] == "dynamic_up_set")
    ].iloc[0]
    main_down = summary.loc[
        (summary["scenario"] == "source_filtered")
        & (summary["set_label"] == "dynamic_down_set")
    ].iloc[0]
    sensitivity_down = summary.loc[
        (summary["scenario"] == "fixed_spot_qc")
        & (summary["set_label"] == "dynamic_down_set")
    ].iloc[0]

    sample_order = [
        "ST5",
        "ST6",
        "ST7",
        "ST8",
        "ST9",
        "ST10",
        "ST11",
        "ST12",
        "ST16",
        "ST17",
        "ST18",
        "ST19",
    ]
    order_lookup = {sample: index for index, sample in enumerate(sample_order)}
    main_patients["order"] = main_patients["sample"].map(order_lookup)
    main_patients = main_patients.sort_values("order")
    qc_lookup = qc_patients.set_index("sample")["spearman_rho"]

    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )
    fig, ax = plt.subplots(figsize=(9.0, 4.8))
    x = np.arange(len(main_patients))
    for index, row in enumerate(main_patients.itertuples(index=False)):
        qc_value = float(qc_lookup.loc[row.sample])
        ax.plot(
            [index, index],
            [row.spearman_rho, qc_value],
            color="#BDBDBD",
            linewidth=1,
            zorder=1,
        )
        ax.scatter(
            index,
            row.spearman_rho,
            color=RESPONSE_COLORS[row.response],
            edgecolor="black",
            linewidth=0.4,
            s=55,
            zorder=3,
        )
        ax.scatter(
            index,
            qc_value,
            facecolor="white",
            edgecolor=RESPONSE_COLORS[row.response],
            linewidth=1.2,
            s=43,
            zorder=2,
        )
    ax.axhline(0, color="#333333", linewidth=0.8)
    ax.axhline(
        main["median_patient_rho"],
        color="#4D4D4D",
        linestyle="--",
        linewidth=1,
        label="Main median",
    )
    for response_label, color in RESPONSE_COLORS.items():
        ax.scatter([], [], color=color, edgecolor="black", s=40, label=response_label)
    ax.scatter([], [], facecolor="white", edgecolor="#555555", s=35, label="Fixed QC")
    ax.set_xticks(x, main_patients["sample"], rotation=45, ha="right")
    ax.set_ylabel("Within-patient Spearman rho")
    ax.set_title(
        "GSE206245: dynamic-up SecAct activity and myeloid-state context\n"
        f"median rho={main['median_patient_rho']:.3f}; "
        f"positive={int(main['positive_patients'])}/{int(main['n_patients'])}; "
        f"Wilcoxon P={main['wilcoxon_p_two_sided']:.4f}"
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, ncol=3, loc="best", fontsize=8)
    fig.tight_layout()
    save_figure(
        fig,
        args.figure_dir / "NPC_SecAct_GSE206245_patient_spatial_correlations_v0.2",
    )

    fig, axes = plt.subplots(4, 6, figsize=(17, 11))
    for sample_index, sample in enumerate(sample_order):
        row_index = sample_index // 3
        pair_index = sample_index % 3
        current = spots.loc[spots["sample"] == sample]
        gsm = current["gsm"].iloc[0]
        response_label = current["response"].iloc[0]
        image, scale = read_lowres(args.data_dir, gsm, sample)
        x_coord = current["pxl_col_in_fullres"].to_numpy(float) * scale
        y_coord = current["pxl_row_in_fullres"].to_numpy(float) * scale
        for offset, feature, cmap, label in (
            (0, "dynamic_up_set", "viridis", "Dynamic-up SecAct"),
            (1, "myeloid_context", "magma", "Myeloid-state context"),
        ):
            ax = axes[row_index, pair_index * 2 + offset]
            ax.imshow(image)
            vmin, vmax = robust_limits(current[feature])
            ax.scatter(
                x_coord,
                y_coord,
                c=current[feature],
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                s=5,
                linewidths=0,
                alpha=0.82,
            )
            suffix = f"{sample} ({response_label})" if offset == 0 else label
            ax.set_title(suffix, color=RESPONSE_COLORS[response_label] if offset == 0 else "black")
            ax.set_axis_off()
    fig.suptitle(
        "All 12 GSE206245 NPC sections (within-section color scaling)",
        fontsize=13,
        y=0.995,
    )
    fig.tight_layout()
    save_figure(
        fig,
        args.figure_dir / "NPC_SecAct_GSE206245_all_section_spatial_maps_v0.2",
    )

    if run["formal_support"]:
        status_text = "met the frozen formal spatial-support rule"
        interpretation = (
            "This provides statistically supported spatial corroboration of a "
            "myeloid-rich context for the bulk treatment-up program."
        )
    elif run["directional_support"]:
        status_text = "met the frozen directional rule but not the formal P < 0.05 rule"
        interpretation = (
            "This is directional spatial corroboration, not statistically "
            "conclusive spatial validation."
        )
    else:
        status_text = "did not meet the frozen spatial-support rule"
        interpretation = (
            "The independent Visium cohort therefore does not corroborate the "
            "proposed positive myeloid spatial context."
        )

    spatial_response = response.loc[
        response["analysis"] == "patient_spatial_rho_by_response"
    ].iloc[0]
    tissue_response = response.loc[
        response["analysis"] == "tissue_dynamic_up_by_response"
    ].iloc[0]
    min_coverage = coverage["signature_coverage_fraction"].min()
    max_coverage = coverage["signature_coverage_fraction"].max()
    total_source_spots = sum(item["source_spots"] for item in run["samples"])
    total_qc_spots = sum(item["fixed_qc_spots"] for item in run["samples"])

    report = f"""# GSE206245 Visium spatial contextualization report v0.2

## Result status

The analysis **{status_text}**. Across {int(main['n_patients'])} independent NPC
patients, the median within-section Spearman correlation between the frozen
dynamic-up SecAct score and the independent myeloid-state context score was
{main['median_patient_rho']:.3f}; {int(main['positive_patients'])}/{int(main['n_patients'])}
patients had a positive correlation. The exact two-sided Wilcoxon P value was
{main['wilcoxon_p_two_sided']:.4f}, and the exact two-sided sign-test P value
was {main['exact_sign_p_two_sided']:.4f}. {interpretation}

## Data and processing

- Twelve baseline FFPE NPC sections from 12 patients (`ST5-ST12`, `ST16-ST19`)
  were reconstructed from the released Matrix Market, barcode, feature, and
  tissue-position files.
- The main analysis retained {total_source_spots:,} source-filtered in-tissue
  spots. Fixed QC retained {total_qc_spots:,} spots with at least 200 detected
  features and 500 total counts.
- ST12 contained 3,668 mutually concordant released matrix columns, barcodes,
  and in-tissue position records, versus 3,669 in Supplementary Data 1. The
  analysis used every observable released spot and did not impute the missing
  publication-level record.
- SecAct signature coverage ranged from {min_coverage:.3%} to
  {max_coverage:.3%} across sections, so the frozen 80% full-inference gate
  {'passed' if min_coverage >= 0.8 else 'failed'}.
- SecAct was run separately within each section. The biological inference unit
  was the patient; spot-level P values were not used.
- The original 14-gene context score in protocol v0.1 was non-computable after
  predictor-overlap removal (2 markers remained versus the required 8). Before
  viewing any GSE206245 activity or spatial endpoint, protocol v0.2 froze a
  20-gene context set from six independent baseline single-cell patients; all
  SecAct signature genes were excluded and marker selection did not use this
  spatial cohort. At least 16 non-zero-variance markers were required per
  scenario. This score represents myeloid state/context, not a calibrated cell
  fraction.

## Frozen endpoint and sensitivity

| Scenario | Patients | Median rho | Positive patients | Minimum LOO median | Wilcoxon P (two-sided) |
|---|---:|---:|---:|---:|---:|
| Source-filtered, dynamic-up | {int(main['n_patients'])} | {main['median_patient_rho']:.3f} | {int(main['positive_patients'])}/{int(main['n_patients'])} | {main['minimum_leave_one_out_median_rho']:.3f} | {main['wilcoxon_p_two_sided']:.4f} |
| Fixed spot QC, dynamic-up | {int(sensitivity['n_patients'])} | {sensitivity['median_patient_rho']:.3f} | {int(sensitivity['positive_patients'])}/{int(sensitivity['n_patients'])} | {sensitivity['minimum_leave_one_out_median_rho']:.3f} | {sensitivity['wilcoxon_p_two_sided']:.4f} |
| Source-filtered, dynamic-down | {int(main_down['n_patients'])} | {main_down['median_patient_rho']:.3f} | {int(main_down['positive_patients'])}/{int(main_down['n_patients'])} | {main_down['minimum_leave_one_out_median_rho']:.3f} | {main_down['wilcoxon_p_two_sided']:.4f} |
| Fixed spot QC, dynamic-down | {int(sensitivity_down['n_patients'])} | {sensitivity_down['median_patient_rho']:.3f} | {int(sensitivity_down['positive_patients'])}/{int(sensitivity_down['n_patients'])} | {sensitivity_down['minimum_leave_one_out_median_rho']:.3f} | {sensitivity_down['wilcoxon_p_two_sided']:.4f} |

The dynamic-down set was also positively associated with the myeloid-state
context in the source-filtered analysis (median rho =
{main_down['median_patient_rho']:.3f}; Wilcoxon P =
{main_down['wilcoxon_p_two_sided']:.4f}), although this attenuated after fixed
spot QC (P = {sensitivity_down['wilcoxon_p_two_sided']:.4f}). The formal
dynamic-up endpoint remains supported under its frozen rules, but the shared
direction limits claims that the spatial association is specific to the
treatment-up program.

## Prespecified secondary response analyses

Subsequent treatment response did not define the primary endpoint. Comparing
the patient spatial correlations between CR/PR responders and SD/PD
non-responders gave an exact rank-permutation P value of
{spatial_response['exact_rank_permutation_p_two_sided']:.4f} (rank-biserial
effect {spatial_response['rank_biserial']:.3f}). The exploratory tissue-level
dynamic-up SecAct comparison gave P =
{tissue_response['exact_rank_permutation_p_two_sided']:.4f} (rank-biserial
effect {tissue_response['rank_biserial']:.3f}). These analyses cannot identify
an ICB-specific predictive biomarker because chemotherapy partners varied and
there was no non-ICB comparator.

## Interpretation limits

GSE206245 is a small retrospective baseline cohort. A positive spatial
association is compatible with a myeloid-rich tissue context but does not
prove myeloid secretion, molecular interaction, treatment-induced change, or
causality. The corresponding dynamic-down association further limits
program-specific interpretation. Spatial autocorrelation makes spots unsuitable as independent
replicates; all inferential language is therefore based on patient-level
summaries. This result is retained regardless of direction and does not erase
the negative GeoMx spatial result.
"""
    args.report_output.write_text(report, encoding="utf-8")
    print(f"REPORT COMPLETE: {args.report_output}", flush=True)


if __name__ == "__main__":
    main()
