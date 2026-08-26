from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from analyze_gse206245_visium_spatial import SAMPLES, read_slice, standardize_rows
from analyze_gse200310_external_spatial import toroidal_shift_null


SEED = 20260823
PROXY_ORDER = ["Myeloid", "Lymphoid", "Epithelial", "Stromal"]


def rank(values: np.ndarray) -> np.ndarray:
    return stats.rankdata(np.asarray(values, dtype=float), method="average")


def residualize(values: np.ndarray, design: np.ndarray) -> np.ndarray:
    return rank(values) - design @ np.linalg.lstsq(design, rank(values), rcond=None)[0]


def vif_values(design_without_intercept: np.ndarray) -> list[float]:
    result = []
    for target in range(design_without_intercept.shape[1]):
        y = design_without_intercept[:, target]
        x = np.delete(design_without_intercept, target, axis=1)
        x = np.column_stack([np.ones(len(x)), x])
        fitted = x @ np.linalg.lstsq(x, y, rcond=None)[0]
        ss_total = np.sum((y - y.mean()) ** 2)
        r_squared = 1 - np.sum((y - fitted) ** 2) / ss_total
        result.append(float(1 / (1 - r_squared)) if r_squared < 1 else np.inf)
    return result


def patient_summary(values: np.ndarray) -> dict:
    values = np.asarray(values, dtype=float)
    loo = np.asarray([np.median(np.delete(values, index)) for index in range(len(values))])
    positive = int((values > 0).sum())
    wilcoxon = stats.wilcoxon(values, alternative="two-sided", method="exact")
    return {
        "n_patients": len(values),
        "median_patient_rho": float(np.median(values)),
        "mean_patient_rho": float(values.mean()),
        "positive_patients": positive,
        "positive_fraction": positive / len(values),
        "minimum_leave_one_out_median_rho": float(loo.min()),
        "all_leave_one_out_medians_positive": bool((loo > 0).all()),
        "exact_wilcoxon_p_two_sided": float(wilcoxon.pvalue),
        "exact_sign_p_two_sided": float(
            stats.binomtest(positive, len(values), 0.5).pvalue
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("data_dir", type=Path)
    parser.add_argument("prior_spot_scores", type=Path)
    parser.add_argument("composition_markers", type=Path)
    parser.add_argument("work_output_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.work_output_dir.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    prior = pd.read_csv(args.prior_spot_scores, sep="\t")
    markers = pd.read_csv(args.composition_markers, sep="\t")
    if markers.groupby("compartment").size().to_dict() != {
        "Epithelial": 30,
        "Lymphoid": 30,
        "Myeloid": 30,
        "Stromal": 30,
    }:
        raise ValueError("Frozen composition proxy marker dimensions changed")

    spot_frames = []
    patient_rows = []
    coverage_rows = []
    null_frames = []
    rng = np.random.default_rng(SEED)

    for metadata in SAMPLES:
        sample = metadata["sample"]
        print(f"ROBUSTNESS LOAD {sample}", flush=True)
        released = read_slice(args.data_dir, metadata)
        current = prior.loc[prior["sample"].eq(sample)].copy()
        if current["barcode"].tolist() != released["barcodes"].tolist():
            raise ValueError(f"Prior-score barcode order mismatch for {sample}")
        genes = released["genes"]
        lookup = pd.Series(np.arange(len(genes)), index=genes)
        selected_genes = sorted(set(markers["gene"].astype(str).str.upper()) & set(genes))
        counts = released["counts"][lookup.loc[selected_genes].to_numpy(int)].toarray()
        expression = pd.DataFrame(
            np.log2(counts / released["total_counts"][None, :] * 1_000_000 + 1),
            index=selected_genes,
            columns=released["barcodes"],
        )
        for compartment in PROXY_ORDER:
            frozen = markers.loc[
                markers["compartment"].eq(compartment), "gene"
            ].astype(str).str.upper().tolist()
            observed = [gene for gene in frozen if gene in expression.index]
            standardized, omitted = standardize_rows(expression.loc[observed])
            if len(standardized) < 10:
                raise ValueError(
                    f"Only {len(standardized)} usable {compartment} proxy markers in {sample}"
                )
            column = f"proxy_{compartment.lower()}"
            current[column] = standardized.mean(axis=0).reindex(current["barcode"]).to_numpy()
            coverage_rows.append(
                {
                    "sample": sample,
                    "compartment": compartment,
                    "frozen_markers": 30,
                    "observed_markers": len(observed),
                    "variable_markers": len(standardized),
                    "omitted_zero_variance_markers": ";".join(omitted),
                }
            )

        covariate_names = [f"proxy_{name.lower()}" for name in PROXY_ORDER] + [
            "log1p_total_umi",
            "detected_genes",
        ]
        covariates = np.column_stack(
            [
                current[f"proxy_{name.lower()}"].to_numpy(float)
                for name in PROXY_ORDER
            ]
            + [
                np.log1p(current["nCount"].to_numpy(float)),
                current["nFeature"].to_numpy(float),
            ]
        )
        ranked_covariates = np.column_stack(
            [rank(covariates[:, index]) for index in range(covariates.shape[1])]
        )
        design = np.column_stack([np.ones(len(current)), ranked_covariates])
        matrix_rank = int(np.linalg.matrix_rank(design))
        if matrix_rank != design.shape[1]:
            raise ValueError(f"Rank-deficient adjustment model for {sample}")
        vifs = vif_values(ranked_covariates)
        residual_up = residualize(current["dynamic_up_set"].to_numpy(float), design)
        residual_context = residualize(
            current["myeloid_context"].to_numpy(float), design
        )
        current["residual_dynamic_up_set"] = residual_up
        current["residual_myeloid_context"] = residual_context
        adjusted_rho = float(stats.spearmanr(residual_up, residual_context).statistic)
        unadjusted_rho = float(
            stats.spearmanr(current["dynamic_up_set"], current["myeloid_context"]).statistic
        )

        rows = current["array_row"].to_numpy(int)
        cols = current["array_col"].to_numpy(int)
        shift_null, overlaps = toroidal_shift_null(
            current["dynamic_up_set"].to_numpy(float),
            current["myeloid_context"].to_numpy(float),
            rows,
            cols,
        )
        if len(shift_null) > 999:
            selected = rng.choice(len(shift_null), size=999, replace=False)
            shift_null = shift_null[selected]
            overlaps = overlaps[selected]
        empirical_p = float(
            (1 + (np.abs(shift_null) >= abs(unadjusted_rho)).sum())
            / (1 + len(shift_null))
        )
        patient_rows.append(
            {
                **metadata,
                "n_spots": len(current),
                "unadjusted_rho": unadjusted_rho,
                "composition_adjusted_rho": adjusted_rho,
                "design_rank": matrix_rank,
                "design_columns": design.shape[1],
                "maximum_vif": max(vifs),
                "n_accepted_toroidal_shifts": len(shift_null),
                "minimum_shift_overlap_spots": int(overlaps.min()),
                "toroidal_shift_p_two_sided": empirical_p,
                **{f"vif_{name}": value for name, value in zip(covariate_names, vifs)},
            }
        )
        null_frames.append(
            pd.DataFrame(
                {
                    "sample": sample,
                    "shift_index": np.arange(len(shift_null)),
                    "spearman_rho_null": shift_null,
                    "overlap_spots": overlaps,
                }
            )
        )
        spot_frames.append(current)
        print(
            f"ROBUSTNESS COMPLETE {sample} unadjusted={unadjusted_rho:.4f} adjusted={adjusted_rho:.4f}",
            flush=True,
        )

    patient = pd.DataFrame(patient_rows)
    nulls = pd.concat(null_frames, ignore_index=True)
    adjusted_summary = patient_summary(patient["composition_adjusted_rho"].to_numpy())
    unadjusted_summary = patient_summary(patient["unadjusted_rho"].to_numpy())
    adjusted_summary["passes_frozen_adjusted_support"] = bool(
        adjusted_summary["median_patient_rho"] > 0
        and adjusted_summary["positive_fraction"] >= 0.75
        and adjusted_summary["all_leave_one_out_medians_positive"]
        and adjusted_summary["exact_wilcoxon_p_two_sided"] < 0.05
    )

    null_by_sample = {
        sample: group["spearman_rho_null"].to_numpy(float)
        for sample, group in nulls.groupby("sample")
    }
    cohort_null = np.empty(999, dtype=float)
    for iteration in range(999):
        draw = [rng.choice(null_by_sample[item["sample"]]) for item in SAMPLES]
        cohort_null[iteration] = np.median(draw)
    observed_median = float(patient["unadjusted_rho"].median())
    cohort_p = float(
        (1 + (np.abs(cohort_null) >= abs(observed_median)).sum()) / 1000
    )

    pd.concat(spot_frames, ignore_index=True).to_csv(
        args.work_output_dir / "GSE206245_composition_adjusted_spot_scores.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    nulls.to_csv(
        args.work_output_dir / "GSE206245_toroidal_shift_nulls.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    pd.DataFrame({"median_rho_null": cohort_null}).to_csv(
        args.work_output_dir / "GSE206245_cohort_toroidal_shift_null.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    patient.to_csv(
        args.output_dir / "NPC_SecAct_GSE206245_composition_spatial_robustness_v0.1.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(coverage_rows).to_csv(
        args.output_dir / "NPC_SecAct_GSE206245_composition_proxy_coverage_v0.1.tsv",
        sep="\t",
        index=False,
    )
    summary = {
        "unadjusted": unadjusted_summary,
        "composition_adjusted": adjusted_summary,
        "cohort_toroidal_shift": {
            "observed_median_rho": observed_median,
            "n_null_replicates": 999,
            "empirical_p_two_sided": cohort_p,
            "null_median": float(np.median(cohort_null)),
            "null_ci95_lower": float(np.quantile(cohort_null, 0.025)),
            "null_ci95_upper": float(np.quantile(cohort_null, 0.975)),
        },
    }
    (args.output_dir / "NPC_SecAct_GSE206245_composition_spatial_robustness_v0.1.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
