from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from importlib.metadata import version
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from secactpy import load_signature, secact_activity_inference
from statsmodels.stats.multitest import multipletests

from analyze_gse206245_visium_spatial import (
    SAMPLES,
    enable_process_gsl_permutation_cache,
    read_slice,
)


SEED = 20260823
LAMBDA = 1e4
N_RAND = 1000


def exact_sign_flip_p(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    signs = np.asarray(list(itertools.product((-1.0, 1.0), repeat=len(values))))
    null = signs @ values / len(values)
    return float((np.abs(null) >= abs(values.mean()) - 1e-12).mean())


def paired_table(activity: pd.DataFrame, pairs: list[tuple[str, str]], modality: str,
                 scenario: str, cell_type: str) -> pd.DataFrame:
    rows = []
    correction = 1 - 3 / (4 * (len(pairs) - 1) - 1)
    for cytokine in activity.index:
        differences = np.asarray([
            float(activity.loc[cytokine, ot] - activity.loc[cytokine, bt])
            for bt, ot in pairs
        ])
        sd = differences.std(ddof=1)
        loo = np.asarray([np.delete(differences, i).mean() for i in range(len(differences))])
        rows.append({
            "modality": modality,
            "scenario": scenario,
            "cell_type": cell_type,
            "cytokine": cytokine,
            "n_patients": len(differences),
            "mean_difference_OT_minus_BT": float(differences.mean()),
            "median_difference_OT_minus_BT": float(np.median(differences)),
            "hedges_gz": float(correction * differences.mean() / sd) if sd else np.nan,
            "positive_patients": int((differences > 0).sum()),
            "negative_patients": int((differences < 0).sum()),
            "all_leave_one_out_same_as_mean": bool(np.all(np.sign(loo) == np.sign(differences.mean()))),
            "exact_sign_flip_p_two_sided": exact_sign_flip_p(differences),
        })
    result = pd.DataFrame(rows)
    result["fdr_bh_43"] = multipletests(result["exact_sign_flip_p_two_sided"], method="fdr_bh")[1]
    return result


def infer(expression: pd.DataFrame, signature: pd.DataFrame) -> pd.DataFrame:
    common = expression.index.intersection(signature.index)
    if len(common) / len(signature) < 0.8:
        raise ValueError(f"CytoSig coverage below 80%: {len(common)}/{len(signature)}")
    result = secact_activity_inference(
        expression.loc[common],
        is_differential=False,
        sig_matrix=signature,
        is_group_sig=False,
        lambda_=LAMBDA,
        n_rand=N_RAND,
        seed=SEED,
        rng_method="gsl",
        use_cache=True,
        sort_genes=True,
        verbose=False,
    )
    return result["zscore"]


def load_bulk(path: Path, signature: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    columns = [c for c in raw.columns if c.startswith("Patient_")]
    expr = raw.loc[raw["gene_symbol"].notna(), ["gene_symbol", *columns]].copy()
    expr["gene_symbol"] = expr["gene_symbol"].astype(str).str.upper()
    expr = expr.set_index("gene_symbol").apply(pd.to_numeric, errors="raise")
    expr = np.log2(expr.groupby(level=0, sort=False).median() + 1)
    activity = infer(expr, signature)
    patients = [12, 13, 19, 26, 28, 29, 33, 35, 39]
    pairs = [(f"Patient_{p}_BT{p}", f"Patient_{p}_OT{p}") for p in patients]
    return activity, paired_table(activity, pairs, "paired_bulk", "main", "all")


def load_scrna(path: Path, signature: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    analyses = [
        ("hard_qc", "B", "main"), ("hard_qc", "Myeloid", "main"),
        ("hard_qc", "T_NK", "main"),
        ("doublet_excluded", "B", "doublet_sensitivity"),
        ("doublet_excluded", "Myeloid", "doublet_sensitivity"),
        ("doublet_excluded", "T_NK", "doublet_sensitivity"),
        ("hard_qc", "PlasmaB", "cell_count_sensitivity"),
        ("doublet_excluded", "PlasmaB", "combined_sensitivity"),
    ]
    score_frames, result_frames = [], []
    for matrix_name, cell_type, scenario in analyses:
        expr = pd.read_csv(path / f"{matrix_name}_log2cpm.tsv.gz", sep="\t", index_col=0)
        expr.index = expr.index.astype(str).str.upper()
        columns = [f"Patient_{p}_{tp}{p}|{cell_type}" for p in range(2, 7) for tp in ("BT", "OT")]
        activity = infer(expr.loc[:, columns], signature)
        activity.index.name = "cytokine"
        score = activity.reset_index().melt(id_vars="cytokine", var_name="sample", value_name="zscore")
        score.insert(0, "cell_type", cell_type)
        score.insert(0, "scenario", scenario)
        score_frames.append(score)
        pairs = [(f"Patient_{p}_BT{p}|{cell_type}", f"Patient_{p}_OT{p}|{cell_type}") for p in range(2, 7)]
        result_frames.append(paired_table(activity, pairs, "scrna", scenario, cell_type))
    return pd.concat(score_frames, ignore_index=True), pd.concat(result_frames, ignore_index=True)


def summarize_spatial(patient: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for scenario in ("source_filtered", "fixed_spot_qc"):
        for cytokine in patient["cytokine"].drop_duplicates():
            values = patient.loc[(patient["scenario"] == scenario) & (patient["cytokine"] == cytokine), "spearman_rho"].to_numpy(float)
            loo = np.asarray([np.median(np.delete(values, i)) for i in range(len(values))])
            positive = int((values > 0).sum())
            rows.append({
                "scenario": scenario, "cytokine": cytokine, "n_patients": len(values),
                "median_patient_rho": float(np.median(values)), "mean_patient_rho": float(values.mean()),
                "positive_patients": positive, "minimum_leave_one_out_median_rho": float(loo.min()),
                "all_leave_one_out_medians_positive": bool(np.all(loo > 0)),
                "wilcoxon_p_two_sided": float(stats.wilcoxon(values, method="exact").pvalue),
                "exact_sign_p_two_sided": float(stats.binomtest(positive, len(values)).pvalue),
            })
    result = pd.DataFrame(rows)
    result["fdr_bh_43"] = np.nan
    for scenario in result["scenario"].unique():
        mask = result["scenario"] == scenario
        result.loc[mask, "fdr_bh_43"] = multipletests(result.loc[mask, "wilcoxon_p_two_sided"], method="fdr_bh")[1]
    main = result.loc[result["scenario"] == "source_filtered"].set_index("cytokine")
    sens = result.loc[result["scenario"] == "fixed_spot_qc"].set_index("cytokine")
    support = (
        (main["fdr_bh_43"] < 0.05) & (main["median_patient_rho"] > 0)
        & (main["positive_patients"] >= 9) & main["all_leave_one_out_medians_positive"]
        & (sens["median_patient_rho"] > 0) & (sens["positive_patients"] >= 9)
    )
    result["passes_frozen_spatial_support"] = result["cytokine"].map(support).fillna(False) & result["scenario"].eq("source_filtered")
    return result


def run_visium(data_dir: Path, prior_spots_path: Path, signature: pd.DataFrame,
                work_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    prior = pd.read_csv(prior_spots_path, sep="\t")
    patient_rows, score_frames, coverage_rows = [], [], []
    for meta in SAMPLES:
        sample = meta["sample"]
        checkpoint = work_dir / f"{sample}_cytosig_activity.tsv.gz"
        released = read_slice(data_dir, meta)
        genes = released["genes"]
        common = signature.index.intersection(genes)
        coverage = len(common) / len(signature)
        if coverage < 0.8:
            raise ValueError(f"{sample} CytoSig coverage below 80%")
        if checkpoint.exists():
            activity = pd.read_csv(checkpoint, sep="\t", index_col=0)
        else:
            lookup = pd.Series(np.arange(len(genes)), index=genes)
            counts = released["counts"][lookup.loc[common].to_numpy(int)].toarray().astype(float)
            expression = pd.DataFrame(
                np.log2(counts / released["total_counts"][None, :] * 1_000_000 + 1),
                index=common, columns=released["barcodes"],
            )
            activity = infer(expression, signature)
            activity.to_csv(checkpoint, sep="\t", compression="gzip")
        current = prior.loc[prior["sample"] == sample].set_index("barcode")
        if activity.columns.tolist() != current.index.tolist():
            activity = activity.loc[:, current.index]
        for scenario, mask in (
            ("source_filtered", np.ones(len(current), dtype=bool)),
            ("fixed_spot_qc", current["passes_fixed_spot_qc"].to_numpy(bool)),
        ):
            bars = current.index[mask]
            for cytokine in activity.index:
                rho = stats.spearmanr(activity.loc[cytokine, bars], current.loc[bars, "myeloid_context"]).statistic
                patient_rows.append({**meta, "scenario": scenario, "cytokine": cytokine,
                                     "n_spots": len(bars), "spearman_rho": float(rho)})
        frame = activity.T.reset_index(names="barcode")
        frame.insert(0, "sample", sample)
        score_frames.append(frame)
        coverage_rows.append({"sample": sample, "observed_signature_genes": len(common),
                              "total_signature_genes": len(signature), "coverage_fraction": coverage})
        print(f"CYTOSIG COMPLETE {sample}", flush=True)
    patient = pd.DataFrame(patient_rows)
    return pd.concat(score_frames, ignore_index=True), patient, pd.DataFrame(coverage_rows)


def plot_results(bulk: pd.DataFrame, scrna: pd.DataFrame, patient: pd.DataFrame,
                 summary: pd.DataFrame, output: Path) -> None:
    b = bulk.set_index("cytokine")["mean_difference_OT_minus_BT"]
    m = scrna.loc[(scrna["scenario"] == "main") & (scrna["cell_type"] == "Myeloid")].set_index("cytokine")["mean_difference_OT_minus_BT"]
    matrix = patient.loc[patient["scenario"] == "source_filtered"].pivot(index="cytokine", columns="sample", values="spearman_rho")
    order = summary.loc[summary["scenario"] == "source_filtered"].sort_values("median_patient_rho", ascending=False)["cytokine"]
    matrix = matrix.loc[order]
    fig, axes = plt.subplots(1, 2, figsize=(13, 9), gridspec_kw={"width_ratios": [1, 2.3]})
    axes[0].scatter(b, m, s=30, color="#3D6B7D")
    rho = stats.spearmanr(b, m).statistic
    axes[0].axhline(0, color="#999999", lw=.8); axes[0].axvline(0, color="#999999", lw=.8)
    axes[0].set_xlabel("Bulk mean OT-BT"); axes[0].set_ylabel("Myeloid mean OT-BT")
    axes[0].set_title(f"A  Paired-effect concordance\nSpearman rho={rho:.3f}", loc="left", fontweight="bold")
    im = axes[1].imshow(matrix.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-.8, vmax=.8)
    axes[1].set_xticks(range(matrix.shape[1]), matrix.columns, rotation=45, ha="right")
    axes[1].set_yticks(range(matrix.shape[0]), matrix.index, fontsize=7)
    axes[1].set_title("B  CytoSig-myeloid spatial correlations", loc="left", fontweight="bold")
    fig.colorbar(im, ax=axes[1], label="Within-patient Spearman rho", shrink=.6)
    fig.tight_layout(); fig.savefig(str(output) + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(str(output) + ".pdf", bbox_inches="tight"); plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("bulk", type=Path); p.add_argument("pseudobulk", type=Path)
    p.add_argument("visium_data", type=Path); p.add_argument("prior_spots", type=Path)
    p.add_argument("work_output", type=Path); p.add_argument("output", type=Path)
    args = p.parse_args(); args.work_output.mkdir(parents=True, exist_ok=True); args.output.mkdir(parents=True, exist_ok=True)
    enable_process_gsl_permutation_cache()
    signature = load_signature("cytosig"); signature.index = signature.index.astype(str).str.upper()
    bulk_scores, bulk_results = load_bulk(args.bulk, signature)
    scrna_scores, scrna_results = load_scrna(args.pseudobulk, signature)
    spot_scores, patient, coverage = run_visium(args.visium_data, args.prior_spots, signature, args.work_output)
    spatial_summary = summarize_spatial(patient)
    b = bulk_results.set_index("cytokine")["mean_difference_OT_minus_BT"]
    m = scrna_results.loc[(scrna_results["scenario"] == "main") & (scrna_results["cell_type"] == "Myeloid")].set_index("cytokine")["mean_difference_OT_minus_BT"]
    concordance = pd.DataFrame([{
        "n_cytokines": len(b), "spearman_rho": float(stats.spearmanr(b, m).statistic),
        "spearman_p": float(stats.spearmanr(b, m).pvalue),
        "same_direction_cytokines": int((np.sign(b) == np.sign(m)).sum()),
    }])
    bulk_scores.to_csv(args.work_output / "bulk_cytosig_zscores.tsv.gz", sep="\t", compression="gzip")
    scrna_scores.to_csv(args.work_output / "scrna_cytosig_zscores.tsv.gz", sep="\t", index=False, compression="gzip")
    spot_scores.to_csv(args.work_output / "visium_cytosig_zscores.tsv.gz", sep="\t", index=False, compression="gzip")
    bulk_results.to_csv(args.output / "NPC_SecAct_CytoSig_bulk_results_v0.1.tsv", sep="\t", index=False)
    scrna_results.to_csv(args.output / "NPC_SecAct_CytoSig_scrna_results_v0.1.tsv", sep="\t", index=False)
    patient.to_csv(args.output / "NPC_SecAct_CytoSig_visium_patient_correlations_v0.1.tsv", sep="\t", index=False)
    spatial_summary.to_csv(args.output / "NPC_SecAct_CytoSig_visium_summary_v0.1.tsv", sep="\t", index=False)
    coverage.to_csv(args.output / "NPC_SecAct_CytoSig_coverage_v0.1.tsv", sep="\t", index=False)
    concordance.to_csv(args.output / "NPC_SecAct_CytoSig_concordance_v0.1.tsv", sep="\t", index=False)
    plot_results(bulk_results, scrna_results, patient, spatial_summary,
                 args.output / "figures" / "NPC_SecAct_CytoSig_benchmark_v0.1")
    run = {"secactpy_version": version("secactpy"), "seed": SEED, "lambda": LAMBDA,
           "n_rand": N_RAND, "rng_method": "gsl", "signature_shape": list(signature.shape),
           "signature_sha256": hashlib.sha256(Path("work/cytosig_audit/official_signature.centroid").read_bytes()).hexdigest(),
           "formal_spatial_support_count": int(spatial_summary["passes_frozen_spatial_support"].sum())}
    (args.work_output / "cytosig_run.json").write_text(json.dumps(run, indent=2), encoding="utf-8")
    print(json.dumps(run, indent=2), flush=True)


if __name__ == "__main__":
    main()
