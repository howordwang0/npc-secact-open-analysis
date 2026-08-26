from __future__ import annotations

import argparse
import itertools
import json
import re
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from secactpy import load_signature, secact_activity_inference
from secactpy.inference import group_signatures
from statsmodels.stats.multitest import multipletests


SEED = 20260823


def parse_sample(column: str) -> tuple[int, str]:
    match = re.fullmatch(r"Patient_(\d+)_(BT|OT)\d+", column)
    if not match:
        raise ValueError(f"Unrecognized sample column: {column}")
    return int(match.group(1)), match.group(2)


def group_level_matrix(expanded_activity: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    grouped_signature = group_signatures(load_signature("secact"), cor_threshold=0.9)
    activity_rows = []
    member_rows = []
    for group_name in grouped_signature.columns:
        members = group_name.split("|")
        representative = members[0]
        activity_rows.append(expanded_activity.loc[representative].rename(group_name))
        for member in members:
            member_rows.append(
                {"signature_group": group_name, "protein": member, "representative": representative}
            )
    return pd.DataFrame(activity_rows), pd.DataFrame(member_rows)


def exact_sign_flip_p(differences: np.ndarray) -> np.ndarray:
    signs = np.asarray(list(itertools.product((-1.0, 1.0), repeat=differences.shape[1])))
    permuted = signs @ differences.T / differences.shape[1]
    observed = np.abs(differences.mean(axis=1))
    return (np.abs(permuted) >= observed[None, :] - 1e-12).mean(axis=0)


def paired_statistics(activity: pd.DataFrame, patient_order: list[int], bootstrap: int) -> pd.DataFrame:
    bt_columns = [f"Patient_{patient}_BT{patient}" for patient in patient_order]
    ot_columns = [f"Patient_{patient}_OT{patient}" for patient in patient_order]
    bt = activity.loc[:, bt_columns].to_numpy(dtype=float)
    ot = activity.loc[:, ot_columns].to_numpy(dtype=float)
    differences = ot - bt
    n = differences.shape[1]
    means = differences.mean(axis=1)
    medians = np.median(differences, axis=1)
    sds = differences.std(axis=1, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        cohen_dz = means / sds
    correction = 1 - 3 / (4 * (n - 1) - 1)
    hedges_gz = correction * cohen_dz

    paired_t = stats.ttest_rel(ot, bt, axis=1, nan_policy="raise")
    wilcoxon_p = []
    for row in differences:
        try:
            wilcoxon_p.append(float(stats.wilcoxon(row, alternative="two-sided", method="auto").pvalue))
        except ValueError:
            wilcoxon_p.append(1.0)

    sign_flip_p = exact_sign_flip_p(differences)
    sign_flip_fdr = multipletests(sign_flip_p, method="fdr_bh")[1]

    rng = np.random.default_rng(SEED)
    sampled_indexes = rng.integers(0, n, size=(bootstrap, n))
    boot_means = differences[:, sampled_indexes].mean(axis=2).T
    observed_sign = np.sign(means)
    same_direction = (np.sign(boot_means) == observed_sign[None, :]).mean(axis=0)

    loo_means = np.stack(
        [np.delete(differences, index, axis=1).mean(axis=1) for index in range(n)], axis=0
    )
    loo_same_direction = (np.sign(loo_means) == observed_sign[None, :]).mean(axis=0)

    result = pd.DataFrame(
        {
            "signature_group": activity.index,
            "n_patients": n,
            "mean_paired_difference_OT_minus_BT": means,
            "median_paired_difference_OT_minus_BT": medians,
            "cohen_dz": cohen_dz,
            "hedges_gz": hedges_gz,
            "paired_t": paired_t.statistic,
            "paired_t_p": paired_t.pvalue,
            "wilcoxon_p": wilcoxon_p,
            "exact_sign_flip_p": sign_flip_p,
            "exact_sign_flip_fdr_bh": sign_flip_fdr,
            "bootstrap_mean_ci_low": np.quantile(boot_means, 0.025, axis=0),
            "bootstrap_mean_ci_high": np.quantile(boot_means, 0.975, axis=0),
            "bootstrap_same_direction_fraction": same_direction,
            "leave_one_patient_out_same_direction_fraction": loo_same_direction,
        }
    )
    result["exploratory_dynamic_group"] = (
        (result["exact_sign_flip_fdr_bh"] < 0.10)
        & (result["hedges_gz"].abs() >= 1.0)
        & (result["bootstrap_same_direction_fraction"] >= 0.90)
    )
    return result.sort_values(
        ["exact_sign_flip_fdr_bh", "exact_sign_flip_p", "hedges_gz"],
        ascending=[True, True, False],
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the frozen Zenodo 13736891 paired SecAct analysis.")
    parser.add_argument("expression_tsv", type=Path)
    parser.add_argument("work_output_dir", type=Path)
    parser.add_argument("result_tsv", type=Path)
    parser.add_argument("--n-rand", type=int, default=1000)
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--rng-method", choices=["srand", "gsl", "numpy"], default="gsl")
    args = parser.parse_args()

    raw = pd.read_csv(args.expression_tsv, sep="\t")
    sample_columns = raw.columns[3:].tolist()
    parsed = {column: parse_sample(column) for column in sample_columns}
    patients = sorted({patient for patient, _ in parsed.values()})
    if len(patients) != 9:
        raise ValueError(f"Expected 9 patients, found {patients}")
    for patient in patients:
        observed = sorted(timepoint for current, timepoint in parsed.values() if current == patient)
        if observed != ["BT", "OT"]:
            raise ValueError(f"Patient {patient} does not have exactly one BT and one OT sample: {observed}")

    signature = load_signature("secact")
    signature_genes = set(signature.index.astype(str))
    expression = raw.loc[raw["gene_symbol"].notna(), ["gene_symbol", *sample_columns]].copy()
    expression["gene_symbol"] = expression["gene_symbol"].astype(str)
    expression = expression.loc[expression["gene_symbol"].isin(signature_genes)].set_index("gene_symbol")
    expression = expression.apply(pd.to_numeric, errors="raise")
    expression = expression.groupby(level=0, sort=False).mean()
    if expression.isna().any().any():
        raise ValueError("Missing expression values found")
    if (expression.to_numpy() < 0).any():
        raise ValueError("Negative expression values found")
    expression = np.log2(expression + 1.0)

    canonical_order = [
        column
        for patient in patients
        for timepoint in ("BT", "OT")
        for column, parsed_value in parsed.items()
        if parsed_value == (patient, timepoint)
    ]
    expression = expression.loc[:, canonical_order]

    secact = secact_activity_inference(
        expression,
        is_differential=False,
        sig_matrix="secact",
        n_rand=args.n_rand,
        seed=SEED,
        rng_method=args.rng_method,
        sort_genes=True,
        verbose=True,
    )
    group_activity, group_members = group_level_matrix(secact["zscore"])
    results = paired_statistics(group_activity, patients, args.bootstrap)

    args.work_output_dir.mkdir(parents=True, exist_ok=True)
    args.result_tsv.parent.mkdir(parents=True, exist_ok=True)
    expression.to_csv(args.work_output_dir / "Zenodo13736891_secact_input_log2.tsv.gz", sep="\t", compression="gzip")
    group_activity.to_csv(
        args.work_output_dir / "Zenodo13736891_secact_group_activity.tsv.gz", sep="\t", compression="gzip"
    )
    group_members.to_csv(args.work_output_dir / "Zenodo13736891_secact_group_members.tsv", sep="\t", index=False)
    for key in ("beta", "se", "zscore", "pvalue"):
        secact[key].to_csv(
            args.work_output_dir / f"Zenodo13736891_secact_{key}.tsv.gz", sep="\t", compression="gzip"
        )
    results.to_csv(args.result_tsv, sep="\t", index=False)

    summary = {
        "dataset": "Zenodo 13736891",
        "doi": "10.5281/zenodo.13736891",
        "patients": len(patients),
        "paired_samples": len(sample_columns),
        "patient_ids": patients,
        "input_signature_genes": int(expression.shape[0]),
        "secactpy_version": version("secactpy"),
        "n_rand": args.n_rand,
        "rng_method": args.rng_method,
        "seed": SEED,
        "independent_signature_groups_tested": int(group_activity.shape[0]),
        "bootstrap_iterations": args.bootstrap,
        "exact_sign_flip_permutations": 2 ** len(patients),
        "minimum_two_sided_exact_p_resolution": 2 / (2 ** len(patients)),
        "exploratory_dynamic_groups": int(results["exploratory_dynamic_group"].sum()),
        "minimum_exact_sign_flip_p": float(results["exact_sign_flip_p"].min()),
        "minimum_exact_sign_flip_fdr_bh": float(results["exact_sign_flip_fdr_bh"].min()),
        "response_labels_used": False,
    }
    (args.work_output_dir / "Zenodo13736891_paired_run.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
