from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from analyze_gse206245_visium_spatial import SAMPLES


SEED = 20260823
N_RANDOM = 1000


def rank_bins(values: pd.Series) -> pd.Series:
    ranks = values.rank(method="first")
    return pd.qcut(ranks, q=10, labels=False, duplicates="drop").astype(int)


def standardize_rows(frame: pd.DataFrame) -> pd.DataFrame:
    sd = frame.std(axis=1, ddof=1)
    if (sd <= 0).any() or sd.isna().any():
        raise ValueError("Invariant full SecAct group activity")
    return frame.sub(frame.mean(axis=1), axis=0).div(sd, axis=0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("full_group_dir", type=Path)
    parser.add_argument("prior_spot_scores", type=Path)
    parser.add_argument("dynamic_sets", type=Path)
    parser.add_argument("work_output_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    args.work_output_dir.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    prior = pd.read_csv(args.prior_spot_scores, sep="\t")
    dynamic = pd.read_csv(args.dynamic_sets, sep="\t")
    up_groups = dynamic.loc[
        dynamic["set_label"].eq("dynamic_up_set"), "signature_group"
    ].astype(str).tolist()
    down_groups = set(
        dynamic.loc[
            dynamic["set_label"].eq("dynamic_down_set"), "signature_group"
        ].astype(str)
    )
    if len(up_groups) != 70 or len(down_groups) != 93:
        raise ValueError("Frozen dynamic group dimensions changed")

    activities = {}
    standardized = {}
    characteristics = []
    observed_rhos = []
    for metadata in SAMPLES:
        sample = metadata["sample"]
        path = args.full_group_dir / f"{sample}_all_group_activity.tsv.gz"
        activity = pd.read_csv(path, sep="\t", index_col=0)
        if activity.shape != (657, metadata["expected_spots"]):
            raise ValueError(f"Invalid full group activity for {sample}: {activity.shape}")
        current = prior.loc[prior["sample"].eq(sample)]
        if activity.columns.tolist() != current["barcode"].tolist():
            raise ValueError(f"Barcode mismatch for {sample}")
        activities[sample] = activity
        standardized[sample] = standardize_rows(activity)
        means = activity.mean(axis=1)
        log_sds = np.log(activity.std(axis=1, ddof=1))
        characteristics.append(
            pd.DataFrame(
                {
                    "group": activity.index,
                    "sample": sample,
                    "mean_activity": means.to_numpy(),
                    "log_sd_activity": log_sds.to_numpy(),
                }
            )
        )
        observed_score = standardized[sample].loc[up_groups].mean(axis=0)
        observed_rhos.append(
            float(stats.spearmanr(observed_score, current["myeloid_context"]).statistic)
        )

    characteristic_long = pd.concat(characteristics, ignore_index=True)
    pooled = (
        characteristic_long.groupby("group", observed=True)[
            ["mean_activity", "log_sd_activity"]
        ]
        .median()
        .rename(
            columns={
                "mean_activity": "median_sample_mean_activity",
                "log_sd_activity": "median_sample_log_sd_activity",
            }
        )
    )
    pooled["mean_decile"] = rank_bins(pooled["median_sample_mean_activity"])
    pooled["sd_decile"] = rank_bins(pooled["median_sample_log_sd_activity"])
    excluded = set(up_groups) | down_groups
    candidates = sorted(set(pooled.index) - excluded)
    if len(candidates) != 494:
        raise ValueError(f"Expected 494 non-dynamic candidate groups, found {len(candidates)}")

    rng = np.random.default_rng(SEED)
    random_rows = []
    membership_rows = []
    contexts = {
        metadata["sample"]: prior.loc[
            prior["sample"].eq(metadata["sample"]), "myeloid_context"
        ].to_numpy(float)
        for metadata in SAMPLES
    }
    for random_index in range(N_RANDOM):
        available = set(candidates)
        selected = []
        up_order = rng.permutation(up_groups)
        for up_group in up_order:
            target_mean = int(pooled.loc[up_group, "mean_decile"])
            target_sd = int(pooled.loc[up_group, "sd_decile"])
            available_list = sorted(available)
            distances = np.asarray(
                [
                    abs(int(pooled.loc[group, "mean_decile"]) - target_mean)
                    + abs(int(pooled.loc[group, "sd_decile"]) - target_sd)
                    for group in available_list
                ]
            )
            closest = [
                group
                for group, distance in zip(available_list, distances)
                if distance == distances.min()
            ]
            chosen = str(rng.choice(closest))
            selected.append(chosen)
            available.remove(chosen)
            membership_rows.append(
                {
                    "random_program": random_index + 1,
                    "frozen_up_group": up_group,
                    "matched_candidate_group": chosen,
                    "decile_manhattan_distance": int(distances.min()),
                }
            )
        if len(selected) != 70 or len(set(selected)) != 70:
            raise ValueError("Random program does not contain 70 distinct groups")
        patient_rhos = []
        positive = 0
        for metadata in SAMPLES:
            sample = metadata["sample"]
            random_score = standardized[sample].loc[selected].mean(axis=0)
            rho = float(stats.spearmanr(random_score, contexts[sample]).statistic)
            patient_rhos.append(rho)
            positive += int(rho > 0)
        random_rows.append(
            {
                "random_program": random_index + 1,
                "median_patient_rho": float(np.median(patient_rhos)),
                "mean_patient_rho": float(np.mean(patient_rhos)),
                "positive_patients": positive,
                "meets_positive_median_75pct_rule": bool(
                    np.median(patient_rhos) > 0 and positive >= 9
                ),
                **{
                    f"rho_{metadata['sample']}": rho
                    for metadata, rho in zip(SAMPLES, patient_rhos)
                },
            }
        )
        if (random_index + 1) % 100 == 0:
            print(f"RANDOM PROGRAMS {random_index + 1}/{N_RANDOM}", flush=True)

    random_table = pd.DataFrame(random_rows)
    observed_median = float(np.median(observed_rhos))
    empirical_p = float(
        (1 + (random_table["median_patient_rho"] >= observed_median).sum())
        / (N_RANDOM + 1)
    )
    percentile = float(
        100 * (random_table["median_patient_rho"] < observed_median).mean()
    )
    summary = {
        "observed_program": {
            "median_patient_rho": observed_median,
            "mean_patient_rho": float(np.mean(observed_rhos)),
            "positive_patients": int((np.asarray(observed_rhos) > 0).sum()),
            "patient_rhos": {
                metadata["sample"]: rho
                for metadata, rho in zip(SAMPLES, observed_rhos)
            },
        },
        "matched_random_programs": {
            "n_programs": N_RANDOM,
            "candidate_groups": len(candidates),
            "program_size": 70,
            "null_median": float(random_table["median_patient_rho"].median()),
            "null_ci95_lower": float(
                random_table["median_patient_rho"].quantile(0.025)
            ),
            "null_ci95_upper": float(
                random_table["median_patient_rho"].quantile(0.975)
            ),
            "observed_percentile": percentile,
            "empirical_p_one_sided": empirical_p,
            "random_programs_meeting_positive_median_75pct_rule": int(
                random_table["meets_positive_median_75pct_rule"].sum()
            ),
        },
        "passes_frozen_specificity_rule": bool(empirical_p < 0.05),
    }
    characteristic_long.to_csv(
        args.work_output_dir / "GSE206245_group_activity_characteristics.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    pooled.reset_index().to_csv(
        args.output_dir / "NPC_SecAct_GSE206245_group_matching_characteristics_v0.1.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(membership_rows).to_csv(
        args.work_output_dir / "GSE206245_matched_random_program_memberships.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    random_table.to_csv(
        args.output_dir / "NPC_SecAct_GSE206245_matched_random_programs_v0.1.tsv",
        sep="\t",
        index=False,
    )
    (args.output_dir / "NPC_SecAct_GSE206245_matched_random_summary_v0.1.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
