from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib
import itertools
import json
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.io import mmread
from secactpy import load_signature, secact_activity_inference
from secactpy.inference import group_signatures


SEED = 20260823
SAMPLES = [
    {"sample": "ST5", "gsm": "GSM6248646", "stage": "IVa", "response": "CR", "expected_spots": 2021},
    {"sample": "ST6", "gsm": "GSM6248647", "stage": "IVb", "response": "CR", "expected_spots": 2614},
    {"sample": "ST7", "gsm": "GSM6248648", "stage": "IVa", "response": "PD", "expected_spots": 1939},
    {"sample": "ST8", "gsm": "GSM6248649", "stage": "IVb", "response": "CR", "expected_spots": 1209},
    {"sample": "ST9", "gsm": "GSM6248650", "stage": "III", "response": "PD", "expected_spots": 2602},
    {"sample": "ST10", "gsm": "GSM6248651", "stage": "IVa", "response": "PR", "expected_spots": 1309},
    {"sample": "ST11", "gsm": "GSM6248652", "stage": "IVb", "response": "CR", "expected_spots": 1272},
    {"sample": "ST12", "gsm": "GSM6248653", "stage": "IVa", "response": "CR", "expected_spots": 3668},
    {"sample": "ST16", "gsm": "GSM6248654", "stage": "IVb", "response": "PD", "expected_spots": 1556},
    {"sample": "ST17", "gsm": "GSM6248655", "stage": "IVb", "response": "SD", "expected_spots": 2067},
    {"sample": "ST18", "gsm": "GSM6248656", "stage": "III", "response": "PR", "expected_spots": 1977},
    {"sample": "ST19", "gsm": "GSM6248657", "stage": "IVb", "response": "PD", "expected_spots": 1086},
]


def released_path(data_dir: Path, gsm: str, sample: str, suffix: str) -> Path:
    path = data_dir / f"{gsm}_NPC_{sample}_{suffix}"
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def read_positions(path: Path) -> pd.DataFrame:
    columns = [
        "barcode",
        "in_tissue",
        "array_row",
        "array_col",
        "pxl_col_in_fullres",
        "pxl_row_in_fullres",
    ]
    table = pd.read_csv(path, header=None, names=columns, compression="gzip")
    if str(table.iloc[0, 0]).lower() == "barcode":
        table = table.iloc[1:].copy()
    table["barcode"] = table["barcode"].astype(str)
    for column in columns[1:]:
        table[column] = pd.to_numeric(table[column], errors="raise")
    if table["barcode"].duplicated().any():
        raise ValueError(f"Duplicate position barcodes in {path}")
    return table.set_index("barcode")


def read_slice(data_dir: Path, metadata: dict) -> dict:
    gsm = metadata["gsm"]
    sample = metadata["sample"]
    features_path = released_path(data_dir, gsm, sample, "features.tsv.gz")
    barcodes_path = released_path(data_dir, gsm, sample, "barcodes.tsv.gz")
    matrix_path = released_path(data_dir, gsm, sample, "matrix.mtx.gz")
    positions_path = released_path(
        data_dir, gsm, sample, "tissue_positions_list.csv.gz"
    )

    features = pd.read_csv(
        features_path, sep="\t", header=None, dtype=str, compression="gzip"
    )
    barcodes = pd.read_csv(
        barcodes_path, sep="\t", header=None, dtype=str, compression="gzip"
    ).iloc[:, 0]
    with gzip.open(matrix_path, "rb") as handle:
        counts = mmread(handle).tocsr()
    if counts.shape != (len(features), len(barcodes)):
        raise ValueError(
            f"Matrix dimensions do not match features/barcodes for {sample}: "
            f"{counts.shape} versus {(len(features), len(barcodes))}"
        )

    if features.shape[1] >= 3:
        gene_expression = features.iloc[:, 2].eq("Gene Expression").to_numpy()
        if gene_expression.any():
            features = features.loc[gene_expression].reset_index(drop=True)
            counts = counts[gene_expression]
    symbols = features.iloc[:, 1].fillna(features.iloc[:, 0]).astype(str).str.upper()
    valid = symbols.str.len().gt(0).to_numpy()
    symbols = symbols.loc[valid].reset_index(drop=True)
    counts = counts[valid]

    unique_symbols = pd.Index(pd.unique(symbols), dtype=str)
    symbol_codes = pd.Categorical(symbols, categories=unique_symbols).codes
    collapse = sparse.csr_matrix(
        (
            np.ones(len(symbol_codes), dtype=np.float32),
            (symbol_codes, np.arange(len(symbol_codes))),
        ),
        shape=(len(unique_symbols), len(symbol_codes)),
    )
    counts = (collapse @ counts).tocsr()

    positions = read_positions(positions_path)
    barcode_values = barcodes.astype(str).to_numpy()
    in_tissue = np.array(
        [
            barcode in positions.index and int(positions.loc[barcode, "in_tissue"]) == 1
            for barcode in barcode_values
        ],
        dtype=bool,
    )
    if not in_tissue.any():
        raise ValueError(f"No in-tissue spots mapped for {sample}")
    counts = counts[:, in_tissue]
    barcode_values = barcode_values[in_tissue]
    positions = positions.loc[barcode_values].copy()

    total_counts = np.asarray(counts.sum(axis=0)).ravel().astype(float)
    detected_features = np.asarray(counts.getnnz(axis=0)).ravel().astype(int)
    if (total_counts <= 0).any():
        raise ValueError(f"Zero-library in-tissue spot in {sample}")
    return {
        "genes": unique_symbols,
        "counts": counts,
        "barcodes": barcode_values,
        "positions": positions,
        "total_counts": total_counts,
        "detected_features": detected_features,
    }


def group_level_matrix(expanded_activity: pd.DataFrame) -> pd.DataFrame:
    grouped = group_signatures(load_signature("secact"), cor_threshold=0.9)
    rows = []
    for group_name in grouped.columns:
        representative = group_name.split("|")[0]
        rows.append(expanded_activity.loc[representative].rename(group_name))
    return pd.DataFrame(rows)


def predictor_genes(grouped_signature: pd.DataFrame, group_name: str) -> set[str]:
    values = pd.to_numeric(grouped_signature[group_name], errors="coerce")
    mask = values.notna() & values.ne(0)
    return set(grouped_signature.index[mask].astype(str).str.upper())


def signature_digest(signature: pd.DataFrame) -> str:
    digest = hashlib.sha256()
    digest.update("\0".join(signature.index.astype(str)).encode("utf-8"))
    digest.update("\0".join(signature.columns.astype(str)).encode("utf-8"))
    digest.update(np.ascontiguousarray(signature.to_numpy(np.float64)).tobytes())
    return digest.hexdigest()


def enable_process_gsl_permutation_cache() -> None:
    ridge_module = importlib.import_module("secactpy.ridge")
    rng_module = importlib.import_module("secactpy.rng")
    cache: dict[tuple[int, int, int], np.ndarray] = {}

    def cached_gsl_table(
        n_genes: int, n_permutations: int, seed: int, verbose: bool = False
    ) -> np.ndarray:
        key = (n_genes, n_permutations, seed)
        if key not in cache:
            if verbose:
                print(f"  Generating process GSL permutation cache: {key}")
            cache[key] = rng_module.GSLRNG(seed).inverse_permutation_table(
                n_genes, n_permutations
            )
        elif verbose:
            print(f"  Reusing process GSL permutation cache: {key}")
        return cache[key]

    ridge_module.get_cached_inverse_perm_table = cached_gsl_table


def standardize_rows(values: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    row_sd = values.std(axis=1, ddof=1)
    retained = row_sd.index[row_sd.gt(0) & row_sd.notna()]
    omitted = sorted(set(values.index) - set(retained))
    standardized = values.loc[retained].sub(
        values.loc[retained].mean(axis=1), axis=0
    ).div(row_sd.loc[retained], axis=0)
    return standardized, omitted


def score_spots(
    group_activity: pd.DataFrame,
    marker_expression: pd.DataFrame,
    set_members: dict[str, set[str]],
    spot_mask: np.ndarray,
) -> tuple[pd.DataFrame, dict]:
    columns = group_activity.columns[spot_mask]
    standardized_activity, omitted_groups = standardize_rows(
        group_activity.loc[:, columns]
    )
    standardized_markers, omitted_markers = standardize_rows(
        marker_expression.loc[:, columns]
    )
    scores = pd.DataFrame(index=columns)
    retained_counts = {}
    for set_label, members in set_members.items():
        retained = sorted(set(standardized_activity.index) & members)
        retained_counts[set_label] = len(retained)
        if set_label == "dynamic_up_set" and len(retained) < int(np.ceil(0.9 * 70)):
            raise ValueError(
                f"Only {len(retained)} dynamic-up groups retained after variance filter"
            )
        scores[set_label] = standardized_activity.loc[retained].mean(axis=0)
    if standardized_markers.shape[0] < 16:
        raise ValueError(
            f"Only {standardized_markers.shape[0]} myeloid markers retained"
        )
    scores["myeloid_context"] = standardized_markers.mean(axis=0)
    details = {
        "retained_group_counts": retained_counts,
        "omitted_zero_variance_groups": omitted_groups,
        "retained_myeloid_markers": list(standardized_markers.index),
        "omitted_zero_variance_markers": omitted_markers,
    }
    return scores, details


def patient_correlations(scores: pd.DataFrame, set_label: str) -> float:
    result = stats.spearmanr(
        scores[set_label].to_numpy(float),
        scores["myeloid_context"].to_numpy(float),
    )
    return float(result.statistic)


def summarize_correlations(
    patient_table: pd.DataFrame, scenario: str, set_label: str
) -> dict:
    values = patient_table.loc[
        (patient_table["scenario"] == scenario)
        & (patient_table["set_label"] == set_label),
        "spearman_rho",
    ].to_numpy(float)
    n = len(values)
    positive = int((values > 0).sum())
    loo_medians = np.array([np.median(np.delete(values, idx)) for idx in range(n)])
    wilcoxon_greater = stats.wilcoxon(
        values, alternative="greater", zero_method="wilcox", method="auto"
    )
    wilcoxon_two_sided = stats.wilcoxon(
        values, alternative="two-sided", zero_method="wilcox", method="auto"
    )
    return {
        "scenario": scenario,
        "set_label": set_label,
        "n_patients": n,
        "median_patient_rho": float(np.median(values)),
        "mean_patient_rho": float(np.mean(values)),
        "min_patient_rho": float(np.min(values)),
        "max_patient_rho": float(np.max(values)),
        "positive_patients": positive,
        "positive_fraction": positive / n,
        "minimum_leave_one_out_median_rho": float(loo_medians.min()),
        "all_leave_one_out_medians_positive": bool((loo_medians > 0).all()),
        "exact_sign_p_greater": float(
            stats.binomtest(positive, n, p=0.5, alternative="greater").pvalue
        ),
        "exact_sign_p_two_sided": float(
            stats.binomtest(positive, n, p=0.5, alternative="two-sided").pvalue
        ),
        "wilcoxon_p_greater": float(wilcoxon_greater.pvalue),
        "wilcoxon_p_two_sided": float(wilcoxon_two_sided.pvalue),
    }


def exact_rank_test(values: np.ndarray, labels: np.ndarray) -> dict:
    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels, dtype=bool)
    ranks = stats.rankdata(values, method="average")
    n_positive = int(labels.sum())
    observed = abs(ranks[labels].mean() - ranks[~labels].mean())
    permutations = []
    for selected in itertools.combinations(range(len(values)), n_positive):
        current = np.zeros(len(values), dtype=bool)
        current[list(selected)] = True
        permutations.append(abs(ranks[current].mean() - ranks[~current].mean()))
    p_value = float(
        np.mean(np.asarray(permutations) >= observed - np.finfo(float).eps * 16)
    )
    u_statistic = float(stats.mannwhitneyu(values[labels], values[~labels]).statistic)
    rank_biserial = 2 * u_statistic / (n_positive * (len(values) - n_positive)) - 1
    return {
        "n_responders": n_positive,
        "n_nonresponders": int(len(values) - n_positive),
        "responder_mean": float(values[labels].mean()),
        "nonresponder_mean": float(values[~labels].mean()),
        "responder_median": float(np.median(values[labels])),
        "nonresponder_median": float(np.median(values[~labels])),
        "rank_biserial": float(rank_biserial),
        "exact_rank_permutation_p_two_sided": p_value,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("data_dir", type=Path)
    parser.add_argument("dynamic_sets_tsv", type=Path)
    parser.add_argument("myeloid_markers_tsv", type=Path)
    parser.add_argument("work_output_dir", type=Path)
    parser.add_argument("spot_scores_output", type=Path)
    parser.add_argument("patient_correlations_output", type=Path)
    parser.add_argument("summary_output", type=Path)
    parser.add_argument("--n-rand", type=int, default=1000)
    parser.add_argument(
        "--rng-method", choices=["srand", "gsl", "numpy"], default="gsl"
    )
    parser.add_argument(
        "--checkpoint-samples",
        nargs="+",
        choices=[sample["sample"] for sample in SAMPLES],
        help="Generate validated SecAct checkpoints for only these samples, then exit",
    )
    args = parser.parse_args()
    use_permutation_cache = args.rng_method == "gsl"
    if use_permutation_cache:
        enable_process_gsl_permutation_cache()
    args.work_output_dir.mkdir(parents=True, exist_ok=True)
    args.spot_scores_output.parent.mkdir(parents=True, exist_ok=True)

    dynamic_sets = pd.read_csv(args.dynamic_sets_tsv, sep="\t")
    set_members = {
        label: set(group["signature_group"].astype(str))
        for label, group in dynamic_sets.groupby("set_label")
    }
    if {key: len(value) for key, value in set_members.items()} != {
        "dynamic_down_set": 93,
        "dynamic_up_set": 70,
    }:
        raise ValueError("Frozen dynamic set sizes changed")
    marker_table = pd.read_csv(args.myeloid_markers_tsv, sep="\t")
    myeloid_markers = marker_table["gene"].astype(str).str.upper().tolist()
    if len(myeloid_markers) != 20 or len(set(myeloid_markers)) != 20:
        raise ValueError("Frozen myeloid marker set must contain 20 unique genes")

    raw_signature = load_signature("secact")
    inference_grouped_signature = group_signatures(raw_signature, cor_threshold=0.9)
    inference_signature_digest = signature_digest(inference_grouped_signature)
    signature = raw_signature.copy()
    signature.index = signature.index.astype(str).str.upper()
    grouped_signature = inference_grouped_signature.copy()
    grouped_signature.index = grouped_signature.index.astype(str).str.upper()
    missing_groups = set().union(*set_members.values()) - set(grouped_signature.columns)
    if missing_groups:
        raise ValueError(f"Frozen groups absent from SecAct signature: {missing_groups}")
    dynamic_up_predictors = set().union(
        *(predictor_genes(grouped_signature, group) for group in set_members["dynamic_up_set"])
    )
    independent_marker_candidates = sorted(
        set(myeloid_markers) - dynamic_up_predictors
    )
    if len(independent_marker_candidates) != 20:
        raise ValueError("Frozen myeloid marker set overlaps SecAct predictor genes")

    spot_frames = []
    patient_rows = []
    coverage_rows = []
    run_samples = []
    aggregate_expression = []
    aggregate_sample_names = []
    samples_to_run = SAMPLES
    if args.checkpoint_samples:
        if len(set(args.checkpoint_samples)) != len(args.checkpoint_samples):
            raise ValueError("Duplicate checkpoint sample requested")
        requested = set(args.checkpoint_samples)
        samples_to_run = [
            metadata for metadata in SAMPLES if metadata["sample"] in requested
        ]

    for metadata in samples_to_run:
        sample = metadata["sample"]
        print(f"LOAD START {sample}", flush=True)
        released = read_slice(args.data_dir, metadata)
        if len(released["barcodes"]) != metadata["expected_spots"]:
            raise ValueError(
                f"{sample} spot count differs from Supplementary Data 1: "
                f"{len(released['barcodes'])} != {metadata['expected_spots']}"
            )
        genes = released["genes"]
        gene_lookup = pd.Series(np.arange(len(genes)), index=genes)
        signature_genes = sorted(set(signature.index) & set(genes))
        coverage_fraction = len(signature_genes) / len(signature.index)
        observed_group_counts = {
            group: len(predictor_genes(grouped_signature, group) & set(genes))
            for group in set_members["dynamic_up_set"]
        }
        unrepresented_groups = sorted(
            group for group, count in observed_group_counts.items() if count == 0
        )
        if coverage_fraction < 0.8 or unrepresented_groups:
            raise ValueError(
                f"Full SecAct coverage gate failed for {sample}: "
                f"coverage={coverage_fraction:.3f}, missing groups={unrepresented_groups}"
            )
        independent_markers = sorted(
            set(independent_marker_candidates) & set(genes)
        )
        if len(independent_markers) < 16:
            raise ValueError(
                f"Only {len(independent_markers)} independent myeloid markers in {sample}"
            )

        needed_genes = sorted(set(signature_genes) | set(independent_markers))
        row_indices = gene_lookup.loc[needed_genes].to_numpy(int)
        dense_counts = released["counts"][row_indices].toarray().astype(np.float64)
        log2cpm = np.log2(
            dense_counts
            / released["total_counts"][None, :]
            * 1_000_000
            + 1
        )
        expression = pd.DataFrame(
            log2cpm, index=needed_genes, columns=released["barcodes"]
        )
        secact_input = expression.loc[signature_genes]
        selected_groups = sorted(set().union(*set_members.values()))
        activity_path = (
            args.work_output_dir / f"{sample}_selected_group_activity.tsv.gz"
        )
        checkpoint_path = args.work_output_dir / f"{sample}_secact_checkpoint.json"
        checkpoint_key = {
            "sample": sample,
            "secactpy_version": version("secactpy"),
            "seed": SEED,
            "n_rand": args.n_rand,
            "rng_method": args.rng_method,
            "permutation_cache_mode": (
                "process_gsl_exact" if use_permutation_cache else "disabled"
            ),
            "n_spots": int(secact_input.shape[1]),
            "input_genes": int(secact_input.shape[0]),
            "inference_common_genes": int(
                len(secact_input.index.intersection(inference_grouped_signature.index))
            ),
            "signature_groups": int(inference_grouped_signature.shape[1]),
            "signature_sha256": inference_signature_digest,
            "matrix_bytes": released_path(
                args.data_dir, metadata["gsm"], sample, "matrix.mtx.gz"
            ).stat().st_size,
        }
        group_activity = None
        if activity_path.exists() and checkpoint_path.exists():
            saved_key = json.loads(checkpoint_path.read_text(encoding="utf-8"))
            if saved_key == checkpoint_key:
                candidate = pd.read_csv(activity_path, sep="\t", index_col=0)
                valid_checkpoint = (
                    candidate.index.is_unique
                    and candidate.columns.tolist() == list(secact_input.columns)
                    and set(selected_groups).issubset(candidate.index)
                    and np.isfinite(candidate.to_numpy(float)).all()
                )
                if valid_checkpoint:
                    group_activity = candidate
                    print(f"SECACT CHECKPOINT {sample}", flush=True)
        if group_activity is None:
            print(
                f"SECACT START {sample} spots={secact_input.shape[1]} "
                f"genes={secact_input.shape[0]}",
                flush=True,
            )
            result = secact_activity_inference(
                secact_input,
                is_differential=False,
                sig_matrix=inference_grouped_signature,
                is_group_sig=False,
                n_rand=args.n_rand,
                seed=SEED,
                rng_method=args.rng_method,
                use_cache=use_permutation_cache,
                sort_genes=True,
                verbose=True,
            )
            group_activity = result["zscore"]
            activity_tmp = activity_path.with_name(activity_path.name + ".tmp")
            group_activity.loc[selected_groups].to_csv(
                activity_tmp, sep="\t", compression="gzip"
            )
            activity_tmp.replace(activity_path)
            checkpoint_tmp = checkpoint_path.with_name(checkpoint_path.name + ".tmp")
            checkpoint_tmp.write_text(
                json.dumps(checkpoint_key, indent=2), encoding="utf-8"
            )
            checkpoint_tmp.replace(checkpoint_path)
        absent_activity = set(selected_groups) - set(group_activity.index)
        if absent_activity:
            raise ValueError(f"SecAct output omitted frozen groups: {absent_activity}")
        selected_activity = group_activity.loc[selected_groups]
        if args.checkpoint_samples:
            print(f"SECACT COMPLETE {sample}", flush=True)
            continue

        marker_expression = expression.loc[independent_markers]
        qc_pass = (
            (released["detected_features"] >= 200)
            & (released["total_counts"] >= 500)
        )
        scenarios = {
            "source_filtered": np.ones(len(released["barcodes"]), dtype=bool),
            "fixed_spot_qc": qc_pass,
        }
        scenario_details = {}
        main_scores = None
        for scenario, mask in scenarios.items():
            if int(mask.sum()) < 500:
                raise ValueError(f"{sample} retains fewer than 500 spots in {scenario}")
            current_scores, details = score_spots(
                selected_activity, marker_expression, set_members, mask
            )
            scenario_details[scenario] = details
            if scenario == "source_filtered":
                main_scores = current_scores
            for set_label in sorted(set_members):
                patient_rows.append(
                    {
                        **metadata,
                        "responder": metadata["response"] in {"CR", "PR"},
                        "scenario": scenario,
                        "set_label": set_label,
                        "n_spots": int(mask.sum()),
                        "spearman_rho": patient_correlations(
                            current_scores, set_label
                        ),
                    }
                )
        assert main_scores is not None

        positions = released["positions"]
        spot_frame = pd.DataFrame(
            {
                "sample": sample,
                "gsm": metadata["gsm"],
                "stage": metadata["stage"],
                "response": metadata["response"],
                "responder": metadata["response"] in {"CR", "PR"},
                "barcode": released["barcodes"],
                "array_row": positions["array_row"].to_numpy(),
                "array_col": positions["array_col"].to_numpy(),
                "pxl_col_in_fullres": positions["pxl_col_in_fullres"].to_numpy(),
                "pxl_row_in_fullres": positions["pxl_row_in_fullres"].to_numpy(),
                "nCount": released["total_counts"],
                "nFeature": released["detected_features"],
                "passes_fixed_spot_qc": qc_pass,
                "dynamic_up_set": main_scores.loc[
                    released["barcodes"], "dynamic_up_set"
                ].to_numpy(),
                "dynamic_down_set": main_scores.loc[
                    released["barcodes"], "dynamic_down_set"
                ].to_numpy(),
                "myeloid_context": main_scores.loc[
                    released["barcodes"], "myeloid_context"
                ].to_numpy(),
            }
        )
        spot_frames.append(spot_frame)

        aggregate_counts = np.asarray(
            released["counts"][gene_lookup.loc[signature_genes].to_numpy(int)].sum(axis=1)
        ).ravel()
        aggregate_log2cpm = np.log2(
            aggregate_counts / released["total_counts"].sum() * 1_000_000 + 1
        )
        aggregate_expression.append(
            pd.Series(aggregate_log2cpm, index=signature_genes, name=sample)
        )
        aggregate_sample_names.append(sample)
        coverage_rows.append(
            {
                "sample": sample,
                "signature_genes_observed": len(signature_genes),
                "signature_genes_total": len(signature.index),
                "signature_coverage_fraction": coverage_fraction,
                "unrepresented_dynamic_up_groups": len(unrepresented_groups),
                "independent_myeloid_markers": len(independent_markers),
            }
        )
        run_samples.append(
            {
                **metadata,
                "source_spots": int(len(released["barcodes"])),
                "fixed_qc_spots": int(qc_pass.sum()),
                "signature_gene_coverage": len(signature_genes),
                "signature_coverage_fraction": coverage_fraction,
                "independent_myeloid_markers": independent_markers,
                "scenario_details": scenario_details,
            }
        )
        print(f"SECACT COMPLETE {sample}", flush=True)

    if args.checkpoint_samples:
        print(
            "CHECKPOINT GENERATION COMPLETE: "
            + ", ".join(metadata["sample"] for metadata in samples_to_run),
            flush=True,
        )
        return

    spot_scores = pd.concat(spot_frames, ignore_index=True)
    expected_total_spots = sum(item["expected_spots"] for item in SAMPLES)
    if len(spot_scores) != expected_total_spots:
        raise ValueError(
            "Released Visium cohort spot count differs from Supplementary Data 1: "
            f"{len(spot_scores)} != {expected_total_spots}"
        )
    patient_table = pd.DataFrame(patient_rows)
    summary_rows = []
    for scenario in ("source_filtered", "fixed_spot_qc"):
        for set_label in sorted(set_members):
            summary_rows.append(
                summarize_correlations(patient_table, scenario, set_label)
            )
    summaries = pd.DataFrame(summary_rows)

    main = summaries.loc[
        (summaries["scenario"] == "source_filtered")
        & (summaries["set_label"] == "dynamic_up_set")
    ].iloc[0]
    sensitivity = summaries.loc[
        (summaries["scenario"] == "fixed_spot_qc")
        & (summaries["set_label"] == "dynamic_up_set")
    ].iloc[0]
    minimum_positive = int(np.ceil(0.75 * main["n_patients"]))
    directional_support = bool(
        main["n_patients"] >= 10
        and main["median_patient_rho"] > 0
        and main["positive_patients"] >= minimum_positive
        and main["all_leave_one_out_medians_positive"]
        and sensitivity["median_patient_rho"] > 0
        and sensitivity["positive_patients"]
        >= int(np.ceil(0.75 * sensitivity["n_patients"]))
    )
    formal_support = bool(
        directional_support and main["wilcoxon_p_two_sided"] < 0.05
    )
    summaries["passes_frozen_directional_support"] = False
    summaries["passes_frozen_formal_support"] = False
    main_index = summaries.index[
        (summaries["scenario"] == "source_filtered")
        & (summaries["set_label"] == "dynamic_up_set")
    ]
    summaries.loc[main_index, "passes_frozen_directional_support"] = directional_support
    summaries.loc[main_index, "passes_frozen_formal_support"] = formal_support

    main_patients = patient_table.loc[
        (patient_table["scenario"] == "source_filtered")
        & (patient_table["set_label"] == "dynamic_up_set")
    ].copy()
    response_spatial = exact_rank_test(
        main_patients["spearman_rho"].to_numpy(float),
        main_patients["responder"].to_numpy(bool),
    )
    response_spatial["analysis"] = "patient_spatial_rho_by_response"

    aggregate = pd.concat(aggregate_expression, axis=1, join="inner")
    aggregate = aggregate.loc[:, aggregate_sample_names]
    print(
        f"TISSUE SECACT START patients={aggregate.shape[1]} genes={aggregate.shape[0]}",
        flush=True,
    )
    aggregate_result = secact_activity_inference(
        aggregate,
        is_differential=False,
        sig_matrix=inference_grouped_signature,
        is_group_sig=False,
        n_rand=args.n_rand,
        seed=SEED,
        rng_method=args.rng_method,
        use_cache=use_permutation_cache,
        sort_genes=True,
        verbose=True,
    )
    aggregate_groups = aggregate_result["zscore"]
    aggregate_standardized, aggregate_omitted = standardize_rows(
        aggregate_groups.loc[sorted(set().union(*set_members.values()))]
    )
    tissue_scores = pd.DataFrame(index=aggregate.columns)
    for set_label, members in set_members.items():
        retained = sorted(set(aggregate_standardized.index) & members)
        tissue_scores[set_label] = aggregate_standardized.loc[retained].mean(axis=0)
    manifest = pd.DataFrame(SAMPLES).set_index("sample")
    tissue_scores = tissue_scores.join(manifest)
    tissue_scores["responder"] = tissue_scores["response"].isin(["CR", "PR"])
    response_tissue = exact_rank_test(
        tissue_scores["dynamic_up_set"].to_numpy(float),
        tissue_scores["responder"].to_numpy(bool),
    )
    response_tissue["analysis"] = "tissue_dynamic_up_by_response"

    spot_scores.to_csv(
        args.spot_scores_output, sep="\t", index=False, compression="gzip"
    )
    patient_table.to_csv(args.patient_correlations_output, sep="\t", index=False)
    summaries.to_csv(args.summary_output, sep="\t", index=False)
    pd.DataFrame(coverage_rows).to_csv(
        args.work_output_dir / "GSE206245_secact_coverage.tsv", sep="\t", index=False
    )
    tissue_scores.reset_index(names="sample").to_csv(
        args.work_output_dir / "GSE206245_tissue_secact_scores.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame([response_spatial, response_tissue]).to_csv(
        args.work_output_dir / "GSE206245_response_secondary.tsv",
        sep="\t",
        index=False,
    )
    run = {
        "dataset": "GSE206245 Visium FFPE NPC",
        "secactpy_version": version("secactpy"),
        "seed": SEED,
        "n_rand": args.n_rand,
        "rng_method": args.rng_method,
        "signature_pre_grouped_once": True,
        "signature_sha256": inference_signature_digest,
        "permutation_cache_mode": (
            "process_gsl_exact" if use_permutation_cache else "disabled"
        ),
        "n_patients": len(SAMPLES),
        "primary_unit": "patient/tumour section",
        "spot_role": "within-section spatial sampling unit",
        "released_vs_supplement_discrepancies": [
            {
                "sample": "ST12",
                "supplementary_data_1_spots": 3669,
                "geo_matrix_barcodes_and_in_tissue_spots": 3668,
                "analysis_value": 3668,
            }
        ],
        "myeloid_marker_candidates_after_predictor_overlap_removal": independent_marker_candidates,
        "myeloid_marker_manifest": str(args.myeloid_markers_tsv),
        "directional_support": directional_support,
        "formal_support": formal_support,
        "response_secondary": [response_spatial, response_tissue],
        "aggregate_zero_variance_groups": aggregate_omitted,
        "samples": run_samples,
    }
    (args.work_output_dir / "GSE206245_visium_run.json").write_text(
        json.dumps(run, indent=2), encoding="utf-8"
    )
    print(json.dumps({"directional_support": directional_support, "formal_support": formal_support}, indent=2), flush=True)


if __name__ == "__main__":
    main()
