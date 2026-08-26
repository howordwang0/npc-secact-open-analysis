from __future__ import annotations

import argparse
import gzip
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


MYELOID_MARKERS = [
    "AIF1",
    "TYROBP",
    "FCER1G",
    "CTSS",
    "CST3",
    "S100A8",
    "S100A9",
    "FCN1",
    "CD14",
    "ITGAM",
    "CSF1R",
    "CD68",
    "CD163",
    "LST1",
]


def ngeomean(values: np.ndarray, axis: int | None = None) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if axis is None:
        shifted = values + 0.5 if np.nanmin(values) < 0.5 else values
        return np.exp(np.nanmean(np.log(shifted)))
    minimum = np.nanmin(values, axis=axis, keepdims=True)
    shifted = values + (minimum < 0.5) * 0.5
    return np.exp(np.nanmean(np.log(shifted), axis=axis))


def parse_dcc(path: Path) -> tuple[dict[str, str], pd.Series]:
    with gzip.open(path, "rt", encoding="utf-8-sig") as handle:
        lines = [line.rstrip("\r\n") for line in handle]
    attributes = {}
    for section in ("NGS_Processing_Attributes", "Scan_Attributes"):
        start = lines.index(f"<{section}>") + 1
        end = lines.index(f"</{section}>")
        for line in lines[start:end]:
            if not line or "," not in line:
                continue
            key, value = line.split(",", 1)
            attributes[key] = value.strip('"')
    start = lines.index("<Code_Summary>") + 1
    end = lines.index("</Code_Summary>")
    counts = {}
    for line in lines[start:end]:
        if not line:
            continue
        rts_id, count = line.split(",", 1)
        counts[rts_id] = int(count)
    return attributes, pd.Series(counts, name=path.name, dtype=float)


def target_aggregate(
    probe_counts: pd.DataFrame,
    probe_table: pd.DataFrame,
    method: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    endogenous = probe_table.loc[probe_table["code_class"] == "Endogenous01"].copy()
    target_rows = []
    qc_rows = []
    for target, block in endogenous.groupby("target", sort=True):
        current = probe_counts.loc[block.index]
        if len(block) > 1:
            target_by_roi = ngeomean(current.to_numpy(), axis=0)
            target_global = float(ngeomean(target_by_roi))
            probe_global = ngeomean(current.to_numpy(), axis=1)
            ratios = probe_global / target_global
            keep = ratios > 0.1
        else:
            ratios = np.ones(1)
            keep = np.ones(1, dtype=bool)
        for rts_id, ratio, retained in zip(block.index, ratios, keep):
            qc_rows.append(
                {
                    "RTS_ID": rts_id,
                    "target": target,
                    "probe_ratio": ratio,
                    "low_probe_ratio": not bool(retained),
                }
            )
        retained_counts = current.iloc[np.flatnonzero(keep)]
        if retained_counts.empty:
            raise ValueError(f"All probes failed the ratio rule for {target}")
        if method == "geomean":
            aggregated = ngeomean(retained_counts.to_numpy(), axis=0)
        elif method == "median":
            aggregated = np.median(retained_counts.to_numpy(), axis=0)
        else:
            raise ValueError(method)
        target_rows.append(pd.Series(aggregated, index=probe_counts.columns, name=target))
    return pd.DataFrame(target_rows), pd.DataFrame(qc_rows)


def q3_normalize(target_counts: pd.DataFrame) -> pd.DataFrame:
    q3 = target_counts.quantile(0.75, axis=0)
    factor = q3 / ngeomean(q3.to_numpy())
    normalized = target_counts.divide(factor, axis=1)
    return np.log2(normalized + 1)


def module_scores(
    expression: pd.DataFrame,
    dynamic_sets: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, list[str]], list[str]]:
    gene_z = expression.sub(expression.mean(axis=1), axis=0).div(
        expression.std(axis=1, ddof=1).replace(0, np.nan), axis=0
    )
    scores = pd.DataFrame(index=expression.columns)
    covered_groups = {}
    up_genes = set()
    for set_label, block in dynamic_sets.groupby("set_label", sort=False):
        group_scores = []
        labels = []
        for group_name in block["signature_group"].astype(str):
            members = sorted(set(group_name.split("|")) & set(expression.index))
            if not members:
                continue
            group_scores.append(gene_z.loc[members].mean(axis=0).rename(group_name))
            labels.append(group_name)
            if set_label == "dynamic_up_set":
                up_genes.update(members)
        covered_groups[set_label] = labels
        scores[set_label] = pd.DataFrame(group_scores).mean(axis=0)

    independent_markers = [
        gene for gene in MYELOID_MARKERS if gene in expression.index and gene not in up_genes
    ]
    if len(independent_markers) < 5:
        raise ValueError(f"Insufficient independent myeloid markers: {independent_markers}")
    scores["independent_myeloid_score"] = gene_z.loc[independent_markers].mean(axis=0)
    return scores, covered_groups, independent_markers


def spearman_exact_permutation(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    x_rank = stats.rankdata(x)
    y_rank = stats.rankdata(y)
    x_centered = x_rank - x_rank.mean()
    y_centered = y_rank - y_rank.mean()
    denominator = np.sqrt(np.sum(x_centered**2) * np.sum(y_centered**2))
    observed = float(np.dot(x_centered, y_centered) / denominator)
    extreme = 0
    total = 0
    for permutation in itertools.permutations(y_centered):
        correlation = float(np.dot(x_centered, permutation) / denominator)
        extreme += abs(correlation) >= abs(observed) - 1e-12
        total += 1
    return observed, extreme / total


def summarize_scenario(
    scores: pd.DataFrame,
    metadata: pd.DataFrame,
    roi_mask: pd.Series,
    aggregation: str,
    roi_scenario: str,
) -> tuple[pd.DataFrame, dict]:
    selected_scores = scores.loc[roi_mask]
    selected_meta = metadata.set_index("roi_id").loc[selected_scores.index]
    combined = selected_scores.join(selected_meta[["PatientID", "Tissue"]])
    patient = (
        combined.groupby(["PatientID", "Tissue"], observed=True)
        .mean(numeric_only=True)
        .reset_index()
    )
    npc = patient.loc[patient["Tissue"] == "NPC"]
    rho, exact_p = spearman_exact_permutation(
        npc["dynamic_up_set"].to_numpy(),
        npc["independent_myeloid_score"].to_numpy(),
    )
    loo = []
    for idx in range(len(npc)):
        reduced = npc.drop(npc.index[idx])
        loo.append(
            stats.spearmanr(
                reduced["dynamic_up_set"], reduced["independent_myeloid_score"]
            ).statistic
        )
    roi_rho = stats.spearmanr(
        combined.loc[combined["Tissue"] == "NPC", "dynamic_up_set"],
        combined.loc[combined["Tissue"] == "NPC", "independent_myeloid_score"],
    ).statistic
    summary = {
        "aggregation": aggregation,
        "roi_scenario": roi_scenario,
        "n_rois": len(combined),
        "n_npc_patients": len(npc),
        "npc_patient_spearman_rho": rho,
        "npc_patient_exact_permutation_p_two_sided": exact_p,
        "npc_leave_one_out_min_rho": float(np.min(loo)),
        "npc_leave_one_out_all_positive": bool(np.all(np.asarray(loo) > 0)),
        "npc_roi_level_descriptive_rho": float(roi_rho),
    }
    patient["aggregation"] = aggregation
    patient["roi_scenario"] = roi_scenario
    return patient, summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("dcc_dir", type=Path)
    parser.add_argument("pkc_gz", type=Path)
    parser.add_argument("roi_metadata_xlsx", type=Path)
    parser.add_argument("dynamic_sets_tsv", type=Path)
    parser.add_argument("work_output_dir", type=Path)
    parser.add_argument("roi_qc_output", type=Path)
    parser.add_argument("patient_scores_output", type=Path)
    parser.add_argument("summary_output", type=Path)
    args = parser.parse_args()
    args.work_output_dir.mkdir(parents=True, exist_ok=True)
    args.roi_qc_output.parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(args.pkc_gz, "rt", encoding="utf-8-sig") as handle:
        pkc = json.load(handle)
    probe_rows = []
    for target in pkc["Targets"]:
        for probe in target["Probes"]:
            probe_rows.append(
                {
                    "RTS_ID": probe["RTS_ID"],
                    "target": target["DisplayName"],
                    "code_class": target["CodeClass"],
                }
            )
    probe_table = pd.DataFrame(probe_rows).set_index("RTS_ID")
    if not probe_table.index.is_unique:
        raise ValueError("PKC RTS_ID values are not unique")

    dcc_rows = []
    count_columns = []
    for path in sorted(args.dcc_dir.glob("*.dcc.gz")):
        attributes, counts = parse_dcc(path)
        extra_ids = set(counts.index) - set(probe_table.index)
        if extra_ids:
            raise ValueError(f"DCC contains RTS_ID values absent from PKC in {path.name}")
        omitted_ids = set(probe_table.index) - set(counts.index)
        gsm = path.name.split("_", 1)[0]
        roi_id = attributes["ID"]
        # DCC files omit a variable number of zero-count probes rather than
        # listing them explicitly. All observed IDs are validated against PKC.
        counts = counts.reindex(probe_table.index, fill_value=0).rename(roi_id)
        count_columns.append(counts)
        dcc_rows.append(
            {
                "roi_id": roi_id,
                "gsm": gsm,
                "dcc_file": path.name,
                "Raw_dcc": int(attributes["Raw"]),
                "Trimmed_dcc": int(attributes["Trimmed"]),
                "Stitched_dcc": int(attributes["Stitched"]),
                "Aligned_dcc": int(attributes["Aligned"]),
                "UMIQ30_dcc": float(attributes["umiQ30"]),
                "RTSQ30_dcc": float(attributes["rtsQ30"]),
                "omitted_zero_probe_ids": len(omitted_ids),
            }
        )
    if len(dcc_rows) != 55:
        raise ValueError(f"Expected 55 DCC files, found {len(dcc_rows)}")
    probe_counts = pd.concat(count_columns, axis=1)
    dcc_meta = pd.DataFrame(dcc_rows)

    roi_meta = pd.read_excel(args.roi_metadata_xlsx)
    if roi_meta["RawReads"].duplicated().any() or dcc_meta["Raw_dcc"].duplicated().any():
        raise ValueError("Raw read counts are not unique for DCC-to-ROI matching")
    metadata = dcc_meta.merge(
        roi_meta,
        left_on="Raw_dcc",
        right_on="RawReads",
        how="left",
        validate="one_to_one",
    )
    if metadata["PatientID"].isna().any():
        raise ValueError("Not all DCC files mapped to ROI metadata")
    metadata = metadata.set_index("roi_id", drop=False)
    probe_counts = probe_counts.loc[:, metadata.index]

    negative_ids = probe_table.index[probe_table["code_class"] == "Negative01"]
    negative_geomean = pd.Series(
        ngeomean(probe_counts.loc[negative_ids].to_numpy(), axis=0),
        index=probe_counts.columns,
    )
    metadata["trimmed_pct"] = 100 * metadata["TrimmedReads"] / metadata["RawReads"]
    metadata["stitched_pct"] = 100 * metadata["StitchedReads"] / metadata["RawReads"]
    metadata["aligned_pct"] = 100 * metadata["AlignedReads"] / metadata["RawReads"]
    metadata["saturation_geomxtools_pct"] = 100 * (
        1 - metadata["DeduplicatedReads"] / metadata["AlignedReads"]
    )
    metadata["negative_probe_geomean"] = negative_geomean
    metadata["passes_geomxtools_default_qc"] = (
        (metadata["RawReads"] >= 1000)
        & (metadata["trimmed_pct"] >= 80)
        & (metadata["stitched_pct"] >= 80)
        & (metadata["aligned_pct"] >= 80)
        & (metadata["saturation_geomxtools_pct"] >= 50)
        & (metadata["negative_probe_geomean"] >= 10)
        & (metadata["AOINucleiCount"] >= 200)
        & (metadata["AOISurfaceArea"] >= 16000)
    )
    metadata.to_csv(args.roi_qc_output, sep="\t", index=False)

    dynamic_sets = pd.read_csv(args.dynamic_sets_tsv, sep="\t")
    all_target_counts = {}
    all_scores = {}
    probe_qc = None
    covered_groups = None
    independent_markers = None
    for method in ("geomean", "median"):
        endogenous_counts, current_probe_qc = target_aggregate(
            probe_counts, probe_table, method
        )
        if method == "geomean":
            probe_qc = current_probe_qc
        # Q3 in GeoMxTools is calculated on the complete target-level matrix.
        other_targets = probe_table.loc[
            probe_table["code_class"] != "Endogenous01"
        ].copy()
        other_rows = []
        for target, block in other_targets.groupby("target", sort=True):
            values = probe_counts.loc[block.index].to_numpy()
            aggregate = ngeomean(values, axis=0) if method == "geomean" else np.median(values, axis=0)
            other_rows.append(pd.Series(aggregate, index=probe_counts.columns, name=target))
        complete_counts = pd.concat([endogenous_counts, pd.DataFrame(other_rows)])
        normalized = q3_normalize(complete_counts).loc[endogenous_counts.index]
        scores, current_covered, current_markers = module_scores(normalized, dynamic_sets)
        all_target_counts[method] = normalized
        all_scores[method] = scores
        if method == "geomean":
            covered_groups = current_covered
            independent_markers = current_markers
        elif current_covered != covered_groups or current_markers != independent_markers:
            raise ValueError("Coverage changed across aggregation methods")
        normalized.to_csv(
            args.work_output_dir / f"GSE289272_{method}_q3_log2.tsv.gz",
            sep="\t",
            compression="gzip",
        )
        scores.to_csv(
            args.work_output_dir / f"GSE289272_{method}_roi_scores.tsv",
            sep="\t",
        )
    assert probe_qc is not None
    assert covered_groups is not None
    assert independent_markers is not None
    probe_qc.to_csv(args.work_output_dir / "GSE289272_probe_qc.tsv", sep="\t", index=False)

    patient_frames = []
    summary_rows = []
    masks = {
        "source_released": pd.Series(True, index=metadata.index),
        "geomxtools_default_qc": metadata["passes_geomxtools_default_qc"],
    }
    for method in ("geomean", "median"):
        for roi_scenario, mask in masks.items():
            patient, scenario_summary = summarize_scenario(
                all_scores[method], metadata, mask, method, roi_scenario
            )
            patient_frames.append(patient)
            summary_rows.append(scenario_summary)
    patient_scores = pd.concat(patient_frames, ignore_index=True)
    summaries = pd.DataFrame(summary_rows)

    main = summaries.loc[
        (summaries["aggregation"] == "geomean")
        & (summaries["roi_scenario"] == "source_released")
    ].iloc[0]
    qc_sensitivity = summaries.loc[
        (summaries["aggregation"] == "geomean")
        & (summaries["roi_scenario"] == "geomxtools_default_qc")
    ].iloc[0]
    median_sensitivity = summaries.loc[
        (summaries["aggregation"] == "median")
        & (summaries["roi_scenario"] == "source_released")
    ].iloc[0]
    passes = bool(
        main["npc_patient_spearman_rho"] > 0
        and main["npc_leave_one_out_all_positive"]
        and qc_sensitivity["npc_patient_spearman_rho"] > 0
        and median_sensitivity["npc_patient_spearman_rho"] > 0
    )
    summaries["passes_frozen_spatial_support"] = False
    summaries.loc[
        (summaries["aggregation"] == "geomean")
        & (summaries["roi_scenario"] == "source_released"),
        "passes_frozen_spatial_support",
    ] = passes
    patient_scores.to_csv(args.patient_scores_output, sep="\t", index=False)
    summaries.to_csv(args.summary_output, sep="\t", index=False)

    run = {
        "dataset": "GSE289272",
        "n_rois": int(len(metadata)),
        "n_patients": int(metadata["PatientID"].nunique()),
        "n_npc_patients": int(metadata.loc[metadata["Tissue"] == "NPC", "PatientID"].nunique()),
        "n_control_patients": int(metadata.loc[metadata["Tissue"] == "CN", "PatientID"].nunique()),
        "n_rois_passing_geomxtools_default_qc": int(metadata["passes_geomxtools_default_qc"].sum()),
        "max_omitted_zero_probe_ids_per_dcc": int(metadata["omitted_zero_probe_ids"].max()),
        "endogenous_targets": int((probe_table["code_class"] == "Endogenous01").groupby(probe_table["target"]).any().sum()),
        "low_ratio_endogenous_probes": int(probe_qc["low_probe_ratio"].sum()),
        "covered_dynamic_groups": {key: len(value) for key, value in covered_groups.items()},
        "independent_myeloid_markers": independent_markers,
        "passes_frozen_spatial_support": passes,
    }
    (args.work_output_dir / "GSE289272_spatial_run.json").write_text(
        json.dumps(run, indent=2), encoding="utf-8"
    )
    print(json.dumps(run, indent=2), flush=True)


if __name__ == "__main__":
    main()
