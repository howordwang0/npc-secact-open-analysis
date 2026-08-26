from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from secactpy import load_signature


BASELINE_SAMPLES = [f"Patient_{patient}_BT{patient}" for patient in range(1, 7)]
TECHNICAL_PATTERN = re.compile(
    r"^(MT-|RP[SL][0-9]|IG[HKL]|TR[ABDG]|MIR|SNOR|SNAR|RNU|AC[0-9]|AL[0-9]|LINC)"
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("counts_tsv", type=Path)
    parser.add_argument("cell_counts_tsv", type=Path)
    parser.add_argument("candidate_output", type=Path)
    parser.add_argument("selected_output", type=Path)
    args = parser.parse_args()
    args.candidate_output.parent.mkdir(parents=True, exist_ok=True)
    args.selected_output.parent.mkdir(parents=True, exist_ok=True)

    counts = pd.read_csv(args.counts_tsv, sep="\t", index_col=0)
    counts.index = counts.index.astype(str)
    cell_counts = pd.read_csv(args.cell_counts_tsv, sep="\t")
    baseline_myeloid = cell_counts.loc[
        cell_counts["scenario"].eq("hard_qc")
        & cell_counts["sample_id"].isin(BASELINE_SAMPLES)
        & cell_counts["broad_cell_type"].eq("Myeloid")
    ]
    if set(baseline_myeloid["sample_id"]) != set(BASELINE_SAMPLES):
        raise ValueError("Not all six frozen baseline samples have a myeloid pseudobulk")
    if not baseline_myeloid["eligible_50_cells"].all():
        raise ValueError("A baseline myeloid pseudobulk contains fewer than 50 cells")

    differences = {}
    myeloid_cpms = {}
    rest_cpms = {}
    for sample in BASELINE_SAMPLES:
        myeloid_column = f"{sample}|Myeloid"
        rest_columns = [
            column
            for column in counts.columns
            if column.startswith(f"{sample}|") and column != myeloid_column
        ]
        if myeloid_column not in counts or len(rest_columns) != 7:
            raise ValueError(f"Incomplete frozen broad-lineage pseudobulks for {sample}")
        myeloid = counts[myeloid_column]
        rest = counts[rest_columns].sum(axis=1)
        myeloid_cpm = myeloid / myeloid.sum() * 1_000_000
        rest_cpm = rest / rest.sum() * 1_000_000
        differences[sample] = np.log2(myeloid_cpm + 1) - np.log2(rest_cpm + 1)
        myeloid_cpms[sample] = myeloid_cpm
        rest_cpms[sample] = rest_cpm

    difference = pd.DataFrame(differences)
    myeloid_cpm = pd.DataFrame(myeloid_cpms)
    rest_cpm = pd.DataFrame(rest_cpms)
    secact_genes = set(load_signature("secact").index.astype(str).str.upper())
    candidates = pd.DataFrame(
        {
            "gene": difference.index,
            "positive_baseline_patients": (difference > 0).sum(axis=1),
            "minimum_log2cpm_difference": difference.min(axis=1),
            "median_log2cpm_difference": difference.median(axis=1),
            "minimum_myeloid_cpm": myeloid_cpm.min(axis=1),
            "median_myeloid_cpm": myeloid_cpm.median(axis=1),
            "maximum_rest_cpm": rest_cpm.max(axis=1),
        }
    )
    candidates["secact_signature_overlap"] = candidates["gene"].str.upper().isin(
        secact_genes
    )
    candidates["technical_gene_exclusion"] = candidates["gene"].str.upper().map(
        lambda gene: bool(TECHNICAL_PATTERN.match(gene))
    )
    candidates["passes_frozen_marker_gate"] = (
        candidates["positive_baseline_patients"].eq(6)
        & candidates["minimum_log2cpm_difference"].gt(0.5)
        & candidates["median_log2cpm_difference"].ge(1.0)
        & candidates["minimum_myeloid_cpm"].ge(1.0)
        & ~candidates["secact_signature_overlap"]
        & ~candidates["technical_gene_exclusion"]
    )
    eligible = candidates.loc[candidates["passes_frozen_marker_gate"]].sort_values(
        [
            "minimum_log2cpm_difference",
            "median_log2cpm_difference",
            "median_myeloid_cpm",
            "gene",
        ],
        ascending=[False, False, False, True],
    )
    eligible = eligible.assign(selection_rank=np.arange(1, len(eligible) + 1))
    candidates = candidates.merge(
        eligible[["gene", "selection_rank"]], on="gene", how="left"
    ).sort_values(
        ["passes_frozen_marker_gate", "selection_rank", "gene"],
        ascending=[False, True, True],
        na_position="last",
    )
    selected = eligible.head(20).copy()
    if len(selected) != 20 or selected["secact_signature_overlap"].any():
        raise ValueError("Unable to freeze 20 non-overlapping myeloid markers")

    candidates.to_csv(args.candidate_output, sep="\t", index=False)
    selected.to_csv(args.selected_output, sep="\t", index=False)
    print(
        f"FROZEN MYELOID MARKERS: eligible={len(eligible)} selected={len(selected)}",
        flush=True,
    )
    print(", ".join(selected["gene"]), flush=True)


if __name__ == "__main__":
    main()
