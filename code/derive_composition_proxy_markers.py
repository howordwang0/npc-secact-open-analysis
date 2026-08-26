from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from secactpy import load_signature


TECHNICAL_PATTERN = re.compile(
    r"^(MT-|RPS|RPL|RP[0-9]|IG[HKL]|TR[ABDG]|MIR|SNOR|SNAR|RNU|AC[0-9]|AL[0-9]|"
    r"LINC|LOC|MALAT1$|NEAT1$|XIST$|TSIX$|RPS4Y|DDX3Y$|KDM5D$|UTY$|EIF1AY$|"
    r"ZFY$|TTTY|USP9Y$)|-AS[0-9]*$"
)
BASELINE_SAMPLES = [f"Patient_{patient}_BT{patient}" for patient in range(1, 7)]


def excluded_genes(
    group_members: pd.DataFrame, myeloid_markers: pd.DataFrame
) -> set[str]:
    genes = set(load_signature("secact").index.astype(str).str.upper())
    genes.update(group_members["protein"].astype(str).str.upper())
    genes.update(myeloid_markers["gene"].astype(str).str.upper())
    return genes


def derive_gse150430(
    expression: pd.DataFrame,
    manifest: pd.DataFrame,
    exclusions: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    compartments = ["Myeloid", "Lymphoid", "Epithelial"]
    medians = {}
    detection = {}
    eligible_counts = {}
    for compartment in compartments:
        selected = manifest.loc[
            manifest["sample_type"].eq("primary_NPC")
            & manifest["compartment"].eq(compartment)
            & manifest["eligible_5_cells"]
        ]
        if len(selected) < 5:
            raise ValueError(f"Fewer than five reference patients for {compartment}")
        values = expression[selected["group_id"]]
        medians[compartment] = values.median(axis=1)
        detection[compartment] = (values > 0).sum(axis=1)
        eligible_counts[compartment] = len(selected)

    median_frame = pd.DataFrame(medians)
    rows = []
    selected_frames = []
    for compartment in compartments:
        others = [name for name in compartments if name != compartment]
        margin = median_frame[compartment] - median_frame[others].max(axis=1)
        table = pd.DataFrame(
            {
                "gene": expression.index,
                "source": "GSE150430",
                "compartment": compartment,
                "eligible_reference_patients": eligible_counts[compartment],
                "positive_reference_patients": detection[compartment].to_numpy(),
                "required_positive_patients": int(
                    np.ceil(eligible_counts[compartment] / 2)
                ),
                "compartment_median": median_frame[compartment].to_numpy(),
                "maximum_other_compartment_median": median_frame[others]
                .max(axis=1)
                .to_numpy(),
                "minimum_reference_margin": margin.to_numpy(),
            }
        )
        table["excluded_signature_or_direct_or_context"] = table["gene"].isin(
            exclusions
        )
        table["technical_gene_exclusion"] = table["gene"].map(
            lambda gene: bool(TECHNICAL_PATTERN.search(gene))
        )
        table["passes_frozen_marker_gate"] = (
            table["positive_reference_patients"].ge(
                table["required_positive_patients"]
            )
            & table["minimum_reference_margin"].gt(0.25)
            & ~table["excluded_signature_or_direct_or_context"]
            & ~table["technical_gene_exclusion"]
        )
        eligible = table.loc[table["passes_frozen_marker_gate"]].sort_values(
            [
                "minimum_reference_margin",
                "positive_reference_patients",
                "compartment_median",
                "gene",
            ],
            ascending=[False, False, False, True],
        )
        eligible = eligible.assign(selection_rank=np.arange(1, len(eligible) + 1))
        selected = eligible.head(30).copy()
        if len(selected) != 30:
            raise ValueError(f"Unable to select 30 {compartment} markers")
        table = table.merge(
            eligible[["gene", "selection_rank"]], on="gene", how="left"
        )
        rows.append(table)
        selected_frames.append(selected)
    return pd.concat(rows, ignore_index=True), pd.concat(selected_frames, ignore_index=True)


def derive_stromal(
    counts: pd.DataFrame,
    cell_counts: pd.DataFrame,
    exclusions: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    differences = {}
    fibro_cpms = {}
    for sample in BASELINE_SAMPLES:
        current = cell_counts.loc[
            cell_counts["scenario"].eq("hard_qc")
            & cell_counts["sample_id"].eq(sample)
        ]
        fibro = current.loc[current["broad_cell_type"].eq("Fibroblasts")]
        if len(fibro) != 1 or int(fibro.iloc[0]["n_cells"]) < 10:
            raise ValueError(f"Frozen fibroblast eligibility failed for {sample}")
        fibro_column = f"{sample}|Fibroblasts"
        rest_columns = [
            column
            for column in counts.columns
            if column.startswith(f"{sample}|") and column != fibro_column
        ]
        if len(rest_columns) != 7:
            raise ValueError(f"Expected seven pooled-rest columns for {sample}")
        fibro_counts = counts[fibro_column]
        rest_counts = counts[rest_columns].sum(axis=1)
        fibro_cpm = fibro_counts / fibro_counts.sum() * 1_000_000
        rest_cpm = rest_counts / rest_counts.sum() * 1_000_000
        differences[sample] = np.log2(fibro_cpm + 1) - np.log2(rest_cpm + 1)
        fibro_cpms[sample] = fibro_cpm
    difference = pd.DataFrame(differences)
    fibro_cpm = pd.DataFrame(fibro_cpms)
    table = pd.DataFrame(
        {
            "gene": difference.index,
            "source": "Zenodo13736891_baseline_scrna",
            "compartment": "Stromal",
            "eligible_reference_patients": 6,
            "positive_reference_patients": (difference > 0).sum(axis=1),
            "required_positive_patients": 6,
            "compartment_median": np.log2(fibro_cpm + 1).median(axis=1),
            "maximum_other_compartment_median": np.nan,
            "minimum_reference_margin": difference.min(axis=1),
            "median_reference_margin": difference.median(axis=1),
            "minimum_compartment_cpm": fibro_cpm.min(axis=1),
        }
    )
    table["excluded_signature_or_direct_or_context"] = table["gene"].isin(exclusions)
    table["technical_gene_exclusion"] = table["gene"].map(
        lambda gene: bool(TECHNICAL_PATTERN.search(gene))
    )
    table["passes_frozen_marker_gate"] = (
        table["positive_reference_patients"].eq(6)
        & table["minimum_reference_margin"].gt(0.5)
        & table["median_reference_margin"].ge(1.0)
        & table["minimum_compartment_cpm"].ge(1.0)
        & ~table["excluded_signature_or_direct_or_context"]
        & ~table["technical_gene_exclusion"]
    )
    eligible = table.loc[table["passes_frozen_marker_gate"]].sort_values(
        [
            "minimum_reference_margin",
            "median_reference_margin",
            "compartment_median",
            "gene",
        ],
        ascending=[False, False, False, True],
    )
    eligible = eligible.assign(selection_rank=np.arange(1, len(eligible) + 1))
    selected = eligible.head(30).copy()
    if len(selected) != 30:
        raise ValueError("Unable to select 30 frozen stromal markers")
    table = table.merge(eligible[["gene", "selection_rank"]], on="gene", how="left")
    return table, selected


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("gse150430_expression", type=Path)
    parser.add_argument("gse150430_manifest", type=Path)
    parser.add_argument("zenodo_counts", type=Path)
    parser.add_argument("zenodo_cell_counts", type=Path)
    parser.add_argument("group_members", type=Path)
    parser.add_argument("myeloid_context_markers", type=Path)
    parser.add_argument("candidate_output", type=Path)
    parser.add_argument("selected_output", type=Path)
    args = parser.parse_args()
    args.candidate_output.parent.mkdir(parents=True, exist_ok=True)
    args.selected_output.parent.mkdir(parents=True, exist_ok=True)

    expression = pd.read_csv(args.gse150430_expression, sep="\t", index_col=0)
    expression.index = expression.index.astype(str).str.upper()
    manifest = pd.read_csv(args.gse150430_manifest, sep="\t")
    group_members = pd.read_csv(args.group_members, sep="\t")
    myeloid_markers = pd.read_csv(args.myeloid_context_markers, sep="\t")
    exclusions = excluded_genes(group_members, myeloid_markers)

    gse_candidates, gse_selected = derive_gse150430(
        expression, manifest, exclusions
    )
    counts = pd.read_csv(args.zenodo_counts, sep="\t", index_col=0)
    counts.index = counts.index.astype(str).str.upper()
    cell_counts = pd.read_csv(args.zenodo_cell_counts, sep="\t")
    exclusions.update(gse_selected["gene"].astype(str).str.upper())
    stromal_candidates, stromal_selected = derive_stromal(
        counts, cell_counts, exclusions
    )
    candidates = pd.concat([gse_candidates, stromal_candidates], ignore_index=True)
    selected = pd.concat([gse_selected, stromal_selected], ignore_index=True)
    if selected.groupby("compartment").size().to_dict() != {
        "Epithelial": 30,
        "Lymphoid": 30,
        "Myeloid": 30,
        "Stromal": 30,
    }:
        raise ValueError("Final compartment marker dimensions changed")
    if selected["gene"].duplicated().any():
        raise ValueError("Compartment marker sets overlap")
    candidates.to_csv(args.candidate_output, sep="\t", index=False)
    selected.to_csv(args.selected_output, sep="\t", index=False)
    for compartment, subset in selected.groupby("compartment"):
        print(f"{compartment}: {', '.join(subset['gene'])}", flush=True)


if __name__ == "__main__":
    main()
