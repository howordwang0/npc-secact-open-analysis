from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


ANNOTATIONS = {
    0: ("T_NK", "Cytotoxic_T", True, "ZNF683, GZMK, CD8A, CCL5"),
    1: ("T_NK", "CD4_T", True, "CD40LG, IL7R, LEF1, TCF7"),
    2: ("T_NK", "Regulatory_T", True, "FOXP3, CCR8, IL2RA, TNFRSF4"),
    3: ("T_NK", "Activated_cytotoxic_T", True, "IFNG, CCL4, FASLG, NKG7"),
    4: ("B", "Naive_B", True, "IGHD, FCER2, TCL1A, CD22"),
    5: ("B", "Memory_activated_B", True, "TNFRSF13B, AIM2, BANK1"),
    6: ("Myeloid", "Monocyte_macrophage", True, "LYZ, FCN1, C1QA, SPP1"),
    7: ("T_NK", "NK", True, "KLRF1, FGFBP2, NCAM1, KLRD1"),
    8: ("PlasmaB", "Plasma_cell", True, "IGHG1, IGHA1, JCHAIN, MZB1"),
    9: ("T_NK", "Tfh_exhausted_T", True, "IL21, CXCL13, TOX2, PDCD1"),
    10: (
        "Ambiguous",
        "Cycling_lymphocyte",
        False,
        "MKI67, TOP2A, RRM2; mixed T/B lineage markers",
    ),
    11: ("Epi_Cancer", "Epithelial_cancer", True, "FGFBP1, KRT8, KRT16, CLDN3"),
    12: ("Fibroblasts", "Fibroblast", True, "DCN, COL1A1, COL3A1, SFRP2"),
    13: ("Endothelial", "Endothelial", True, "ACKR1, VWF, EMCN, CLDN5"),
    14: (
        "Ambiguous",
        "T_B_doublet",
        False,
        "Concurrent T/B markers; 80.48% in fixed top-3.3% doublet set",
    ),
    15: ("Mast", "Mast", True, "TPSB2, CPA3, TPSAB1, MS4A2"),
    16: ("Myeloid", "Plasmacytoid_DC", True, "SCT, CLEC4C, LILRA4, PTCRA"),
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("cell_metadata", type=Path)
    parser.add_argument("composition", type=Path)
    parser.add_argument("output_tsv", type=Path)
    parser.add_argument("output_note", type=Path)
    args = parser.parse_args()

    metadata = pd.read_csv(args.cell_metadata, sep="\t", index_col=0)
    composition = pd.read_csv(args.composition, sep="\t")
    observed = sorted(metadata["leiden_0p5"].astype(int).unique())
    if observed != sorted(ANNOTATIONS):
        raise ValueError(f"Annotation keys do not match observed clusters: {observed}")

    rows = []
    for cluster in observed:
        cluster_meta = metadata.loc[metadata["leiden_0p5"].astype(int) == cluster]
        annotation = ANNOTATIONS[cluster]
        cluster_comp = composition.loc[composition["leiden_0p5"].astype(int) == cluster]
        rows.append(
            {
                "cluster": cluster,
                "broad_cell_type": annotation[0],
                "descriptive_subtype": annotation[1],
                "include_in_lineage_secact": annotation[2],
                "n_cells": len(cluster_meta),
                "n_samples_with_cells": int((cluster_comp["n_cells"] > 0).sum()),
                "fixed_top_3p3pct_doublet_n": int(
                    cluster_meta["fixed_top_3p3pct_doublet"].sum()
                ),
                "fixed_top_3p3pct_doublet_fraction": float(
                    cluster_meta["fixed_top_3p3pct_doublet"].mean()
                ),
                "supporting_markers": annotation[3],
            }
        )

    result = pd.DataFrame(rows)
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output_tsv, sep="\t", index=False)

    included = result.loc[result["include_in_lineage_secact"]]
    excluded = result.loc[~result["include_in_lineage_secact"]]
    note = f"""# Frozen scRNA-seq cluster annotation v0.1

## Scope

This file freezes the broad lineage mapping for the 90,315 hard-QC cells from
Zenodo record 13736891 before SecAct is run. Harmony coordinates were used only
for clustering and visualization. Raw UMI counts will be used for aggregation
and SecAct.

## Evidence used

- Canonical marker means and detection fractions across Leiden 0.5 clusters.
- Descriptive cluster-versus-rest markers calculated directly from raw counts.
- Sample composition and fixed top-3.3% Scrublet-score enrichment.
- The source study supplementary marker table, including the
  `DC_C4_LILRA4` correspondence for cluster 16.

## Frozen decisions

- Included clusters: {', '.join(map(str, included['cluster']))}.
- Excluded ambiguous clusters: {', '.join(map(str, excluded['cluster']))}.
- Cluster 10 is excluded because cell-cycle markers dominate and T/B lineage
  assignment is unresolved.
- Cluster 14 is excluded because it co-expresses T and B markers and 80.48% of
  its cells fall in the independently frozen doublet-sensitivity set.
- Cluster 16 is assigned to Myeloid / plasmacytoid DC based on concordant
  `SCT`, `CLEC4C`, `LILRA4`, and `PTCRA` markers and the source supplement.

## Downstream unit of analysis

Raw counts will be summed by patient, timepoint, and broad cell type. A group
requires at least 50 cells in the main analysis; 20 cells is the frozen
sensitivity threshold. Cells and clusters are not independent replicates.
"""
    args.output_note.write_text(note, encoding="utf-8")
    print(f"FROZEN ANNOTATION WRITTEN: {args.output_tsv}", flush=True)


if __name__ == "__main__":
    main()
