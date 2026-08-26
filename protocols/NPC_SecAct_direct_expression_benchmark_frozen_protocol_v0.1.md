# Direct ligand-expression benchmark frozen protocol v0.1

**Freeze date:** 2026-08-25  
**Freeze point:** All SecAct, single-cell, GeoMx, and Visium primary analyses
were complete. No direct ligand-program score described below had been
calculated or viewed.

## Purpose

Test whether the frozen SecAct treatment programs provide information beyond a
simpler score formed from the expression of the corresponding secreted-protein
coding genes. This benchmark cannot alter the 70-group treatment-up set, the
93-group treatment-down set, cohort eligibility, or any primary SecAct result.

## Frozen direct-expression score

1. Use `Zenodo13736891_secact_group_members.tsv` to map each already frozen
   correlated SecAct group to its secreted-protein gene symbol(s).
2. Use only coding genes observed in the relevant expression matrix. Do not
   substitute aliases, related family members, receptors, downstream genes, or
   hand-selected markers.
3. Standardize each observed coding gene with sample SD within the relevant
   analysis domain: all 18 paired bulk samples; the 10 paired pseudobulks for a
   given cell type and QC scenario; or spots within one Visium section and QC
   scenario.
4. Average observed genes within each correlated group, then average groups
   within the 70-group or 93-group set. Groups, not the number of measured
   member genes, receive equal weight.
5. Record the number and fraction of frozen groups represented. No missing
   group is imputed. Any coding gene overlapping the 20-gene independent
   myeloid-context score is excluded from the Visium direct score before group
   coverage is assessed, preventing a shared-gene correlation.
6. For the complete bulk table, collapse duplicate non-empty gene symbols by
   the across-row median and transform abundance as `log2(x+1)`. Single-cell
   and Visium raw counts are summed for duplicate symbols before
   `log2(CPM+1)` normalization.

## Paired bulk benchmark

- Use the same nine BT/OT patient pairs as the SecAct dynamic analysis.
- Report OT-minus-BT mean and median, Hedges gz, patient directions, all
  leave-one-patient-out mean directions, and the two-sided exact sign-flip P.
- This analysis is descriptive because the 70/93 groups were selected from
  SecAct behavior in the same bulk cohort. It cannot establish independent
  incremental value.

## Single-cell benchmark

- Use the same five paired patients and eligible broad lineages as the frozen
  SecAct analysis: B, Myeloid, and T/NK at the 50-cell threshold, with hard-QC
  main and fixed doublet-excluded sensitivity matrices.
- Apply the same expected directions, patient-direction counts, leave-one-out
  checks, and exact sign flips.
- Plasma cells remain a 20-cell sensitivity only. No ineligible lineage is
  rescued through the benchmark.
- Report direct-expression and SecAct results side by side; do not declare one
  superior solely because one crosses an arbitrary P-value threshold.

## Visium benchmark and frozen incremental test

- Reconstruct `log2(CPM+1)` for the secreted-protein coding genes from the same
  released matrices used in the SecAct analysis.
- In each section, calculate the within-section Spearman correlation between
  the frozen direct treatment-up expression score and the independent
  myeloid-context score. Repeat after fixed spot QC.
- Summarize the 12 patient correlations using the same median, sign count,
  leave-one-out medians, exact sign test, and exact one-sample Wilcoxon test.
- The prespecified incremental comparison is the paired patient difference
  `SecAct rho - direct-expression rho` in the source-filtered scenario.
  Evidence that SecAct adds spatial association requires a positive median
  paired difference and a two-sided exact paired Wilcoxon P<0.05. Fixed spot QC
  is a required directional sensitivity.
- If this rule is not met, the manuscript must state that incremental value
  over direct ligand expression was not demonstrated, regardless of the
  primary SecAct spatial result.

## Interpretation boundaries

Direct coding-gene expression is not a protein secretion assay and SecAct is
not a direct activity assay. A stronger spatial association for either score
does not prove secretion, source, receptor engagement, treatment causality, or
clinical prediction. This benchmark tests relative transcriptomic association
under one frozen construction only.
