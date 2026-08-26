# CytoSig comparative benchmark frozen protocol v0.1

**Freeze date:** 2026-08-25  
**Freeze point:** The CytoSig software and signature audit was complete, but no
CytoSig activity score had been calculated for any NPC dataset.

## Purpose and fixed implementation

This is the only alternative signaling-inference benchmark selected for the
current manuscript. No second method will be screened and retained according
to its results. Use SecActPy 0.3.1 with `sig_matrix="cytosig"`, no signature
grouping, ridge penalty 10,000 (the standalone CytoSig default), 1,000
permutations, seed 20260823, GSL-compatible randomization, sorted genes, and
the activity z score as the reported score.

The bundled `CytoSig.tsv.gz` was audited before this freeze against the
standalone CytoSig default `signature.centroid` at repository commit
`ce53cc7f496f996dc0f759652e21afcd0f20bb8e`. Both contain 4,881 genes and 43
cytokines; indices, columns, and every numeric value are identical. The
decompressed matrix SHA-256 is
`60a3152afe429a7dbb41fa0597ab9730141ce858af16a35f7063e64bfc239bde`.
CytoSig is GPLv3; SecActPy 0.3.1 is MIT licensed.

## Complete-reporting rule

All 43 default, high-confidence CytoSig cytokines are analyzed and retained.
The expanded 51-cytokine and beta signatures are prohibited. No cytokine,
direction, cell type, QC scenario, or endpoint may be selected after viewing
results. Nominal P values and BH FDR values are both reported. A nominally
positive cytokine cannot be presented as multiplicity-adjusted evidence.

## Paired bulk analysis

- Use the same nine BT/OT patients and released expression matrix as the
  frozen SecAct paired analysis.
- Collapse duplicate non-empty symbols by median, transform as `log2(x+1)`,
  and infer all 18 samples jointly after CytoSig's across-sample row centering.
- For every cytokine report OT-minus-BT mean and median, Hedges gz, patient
  direction count, leave-one-patient-out direction, two-sided exact sign-flip
  P, and BH FDR across all 43 cytokines.
- This is an agnostic benchmark; no CytoSig result can redefine the frozen
  SecAct 70-group up or 93-group down programs.

## Paired single-cell pseudobulk analysis

- Statistical units are the same five paired patients, never cells.
- Main analyses use hard-QC B, Myeloid, and T/NK pseudobulks with at least 50
  cells per patient-timepoint lineage. Fixed 3.3% doublet-excluded matrices are
  required sensitivities. Plasma/B remains a 20-cell sensitivity only.
- Infer the 10 patient-timepoint pseudobulks jointly within each lineage and
  QC scenario. For every cytokine report the same paired statistics and BH FDR
  across 43 cytokines within that analysis.
- With five pairs, the minimum attainable two-sided exact sign-flip P is
  0.0625. Single-cell results are therefore directional/effect-size evidence,
  not stand-alone significance claims.
- Across all 43 cytokines, report the Spearman correlation and sign concordance
  between paired-bulk and main hard-QC Myeloid mean OT-minus-BT effects. This
  global concordance is descriptive because cytokine signatures are
  correlated.

## Visium spatial analysis

- Use the same 12 source-filtered GSE206245 sections and the same fixed spot-QC
  sensitivity masks as the frozen SecAct analysis. ST12 remains 3,668 released
  spots; no spot is imputed.
- Reconstruct `log2(CPM+1)` from released raw counts, collapse duplicate gene
  symbols by summed counts, and infer all 43 cytokines separately within each
  section. At least 80% of the 4,881 signature genes must be observed in every
  section.
- Use the already frozen 20-gene independent myeloid-context score. For each
  cytokine and section calculate the within-section Spearman correlation
  between CytoSig activity and myeloid context. Spots are spatial sampling
  units only; inferential N is 12 patients.
- For every cytokine summarize the 12 patient correlations by median, positive
  patient count, leave-one-patient-out medians, two-sided exact one-sample
  Wilcoxon P, exact sign P, and BH FDR across 43 cytokines. Repeat correlations
  on the fixed-QC spot subset without refitting CytoSig.
- Formal spatial support for an individual cytokine requires source-filtered
  BH FDR <0.05, positive median correlation, at least 9/12 positive patients,
  all leave-one-patient-out medians positive, and preservation of a positive
  median with at least 9/12 positive patients under fixed spot QC.

## Interpretation boundaries

CytoSig infers downstream transcriptional response to cytokine signaling. It
does not measure cytokine abundance, secretion, producing-cell identity,
receptor engagement, protein activity, treatment causality, or clinical
prediction. Agreement with SecAct supports convergent transcriptomic context;
disagreement is a reportable method difference and cannot be hidden or used
to revise either method's frozen rules.
