# GSE206245 Visium spatial contextualization frozen protocol v0.2

**Freeze point:** 2026-08-24, after technical reconstruction and coverage QC
of ST5, but before any GSE206245 SecAct activity or spatial endpoint was
calculated or viewed.

## Reason for the v0.2 amendment

Protocol v0.1 carried forward 14 canonical myeloid genes and required removal
of every gene occurring in a frozen dynamic-up SecAct predictor signature.
Whole-transcriptome feasibility auditing showed that only two of the 14 genes
remained, below the prespecified minimum of eight. The v0.1 spatial endpoint
was therefore non-computable and was retired without viewing an activity score,
correlation, response contrast, or outcome-stratified result. Its failure is
retained as a protocol deviation/feasibility result rather than hidden.

Version 0.2 replaces only the independent myeloid-context marker definition.
The cohort, patient manifest, SecAct sets, normalization, QC, primary
correlation, patient-level inference, direction rule, sensitivity analysis,
and interpretation limits remain as specified in v0.1.

## Released-input discrepancy

During byte-validated reconstruction, the released ST12 Matrix Market header,
barcode file, and `in_tissue=1` position records independently agreed on 3,668
spots. Supplementary Data 1 reports 3,669. No unmatched or duplicated released
barcode exists. The analysis therefore uses the 3,668 observable GEO spots and
records the one-spot publication/release discrepancy; it does not impute an
unreleased spot. This input correction was made before any ST12 activity was
calculated.

## Independent marker derivation

Markers are derived from the already frozen broad-lineage annotations of the
independent Zenodo 13736891 single-cell cohort. No GSE206245 activity or spatial
association is used for selection.

1. Use hard-QC baseline pseudobulks from patients P1-P6.
2. Within each patient, compare the Myeloid pseudobulk with the raw-count sum
   of the other seven included broad lineages. Normalize each side as
   `log2(CPM + 1)` and calculate the within-patient difference.
3. Exclude every gene in the 7,919-gene SecAct signature, preventing direct
   predictor overlap with any dynamic-up group.
4. Exclude mitochondrial, ribosomal, immunoglobulin/TCR, common small-RNA, and
   common uncharacterized/non-coding prefixes using the frozen regular
   expression in `work/derive_scrna_myeloid_context_markers.py`.
5. A candidate must be higher in Myeloid in all six baseline patients, have a
   minimum within-patient log2-CPM difference greater than 0.5, a median
   difference at least 1.0, and Myeloid CPM at least 1 in every patient.
6. Rank eligible genes by descending minimum patient difference, then median
   difference, then median Myeloid CPM, with gene symbol as the final
   deterministic tie-breaker. Freeze the first 20 genes.

This produced 539 eligible genes and the following frozen set:
`MS4A6A`, `OLR1`, `FCN1`, `S100A8`, `IL1B`, `LGALS2`, `IGSF6`, `FPR1`,
`FPR3`, `P2RY13`, `C15orf48`, `CSF2RA`, `LILRB2`, `CXCL8`, `LILRB4`,
`PLA2G7`, `CLEC4E`, `CLEC10A`, `CD300E`, and `ADAMDEC1`.

The reproducible candidate audit is
`outputs/NPC_SecAct_Zenodo13736891_myeloid_marker_candidates_v0.1.tsv`; the
frozen 20-gene manifest is
`outputs/NPC_SecAct_Zenodo13736891_myeloid_context_markers_v0.1.tsv`.

## Spatial scoring and eligibility

For each section, use `log2(CPM + 1)` expression, standardize every retained
marker across source-filtered spots with sample SD, and average markers with
equal weight. At least 16 of 20 markers must be present and retain non-zero
variance in both the source-filtered and fixed-QC scenarios. This score is a
myeloid-state context score; it must not be interpreted as a calibrated cell
fraction.

The primary patient statistic remains the within-section Spearman correlation
between the frozen 70-group dynamic-up SecAct score and this independent
myeloid-context score. Inference remains based on the 12 patient correlations,
never on spot-level P values.

Directional support still requires a positive patient median, at least 75%
positive patients, all leave-one-patient-out medians positive, and preservation
of the positive median and 75% sign threshold under fixed spot QC. Formal
support additionally requires a two-sided exact one-sample Wilcoxon
`P < 0.05`.

## Remaining interpretation limits

The analysis tests spatial co-occurrence, not cellular secretion, molecular
interaction, causality, or a calibrated myeloid abundance. The cohort remains
small, retrospective, baseline-only, and treated with heterogeneous
toripalimab-plus-chemotherapy regimens. Response comparisons remain secondary
and cannot rescue a failed primary spatial endpoint.
