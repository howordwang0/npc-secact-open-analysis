# GSE299775 frozen Original Article extension: results report v0.1

**Analysis date:** 2026-08-27  
**Frozen protocol:** `outputs/NPC_SecAct_Original_Article_GSE299775_frozen_protocol_v0.1.md`  
**Dataset:** GSE299775, PanCK-positive GeoMx WTA compartments

## Cohort and endpoint definition

The GEO series matrix states that patients who developed distant metastasis within three years after treatment were classified as metastatic NPC (MET-NPC), whereas those without early metastasis were classified as non-metastatic NPC (nonMET-NPC). The primary analysis included one PanCK-positive ROI from each of 24 independent primary NPC patients: 14 MET and 10 nonMET. Four normal ROIs from three donors were descriptive only.

All 7,807 of the 7,919 SecAct predictor genes previously documented as covered were retained in the formal run. SecActPy 0.3.1 returned all expected activity outputs, and all 70 independently frozen treatment-up signature groups had nonzero variance in both the SecAct and direct-expression implementations. There were no sample exclusions or missing composite scores.

## Frozen primary result

The frozen SecAct treatment-up composite was lower in MET than nonMET tumors (mean -0.249 versus 0.349; mean difference -0.598). The effect magnitude was moderate to large (Hedges g -0.738), but the prespecified statistical-support rule was not met:

- exact two-sided label-permutation P = 0.07886 (all 1,961,256 allocations);
- Welch two-sided P = 0.07275;
- Mann-Whitney two-sided P = 0.05704;
- 2,000-bootstrap Hedges g 95% interval, -1.783 to 0.014.

This is a negative primary result under the frozen protocol. The point estimate is potentially meaningful, but its uncertainty includes no difference. No individual SecAct group was screened and no sample, phenotype, direction or threshold was changed.

## Prespecified direct-expression benchmark

The same 70-group program measured by direct ligand expression was also lower in MET than nonMET tumors (mean -0.191 versus 0.267; mean difference -0.458). This benchmark met its association criteria:

- Hedges g -1.023; 2,000-bootstrap 95% interval, -2.244 to -0.200;
- exact two-sided label-permutation P = 0.01884;
- Welch two-sided P = 0.03208;
- Mann-Whitney two-sided P = 0.02811.

SecAct and direct-expression patient scores were correlated (post-result descriptive Spearman rho 0.716; P = 0.0000844). The direct-expression effect was larger in absolute value by 0.285 at the point estimate. However, the paired-bootstrap 95% interval for `abs(g SecAct) - abs(g direct)` was -1.242 to 0.598. Therefore, this dataset does not statistically establish superiority of direct expression; it establishes that direct expression was associated with the endpoint while SecAct did not show supported association or incremental separation.

Apparent AUC values were 0.264 for SecAct and 0.229 for direct expression with MET fixed as the positive class. Their orientation reflects lower scores in MET. They are descriptive only and must not be reported as validated prediction performance or converted into optimized cutoffs.

## Order and confounding audit

MET status was not completely ordered by acquisition position. Plate column (P = 0.226), acquisition order (P = 0.141) and the overall plate-row table (descriptive chi-square P = 0.201) did not show a statistically clear phenotype association, although plate row B contained six MET and no nonMET samples. Expression PC1 differed by phenotype (Hedges g 1.224; Welch P = 0.0109), whereas PC2 and PC3 did not. PC1 may contain the biological metastasis signal, residual technical structure, or both; it is not an independently measured batch variable. Consistent with the frozen protocol, no outcome-selected expression-PC adjustment was performed. GEO reports upstream QC, batch-effect removal and Q3 normalization.

## Interpretation for the full manuscript

GSE299775 should be retained as an independent, clinically annotated epithelial-compartment test, but not called validation of a metastasis predictor. The result strengthens the manuscript's method-boundary argument:

1. the treatment-up SecAct program has reproducible myeloid-rich localization in the existing single-cell and Visium analyses;
2. in an independent PanCK-positive clinical spatial cohort, its phenotype association remained uncertain despite a moderate-to-large point estimate;
3. the matched direct-expression program showed endpoint association, while SecAct showed no demonstrated incremental advantage;
4. this pattern is concordant with the previous spatial benchmark in which direct expression was at least as strong as SecAct.

The appropriate Original Article claim is that multimodal benchmarking defines where SecAct localization is reproducible but incremental clinical or spatial value over direct ligand expression is not demonstrated. The data do not support a secretion-measurement, metastatic mechanism, prospective prediction or clinical-utility claim.

## Original Article decision

Proceed with the Original Article upgrade as a rigorous multimodal benchmarking study, not as a biomarker paper. GSE299775 adds 24 independent endpoint-annotated patients and a distinct whole-transcriptome spatial epithelial compartment, making the negative incremental-value conclusion more informative. It does not by itself guarantee a Q2 acceptance; the next upgrade should integrate this result with the existing bulk, single-cell, Visium, GeoMx, CytoSig, direct-expression and clinical negative cohorts in a single prespecified evidence hierarchy.

## Reproducible outputs

- Formal run record: `work/gse299775_original_article/GSE299775_original_article_run.json`
- Patient-level scores: `work/gse299775_original_article/GSE299775_patient_level_program_scores.tsv`
- Primary and benchmark statistics: `work/gse299775_original_article/GSE299775_primary_and_benchmark_statistics.tsv`
- Paired bootstrap draws: `work/gse299775_original_article/GSE299775_paired_bootstrap_effects.tsv.gz`
- Batch/order audit: `work/gse299775_original_article/GSE299775_batch_order_audit.tsv`
- Full SecAct matrices and group activity: `work/gse299775_original_article/GSE299775_secact_*.tsv.gz`

No deviations from frozen sample selection, scoring, inference or primary testing were made.
