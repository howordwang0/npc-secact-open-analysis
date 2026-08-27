# GSE299775 Original Article extension: frozen analysis protocol v0.1

**Frozen:** 2026-08-27 (Asia/Shanghai), before examining any GSE299775 activity-phenotype association  
**Status:** prospective extension of the completed public-data benchmark  
**Rationale:** this is not a validation analysis of GSE103611 candidates. The earlier GSE103611 screen produced no candidate and its stopping rule remains intact. This protocol asks a different, independently anchored question using the 70-group treatment-up program frozen in the nine-patient paired treatment dataset.

## Research question

In PanCK-positive whole-transcriptome spatial compartments from independent nasopharyngeal carcinoma patients, is the independently frozen treatment-up SecAct program associated with the reported MET versus nonMET phenotype, and does SecAct show greater phenotype separation than the same program measured by direct ligand expression?

## Analysis population and unit

- Primary population: 24 NPC patients in GSE299775, one PanCK-positive GeoMx WTA ROI per patient (MET 14, nonMET 10).
- Statistical unit: patient. No ROI, gene, protein or signature-group pseudoreplication is permitted.
- The four normal PanCK-positive ROIs are excluded from the primary test. The two ROIs from normal donor NP-1 will be averaged before any descriptive normal comparison, yielding three normal donors.
- No patient or outlier will be removed after inspecting program scores. A sample is excluded only for a prespecified technical failure: missing phenotype, non-NPC status for the primary analysis, non-finite composite score, or fewer than 63 of the 70 frozen groups represented.

## Locked inputs

| Input | SHA-256 |
|---|---|
| `work/data_audit/GSE299775_count_matrix.txt.gz` | `F9573B58267879E87CEE214943A2E178B6F65BC40F2D89F78F1616DA3F4037D3` |
| `work/data_audit/GSE299775_samples.csv` | `103C07D80C2F53109DFC0B62862270459200F8F0266615179F14D461256518BF` |
| `outputs/NPC_SecAct_Zenodo13736891_dynamic_sets_v0.1.tsv` | `18C1C24F687AFEBD752F07F9028FFA2FA1D84C6C63553197435FDC3F22765635` |
| `work/zenodo13736891_bulk/secact_paired_v0.1/Zenodo13736891_secact_group_members.tsv` | `2E10B85F58FADAFEA78ABD9670B064CEC54F29C9D619202686E7433C93D4555E` |

The frozen program contains exactly 70 correlation-collapsed SecAct signature groups. Individual groups will not be screened against phenotype and no post hoc subset will be created.

## Preprocessing and SecAct inference

1. Preserve the GEO matrix sample order and aggregate duplicate gene symbols by their median if any occur.
2. Transform the non-negative Q3-normalized count matrix with `log2(x + 1)`.
3. Run SecActPy 0.3.1 once using `is_differential=False`, the SecAct signature, 1,000 randomizations, GSL random-number generation, seed 20260823 and sorted genes.
4. Collapse correlated output proteins according to the locked membership table by retaining the recorded representative.
5. Using the 24 NPC patients only, standardize each represented SecAct group across patients and average the 70 standardized group values with equal weight. Normal samples are projected using the NPC means and standard deviations.
6. Construct the direct-expression comparator from the same 70 groups: standardize measured ligand genes across the 24 NPC patients, average genes within each signature group, then average groups with equal weight. Normal samples are projected using the corresponding NPC parameters.

## Prespecified QC and order/batch audit

Before interpreting phenotype association, record the plate row, well column and acquisition order parsed from the ROI identifier. Report their distributions by phenotype, Welch/Mann-Whitney tests for numeric order variables, and phenotype associations of expression PC1-PC3. GEO states that QC, batch-effect removal and Q3 normalization were performed. No additional batch or expression-PC adjustment will be selected after seeing the outcome. Complete or near-complete phenotype-order confounding will downgrade all association findings to exploratory.

## Primary and benchmark analyses

The sole primary endpoint is the frozen 70-group SecAct composite, with MET minus nonMET as the fixed contrast. Report group means, mean difference, Hedges g, 2,000-stratified-bootstrap 95% CI for Hedges g, Welch two-sided P, Mann-Whitney two-sided P, and an exact two-sided label-permutation P for the mean difference. With one primary endpoint, no multiplicity correction is applied.

The direct-expression composite is a prespecified benchmark, analyzed identically. Report the paired bootstrap distribution and 95% interval for `abs(Hedges g SecAct) - abs(Hedges g direct expression)`. Apparent AUC values are descriptive only; no cutoff, internally optimized direction, cross-validation claim or prediction claim is allowed.

## Interpretation and go/no-go rules

- **Statistically supported association:** exact two-sided permutation P < 0.05 and the bootstrap 95% CI for Hedges g excludes zero.
- **Potentially meaningful magnitude:** absolute Hedges g >= 0.50. This magnitude criterion cannot override a non-significant primary test.
- **SecAct incremental separation:** the paired-bootstrap 95% interval for the absolute-effect difference is entirely above zero. A larger point estimate alone is insufficient.
- **Original Article upgrade go:** retain the full multimodal study and add GSE299775 as an independent phenotype-linked epithelial-compartment analysis if the association is supported, or if it provides a clear and informative negative benchmark that materially strengthens the method-boundary conclusion. It must never be described as external validation of a predictor.
- **No-go for biomarker language:** if the primary result is unsupported, the cohort remains in the paper as a complete negative result; do not screen the 70 groups, change the phenotype, remove samples, optimize a cutoff or substitute direct expression as a new primary endpoint.

## Reporting boundaries

This retrospective secondary analysis can establish association and method benchmarking only. `MET` is used exactly as supplied by GEO and must be defined from the source publication before manuscript submission. The analysis does not measure secreted protein, prove a mechanism, establish prospective prediction, or demonstrate clinical utility. All deviations from this protocol must be dated and reported.
