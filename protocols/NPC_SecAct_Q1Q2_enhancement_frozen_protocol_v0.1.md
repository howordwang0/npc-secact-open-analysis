# NPC SecAct Q1/Q2 enhancement frozen protocol v0.1

**Freeze date:** 2026-08-25  
**Freeze point:** Before downloading or inspecting the expression values, cell annotations, SecAct scores, spatial scores, or endpoints of GSE150430 and GSE200310, and before running the new GSE206245 robustness analyses. Public GEO landing-page metadata and file names were checked only to establish feasibility.  
**Purpose:** Strengthen the fully open-data manuscript through independent localization, external spatial replication, and reviewer-facing robustness analyses. This protocol does not replace or suppress any completed positive, negative, or non-computable result.

## 1. Immutable scientific boundary

The frozen 70-group treatment-up and 93-group treatment-down SecAct programs from the paired bulk cohort are the only SecAct programs tested. No gene or signature group will be selected, removed, reweighted, or renamed using GSE150430, GSE200310, or the new GSE206245 results.

The analyses can support treatment-associated inferred secreted-activity remodeling and localization to a myeloid-rich spatial context. They cannot establish secretion, a cellular source, a ligand-receptor mechanism, treatment causality, clinical prediction, or superiority of SecAct over direct expression. Direct-expression and CytoSig comparisons and all null findings will be retained.

## 2. Open inputs and immutable roles

### GSE150430: independent single-cell localization

- Public processed human scRNA-seq matrix: `GSE150430_npc_scRNA_hg19_processed_data.txt.gz`.
- Cohort role: independent baseline localization only; no treatment-effect replication and no candidate discovery.
- Primary population: primary NPC tumor samples. The normal sample is excluded from the primary endpoint and may be shown descriptively.
- Statistical unit: patient, never cell.

### GSE200310: external spatial replication

- Public processed/raw spatial archive: `GSE200310_RAW.tar`.
- Cohort role: locked, descriptive external replication in two primary NPC specimens.
- Statistical unit: specimen. With only two specimens, no formal population-level claim or inferential P value will be used to declare replication.

### GSE206245: robustness of the existing 12-specimen spatial result

- Existing byte-validated Visium matrices, spot coordinates, and QC rules remain unchanged.
- Cohort role: composition adjustment, spatial-autocorrelation-aware inference, and matched-program negative controls.
- Statistical unit: patient/specimen for population inference. Spot-level calculations are intermediate quantities only.

## 3. Acquisition and source audit

Before biological analysis, record for every source file: exact URL, access date, byte size, SHA-256, compression/archive integrity, extracted member list, matrix dimensions, feature identifier type, barcode count, sample/patient naming structure, and annotation fields. Duplicate feature and barcode handling must be deterministic and reported. No undocumented repair of released data is allowed.

If patient identity cannot be recovered without ambiguity, GSE150430 will be restricted to descriptive cell-level visualization and will not contribute to patient-level support. If author annotations are absent, broad lineages will be assigned using a frozen marker-based workflow without looking at SecAct scores.

## 4. GSE150430 independent single-cell endpoint

### 4.1 Annotation and pseudobulk

Prefer author-supplied cell annotations when embedded in the released matrix. Map them once to the following broad lineages using label semantics and a marker audit: myeloid, T/NK, B, plasma B, malignant/epithelial, fibroblast, endothelial, and mast. Ambiguous, cycling-only, doublet, low-quality, or unclassified labels are excluded rather than reassigned to maximize an effect.

For each patient and broad lineage, sum released expression/count values at the gene level and normalize as `log2(CPM + 1)`. A patient-lineage pseudobulk is eligible with at least 50 cells. The primary contrast requires at least 50 myeloid cells and at least 50 pooled eligible non-myeloid cells in the same patient. The pooled non-myeloid pseudobulk is formed by summing all eligible non-myeloid cells before normalization; it is not an unweighted mean of lineage scores.

### 4.2 SecAct scoring

Use SecActpy with the locked SecAct signature, 1,000 randomizations, seed 20260823, GSL random-number generator, and sorted genes. Infer activities jointly across all eligible patient-compartment pseudobulks. Collapse correlated signatures with the existing `cor_threshold=0.9` grouping. Standardize each group across all eligible pseudobulks using sample SD, then average the frozen 70 treatment-up groups and 93 treatment-down groups separately with equal group weight.

### 4.3 Primary and secondary tests

The primary endpoint is the within-patient difference in the frozen treatment-up score, myeloid minus pooled non-myeloid. Report every patient, median difference, proportion positive, Hodges-Lehmann shift with exact-compatible confidence interval where computable, two-sided exact Wilcoxon signed-rank P value, and all leave-one-patient-out medians.

Independent localization support requires all of the following:

1. at least eight eligible tumor patients;
2. positive median myeloid-minus-non-myeloid difference;
3. at least 75% of patients positive;
4. every leave-one-patient-out median positive; and
5. two-sided patient-level P < 0.05.

The frozen treatment-down program is a negative/specificity comparator and is reported under the same calculations without an expected direction. Secondary analyses rank all eligible broad lineages within patient and report direct ligand-expression and CytoSig scores using the same patient-level structure. These secondary results cannot rescue a failed primary endpoint.

## 5. GSE200310 descriptive spatial replication

Reconstruct each released spatial specimen without changing source spot inclusion. Apply the existing GSE206245 normalization and fixed scores:

- frozen 70-group treatment-up and 93-group treatment-down SecAct programs;
- the frozen 20-gene independent myeloid-context score;
- matched direct ligand-expression program scores; and
- CytoSig scores using the already locked implementation.

For each specimen, report spot count, score coverage, Spearman rho between each program and myeloid context, bootstrap confidence interval using spatial blocks, and a spatial-shift empirical P value as described below. The direction of all measures is shown even when unfavorable.

Descriptive concordance is defined in advance as a positive treatment-up SecAct/myeloid-context rho in both specimens. This label does not constitute formal external validation. The treatment-down program, direct expression, and CytoSig results determine specificity and comparative interpretation; SecAct is not called superior unless a prespecified paired comparison supports it, which is not expected from n=2.

## 6. GSE206245 composition adjustment

### 6.1 Reference construction

Construct a four-compartment reference from GSE150430 using only its primary NPC tumors: myeloid, lymphoid (T/NK, B, and plasma B), malignant/epithelial, and stromal (fibroblast, endothelial, and mast). For each compartment, calculate the median patient pseudobulk profile on the linear CPM scale.

Select deconvolution markers using GSE150430 alone. A marker must be detected at CPM >= 1 in at least half of eligible reference patients for its compartment, have median `log2(CPM + 1)` at least 1 higher than the maximum of the other compartments, and rank among the top 50 by this difference. Exclude mitochondrial, ribosomal, immunoglobulin/TCR, sex-chromosome, common small-RNA, all 7,919 SecAct predictor genes, all frozen direct-expression genes, and all 20 myeloid-context genes. The union of retained markers is fixed before fitting any spatial association.

Estimate non-negative compartment coefficients for each Visium spot by NNLS on shared marker-gene CPM values and normalize coefficients to sum to one. These are reference-based composition estimates, not calibrated cell fractions. Specimens with fewer than 20 total usable markers or fewer than five markers for any compartment are non-computable and retained as such.

### 6.2 Adjusted association

Within each specimen, rank-transform continuous variables and residualize both the treatment-up SecAct score and the frozen myeloid-context score against estimated myeloid proportion, estimated malignant proportion, `log1p(total UMI)`, and detected-gene count using ordinary least squares. The adjusted specimen statistic is Spearman rho between the two residual vectors. Variance-inflation factors, coefficient rank, and residual distributions are audited; a rank-deficient model is non-computable rather than simplified after viewing results.

Report all 12 adjusted rhos, median, proportion positive, leave-one-patient-out medians, and a two-sided exact one-sample Wilcoxon test against zero. Robust adjusted support uses the existing spatial rule: positive median, at least 75% positive, all leave-one-out medians positive, and P < 0.05. Unadjusted and adjusted estimates are displayed together.

## 7. Spatial-autocorrelation-aware null

For GSE206245 and GSE200310, map each spot to released Visium array-row/array-column coordinates. Hold the myeloid-context map fixed and toroidally shift the treatment-up score map over the specimen bounding lattice using deterministic seed 20260823. Exclude the zero shift and accept a shift only when at least 70% of observed spots map to an observed shifted coordinate. Generate 999 accepted shifts per specimen; if fewer than 999 unique accepted shifts exist, enumerate all accepted shifts and report the exact count.

For each shift, recompute rho on overlapping spots. The specimen empirical two-sided P value is `(1 + number(|rho_null| >= |rho_observed|)) / (1 + B)`. For GSE206245, create 999 cohort null replicates by independently selecting one accepted shift per specimen and recomputing the median of 12 specimen rhos. Compare the observed median with this null using the same plus-one rule. This analysis supplements rather than replaces patient-level inference and must be labeled as an autocorrelation-aware sensitivity analysis, because toroidal wrap-around is an approximation for irregular tissue boundaries.

## 8. Matched random-program negative controls

Use all non-dynamic SecAct signature groups as the candidate pool; exclude every group in the frozen 70-up and 93-down programs. Within each GSE206245 specimen, characterize every candidate group's standardized activity by mean abundance and SD across eligible spots. Pool these characteristics across specimens and assign decile bins for mean and SD.

Generate 1,000 random 70-group programs with seed 20260823. For each frozen treatment-up group, sample one distinct candidate from the same two-dimensional bin; when no unused candidate exists, expand by one adjacent decile in a fixed Manhattan-distance order. Reject a draw containing duplicates. Average sampled standardized group activities with equal weight, calculate each specimen's rho with the frozen myeloid-context score, and record the 12-specimen median rho. The empirical one-sided specificity P value is `(1 + number(random median rho >= observed median rho)) / 1001`.

Also report the percentile of the observed program, the full null distribution, and the proportion of random programs meeting the same positive-median/75%-positive rule. Failure to exceed matched random programs limits biological specificity and must remain in the abstract or Discussion if it materially changes the interpretation.

## 9. Missingness, multiplicity, and deviations

- No endpoint-driven imputation is permitted.
- All exclusions and non-computable specimens are listed with reasons.
- The GSE150430 treatment-up contrast is the only new confirmatory endpoint. The GSE206245 adjusted analysis and cohort spatial-shift test are prespecified robustness endpoints. GSE200310 is descriptive.
- Secondary lineage, treatment-down, direct-expression, CytoSig, and individual-cytokine results are multiplicity-aware and cannot replace a failed primary endpoint.
- Any implementation change required by released file structure will be written as a dated protocol amendment before examining the affected biological result.

## 10. Journal escalation rule

A CAS Zone 1 submission is scientifically defensible only if the independent GSE150430 primary endpoint supports myeloid localization, both GSE200310 specimens are directionally concordant, GSE206245 remains supported after composition adjustment and spatial-null testing, and the frozen program exceeds matched random programs. Even then, the absence of protein/perturbation validation and the existing direct-expression disadvantage remain explicit limitations.

If one or more of these conditions fail but the main cross-platform observations remain coherent, the manuscript will be positioned for a CAS Zone 2 journal with an atlas/methodological-contextualization claim. Journal choice will follow the completed evidence, not determine which results are retained.
