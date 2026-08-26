# NPC SecAct GSE102349 survival-analysis frozen protocol v0.1

**Amendment date:** 2026-08-26  
**Purpose:** Pre-specified feasibility audit and, only if all gates pass, a
patient-level external public-cohort prognostic association analysis.  This
document is frozen before any SecAct score is compared with the clinical
endpoint.

## 1. Data and provenance

- Dataset: GEO GSE102349, treatment-naive NPC tumor transcriptomes.
- Expression input: the released gene-level FPKM matrix
  `GSE102349_NPC_mRNA_processed.txt.gz`.
- Metadata input: the released sample annotation table audited from GEO.
- Source publication: Zhang et al., *Molecular Cancer Research* (2017), PMID
  28851814, DOI 10.1158/1541-7786.MCR-17-0134.  The original publication
  will be checked for the exact time-to-event definition and unit before the
  result is described as PFS or another named endpoint.
- No outcome values will be imputed.  No controlled-access or hospital data
  are required for this analysis.

## 2. Statistical unit and eligibility

The statistical unit is the patient, never an array column, technical replicate,
spot, ROI, or cell.  A sample is eligible only when its GEO accession maps to a
single patient identifier, its expression column is present in the released
matrix, and its event/time fields can be audited.  Patient identifiers are
parsed from the metadata title/description and cross-checked against GEO
accessions.

If multiple samples map to one patient, they will be handled before outcome
association as follows: samples with conflicting endpoint values or
non-equivalent clinical material will be excluded from the primary analysis;
otherwise expression values will be averaged on the log2 scale and the patient
will contribute once.  The number and reason for every exclusion will be
reported.  The paired/unpaired field is descriptive metadata and does not create
additional patients or justify paired survival statistics.

## 3. Primary endpoint and censoring

- Primary endpoint: the released patient-level numeric `time_to_event` value,
  with the time unit reported exactly as verified from the source publication.
- Event coding: `event = Disease progression` is coded 1 (event); `event =
  Last follow-up` is coded 0 (right-censored).
- Other labels, including `N/A`, `Not evaluated`, blank, or nonnumeric time
  values, are treated as missing and excluded without imputation.
- The feasibility audit must reproduce 113 expression samples, 88 numeric
  time values, 16 progression events, 72 censored observations, and 25
  missing endpoint records.  Any discrepancy stops the inferential analysis
  until resolved.

## 4. Frozen molecular exposures

The exposures are fixed before endpoint association:

1. The 70 treatment-up SecAct groups and 93 treatment-down SecAct groups in
   `outputs/NPC_SecAct_Zenodo13736891_dynamic_sets_v0.1.tsv`.
2. SecAct activity is inferred with the project SecActPy environment, using
   the released SecAct signature, log2(FPKM + 1), `n_rand=100`, NumPy RNG,
   seed 20260823, and no phenotype information.
3. Group activity is the activity of the frozen representative protein for
   each precomputed correlated group.  Program scores are the equal-weight
   mean of the corresponding group activities, standardized within the
   analyzed patient set only after the activity matrix has been generated.
4. A direct ligand-expression comparator is calculated from the same frozen
   group members: within each group, member genes are averaged after log2(FPKM
   + 1), then the 70-group or 93-group equal-weight mean is formed.  This is a
   comparator, not a replacement endpoint and not a new discovery screen.

No protein, group, gene, cutoff, or direction may be selected using GSE102349
outcomes.  No single-protein claim will be made from this analysis.

## 5. Primary and secondary models

- Primary association: continuous standardized treatment-up program score in a
  Cox proportional-hazards model with the released time/event fields.
- Secondary associations: continuous treatment-down score and the direct
  expression comparators for up and down programs.
- Effect measure: hazard ratio per one standard-deviation increase, 95% CI,
  Wald P value, and proportional-hazards diagnostic where computable.
- Clinical baseline: stage (I-IV/Unknown) may be included as one prespecified
  categorical covariate only if its distribution and missingness are stable in
  the audit.  No model will contain more than one clinical covariate because
  there are only 16 events.  If the stage model is unstable or singular, the
  unadjusted continuous association remains descriptive only.
- With 16 events, no machine-learning prediction, LASSO, best-cutoff search,
  time-dependent ROC claim, nomogram, or complex multivariable model is
  permitted.

## 6. Multiplicity and interpretation

The four prespecified molecular associations (up SecAct, down SecAct, up direct
expression, down direct expression) are the multiplicity family.  Benjamini-
Hochberg FDR is reported for this family; two-sided FDR < 0.05 is the nominal
threshold for a statistically supported association.  Results are described as
prognostic associations, not validated clinical prediction, causal effects,
secreted-protein proof, or evidence that SecAct is superior to expression.

## 7. Feasibility and stop gates

The clinical extension proceeds only if all gates pass:

1. One-to-one expression-column/sample mapping and a defensible patient-level
   collapse decision.
2. Exact reproduction of endpoint counts and explicit missingness handling.
3. No unresolved duplicate patients or conflicting endpoint records.
4. At least 10 progression events after eligibility filtering and a finite,
   computable primary model.
5. Source publication verification of the endpoint definition/unit, or an
   explicit limitation stating that only the released time-to-event field is
   analyzed.

If any gate fails, the survival route is stopped and the existing Q1/Q2
manuscript is retained unchanged apart from a transparent feasibility/no-go
note.  A positive-looking result cannot override a failed gate.

## 8. Reproducibility

The audit script, frozen analysis script, input hashes, model outputs, and
number audit will be saved under `work/` and `outputs/`.  All analyses will be
rerun from the frozen files before any manuscript amendment.
