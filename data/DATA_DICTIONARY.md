# Derived-data dictionary

All identifiers are de-identified identifiers already released by the source
repositories. No names, contact details, dates of birth, addresses, or locally
held hospital identifiers are included.

- `GSE102349_patient_scores`: four continuous standardized and raw program
  scores per released sample. Endpoint modeling results are separate.
- `Zenodo13736891_bulk_program_scores`: one treatment-up and one
  treatment-down frozen-program score per paired bulk patient-timepoint.
- `GSE150430_primary_patient_differences`: paired myeloid and pooled
  non-myeloid scores and within-patient differences.
- `GSE206245_patient_correlations`: one within-section correlation per
  specimen, scenario, and frozen program.
- `GSE206245_composition_spatial_robustness`: unadjusted and
  composition-adjusted specimen correlations plus spatial-null diagnostics.
- `GSE200310_external_spatial_summary`: descriptive correlations and spatial
  sensitivity results for two external specimens.
- `GSE289272_patient_scores`: patient-aggregated GeoMx scores; ROI-level values
  are not used as independent patients.
- `GSE299775_patient_level_program_scores`: one PanCK-positive whole-
  transcriptome ROI per NPC patient plus descriptive normal-donor scores; the
  24 NPC patients are the primary inferential units.
- `GSE299775_primary_and_benchmark_statistics`: frozen SecAct primary endpoint
  and matched direct-expression benchmark, including exact permutation tests
  and bootstrap effect intervals.
- `GSE299775_paired_bootstrap_effects`: paired bootstrap draws used to compare
  absolute SecAct and direct-expression effect magnitudes.
- `Zenodo13736891_scrna_dynamic_scores`: patient, timepoint, lineage, scenario,
  and frozen-program score.
- `GSE132112`, `GSE103611`, `Zenodo13736891_paired_results`, direct-expression,
  and CytoSig files: complete group/signature-level results retained to expose
  multiplicity and negative findings.

Column names are intentionally preserved from the analysis outputs. Statistical
definitions and eligibility rules are given in the corresponding frozen
protocols.
