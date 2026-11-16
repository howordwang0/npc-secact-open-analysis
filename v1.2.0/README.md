# npc-secact-open-analysis — revision audit v1.2.0

Audit code and derived outputs supporting the revision of:

> Exploratory multimodal prioritization of secreted-protein activity signals associated with
> treatment response in nasopharyngeal carcinoma (Cancer Investigation, submission 261212126)

This release accompanies the response to reviewers. It adds the candidate-specific analyses that
were not part of the earlier public v1.1.0 archive.

## Contents

`scripts/`
- `audit_GSE153125_full_screen.py` — full 1,170-program screen (Welch and exact two-sided Mann-Whitney P values, separate BH q values, Hedges g)
- `audit_candidate_bootstrap.py` — patient-stratified bootstrap 95% CI and 5,000-resample rank stability
- `audit_spatial_depth.py` — section-level sensitivity adjustment for sequencing depth and detected-gene count
- `GSE150430_candidate_expression_analysis.py` — cell-type and myeloid-subtype detection summaries for WNT7A, VSTM2L and WNT receptor-context genes

`audit_outputs/`
- `GSE153125_all_program_results.tsv` — all 1,170 programs, group means/SD, Cohen d, corrected Hedges g, Welch P/q, exact Mann-Whitney P/q
- `candidate_response_bootstrap_audit.json`, `candidate_spatial_depth_adjustment.tsv/.json`
- `gse62328_candidate_recurrence_audit_v0.1.tsv/.json`, `expanded_candidate_shortlist_v0.2.tsv`

`audit_outputs/GSE150430/`
- `celltype_expression_summary.tsv`, `myeloid_subtype_expression_summary.tsv`
- `patient_level_celltype_expression.tsv` (patient × cell type × gene detection counts)
- `Supplementary_Table_S1_GSE150430_candidate_cell_expression.tsv`
- `analysis_metadata.json` (patient set P02–P15, 48,584 cells in matrix, gene panel, ≥20-cell threshold)

Input matrices are not redistributed here; they derive from public GEO accessions (GSE153125,
GSE150430, GSE206245) and are documented in the manuscript.

## Reproduction

Requires Python with `numpy`, `pandas`, `scipy`.

```
python scripts/audit_GSE153125_full_screen.py --zscore GSE153125_secact_zscore.tsv.gz \
    --manifest GSE153125_patient_manifest.tsv --out GSE153125_all_program_results.tsv
python scripts/audit_candidate_bootstrap.py --zscore GSE153125_secact_zscore.tsv.gz \
    --manifest GSE153125_patient_manifest.tsv --out candidate_response_bootstrap_audit.json
python scripts/audit_spatial_depth.py --spatial-dir GSE206245 --out-dir .
python scripts/GSE150430_candidate_expression_analysis.py \
    --matrix GSE150430_npc_scRNA_hg19_processed_data.txt.gz --out-dir GSE150430
```

## Scope

The scripts reproduce downstream statistical summaries from supplied derived matrices. They do not
rerun SecAct from raw expression, and they do not include a cell-fraction deconvolution estimate.
GSE150430 values are author-processed rather than raw counts, so the cell-type analysis is
descriptive and no count-based differential-expression test is reported.
