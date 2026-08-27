# NPC SecAct open-data analysis

Reproducibility package for the Original Article **Patient-level multimodal
benchmarking reveals reproducible myeloid localization but limited incremental
value of secreted-protein activity inference in nasopharyngeal carcinoma**.

Archived release DOI: https://doi.org/10.5281/zenodo.22110691

## Scope

This repository contains analysis code, frozen protocols, de-identified
patient/specimen-level derived results, source-data audits, software manifests,
and submitted vector figures. It does not redistribute large source matrices,
raw sequencing data, controlled-access data, or the SecAct software package.

The analyses use fully open data from GEO accessions GSE102349, GSE132112,
GSE103611, GSE150430, GSE200310, GSE206245, GSE289272, and GSE299775, plus
Zenodo record 13736891. Source files should be downloaded from the repositories named in
`metadata/NPC_SecAct_dataset_manifest_v0.3.csv`.

## Directory structure

- `code/`: analysis, audit, and figure-generation scripts.
- `data/derived/`: machine-readable patient/specimen-level results and locked
  signature-group outputs used in the manuscript.
- `protocols/`: frozen protocols and dated amendments.
- `metadata/`: accession roles, source audits, hashes, and software versions.
- `figures/`: vector PDF manuscript figures.
- `environment/`: exact Python package snapshots.

## Reproducibility levels

1. **Result verification:** inspect the patient/specimen-level tables and run
   the manuscript audit using the retained derived outputs.
2. **Analysis rerun:** download accession-specific source files, arrange them
   as described in `REPRODUCIBILITY.md`, create the matching environment, and
   run the scripts in the frozen order.

Cells, spots, and regions are never treated as independent patients for the
manuscript's population-level claims. The repository does not support claims
of measured protein secretion, a proven cellular source, treatment causality,
SecAct superiority, incremental clinical value, or a validated clinical
predictor. Version 1.1.0 adds the independently frozen GSE299775 PanCK-positive
metastasis extension and its matched direct-expression benchmark.

## Citation and license

Citation metadata are provided in `CITATION.cff`. Cite the archived v1.0.0
release using https://doi.org/10.5281/zenodo.22110691. Code is released under
the MIT License. Source datasets and third-party software retain their original
licenses and terms.
