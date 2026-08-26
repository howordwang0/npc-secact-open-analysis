# Reproducibility instructions

## 1. Environment

Python 3.12 was used. Create separate environments when necessary because the
single-cell and core analysis snapshots are recorded independently:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r environment\requirements-secact.txt
```

SecActpy 0.3.1 and its bundled signature are required. CytoSig files must be
obtained under their source terms and placed as expected by
`code/run_cytosig_benchmark.py`.

## 2. Source data

Download source files from the GEO accessions and Zenodo record listed in
`metadata/NPC_SecAct_dataset_manifest_v0.3.csv`. Do not substitute updated or
reprocessed matrices without recording a protocol deviation. Verify source
hashes where supplied in `metadata/NPC_SecAct_software_and_hash_manifest_v0.3.json`.

The scripts were written with the repository root containing sibling `code/`,
`outputs/`, and `work/` directories. Source matrices belong under `work/` in
the accession-specific paths referenced near the top of each script. Derived
tables are written to `outputs/`. The released `data/derived/` files allow
result-level verification without redistributing large source matrices.

## 3. Frozen analysis sequence

1. Run source and mapping audits before phenotype analysis.
2. Run GSE132112 and GSE103611 clinical screens.
3. Derive the treatment-up and treatment-down groups in the nine-patient
   paired bulk cohort.
4. Run paired single-cell pseudobulk analysis and independent GSE150430
   localization.
5. Run GeoMx and GSE206245 Visium analyses.
6. Run composition, toroidal-shift, matched-program, direct-expression, and
   CytoSig controls.
7. Run GSE200310 descriptive spatial replication.
8. Audit and run the GSE102349 time-to-event extension.
9. Generate figures and run the manuscript number audit.

Exact gates, seeds, eligibility rules, amendments, and stopping rules are in
`protocols/`. Negative and non-computable branches must be retained.

## 4. Expected principal checks

- GSE150430: 14 eligible NPC patients; treatment-up myeloid-minus-non-myeloid
  difference positive in 14/14.
- GSE206245: 12 specimens; treatment-up median rho 0.343 and adjusted median
  rho 0.243.
- Matched programs: observed treatment-up program at the 98.9th percentile of
  1,000 controls.
- GSE102349: 88 evaluable records and 16 progression events; none of four
  associations passes BH-FDR < 0.05.

Use `SHA256SUMS.txt` to verify the released package before rerunning analyses.
