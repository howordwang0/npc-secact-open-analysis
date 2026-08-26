# GSE102349 protocol audit report v0.1

Audit performed before any molecular-to-endpoint association.

- Metadata rows: 113
- Expression columns: 113
- Unique patient IDs: 113
- Numeric time values: 88
- Progression events: 16
- Censored last-follow-up records: 72
- Missing endpoint records: 25
- Clinical stage distribution: {"III": 41, "Unknown": 40, "IV": 25, "I": 5, "II": 2}
- Cohort field distribution: {"Paired": 55, "Unpaired": 50, "N/A": 8}
- Paired-field distribution: {"Paired": 55, "Unpaired": 50, "N/A": 8}

## Gate status

- One-to-one expression/patient mapping: True
- Exact endpoint-count reproduction: True
- No duplicate patient IDs: True
- At least 10 progression events: True

The GEO metadata do not prove that underlying patients or institutions are independent of every prior cohort; this remains an explicit limitation.
