# GSE206245 publication Source Data audit v0.1

**Audit date:** 2026-08-24, before any GSE206245 SecAct score was calculated.

## Verified sources

- GEO series GSE206245 and its released file manifest.
- Liu et al., *Nature Communications* (2024), PMID 39231979, PMCID
  PMC11375053, DOI 10.1038/s41467-024-52153-4.
- Publication Supplementary Data 1 (`ST cohort` sheet).
- Publication Source Data (`sFig11c` sheet).

## Spot-count reconciliation

Supplementary Data 1 reports 23,321 released Visium spots across the 12 FFPE
NPC sections. The `sFig11c` Source Data sheet contains 21,753 rows with a
sample and region label. The difference is confined to four sections.

| Sample | Supplementary Data 1 spots | `sFig11c` region rows | Difference |
|---|---:|---:|---:|
| ST5 | 2,021 | 2,021 | 0 |
| ST6 | 2,614 | 2,614 | 0 |
| ST7 | 1,939 | 1,939 | 0 |
| ST8 | 1,209 | 1,209 | 0 |
| ST9 | 2,602 | 2,069 | -533 |
| ST10 | 1,309 | 1,309 | 0 |
| ST11 | 1,272 | 1,272 | 0 |
| ST12 | 3,669 | 3,521 | -148 |
| ST16 | 1,556 | 1,556 | 0 |
| ST17 | 2,067 | 2,067 | 0 |
| ST18 | 1,977 | 1,245 | -732 |
| ST19 | 1,086 | 931 | -155 |
| **Total** | **23,321** | **21,753** | **-1,568** |

## Decision

The Source Data rows do not include spot barcodes, and the file does not state
the exact exclusion mechanism for the 1,568 absent rows. It is therefore not
possible to map the region-labeled subset back to the released matrices
without an unsupported row-order assumption.

The frozen independent analysis retains all GEO matrix barcodes that map to a
source-designated in-tissue position (expected total 23,321), with the fixed
count/feature QC sensitivity specified before scoring. Publication region
labels are not used to select spots, define the primary endpoint, or alter the
myeloid marker module. This preserves reproducibility and prevents an
undocumented author-side filter from becoming a post-hoc selection rule.
