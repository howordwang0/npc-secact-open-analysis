# GSE289272 GeoMx spatial contextualization frozen protocol v0.1

**Freeze point:** DCC structure, PKC structure, ROI metadata, patient counts,
and panel coverage have been audited. No dynamic-set expression score,
myeloid score, correlation, or tissue contrast has been calculated.

## Purpose

Determine whether the bulk-frozen treatment-up set has a patient-level spatial
expression context consistent with the directional myeloid localization seen
in the paired scRNA-seq analysis. This dataset has no paired treatment
timepoints and cannot validate treatment dynamics.

## Units and eligibility

- The dataset contains 55 TLS ROIs: 40 from nine NPC patients and 15 from
  three inflammatory-control patients.
- The ROI is a technical/spatial sampling unit. The patient is the biological
  unit for inference.
- Main analysis retains all source-released ROIs because all source `QCFlags`
  are blank.
- Sensitivity analysis excludes ROIs failing current GeoMxTools default
  thresholds: raw reads >=1,000; trimmed, stitched, and aligned fractions
  >=80%; sequencing saturation >=50%; negative-probe geometric mean >=10;
  nuclei >=200; area >=16,000 square micrometers.

## Probe processing and normalization

1. Map each DCC `RTS_ID` through the supplied CTA v1.0 PKC.
2. Flag globally low-ratio probes using the GeoMxTools rule: probe geometric
   mean across ROIs divided by the corresponding target geometric mean <=0.1.
3. Exclude low-ratio probes and aggregate remaining probes per target using
   the GeoMxTools geometric mean convention (values below 0.5 are shifted by
   0.5 before the geometric mean).
4. Normalize endogenous target counts by Q3: divide each ROI by its target
   75th percentile relative to the geometric mean Q3 factor, then use
   `log2(Q3-normalized count + 1)`.
5. Median probe aggregation is retained as a technical sensitivity analysis.

## SecAct coverage gate and fallback

- CTA contains 1,779 endogenous targets and covers only 745/7,919 SecAct
  signature genes (9.408%). Full SecAct inference is prohibited.
- Direct-expression fallback covers 35/70 frozen bulk-up signature groups and
  25/93 frozen bulk-down groups.
- For each covered signature group, standardize each observed member gene
  across ROIs, average members within group, then average groups within the
  frozen set. Groups, not individual panel genes, receive equal weight.
- Scores are explicitly partial direct-expression modules, not SecAct
  activity estimates.

## Frozen spatial endpoint

The primary endpoint is the Spearman association, across the nine NPC
patients, between the patient-mean partial bulk-up direct-expression score and
an independent myeloid marker score. The marker set is frozen as `AIF1`,
`TYROBP`, `FCER1G`, `CTSS`, `CST3`, `S100A8`, `S100A9`, `FCN1`, `CD14`,
`ITGAM`, `CSF1R`, `CD68`, `CD163`, and `LST1`, after removing any overlap with
the covered bulk-up genes.

Directional spatial support requires all of the following:

- Spearman rho >0 in NPC patient means;
- all leave-one-patient-out rho values remain >0;
- rho remains >0 after default GeoMx ROI QC;
- rho remains >0 with median rather than geometric-mean probe aggregation.

Exact or permutation P values and NPC-versus-inflammatory-control contrasts
are secondary context. They do not override the directional rule. ROI-level
correlations are descriptive and cannot be used as inferential sample sizes.

## Interpretation limits

- All ROIs are TLS-selected, so findings do not represent the whole tumor.
- The panel is targeted and covers only half of the frozen groups.
- There are only nine NPC patients and three inflammatory controls.
- Positive association would support a myeloid-rich TLS context, not prove
  myeloid secretion, treatment response, or causal mechanism.
