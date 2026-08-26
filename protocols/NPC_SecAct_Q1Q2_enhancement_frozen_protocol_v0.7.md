# NPC SecAct Q1/Q2 enhancement frozen protocol v0.7 clarification

**Amendment date:** 2026-08-25  
**Timing:** Written before exporting the first complete 657-group GSE206245 activity matrix and before calculating any matched-random-program result.  
**Inheritance:** The matched-program purpose, 1,000 draws, seed, set size, candidate exclusions, distinct-group rule, endpoint, and empirical P rule from v0.1 remain unchanged.

## Activity-distribution matching clarification

In v0.1, “standardized activity” refers to SecAct's inferred group-level z-score activity before the additional within-section centering/scaling used to construct program scores. It does not mean a vector already forced to mean zero and SD one within each specimen, which would make mean/SD matching non-informative.

For each of all 657 groups, calculate within each GSE206245 specimen the mean and sample SD of the exported SecAct group-level z-score activity. Define the pooled matching characteristics as the median of the 12 specimen means and the median of the 12 specimen log-SDs. Assign separate deciles to those two pooled characteristics using deterministic rank-based `qcut` with duplicate boundaries dropped.

For each frozen treatment-up group and each random draw, select without replacement from non-dynamic candidate groups. Candidate distance is the Manhattan distance between the two decile coordinates; sample uniformly from the closest distance stratum that still contains unused candidates. This implements v0.1's adjacent-bin expansion deterministically while remaining computable when exact two-dimensional bins are sparse.

After selecting 70 distinct groups, standardize each selected group's activity within each specimen using sample SD, average the 70 groups with equal weight, and correlate the random score with the frozen 20-gene myeloid-context score. The 12-specimen median rho is the random-program statistic. No observed correlation is used in matching or sampling.
