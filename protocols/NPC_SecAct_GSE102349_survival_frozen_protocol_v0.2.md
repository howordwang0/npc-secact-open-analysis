# NPC SecAct GSE102349 survival-analysis frozen protocol v0.2 amendment

**Amendment date:** 2026-08-26  
**Reason:** Input-level audit of the direct-expression comparator identified
two frozen groups with no represented member gene in the released GSE102349
FPKM matrix. This amendment is based only on the expression/signature mapping,
before reviewing or selecting any endpoint association.

## Direct-expression coverage rule

For the direct ligand-expression comparator, a group is computable when at least
one member gene in its frozen `signature_group` is present in the GSE102349 gene
matrix. A group with zero represented members is excluded from that comparator
before patient-level scores are calculated; it is not imputed and it does not
change the SecAct program score. Member genes present in a multi-member group
are averaged on the log2(FPKM + 1) scale. The retained group count and every
excluded group are reported.

The coverage gate is 90% of frozen groups represented by at least one member for
each comparator. The observed input audit must therefore be recorded as 68/70
for treatment-up and 91/93 for treatment-down, if those counts reproduce.
The four-association multiplicity family remains unchanged: the two SecAct
program scores and the two direct-expression comparators. This amendment does
not permit outcome-informed group removal, threshold tuning, or candidate
rediscovery.
