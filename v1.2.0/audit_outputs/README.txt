Supplementary audit files for the Cancer Investigation revision

These files provide the complete GSE153125 response-screen results and reproducibility materials for the candidate and spatial sensitivity analyses.

The supplied derived GSE153125 activity matrix and patient manifest are in inputs/GSE153125/. The raw public data are available from GEO as cited in the manuscript. The supplied spatial spot summaries and twelve selected activity matrices are in inputs/GSE206245/.

Contents:
- GSE153125_all_program_results.tsv: all 1,170 programs, group means/SDs, Cohen d, corrected Hedges g, Welch P/q, and exact two-sided Mann–Whitney P/q.
- candidate_response_bootstrap_audit.json and audit_candidate_bootstrap.py: patient-stratified confidence intervals and rank stability.
- candidate_spatial_depth_adjustment.tsv/.json and audit_spatial_depth.py: section-level technical sensitivity adjustment for sequencing depth and detected-gene count. This is not cell-fraction deconvolution.
- GSE62328 recurrence audit files and expanded_candidate_shortlist_v0.2.tsv: exploratory context outputs.

The full-screen audit reproduced 67 Welch BH q<0.05 programs, 12 exact Mann–Whitney BH q<0.05 programs, and 55 Welch-only findings.

Reproduction (requires Python with numpy, pandas, and scipy):
1. python audit_GSE153125_full_screen.py --zscore inputs/GSE153125/GSE153125_secact_zscore.tsv.gz --manifest inputs/GSE153125/GSE153125_patient_manifest.tsv --out GSE153125_all_program_results.tsv
2. python audit_candidate_bootstrap.py --zscore inputs/GSE153125/GSE153125_secact_zscore.tsv.gz --manifest inputs/GSE153125/GSE153125_patient_manifest.tsv --out candidate_response_bootstrap_audit.json
3. python audit_spatial_depth.py --spatial-dir inputs/GSE206245 --out-dir .

The scripts reproduce downstream statistical summaries from the supplied derived score matrices; they do not rerun SecAct from raw expression. These revision analyses and outputs are not part of the cited public v1.1.0 archive.

File checksums (SHA-256):
fb6becfecf83ae44c4b3ac08646b4681be4e0ef5404d935e2fbd9e393b65fad7  audit_candidate_bootstrap.py
fcdc51cf8a68e53a24c6ca966f239455ee44af918975f9ee933ddf07e025edc7  audit_GSE153125_full_screen.py
a129ef11dbaffc2b73df3f9749d4f38c32be80cc60f858756696f04aa82e00e8  audit_spatial_depth.py
89a763536b5138b389c296d86e0b03ba87c83f41d4701fe872d78eda7da13267  candidate_response_bootstrap_audit.json
7120e4b514541d835e205d93f53707f9abbf7edad0e3a89793ee8acb0af40460  candidate_spatial_depth_adjustment.json
3f08cd834c98dd6a01232fed556f7399edd07978c2963a39990a0899cb604f78  candidate_spatial_depth_adjustment.tsv
127979d001711a629ed060c4d90feafa79cbc8fc4b3b9aff10e2b52211856bc8  expanded_candidate_shortlist_v0.2.tsv
3e165c30b28f0a10050686efb9842b9632d15446183aded6472d2f3a8ed954df  GSE153125_all_program_results.tsv
b4858ddfc5e29533a8da0594618f3af2507b3e4a8bccb32331e4ba76de3f5d92  gse62328_candidate_recurrence_audit_v0.1.json
fc078830e790b28abb309f3c390c4b751ac846d054e462224af1207e5d62a518  gse62328_candidate_recurrence_audit_v0.1.tsv
fd2b3085e2a1e394d3a4b473b4915b89b5d96edaaac275d7963b5aa991315897  inputs/GSE153125/GSE153125_patient_manifest.tsv
3a8a69ff6a4073ff02f96f2d468dcc25310e65111fc37a457114d3c8edd41e57  inputs/GSE153125/GSE153125_secact_zscore.tsv.gz
5e82b2def2ba2fd2ad8a56483e965d14c1d36860bcad6db5272faab9f095745e  inputs/GSE206245/GSE206245_spot_scores.tsv.gz
a07930ecba71f92c41ee5b0e431ca3aff02905d3d23d082d9238e8150dc8e76d  inputs/GSE206245/ST10_selected_group_activity.tsv.gz
36a93972c201e07f0db1313cd046febf09fd79989a27ec8ad45fa9c95b55a907  inputs/GSE206245/ST11_selected_group_activity.tsv.gz
7ee9308785299e1d55eb77a4ebd26dd2562f57b963741553c91f194c9bcfb74a  inputs/GSE206245/ST12_selected_group_activity.tsv.gz
0963ae5b6cd9751a8bb50b1433515c27950f532af5b1e49028b9dc4f9fc7183b  inputs/GSE206245/ST16_selected_group_activity.tsv.gz
4acb55134dfaa4cca9df565748831f0c88dc488b41357cf492b5a08724de1921  inputs/GSE206245/ST17_selected_group_activity.tsv.gz
2172ccda66d0973dc2f38f881a8475cadfb5ac8572fe83b2ec7c2195390713c0  inputs/GSE206245/ST18_selected_group_activity.tsv.gz
08efde6335a4eb1dfd33a7dec5c2b62472339152869a8b55b988427b1c6c129c  inputs/GSE206245/ST19_selected_group_activity.tsv.gz
93f6dac2f94d3515b2f117aaf1a3b5089a9d3e169abd2d4c9f32ac3289f20cc2  inputs/GSE206245/ST5_selected_group_activity.tsv.gz
d385604fa570f5f117dc0aac19f6fc97951d7622a7060bea891a5cededde8df5  inputs/GSE206245/ST6_selected_group_activity.tsv.gz
00345b68f0f10921cb14a8da8ea65d8e3372ab5882fb99c46cb95e726e4f0648  inputs/GSE206245/ST7_selected_group_activity.tsv.gz
56f45102e379569de05ec6ff15d8c20cb300e9723eff70bf79f805934e5642ad  inputs/GSE206245/ST8_selected_group_activity.tsv.gz
c81e67da7b8bcf0188018c68605956ef60c6ee5e512a5fdd67d194712c5c50e7  inputs/GSE206245/ST9_selected_group_activity.tsv.gz


GSE150430 candidate cell-type expression context
- GSE150430_candidate_expression_analysis.py: reproducible script for the descriptive candidate/receptor expression summaries.
- inputs/GSE150430/patient_level_celltype_expression.tsv: patient-by-original-cell-type summaries for WNT7A, VSTM2L, FZD1-FZD10, LRP5/6, and ROR1/2.
- inputs/GSE150430/celltype_expression_summary.tsv: median summaries across patient-cell-type groups with at least 20 cells. This file is also provided as Supplementary_Table_S1_GSE150430_candidate_cell_expression.tsv.
- inputs/GSE150430/myeloid_subtype_expression_summary.tsv: corresponding summaries for macrophages, monocytes, DC1, DC2, DC3, and pDC.
- inputs/GSE150430/analysis_metadata.json: input, patient selection, thresholds, and interpretation limits.

Reproduction:
python GSE150430_candidate_expression_analysis.py --matrix GSE150430_npc_scRNA_hg19_processed_data.txt.gz --outdir inputs/GSE150430

The source matrix is the author-processed GSE150430 expression matrix and is not redistributed in this package. It remains available from the local project copy and the public accession. These summaries are descriptive transcript context only; they do not perform raw-count differential expression, de novo reclustering, protein quantification, or ligand-receptor interaction inference.
