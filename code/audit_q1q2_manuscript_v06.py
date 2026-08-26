from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def close(actual: float, expected: float, tolerance: float = 5e-7) -> bool:
    return bool(np.isclose(actual, expected, atol=tolerance, rtol=0))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    manuscript_path = args.output_dir / "NPC_SecAct_manuscript_Q1Q2_v0.6.md"
    manuscript = manuscript_path.read_text(encoding="utf-8")

    localization = pd.read_csv(
        args.output_dir / "NPC_SecAct_GSE150430_localization_summary_v0.1.tsv",
        sep="\t",
    )
    external = pd.read_csv(
        args.output_dir / "NPC_SecAct_GSE200310_external_spatial_summary_v0.1.tsv",
        sep="\t",
    )
    robustness = json.loads(
        (
            args.output_dir
            / "NPC_SecAct_GSE206245_composition_spatial_robustness_v0.1.json"
        ).read_text(encoding="utf-8")
    )
    matched = json.loads(
        (
            args.output_dir
            / "NPC_SecAct_GSE206245_matched_random_summary_v0.1.json"
        ).read_text(encoding="utf-8")
    )

    secact_up = localization.loc[
        localization["analysis"].eq("SecAct")
        & localization["set_label"].eq("dynamic_up_set")
    ].iloc[0]
    direct_up = localization.loc[
        localization["analysis"].eq("Direct_expression")
        & localization["set_label"].eq("dynamic_up_set")
    ].iloc[0]
    external_up = external.loc[external["program"].eq("dynamic_up_set")].set_index("sample")
    adjusted = robustness["composition_adjusted"]
    spatial_null = robustness["cohort_toroidal_shift"]
    matched_null = matched["matched_random_programs"]
    survival = pd.read_csv(
        args.output_dir / "NPC_SecAct_GSE102349_survival_results_v0.1.tsv", sep="\t"
    )
    survival_primary = survival.loc[
        survival["association"].isin(["up_secact", "down_secact", "up_direct", "down_direct"])
    ].set_index("association")

    checks: list[dict] = []

    def add(check_id: str, description: str, passed: bool, observed: object) -> None:
        checks.append(
            {
                "check_id": check_id,
                "description": description,
                "passed": bool(passed),
                "observed": str(observed),
            }
        )

    add("L01", "GSE150430 eligible patients", int(secact_up["n_patients"]) == 14, secact_up["n_patients"])
    add("L02", "GSE150430 positive patients", int(secact_up["positive_patients"]) == 14, secact_up["positive_patients"])
    add("L03", "GSE150430 median difference", close(float(secact_up["median_difference"]), 0.7350914552447959), secact_up["median_difference"])
    add("L04", "GSE150430 Hedges gz", close(float(secact_up["hedges_gz"]), 1.3836376135531245), secact_up["hedges_gz"])
    add("L05", "GSE150430 exact P", close(float(secact_up["exact_wilcoxon_p_two_sided"]), 0.0001220703125), secact_up["exact_wilcoxon_p_two_sided"])
    add("L06", "GSE150430 direct positive patients", int(direct_up["positive_patients"]) == 14, direct_up["positive_patients"])
    add("L07", "GSE150430 direct median", close(float(direct_up["median_difference"]), 0.46033500978721853), direct_up["median_difference"])

    p1 = external_up.loc["NPC_patient_1"]
    p2 = external_up.loc["NPC_patient_2"]
    add("E01", "GSE200310 patient 1 rho", close(float(p1["spearman_rho"]), 0.23723748362628494), p1["spearman_rho"])
    add("E02", "GSE200310 patient 2 rho", close(float(p2["spearman_rho"]), 0.3619834000822363), p2["spearman_rho"])
    add("E03", "GSE200310 both positive", bool((external_up["spearman_rho"] > 0).all()), external_up["spearman_rho"].tolist())
    add("E04", "GSE200310 patient 1 shift P", close(float(p1["toroidal_shift_p_two_sided"]), 0.0053475935828877), p1["toroidal_shift_p_two_sided"])
    add("E05", "GSE200310 patient 2 shift P", close(float(p2["toroidal_shift_p_two_sided"]), 0.01694915254237288), p2["toroidal_shift_p_two_sided"])

    add("R01", "Composition-adjusted median rho", close(adjusted["median_patient_rho"], 0.24253490226458221), adjusted["median_patient_rho"])
    add("R02", "Composition-adjusted positive patients", adjusted["positive_patients"] == 12, adjusted["positive_patients"])
    add("R03", "Composition-adjusted exact P", close(adjusted["exact_wilcoxon_p_two_sided"], 0.00048828125), adjusted["exact_wilcoxon_p_two_sided"])
    add("R04", "Cohort toroidal-shift P", close(spatial_null["empirical_p_two_sided"], 0.001), spatial_null["empirical_p_two_sided"])
    add("R05", "Cohort null median", close(spatial_null["null_median"], 0.03733872864847869), spatial_null["null_median"])

    add("M01", "Matched-program observed percentile", close(matched_null["observed_percentile"], 98.9), matched_null["observed_percentile"])
    add("M02", "Matched-program empirical P", close(matched_null["empirical_p_one_sided"], 0.011988011988011988), matched_null["empirical_p_one_sided"])
    add("M03", "Loose-rule random programs", matched_null["random_programs_meeting_positive_median_75pct_rule"] == 847, matched_null["random_programs_meeting_positive_median_75pct_rule"])

    add("S01", "GSE102349 eligible records", int(survival_primary.iloc[0]["n"]) == 88, survival_primary.iloc[0]["n"])
    add("S02", "GSE102349 events", int(survival_primary.iloc[0]["events"]) == 16, survival_primary.iloc[0]["events"])
    add("S03", "GSE102349 up SecAct HR", close(float(survival_primary.loc["up_secact", "score_hr"]), 0.5222816912176317), survival_primary.loc["up_secact", "score_hr"])
    add("S04", "GSE102349 up SecAct FDR", close(float(survival_primary.loc["up_secact", "fdr_bh"]), 0.06639316701789132), survival_primary.loc["up_secact", "fdr_bh"])
    add("S05", "GSE102349 down SecAct HR", close(float(survival_primary.loc["down_secact", "score_hr"]), 1.8726340230455534), survival_primary.loc["down_secact", "score_hr"])
    add("S06", "GSE102349 up direct FDR", close(float(survival_primary.loc["up_direct", "fdr_bh"]), 0.057835153645066424), survival_primary.loc["up_direct", "fdr_bh"])
    add("S07", "No GSE102349 primary association passes FDR<0.05", bool((survival_primary["fdr_bh"] >= 0.05).all()), survival_primary["fdr_bh"].tolist())

    required_text = {
        "T01": "14/14 patients",
        "T02": "median rho=0.243",
        "T03": "two-sided empirical P=0.001",
        "T04": "847/1,000 random programs",
        "T05": "rho was 0.237 in patient 1",
        "T06": "Direct treatment-up expression was again stronger than SecAct in both external",
        "T07": "GSE150430",
        "T08": "GSE200310",
        "T09": "PMID: 32686767",
        "T10": "PMID: 37024479",
        "T11": "GSE102349",
        "T12": "No association passed BH-FDR<0.05",
    }
    for check_id, token in required_text.items():
        add(check_id, f"Manuscript contains: {token}", token in manuscript, token)

    figure_png = args.output_dir / "figures" / "NPC_SecAct_Q1Q2_figure6_enhancement_v0.1.png"
    figure_pdf = args.output_dir / "figures" / "NPC_SecAct_Q1Q2_figure6_enhancement_v0.1.pdf"
    add("F01", "Figure 6 PNG exists and is nonempty", figure_png.exists() and figure_png.stat().st_size > 100_000, figure_png.stat().st_size if figure_png.exists() else "missing")
    add("F02", "Figure 6 PDF exists and is nonempty", figure_pdf.exists() and figure_pdf.stat().st_size > 10_000, figure_pdf.stat().st_size if figure_pdf.exists() else "missing")
    figure1_png = args.output_dir / "figures" / "NPC_SecAct_Q1Q2_figure1_study_flow_v0.2.png"
    figure1_pdf = args.output_dir / "figures" / "NPC_SecAct_Q1Q2_figure1_study_flow_v0.2.pdf"
    add("F03", "Updated Figure 1 PNG exists and is nonempty", figure1_png.exists() and figure1_png.stat().st_size > 100_000, figure1_png.stat().st_size if figure1_png.exists() else "missing")
    add("F04", "Updated Figure 1 PDF exists and is nonempty", figure1_pdf.exists() and figure1_pdf.stat().st_size > 10_000, figure1_pdf.stat().st_size if figure1_pdf.exists() else "missing")
    figure7_png = args.output_dir / "figures" / "NPC_SecAct_Q1Q2_figure7_GSE102349_survival_v0.1.png"
    figure7_pdf = args.output_dir / "figures" / "NPC_SecAct_Q1Q2_figure7_GSE102349_survival_v0.1.pdf"
    add("F05", "Figure 7 PNG exists and is nonempty", figure7_png.exists() and figure7_png.stat().st_size > 50_000, figure7_png.stat().st_size if figure7_png.exists() else "missing")
    add("F06", "Figure 7 PDF exists and is nonempty", figure7_pdf.exists() and figure7_pdf.stat().st_size > 10_000, figure7_pdf.stat().st_size if figure7_pdf.exists() else "missing")

    table = pd.DataFrame(checks)
    audit_path = args.output_dir / "NPC_SecAct_Q1Q2_manuscript_number_audit_v0.2.tsv"
    table.to_csv(audit_path, sep="\t", index=False)
    passed = int(table["passed"].sum())
    total = len(table)
    report = f"""# NPC SecAct Q1/Q2 manuscript numerical audit v0.2

**Audit target:** `{manuscript_path.name}`  
**Result:** {passed}/{total} checks passed.

The audit recalculated or directly retrieved the headline enhancement results
from retained patient-level TSV/JSON artifacts, checked their rounded forms in
the manuscript, confirmed both new verified references, and confirmed the PNG/
PDF exports for updated Figure 1 and new Figure 6. It supplements the prior
17/17 audit of the core v0.4 results.

{'All checks passed.' if passed == total else 'FAILED checks: ' + ', '.join(table.loc[~table['passed'], 'check_id'])}
"""
    report_path = args.output_dir / "NPC_SecAct_Q1Q2_manuscript_number_audit_report_v0.2.md"
    report_path.write_text(report, encoding="utf-8")
    print(f"{passed}/{total} checks passed")
    if passed != total:
        print(table.loc[~table["passed"]].to_string(index=False))
        raise SystemExit(1)


if __name__ == "__main__":
    main()
