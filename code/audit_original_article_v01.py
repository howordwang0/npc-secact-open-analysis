from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
ANALYSIS = ROOT / "work" / "gse299775_original_article"
MANUSCRIPT = OUT / "NPC_SecAct_Original_Article_v0.1.md"


def close(left: float, right: float, tolerance: float = 5e-4) -> bool:
    return bool(np.isclose(float(left), float(right), atol=tolerance, rtol=0))


def main() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8")
    stats = pd.read_csv(
        ANALYSIS / "GSE299775_primary_and_benchmark_statistics.tsv",
        sep="\t",
        index_col="program",
    )
    scores = pd.read_csv(
        ANALYSIS / "GSE299775_patient_level_program_scores.tsv", sep="\t"
    )
    run = json.loads(
        (ANALYSIS / "GSE299775_original_article_run.json").read_text(encoding="utf-8")
    )
    secact = stats.loc["SecAct_treatment_up_composite"]
    direct = stats.loc["direct_expression_treatment_up_composite"]

    rows: list[dict] = []

    def add(check_id: str, item: str, passed: bool, observed: object) -> None:
        rows.append({"check_id": check_id, "item": item, "passed": bool(passed), "observed": observed})

    npc = scores.loc[scores["is_npc"].astype(str).str.lower().eq("true")]
    add("N01", "24 independent NPC samples", len(npc) == 24 and npc["patient_id"].nunique() == 24, len(npc))
    add("N02", "MET/nonMET counts 14/10", npc["cell_type"].value_counts().to_dict() == {"MET": 14, "nonMET": 10}, npc["cell_type"].value_counts().to_dict())
    add("N03", "No NPC score missing", not npc[["secact_treatment_up_composite", "direct_treatment_up_composite"]].isna().any().any(), int(npc.isna().sum().sum()))
    add("N04", "70 SecAct groups retained", run["secact_groups_retained"] == 70, run["secact_groups_retained"])
    add("N05", "70 direct groups retained", run["direct_groups_retained"] == 70, run["direct_groups_retained"])
    add("N06", "Exact allocation count", int(secact["exact_permutations"]) == 1961256, secact["exact_permutations"])

    numeric_checks = [
        ("S01", "SecAct Hedges g", secact["hedges_g"], -0.7376501, "g=-0.738"),
        ("S02", "SecAct exact P", secact["exact_permutation_p_two_sided"], 0.0788592, "P=0.07886"),
        ("S03", "SecAct bootstrap lower", secact["hedges_g_boot_ci_low"], -1.7829607, "-1.783"),
        ("S04", "SecAct bootstrap upper", secact["hedges_g_boot_ci_high"], 0.0136077, "0.014"),
        ("D01", "Direct Hedges g", direct["hedges_g"], -1.0225519, "g=-1.023"),
        ("D02", "Direct exact P", direct["exact_permutation_p_two_sided"], 0.0188359, "P=0.01884"),
        ("D03", "Direct bootstrap lower", direct["hedges_g_boot_ci_low"], -2.2439124, "-2.244"),
        ("D04", "Direct bootstrap upper", direct["hedges_g_boot_ci_high"], -0.1996072, "-0.200"),
    ]
    for check_id, item, observed, expected, manuscript_token in numeric_checks:
        add(check_id, item, close(observed, expected) and manuscript_token in text, observed)

    delta_ci = run["absolute_effect_difference_boot_ci"]
    add("I01", "Incremental interval lower", close(delta_ci[0], -1.2415106) and "-1.242" in text, delta_ci[0])
    add("I02", "Incremental interval upper", close(delta_ci[1], 0.5982428) and "0.598" in text, delta_ci[1])
    add("I03", "Primary association not supported", run["primary_association_supported"] is False and "negative" in text.lower(), run["primary_association_supported"])
    add("I04", "SecAct increment not supported", run["secact_incremental_separation_supported"] is False and "No supported SecAct association or incremental separation" in text, run["secact_incremental_separation_supported"])

    before_references, reference_block = text.split("## Verified references", 1)
    cited: set[int] = set()
    for block in re.findall(r"\[([0-9,\- ]+)\]", before_references):
        cited.update(int(value) for value in re.findall(r"\d+", block))
    references = [int(value) for value in re.findall(r"^(\d+)\. ", reference_block, re.MULTILINE)]
    add("R01", "References numbered 1-17", references == list(range(1, 18)), references)
    add("R02", "All references cited", cited == set(range(1, 18)), sorted(cited))
    add("R03", "GSE299775 source PMID and DOI verified in manuscript", "PMID: 40759888" in text and "10.1038/s41467-025-61597-1" in text, "40759888; 10.1038/s41467-025-61597-1")

    forbidden = ["TO BE COMPLETED", "AUTHOR CONFIRMATION", "requires approval by all authors", "manuscript structuring"]
    add("T01", "No unresolved declaration placeholders", not any(token in text for token in forbidden), [token for token in forbidden if token in text])
    add("T02", "Confirmed concise AI statement", "OpenAI Codex assisted with code development and English-language editing." in text, "present" if "OpenAI Codex assisted with code development and English-language editing." in text else "missing")
    add("T03", "GSE299775 in availability statement", "`GSE299775`" in text and "next versioned repository release" in text, "present")

    required_files = [
        OUT / "NPC_SecAct_Original_Article_GSE299775_frozen_protocol_v0.1.md",
        OUT / "NPC_SecAct_GSE299775_Original_Article_report_v0.1.md",
        OUT / "figures" / "NPC_SecAct_Original_Article_figure1_study_flow_v0.1.png",
        OUT / "figures" / "NPC_SecAct_Original_Article_figure1_study_flow_v0.1.pdf",
        OUT / "figures" / "NPC_SecAct_GSE299775_Original_Article_figure_v0.1.png",
        OUT / "figures" / "NPC_SecAct_GSE299775_Original_Article_figure_v0.1.pdf",
    ]
    for index, path in enumerate(required_files, start=1):
        add(f"F{index:02d}", f"Required artifact {path.name}", path.exists() and path.stat().st_size > 0, path.stat().st_size if path.exists() else "missing")

    abstract = text.split("## Abstract", 1)[1].split("## Introduction", 1)[0]
    main_text = text.split("## Introduction", 1)[1].split("## Verified references", 1)[0]
    word_pattern = re.compile(r"\b[A-Za-z0-9][A-Za-z0-9_.-]*\b")
    abstract_words = len(word_pattern.findall(abstract))
    main_words = len(word_pattern.findall(main_text))
    add("W01", "Abstract word count", 200 <= abstract_words <= 350, abstract_words)
    add("W02", "Main-text word count", 4000 <= main_words <= 9000, main_words)

    result = pd.DataFrame(rows)
    result.to_csv(OUT / "NPC_SecAct_Original_Article_number_audit_v0.1.tsv", sep="\t", index=False)
    failed = result.loc[~result["passed"]]
    report = [
        "# Original Article v0.1 numerical and submission audit",
        "",
        f"- Checks: {len(result)}",
        f"- Passed: {int(result['passed'].sum())}",
        f"- Failed: {len(failed)}",
        f"- Abstract words: {abstract_words}",
        f"- Main-text words (Introduction through declarations/tables): {main_words}",
        "",
    ]
    if failed.empty:
        report.append("All prespecified checks passed.")
    else:
        report.extend(["## Failed checks", ""])
        for row in failed.itertuples(index=False):
            report.append(f"- {row.check_id}: {row.item}; observed={row.observed}")
    (OUT / "NPC_SecAct_Original_Article_number_audit_report_v0.1.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "checks": len(result),
        "passed": int(result["passed"].sum()),
        "failed": int(len(failed)),
        "abstract_words": abstract_words,
        "main_words": main_words,
        "failed_ids": failed["check_id"].tolist(),
    }, indent=2))
    if not failed.empty:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
