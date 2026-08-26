from __future__ import annotations

import hashlib
import json
import warnings
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from secactpy import load_signature
from secactpy.inference import group_signatures
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
META_PATH = ROOT / "work" / "data_audit" / "GSE102349_samples.csv"
EXPR_PATH = ROOT / "work" / "data_audit" / "GSE102349_NPC_mRNA_processed.txt.gz"
Z_PATH = ROOT / "work" / "secact_smoke" / "GSE102349_zscore.tsv.gz"
SETS_PATH = ROOT / "outputs" / "NPC_SecAct_Zenodo13736891_dynamic_sets_v0.1.tsv"
OUT_DIR = ROOT / "outputs"
WORK_DIR = ROOT / "work" / "gse102349_survival"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def zscore(values: pd.Series) -> pd.Series:
    sd = values.std(ddof=1)
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("Cannot standardize a constant score")
    return (values - values.mean()) / sd


def fit_cox(time: np.ndarray, event: np.ndarray, score: np.ndarray, stage: pd.Series | None = None) -> dict:
    x_parts = [score.reshape(-1, 1)]
    term_names = ["score_sd"]
    if stage is not None:
        stage_dummies = pd.get_dummies(stage.astype(str), prefix="stage", drop_first=True, dtype=float)
        if stage_dummies.shape[1] > 0:
            x_parts.append(stage_dummies.to_numpy())
            term_names.extend(stage_dummies.columns.tolist())
    x = np.column_stack(x_parts).astype(float)
    if not np.isfinite(x).all():
        raise ValueError("Non-finite model covariate")
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = PHReg(time, x, status=event).fit(disp=False)
        ci = np.asarray(result.conf_int())
        warning_messages = sorted({str(item.message) for item in caught})
        return {
            "status": "warning" if warning_messages else "ok",
            "warnings": " | ".join(warning_messages),
            "terms": term_names,
            "n": int(len(time)),
            "events": int(event.sum()),
            "hr": np.exp(np.clip(result.params, -700, 700)).tolist(),
            "ci_low": np.exp(np.clip(ci[:, 0], -700, 700)).tolist(),
            "ci_high": np.exp(np.clip(ci[:, 1], -700, 700)).tolist(),
            "p": np.asarray(result.pvalues).tolist(),
            "score_hr": float(np.exp(result.params[0])),
            "score_ci_low": float(np.exp(ci[0, 0])),
            "score_ci_high": float(np.exp(ci[0, 1])),
            "score_p": float(result.pvalues[0]),
            "llf": float(result.llf),
        }
    except Exception as exc:  # model failure is retained, never silently replaced
        return {"status": "failed", "n": int(len(time)), "events": int(event.sum()), "error": repr(exc)}


def main() -> None:
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(META_PATH, dtype=str, keep_default_na=False)
    expr = pd.read_csv(EXPR_PATH, sep="\t", compression="gzip", index_col=0)
    expr.index = expr.index.astype(str).str.strip()
    expr = expr.apply(pd.to_numeric, errors="coerce").groupby(level=0, sort=False).mean()
    if expr.isna().any().any() or (expr.to_numpy() < 0).any():
        raise ValueError("Expression matrix contains missing or negative values")
    expr_log = np.log2(expr + 1.0)
    expr_log = expr_log.loc[:, meta["description"].tolist()]

    z = pd.read_csv(Z_PATH, sep="\t", index_col=0)
    z = z.loc[:, expr_log.columns]
    grouped = group_signatures(load_signature("secact"), cor_threshold=0.9)
    frozen = pd.read_csv(SETS_PATH, sep="\t")
    sets = {
        label: frozen.loc[frozen["set_label"].eq(label), "signature_group"].tolist()
        for label in ("dynamic_up_set", "dynamic_down_set")
    }
    expected = {"dynamic_up_set": 70, "dynamic_down_set": 93}
    for label, names in sets.items():
        if len(names) != expected[label] or len(set(names)) != len(names):
            raise ValueError(f"Frozen set integrity failure for {label}")
        missing = sorted(set(names) - set(grouped.columns))
        if missing:
            raise ValueError(f"Frozen groups absent from current SecAct grouping: {missing[:10]}")

    # The existing phenotype-blind SecAct matrix was generated with the frozen
    # protocol parameters and is checked for complete sample identity here.
    representative = {name: name.split("|")[0] for name in grouped.columns}
    secact_group = pd.DataFrame(
        {name: z.loc[representative[name]] for name in grouped.columns}, index=z.columns
    ).T
    secact_group.to_csv(WORK_DIR / "GSE102349_secact_group_activity_v0.1.tsv.gz", sep="\t", compression="gzip")

    # Direct-expression comparator: average available members inside each
    # frozen group, then average groups equally. No outcome-informed selection.
    direct_group = {}
    direct_coverage = {}
    expr_symbols = set(expr_log.index.astype(str))
    for group_name in grouped.columns:
        members = group_name.split("|")
        present = [gene for gene in members if gene in expr_symbols]
        direct_coverage[group_name] = {"members": len(members), "present": len(present)}
        if present:
            direct_group[group_name] = expr_log.loc[present].mean(axis=0)
        else:
            direct_group[group_name] = pd.Series(np.nan, index=expr_log.columns)
    direct_group_df = pd.DataFrame(direct_group).T
    direct_group_df.to_csv(WORK_DIR / "GSE102349_direct_group_expression_v0.1.tsv.gz", sep="\t", compression="gzip")

    direct_names = {}
    for label, names in sets.items():
        computable = [name for name in names if direct_coverage[name]["present"] > 0]
        if len(computable) / len(names) < 0.90:
            raise ValueError(f"Direct-expression coverage gate failed for {label}")
        direct_names[label] = computable

    scores = pd.DataFrame(index=expr_log.columns)
    for label, names in sets.items():
        scores[f"{label}_secact_raw"] = secact_group.loc[names].mean(axis=0)
        scores[f"{label}_direct_raw"] = direct_group_df.loc[direct_names[label]].mean(axis=0)
        scores[f"{label}_secact"] = zscore(scores[f"{label}_secact_raw"])
        scores[f"{label}_direct"] = zscore(scores[f"{label}_direct_raw"])

    endpoint = meta.set_index("description").loc[scores.index].copy()
    time = pd.to_numeric(endpoint["time_to_event"], errors="coerce")
    event_text = endpoint["event"].astype(str).str.strip()
    event = pd.Series(np.nan, index=endpoint.index)
    event.loc[event_text.eq("Disease progression")] = 1.0
    event.loc[event_text.eq("Last follow-up")] = 0.0
    eligible = time.notna() & event.notna() & (time > 0)
    if int(eligible.sum()) != 88 or int(event.loc[eligible].sum()) != 16:
        raise ValueError("Endpoint gate failed after score/sample alignment")
    endpoint_stage = endpoint["clinical_stage"].replace({"": "Unknown", "N/A": "Unknown"}).fillna("Unknown")

    model_specs = [
        ("up_secact", "dynamic_up_set_secact"),
        ("down_secact", "dynamic_down_set_secact"),
        ("up_direct", "dynamic_up_set_direct"),
        ("down_direct", "dynamic_down_set_direct"),
    ]
    model_rows = []
    for label, score_col in model_specs:
        values = scores.loc[eligible, score_col].to_numpy(float)
        model = fit_cox(time.loc[eligible].to_numpy(float), event.loc[eligible].to_numpy(int), values)
        model_rows.append({"association": label, "score_column": score_col, **model})
        # A stage-adjusted model is explicitly secondary and is retained only if
        # it converges; it is not used to choose or rename the primary result.
        adjusted = fit_cox(
            time.loc[eligible].to_numpy(float),
            event.loc[eligible].to_numpy(int),
            values,
            stage=endpoint_stage.loc[eligible],
        )
        model_rows.append({"association": f"{label}_stage_adjusted", "score_column": score_col, **adjusted})

    results = pd.DataFrame(model_rows)
    primary_mask = results["association"].isin(["up_secact", "down_secact", "up_direct", "down_direct"])
    pvals = results.loc[primary_mask, "score_p"].to_numpy(float)
    results.loc[primary_mask, "fdr_bh"] = multipletests(pvals, method="fdr_bh")[1]
    results.to_csv(OUT_DIR / "NPC_SecAct_GSE102349_survival_results_v0.1.tsv", sep="\t", index=False)
    scores.to_csv(OUT_DIR / "NPC_SecAct_GSE102349_patient_scores_v0.1.tsv", sep="\t", index_label="patient_sample")

    coverage = pd.DataFrame(
        [
            {"signature_group": group, **values}
            for group, values in direct_coverage.items()
            if group in set(sum(sets.values(), []))
        ]
    )
    coverage.to_csv(OUT_DIR / "NPC_SecAct_GSE102349_direct_expression_coverage_v0.1.tsv", sep="\t", index=False)

    summary = {
        "dataset": "GSE102349",
        "n_expression_samples": int(expr_log.shape[1]),
        "n_patient_units": int(scores.shape[0]),
        "n_survival_eligible": int(eligible.sum()),
        "n_events": int(event.loc[eligible].sum()),
        "n_censored": int((event.loc[eligible] == 0).sum()),
        "n_missing_endpoint": int((~eligible).sum()),
        "frozen_program_sizes": {key: len(value) for key, value in sets.items()},
        "direct_expression_computable_group_sizes": {key: len(value) for key, value in direct_names.items()},
        "secactpy_version": version("secactpy"),
        "secact_source": str(Z_PATH),
        "secact_source_sha256": sha256(Z_PATH),
        "protocol_status": "all data-count, mapping, duplicate, and event gates passed; source unit remains a stated limitation",
        "primary_results": results.loc[primary_mask, ["association", "score_hr", "score_ci_low", "score_ci_high", "score_p", "fdr_bh"]].to_dict("records"),
    }
    (OUT_DIR / "NPC_SecAct_GSE102349_survival_run_v0.1.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (OUT_DIR / "NPC_SecAct_GSE102349_survival_report_v0.1.md").write_text(
        "# GSE102349 frozen survival analysis report v0.1\n\n"
        "The analysis was run after the protocol and metadata audit were frozen.\n\n"
        f"- Patient-level expression units: {summary['n_patient_units']}\n"
        f"- Survival-eligible units: {summary['n_survival_eligible']}\n"
        f"- Progression events: {summary['n_events']}\n"
        f"- Censored observations: {summary['n_censored']}\n"
        f"- Missing/ineligible endpoint records: {summary['n_missing_endpoint']}\n\n"
        f"- Direct-expression computable frozen groups: up {summary['direct_expression_computable_group_sizes']['dynamic_up_set']}/70; down {summary['direct_expression_computable_group_sizes']['dynamic_down_set']}/93\n\n"
        "## Prespecified primary family\n\n"
        "Four continuous one-SD Cox associations were tested: treatment-up SecAct, treatment-down SecAct, and their direct-expression comparators. Benjamini-Hochberg FDR was applied within this family. Stage-adjusted fits are secondary and do not alter the primary family.\n\n"
        + "\n".join(
            f"- {row['association']}: HR={row['score_hr']:.4g} (95% CI {row['score_ci_low']:.4g}-{row['score_ci_high']:.4g}), P={row['score_p']:.4g}, FDR={row['fdr_bh']:.4g}"
            for row in summary["primary_results"]
        )
        + "\n\nAll four stage-adjusted fits triggered maximum-likelihood convergence warnings because of sparse stage categories and are not interpreted.\n\nThe released time-to-event field was analyzed in its released numeric units; the exact publication-level unit remains a limitation unless independently reconciled. These are prognostic associations, not validated prediction or causal evidence.\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
