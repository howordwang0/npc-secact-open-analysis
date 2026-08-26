from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from secactpy import secact_activity_inference

from run_gse132112_discovery import compare_groups, group_level_matrix, ligand_expression_comparator


SEED = 20260823


def main() -> None:
    parser = argparse.ArgumentParser(description="Run pre-specified GSE103611 metastasis discovery.")
    parser.add_argument("expression_tsv_gz", type=Path)
    parser.add_argument("manifest_csv", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--n-rand", type=int, default=1000)
    parser.add_argument("--bootstrap", type=int, default=2000)
    args = parser.parse_args()

    expression = pd.read_csv(args.expression_tsv_gz, sep="\t", index_col=0)
    manifest = pd.read_csv(args.manifest_csv)
    sample_order = manifest["expression_sample_id"].tolist()
    if set(sample_order) != set(expression.columns):
        raise ValueError("Manifest and expression samples differ")
    expression = expression.loc[:, sample_order]
    result = secact_activity_inference(
        expression,
        is_differential=False,
        sig_matrix="secact",
        n_rand=args.n_rand,
        seed=SEED,
        rng_method="gsl",
        sort_genes=True,
        verbose=True,
    )
    activity = result["zscore"]
    group_activity, group_members = group_level_matrix(activity)
    labels = manifest.set_index("expression_sample_id")["metastasis_binary"].astype(int)
    stats = compare_groups(group_activity, labels, "metastasis", "non-metastasis", args.bootstrap, 10)
    stats["candidate_rule_pass"] = (
        (stats["welch_fdr_bh"] < 0.10)
        & (stats["hedges_g"].abs() >= 0.80)
        & (stats["bootstrap_same_direction_fraction"] >= 0.90)
    )
    passing = stats.index[stats["candidate_rule_pass"]].tolist()
    frozen = stats.loc[passing].sort_values(["welch_fdr_bh", "welch_p"]).head(10).index.tolist()
    stats["frozen_for_gse299775"] = stats.index.isin(frozen)

    ligand = ligand_expression_comparator(expression, group_members, sample_order)
    ligand_stats = compare_groups(ligand, labels, "metastasis", "non-metastasis", args.bootstrap, 11)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name in ("beta", "se", "zscore", "pvalue"):
        result[name].to_csv(args.output_dir / f"GSE103611_secact_{name}.tsv.gz", sep="\t", compression="gzip")
    group_activity.to_csv(args.output_dir / "GSE103611_secact_group_activity.tsv.gz", sep="\t", compression="gzip")
    group_members.to_csv(args.output_dir / "GSE103611_secact_group_members.tsv", sep="\t", index=False)
    stats.to_csv(args.output_dir / "GSE103611_metastasis_discovery.tsv", sep="\t")
    ligand_stats.to_csv(args.output_dir / "GSE103611_ligand_expression.tsv", sep="\t")
    pd.DataFrame({"signature_group": frozen}).to_csv(args.output_dir / "GSE103611_frozen_candidates.tsv", sep="\t", index=False)

    summary = {
        "dataset": "GSE103611",
        "patients": 48,
        "metastasis": 24,
        "non_metastasis": 24,
        "signature_groups_tested": int(len(stats)),
        "n_rand": args.n_rand,
        "rng_method": "gsl",
        "seed": SEED,
        "bootstrap_iterations": args.bootstrap,
        "candidate_rule": "BH FDR < 0.10, absolute Hedges g >= 0.80, bootstrap same-direction fraction >= 0.90",
        "candidate_rule_passed": int(stats["candidate_rule_pass"].sum()),
        "frozen_candidate_groups": frozen,
        "sample_order_warning": "Metastasis labels occupy columns 1-24 and non-metastasis labels columns 25-48; GEO has no usable batch field",
    }
    (args.output_dir / "GSE103611_metastasis_run.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
