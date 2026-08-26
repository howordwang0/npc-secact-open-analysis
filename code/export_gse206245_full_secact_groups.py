from __future__ import annotations

import argparse
import gc
import json
from importlib.metadata import version
from pathlib import Path

import pandas as pd
from secactpy import load_signature, secact_activity_inference

from analyze_gse206245_visium_spatial import (
    SAMPLES,
    enable_process_gsl_permutation_cache,
    group_level_matrix,
    read_slice,
)


SEED = 20260823


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("data_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument(
        "samples",
        nargs="+",
        choices=[metadata["sample"] for metadata in SAMPLES],
    )
    parser.add_argument("--n-rand", type=int, default=1000)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    enable_process_gsl_permutation_cache()
    signature_genes = set(load_signature("secact").index.astype(str).str.upper())
    for sample in args.samples:
        metadata = next(item for item in SAMPLES if item["sample"] == sample)
        output = args.output_dir / f"{sample}_all_group_activity.tsv.gz"
        run_output = args.output_dir / f"{sample}_all_group_activity_run.json"
        if output.exists() and run_output.exists():
            activity = pd.read_csv(output, sep="\t", index_col=0)
            if activity.shape == (657, metadata["expected_spots"]):
                print(f"VALID CHECKPOINT {sample} {activity.shape}", flush=True)
                continue
            raise ValueError(f"Invalid existing checkpoint for {sample}: {activity.shape}")

        released = read_slice(args.data_dir, metadata)
        observed = [gene for gene in released["genes"] if gene in signature_genes]
        lookup = pd.Series(range(len(released["genes"])), index=released["genes"])
        counts = released["counts"][lookup.loc[observed].to_numpy(int)].toarray()
        expression = pd.DataFrame(
            __import__("numpy").log2(
                counts / released["total_counts"][None, :] * 1_000_000 + 1
            ),
            index=observed,
            columns=released["barcodes"],
        )
        print(
            f"FULL SECACT START {sample} genes={len(observed)} spots={expression.shape[1]}",
            flush=True,
        )
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
        activity = group_level_matrix(result["zscore"])
        if activity.shape != (657, metadata["expected_spots"]):
            raise ValueError(f"Unexpected full group activity shape: {activity.shape}")
        activity.to_csv(output, sep="\t", compression="gzip")
        run = {
            "dataset": "GSE206245",
            "sample": sample,
            "secactpy_version": version("secactpy"),
            "seed": SEED,
            "n_rand": args.n_rand,
            "rng_method": "gsl",
            "n_input_genes": len(observed),
            "n_spots": expression.shape[1],
            "n_groups": activity.shape[0],
        }
        run_output.write_text(json.dumps(run, indent=2), encoding="utf-8")
        print(json.dumps(run, indent=2), flush=True)
        del released, counts, expression, result, activity
        gc.collect()


if __name__ == "__main__":
    main()
