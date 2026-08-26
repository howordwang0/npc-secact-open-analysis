from __future__ import annotations

import argparse
import csv
import json
import sqlite3
from pathlib import Path

from secactpy import load_signature


def main() -> None:
    parser = argparse.ArgumentParser(description="Build probe-to-symbol mapping from Bioconductor SQLite databases.")
    parser.add_argument("chip_database", type=Path)
    parser.add_argument("organism_database", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()

    chip = sqlite3.connect(args.chip_database)
    organism = sqlite3.connect(args.organism_database)
    gene_symbols = {
        gene_id: symbol
        for gene_id, symbol in organism.execute(
            "SELECT genes.gene_id, gene_info.symbol "
            "FROM genes JOIN gene_info USING (_id) WHERE gene_info.symbol IS NOT NULL"
        )
    }
    mappings = []
    for probe_id, gene_id, is_multiple in chip.execute(
        "SELECT probe_id, gene_id, is_multiple FROM probes WHERE gene_id IS NOT NULL"
    ):
        symbol = gene_symbols.get(str(gene_id))
        if symbol:
            mappings.append((str(probe_id), str(gene_id), symbol, int(is_multiple)))

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["probe_id", "entrez_id", "symbol", "is_multiple"])
        writer.writerows(mappings)

    signature_genes = set(load_signature("secact").index.astype(str))
    mapped_symbols = {row[2] for row in mappings}
    result = {
        "chip_database": str(args.chip_database),
        "organism_database": str(args.organism_database),
        "mapping_rows": len(mappings),
        "unique_probes": len({row[0] for row in mappings}),
        "unique_symbols": len(mapped_symbols),
        "secact_signature_symbols_mappable": len(mapped_symbols & signature_genes),
        "secact_signature_coverage_percent": round(
            100 * len(mapped_symbols & signature_genes) / len(signature_genes), 3
        ),
        "probe_examples": [row[0] for row in mappings[:5]],
        "output_csv": str(args.output_csv),
    }
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
