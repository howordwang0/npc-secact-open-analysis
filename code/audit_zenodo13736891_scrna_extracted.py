from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
from pathlib import Path


def count_lines(path: Path) -> int:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return sum(1 for _ in handle)


def matrix_dimensions(path: Path) -> tuple[int, int, int]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = next(handle).strip()
        if header != "%%MatrixMarket matrix coordinate integer general":
            raise ValueError(f"Unexpected Matrix Market header in {path}: {header}")
        for line in handle:
            if not line.startswith("%"):
                rows, columns, nonzero = map(int, line.split())
                return rows, columns, nonzero
    raise ValueError(f"Missing matrix dimensions in {path}")


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("expected_counts", type=Path)
    parser.add_argument("archive_dir", type=Path)
    parser.add_argument("extracted_dir", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    with args.expected_counts.open(encoding="utf-8-sig", newline="") as handle:
        expected = {
            row["archive_sample_id"]: int(row["total_cells"])
            for row in csv.DictReader(handle)
        }

    rows = []
    with args.manifest.open(encoding="utf-8-sig", newline="") as handle:
        for item in csv.DictReader(handle):
            sample_id = item["sample_id"]
            archive = args.archive_dir / item["filename"]
            sample_dir = args.extracted_dir / sample_id
            barcodes = sample_dir / "barcodes.tsv.gz"
            features = sample_dir / "features.tsv.gz"
            matrix = sample_dir / "matrix.mtx.gz"
            missing = [str(path) for path in (barcodes, features, matrix) if not path.exists()]
            if missing:
                raise FileNotFoundError(f"Missing files for {sample_id}: {missing}")

            barcode_count = count_lines(barcodes)
            feature_count = count_lines(features)
            matrix_rows, matrix_columns, matrix_nonzero = matrix_dimensions(matrix)
            if matrix_rows != feature_count or matrix_columns != barcode_count:
                raise ValueError(
                    f"Dimension mismatch for {sample_id}: "
                    f"matrix={matrix_rows}x{matrix_columns}, "
                    f"features={feature_count}, barcodes={barcode_count}"
                )

            retained = expected[sample_id]
            archive_md5 = md5(archive)
            rows.append(
                {
                    "patient_id": item["patient_id"],
                    "timepoint": item["timepoint"],
                    "sample_id": sample_id,
                    "archive_filename": item["filename"],
                    "archive_size_bytes": archive.stat().st_size,
                    "archive_md5": archive_md5,
                    "archive_integrity_ok": str(
                        archive.stat().st_size == int(item["size_bytes"])
                        and archive_md5 == item["md5"]
                    ).lower(),
                    "raw_barcodes": barcode_count,
                    "features": feature_count,
                    "matrix_nonzero": matrix_nonzero,
                    "author_retained_cells": retained,
                    "author_retained_fraction_of_raw": f"{retained / barcode_count:.6f}",
                }
            )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
