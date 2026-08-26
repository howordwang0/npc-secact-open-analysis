from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_manifest(root: Path) -> None:
    files = sorted(
        path
        for path in root.rglob("*")
        if path.is_file()
        and path.name != "SHA256SUMS.txt"
        and ".git" not in path.relative_to(root).parts
        and "__pycache__" not in path.relative_to(root).parts
    )
    lines = [f"{sha256(path)}  {path.relative_to(root).as_posix()}" for path in files]
    (root / "SHA256SUMS.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"WROTE {len(lines)} checksums")


def verify_manifest(root: Path) -> None:
    manifest = root / "SHA256SUMS.txt"
    failures = []
    rows = [line for line in manifest.read_text(encoding="utf-8").splitlines() if line]
    for row in rows:
        expected, relative = row.split("  ", 1)
        path = root / Path(relative)
        observed = sha256(path) if path.is_file() else "MISSING"
        if observed != expected:
            failures.append((relative, expected, observed))
    if failures:
        for failure in failures:
            print("FAIL", *failure)
        raise SystemExit(1)
    print(f"PASS {len(rows)}/{len(rows)} checksums")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    root = args.root.resolve()
    if not root.is_dir():
        raise FileNotFoundError(root)
    if args.write:
        write_manifest(root)
    verify_manifest(root)


if __name__ == "__main__":
    main()
