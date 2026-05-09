#!/usr/bin/env python3
"""
Script: 00_check_project_inputs.py

Purpose:
    Checks that the local input files required for the public-facing
    Gpnmb SCI analysis repository are present.

Inputs:
    - config/paths.yml
    - raw count table
    - raw sample metadata
    - optional reference folders/files listed in config/paths.yml

Outputs:
    - results/project_setup/input_manifest.tsv

Notes:
    This script does not inspect biological results or run analysis.
    It only checks that expected files and folders exist locally.
"""

from pathlib import Path
import sys
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))

from python.io_helpers import load_paths, project_path, ensure_dir


def describe_path(label: str, path: Path, required: bool = True) -> dict:
    exists = path.exists()

    if exists and path.is_file():
        path_type = "file"
        size_mb = path.stat().st_size / 1024**2
    elif exists and path.is_dir():
        path_type = "directory"
        size_mb = None
    else:
        path_type = "missing"
        size_mb = None

    return {
        "label": label,
        "path": str(path),
        "required": required,
        "exists": exists,
        "type": path_type,
        "size_mb": round(size_mb, 3) if size_mb is not None else "",
    }


def main() -> None:
    paths = load_paths()

    checks = [
        describe_path(
            "raw_counts",
            project_path(paths["inputs"]["raw_counts"], paths),
            required=True,
        ),
        describe_path(
            "raw_metadata",
            project_path(paths["inputs"]["raw_metadata"], paths),
            required=True,
        ),
        describe_path(
            "msigdb_dir",
            project_path(paths["references"]["msigdb_dir"], paths),
            required=False,
        ),
        describe_path(
            "marker_panels",
            project_path(paths["references"]["marker_panels"], paths),
            required=False,
        ),
        describe_path(
            "tabpar_raw_dir",
            project_path(paths["references"]["tabpar_raw_dir"], paths),
            required=False,
        ),
    ]

    manifest = pd.DataFrame(checks)

    out_dir = ensure_dir(project_path("results/project_setup", paths))
    out_file = out_dir / "input_manifest.tsv"
    manifest.to_csv(out_file, sep="\t", index=False)

    missing_required = manifest.query("required == True and exists == False")

    print(manifest.to_string(index=False))
    print(f"\nWrote manifest: {out_file}")

    if not missing_required.empty:
        missing = ", ".join(missing_required["label"].tolist())
        raise SystemExit(f"\nMissing required input(s): {missing}")

    print("\nInput check complete.")


if __name__ == "__main__":
    main()
