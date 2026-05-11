#!/usr/bin/env python3
"""
Script: 09_build_tabpar_reference.py

Purpose:
    Builds a filtered Tabula Paralytica reference from the public
    GSE234774 single-cell RNA-seq matrix for downstream deconvolution
    and marker/reference support.

Inputs:
    - data/references/tabpar/raw_or_private_files/GSE234774_rnaseq_meta.txt.gz
    - data/references/tabpar/raw_or_private_files/GSE234774_rnaseq_barcodes.txt.gz
    - data/references/tabpar/raw_or_private_files/GSE234774_rnaseq_features.txt.gz
    - data/references/tabpar/raw_or_private_files/GSE234774_rnaseq_filtered_scRNA.mtx.gz

Outputs:
    - data/references/tabpar/processed/tabpar_timecourse_7d14d_metadata.tsv.gz
    - data/references/tabpar/processed/tabpar_timecourse_7d14d_barcodes.tsv.gz
    - data/references/tabpar/processed/tabpar_timecourse_7d14d_features.tsv.gz
    - data/references/tabpar/processed/tabpar_timecourse_7d14d_counts.mtx.gz
    - data/references/tabpar/processed/tabpar_timecourse_7d14d_pseudobulk_counts.tsv.gz
    - data/references/tabpar/processed/tabpar_timecourse_7d14d_pseudobulk_metadata.tsv
    - results/tabpar_reference/celltype_counts.tsv
    - results/tabpar_reference/pseudobulk_group_counts.tsv

Notes:
    The script keeps the matrix sparse throughout. It filters to the
    Time course 7d and 14d samples, maps detailed TabPar labels to broad
    reference labels, removes lymphocytes/pericyte-mural cells/unmapped
    cells, and builds pseudobulk profiles by library_clean × broad label.
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread, mmwrite

sys.path.append(str(Path(__file__).resolve().parents[1]))

from python.io_helpers import load_paths, project_path, ensure_dir, require_file


def assign_broad_labels(meta: pd.DataFrame) -> pd.Series:
    layer4 = meta["layer4"].astype(str)
    layer1 = meta["layer1"].astype(str)

    broad = np.full(meta.shape[0], None, dtype=object)

    broad[np.isin(layer4, [
        "Astrocytes, fibrous",
        "Astrocytes, protoplasmic",
        "Astrocytes, reactive",
        "Astrocytes, cycling",
    ])] = "Astrocytes"

    broad[np.isin(layer4, [
        "Ependymal cells",
        "Neural stem cells",
    ])] = "Ependymal_NS"

    broad[np.isin(layer4, [
        "Microglia",
        "Microglia, activated",
        "Microglia, homeostatic",
        "Microglia, homeostatic-injury susceptible",
        "Myeloid, dividing",
    ])] = "Microglia"

    broad[np.isin(layer4, [
        "Macrophages, inflammatory",
        "Macrophages, chemotaxis-inducing",
        "Infiltrating, bordering leukocytes",
    ])] = "Macrophages"

    broad[np.isin(layer4, [
        "NK, T cells",
        "B cells",
    ])] = "Lymphocytes"

    broad[np.isin(layer4, [
        "Differentiation-committed oligodendrocyte precursors (COPs)",
        "Newly formed oligodendrocytes (NFOL)",
        "Oligodendrocytes precursor cells (OPC)",
    ])] = "Oligo_precursor"

    broad[np.isin(layer4, [
        "Mature oligodendrocytes (MOL)",
        "Mature oligodendrocytes (MOL), ischemic",
        "Myelin forming oligodendrocytes (MFOL)",
    ])] = "Oligo_mature"

    broad[np.isin(layer4, [
        "Vascular endothelial cells",
        "Vascular endothelial cells, arterial",
        "Vascular endothelial cells, venous",
    ])] = "Endothelial"

    broad[np.isin(layer4, [
        "Pericytes",
        "Mural cells",
    ])] = "Pericyte_mural"

    broad[np.isin(layer4, [
        "Vascular leptomeningeal cells",
    ])] = "VLMC_stromal"

    broad[(pd.isna(broad)) & (layer1 == "Neurons")] = "Neurons"

    return pd.Series(broad, index=meta.index, name="ref_label_broad")


def main() -> None:
    paths = load_paths()

    raw_dir = project_path(paths["references"]["tabpar_raw_dir"], paths)
    processed_dir = ensure_dir(project_path(paths["references"]["tabpar_processed_dir"], paths))
    qc_dir = ensure_dir(project_path("results/tabpar_reference", paths))

    meta_file = require_file(raw_dir / "GSE234774_rnaseq_meta.txt.gz", "TabPar metadata")
    barcode_file = require_file(raw_dir / "GSE234774_rnaseq_barcodes.txt.gz", "TabPar barcodes")
    feature_file = require_file(raw_dir / "GSE234774_rnaseq_features.txt.gz", "TabPar features")
    matrix_file = require_file(raw_dir / "GSE234774_rnaseq_filtered_scRNA.mtx.gz", "TabPar count matrix")

    print("Loading metadata, barcodes, and features...")
    meta = pd.read_csv(meta_file, sep="\t", compression="gzip", dtype=str)
    barcodes = pd.read_csv(barcode_file, sep="\t", header=None, names=["barcode"], compression="gzip")
    features = pd.read_csv(feature_file, sep="\t", header=None, names=["gene_symbol"], compression="gzip")

    required_meta_cols = [
        "barcode",
        "library_clean",
        "experiment_clean",
        "label_clean",
        "replicate_clean",
        "layer1",
        "layer4",
    ]
    missing_cols = [c for c in required_meta_cols if c not in meta.columns]
    if missing_cols:
        raise ValueError(f"Metadata missing required columns: {missing_cols}")

    if len(meta) != len(barcodes):
        raise ValueError(f"Metadata rows ({len(meta)}) do not match barcodes ({len(barcodes)}).")

    if not (meta["barcode"].to_numpy() == barcodes["barcode"].to_numpy()).all():
        raise ValueError("Metadata barcode order does not match barcode file order.")

    meta["ref_label_broad"] = assign_broad_labels(meta)

    selected = (
        (meta["experiment_clean"] == "Time course")
        & (meta["label_clean"].isin(["7d", "14d"]))
        & (meta["ref_label_broad"].notna())
        & (~meta["ref_label_broad"].isin(["Lymphocytes", "Pericyte_mural"]))
    )

    selected_idx = np.flatnonzero(selected.to_numpy())

    ref_meta = meta.loc[selected].copy()
    ref_meta.insert(0, "cell_index_original", selected_idx)
    ref_barcodes = barcodes.iloc[selected_idx].copy()

    if len(ref_meta) == 0:
        raise ValueError("No cells remain after TabPar filtering.")

    print(f"Selected cells: {len(ref_meta):,}")
    print("Cell counts by broad label:")
    print(ref_meta["ref_label_broad"].value_counts().to_string())

    celltype_counts = (
        ref_meta
        .groupby(["label_clean", "ref_label_broad"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["label_clean", "ref_label_broad"])
    )
    celltype_counts.to_csv(qc_dir / "celltype_counts.tsv", sep="\t", index=False)

    pseudobulk_meta = (
        ref_meta
        .groupby(["library_clean", "label_clean", "ref_label_broad"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["library_clean", "ref_label_broad"])
    )
    pseudobulk_meta["pseudobulk_sample"] = (
        pseudobulk_meta["library_clean"] + "__" + pseudobulk_meta["ref_label_broad"]
    )

    pseudobulk_meta.to_csv(
        processed_dir / "tabpar_timecourse_7d14d_pseudobulk_metadata.tsv",
        sep="\t",
        index=False,
    )
    pseudobulk_meta.to_csv(
        qc_dir / "pseudobulk_group_counts.tsv",
        sep="\t",
        index=False,
    )

    print("Loading sparse matrix. This may take a while...")
    counts = mmread(matrix_file).tocsr()

    expected_shape = (len(features), len(barcodes))
    if counts.shape != expected_shape:
        raise ValueError(f"Matrix shape {counts.shape} does not match features/barcodes {expected_shape}.")

    print("Subsetting sparse matrix...")
    ref_counts = counts[:, selected_idx].tocsc()

    print("Writing filtered reference matrix and metadata...")
    ref_meta.to_csv(
        processed_dir / "tabpar_timecourse_7d14d_metadata.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    ref_barcodes.to_csv(
        processed_dir / "tabpar_timecourse_7d14d_barcodes.tsv.gz",
        sep="\t",
        index=False,
        header=False,
        compression="gzip",
    )
    features.to_csv(
        processed_dir / "tabpar_timecourse_7d14d_features.tsv.gz",
        sep="\t",
        index=False,
        header=False,
        compression="gzip",
    )
    mmwrite(
        str(processed_dir / "tabpar_timecourse_7d14d_counts.mtx"),
        ref_counts,
    )

    # scipy mmwrite does not write gzip directly.
    import gzip
    import shutil

    with open(processed_dir / "tabpar_timecourse_7d14d_counts.mtx", "rb") as src:
        with gzip.open(processed_dir / "tabpar_timecourse_7d14d_counts.mtx.gz", "wb") as dst:
            shutil.copyfileobj(src, dst)
    Path(processed_dir / "tabpar_timecourse_7d14d_counts.mtx").unlink()

    print("Building pseudobulk matrix...")
    pb_counts = []
    pb_names = []

    for _, row in pseudobulk_meta.iterrows():
        mask = (
            (ref_meta["library_clean"] == row["library_clean"])
            & (ref_meta["ref_label_broad"] == row["ref_label_broad"])
        ).to_numpy()

        cell_positions = np.flatnonzero(mask)
        summed = np.asarray(ref_counts[:, cell_positions].sum(axis=1)).ravel()
        pb_counts.append(summed)
        pb_names.append(row["pseudobulk_sample"])

    pb_matrix = pd.DataFrame(
        np.vstack(pb_counts).T,
        index=features["gene_symbol"],
        columns=pb_names,
    )
    pb_matrix.index.name = "gene_symbol"

    pb_matrix.to_csv(
        processed_dir / "tabpar_timecourse_7d14d_pseudobulk_counts.tsv.gz",
        sep="\t",
        compression="gzip",
    )

    print("TabPar reference construction complete.")
    print(f"Processed outputs: {processed_dir}")
    print(f"QC outputs: {qc_dir}")


if __name__ == "__main__":
    main()
