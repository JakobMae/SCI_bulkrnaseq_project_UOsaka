#!/usr/bin/env python3
"""
Script: 02_pca_sample_identity.py

Purpose:
    Performs PCA-based sample identity analysis using normalized logCPM
    values exported by the limma-voom pipeline.

Inputs:
    - results/limma_voom/objects/logcpm_annotated.tsv
    - results/limma_voom/objects/sample_metadata.tsv

Outputs:
    - results/tables/pca/pca_scores.tsv
    - results/tables/pca/pca_variance_explained.tsv
    - results/figures/pca/pca_all_samples.pdf
    - results/figures/pca/pca_all_samples.png

Notes:
    PCA is performed on the 2,000 most variable genes across samples.
    This script is for sample-level QC and figure generation only.
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

sys.path.append(str(Path(__file__).resolve().parents[1]))

from python.io_helpers import load_paths, project_path, ensure_dir, require_file
from python.plotting_theme import GROUP_COLORS, GROUP_LABELS, GROUP_MARKERS


N_TOP_VARIABLE_GENES = 2000


def main() -> None:
    paths = load_paths()

    logcpm_file = require_file(
        project_path("results/limma_voom/objects/logcpm_annotated.tsv", paths),
        "logCPM expression matrix",
    )
    metadata_file = require_file(
        project_path("results/limma_voom/objects/sample_metadata.tsv", paths),
        "sample metadata",
    )

    table_dir = ensure_dir(project_path("results/tables/pca", paths))
    figure_dir = ensure_dir(project_path("results/figures/pca", paths))

    logcpm = pd.read_csv(logcpm_file, sep="\t")
    metadata = pd.read_csv(metadata_file, sep="\t")

    annotation_cols = ["gene_id", "gene_symbol", "gene_name"]
    sample_cols = [c for c in logcpm.columns if c not in annotation_cols]

    missing_samples = set(metadata["sample_id"]) - set(sample_cols)
    if missing_samples:
        raise ValueError(f"Samples in metadata missing from logCPM table: {sorted(missing_samples)}")

    expression = logcpm.set_index("gene_id")[metadata["sample_id"]]

    gene_variance = expression.var(axis=1).sort_values(ascending=False)
    top_genes = gene_variance.head(N_TOP_VARIABLE_GENES).index
    expression_top = expression.loc[top_genes].T

    scaled = StandardScaler().fit_transform(expression_top)
    pca = PCA(n_components=5)
    scores = pca.fit_transform(scaled)

    score_table = pd.DataFrame(
        scores[:, :5],
        columns=[f"PC{i}" for i in range(1, 6)],
    )
    score_table.insert(0, "sample_id", expression_top.index)
    score_table = score_table.merge(metadata, on="sample_id", how="left")

    variance_table = pd.DataFrame({
        "PC": [f"PC{i}" for i in range(1, 6)],
        "variance_explained": pca.explained_variance_ratio_[:5],
        "percent_variance_explained": pca.explained_variance_ratio_[:5] * 100,
    })

    score_table.to_csv(table_dir / "pca_scores.tsv", sep="\t", index=False)
    variance_table.to_csv(table_dir / "pca_variance_explained.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(4.2, 3.6))

    for group, group_df in score_table.groupby("group"):
        ax.scatter(
            group_df["PC1"],
            group_df["PC2"],
            label=GROUP_LABELS.get(group, group),
            color=GROUP_COLORS.get(group, "black"),
            marker=GROUP_MARKERS.get(group, "o"),
            s=55,
            edgecolor="black",
            linewidth=0.4,
        )

    pc1 = variance_table.loc[0, "percent_variance_explained"]
    pc2 = variance_table.loc[1, "percent_variance_explained"]

    ax.set_xlabel(f"PC1 ({pc1:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({pc2:.1f}% variance)")
    ax.set_title("PCA of bulk RNA-seq samples")
    ax.legend(frameon=False, fontsize=7, loc="best")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(figure_dir / "pca_all_samples.pdf")
    fig.savefig(figure_dir / "pca_all_samples.png", dpi=300)
    plt.close(fig)

    print("PCA analysis complete.")
    print(f"Wrote tables to: {table_dir}")
    print(f"Wrote figures to: {figure_dir}")


if __name__ == "__main__":
    main()
