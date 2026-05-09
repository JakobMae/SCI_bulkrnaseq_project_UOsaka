#!/usr/bin/env python3
"""
Script: 03_deg_volcano_concordance.py

Purpose:
    Generates DEG summary tables, ranked gene lists, volcano plots,
    and the WT-vs-KO injury-response concordance plot.

Inputs:
    - results/limma_voom/deg/deg_WT_SCI_vs_SHAM_dpi7_annotated.tsv
    - results/limma_voom/deg/deg_KO_SCI_vs_SHAM_dpi7_annotated.tsv
    - results/limma_voom/deg/deg_DTR_DTX_vs_PBS_dpi14_annotated.tsv
    - results/limma_voom/deg/deg_Injury_interaction_WTvsKO_annotated.tsv

Outputs:
    - results/tables/deg_summary/deg_counts_by_contrast.tsv
    - results/tables/ranked_gene_lists/*.tsv
    - results/figures/volcano/*.pdf and *.png
    - results/figures/concordance/wt_ko_injury_logfc_concordance.pdf and *.png

Notes:
    DEG thresholds are FDR < 0.05 and |log2FC| >= 1 for standard
    pairwise contrasts. The interaction contrast is summarized at
    FDR < 0.05 without an additional logFC threshold.
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parents[1]))

from python.io_helpers import load_paths, project_path, ensure_dir, require_file


STANDARD_FDR = 0.05
STANDARD_LOGFC = 1.0
INTERACTION_FDR = 0.05

VOLCANO_COLORS = {
    "up": "#D73027",
    "down": "#4575B4",
    "ns": "#BDBDBD",
}

CONTRASTS = {
    "WT_SCI_vs_SHAM_dpi7": {
        "file": "deg_WT_SCI_vs_SHAM_dpi7_annotated.tsv",
        "title": "WT SCI vs sham, 7 dpi",
        "threshold_type": "standard",
    },
    "KO_SCI_vs_SHAM_dpi7": {
        "file": "deg_KO_SCI_vs_SHAM_dpi7_annotated.tsv",
        "title": "Gpnmb KO SCI vs sham, 7 dpi",
        "threshold_type": "standard",
    },
    "DTR_DTX_vs_PBS_dpi14": {
        "file": "deg_DTR_DTX_vs_PBS_dpi14_annotated.tsv",
        "title": "Gpnmb-DTR DTX vs PBS, 14 dpi",
        "threshold_type": "standard",
    },
    "Injury_interaction_WTvsKO": {
        "file": "deg_Injury_interaction_WTvsKO_annotated.tsv",
        "title": "Injury interaction: WT vs Gpnmb KO",
        "threshold_type": "interaction",
    },
}


REQUIRED_COLUMNS = [
    "gene_id",
    "gene_symbol",
    "gene_name",
    "logFC",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.Val",
    "B",
    "gene_id_stripped",
]


def read_deg_table(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t")
    missing = [c for c in REQUIRED_COLUMNS if c not in table.columns]
    if missing:
        raise ValueError(f"{path.name} is missing required columns: {missing}")
    return table


def add_deg_status(table: pd.DataFrame, threshold_type: str) -> pd.DataFrame:
    table = table.copy()

    if threshold_type == "standard":
        is_sig = (table["adj.P.Val"] < STANDARD_FDR) & (table["logFC"].abs() >= STANDARD_LOGFC)
    elif threshold_type == "interaction":
        is_sig = table["adj.P.Val"] < INTERACTION_FDR
    else:
        raise ValueError(f"Unknown threshold type: {threshold_type}")

    table["deg_status"] = "ns"
    table.loc[is_sig & (table["logFC"] > 0), "deg_status"] = "up"
    table.loc[is_sig & (table["logFC"] < 0), "deg_status"] = "down"

    return table


def summarize_deg(table: pd.DataFrame, contrast: str, threshold_type: str) -> dict:
    table = add_deg_status(table, threshold_type)

    return {
        "contrast": contrast,
        "threshold_type": threshold_type,
        "n_total_genes": len(table),
        "n_significant": int((table["deg_status"] != "ns").sum()),
        "n_up": int((table["deg_status"] == "up").sum()),
        "n_down": int((table["deg_status"] == "down").sum()),
    }


def export_ranked_list(table: pd.DataFrame, out_file: Path) -> None:
    ranked = table.sort_values("t", ascending=False)
    ranked.to_csv(out_file, sep="\t", index=False)


def plot_volcano(table: pd.DataFrame, contrast: str, title: str, threshold_type: str, out_dir: Path) -> None:
    table = add_deg_status(table, threshold_type)
    table["minus_log10_fdr"] = -np.log10(table["adj.P.Val"].clip(lower=np.finfo(float).tiny))

    fig, ax = plt.subplots(figsize=(4.0, 3.6))

    for status in ["ns", "down", "up"]:
        subset = table[table["deg_status"] == status]
        ax.scatter(
            subset["logFC"],
            subset["minus_log10_fdr"],
            s=8,
            color=VOLCANO_COLORS[status],
            alpha=0.75,
            linewidth=0,
            label=status,
        )

    ax.axhline(-np.log10(STANDARD_FDR), color="black", linewidth=0.6, linestyle="--")

    if threshold_type == "standard":
        ax.axvline(-STANDARD_LOGFC, color="black", linewidth=0.6, linestyle="--")
        ax.axvline(STANDARD_LOGFC, color="black", linewidth=0.6, linestyle="--")

    ax.set_xlabel("log2 fold-change")
    ax.set_ylabel("-log10 FDR")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=7, markerscale=1.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(out_dir / f"volcano_{contrast}.pdf")
    fig.savefig(out_dir / f"volcano_{contrast}.png", dpi=300)
    plt.close(fig)


def plot_wt_ko_concordance(wt: pd.DataFrame, ko: pd.DataFrame, out_dir: Path) -> None:
    merged = wt[["gene_id", "gene_symbol", "logFC"]].rename(columns={"logFC": "logFC_WT"})
    merged = merged.merge(
        ko[["gene_id", "logFC"]].rename(columns={"logFC": "logFC_KO"}),
        on="gene_id",
        how="inner",
    )

    r = merged[["logFC_WT", "logFC_KO"]].corr(method="pearson").iloc[0, 1]

    fig, ax = plt.subplots(figsize=(4.0, 3.8))
    ax.scatter(
        merged["logFC_WT"],
        merged["logFC_KO"],
        s=8,
        color="#4D4D4D",
        alpha=0.45,
        linewidth=0,
    )

    lim = max(abs(merged["logFC_WT"]).max(), abs(merged["logFC_KO"]).max())
    ax.plot([-lim, lim], [-lim, lim], color="black", linewidth=0.8, linestyle="--")

    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("WT SCI vs sham log2FC")
    ax.set_ylabel("KO SCI vs sham log2FC")
    ax.set_title(f"WT and KO injury-response concordance\nPearson r = {r:.3f}")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(out_dir / "wt_ko_injury_logfc_concordance.pdf")
    fig.savefig(out_dir / "wt_ko_injury_logfc_concordance.png", dpi=300)
    plt.close(fig)

    merged.to_csv(out_dir / "wt_ko_injury_logfc_concordance_points.tsv", sep="\t", index=False)


def main() -> None:
    paths = load_paths()

    deg_dir = project_path("results/limma_voom/deg", paths)
    summary_dir = ensure_dir(project_path("results/tables/deg_summary", paths))
    ranked_dir = ensure_dir(project_path("results/tables/ranked_gene_lists", paths))
    volcano_dir = ensure_dir(project_path("results/figures/volcano", paths))
    concordance_dir = ensure_dir(project_path("results/figures/concordance", paths))

    deg_tables = {}
    summaries = []

    for contrast, spec in CONTRASTS.items():
        path = require_file(deg_dir / spec["file"], f"DEG table for {contrast}")
        table = read_deg_table(path)
        deg_tables[contrast] = table

        summaries.append(summarize_deg(table, contrast, spec["threshold_type"]))
        export_ranked_list(table, ranked_dir / f"ranked_{contrast}.tsv")
        plot_volcano(table, contrast, spec["title"], spec["threshold_type"], volcano_dir)

    summary_table = pd.DataFrame(summaries)
    summary_table.to_csv(summary_dir / "deg_counts_by_contrast.tsv", sep="\t", index=False)

    plot_wt_ko_concordance(
        deg_tables["WT_SCI_vs_SHAM_dpi7"],
        deg_tables["KO_SCI_vs_SHAM_dpi7"],
        concordance_dir,
    )

    print("DEG figure and summary export complete.")
    print(f"Wrote DEG summaries to: {summary_dir}")
    print(f"Wrote ranked lists to: {ranked_dir}")
    print(f"Wrote volcano plots to: {volcano_dir}")
    print(f"Wrote concordance plot to: {concordance_dir}")


if __name__ == "__main__":
    main()
