#!/usr/bin/env Rscript

# ============================================================
# Script: 05_leading_edge_driver_genes.R
#
# Purpose:
#   Summarizes recurrent leading-edge genes from pathway enrichment
#   results. This identifies genes that repeatedly contribute to
#   significant enriched pathways within each contrast.
#
# Inputs:
#   - results/enrichment/leading_edge/driver_genes_leadingedge_<contrast>.tsv
#   - results/limma_voom/deg/deg_<contrast>_annotated.tsv
#
# Outputs:
#   - results/leading_edge/driver_gene_tables/recurrent_leading_edge_<contrast>.tsv
#   - results/leading_edge/driver_gene_tables/recurrent_reactome_leading_edge_<contrast>.tsv
#
# Notes:
#   This script does not perform pathway statistics. It summarizes
#   outputs from 04_pathway_enrichment_camera_fgsea.R.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(yaml)
})

CONTRASTS <- c(
  "WT_SCI_vs_SHAM_dpi7",
  "KO_SCI_vs_SHAM_dpi7",
  "DTR_DTX_vs_PBS_dpi14",
  "Injury_interaction_WTvsKO"
)

load_paths <- function(config_file = "config/paths.yml") {
  if (!file.exists(config_file)) {
    stop("Missing config/paths.yml.")
  }
  yaml::read_yaml(config_file)
}

as_project_path <- function(path, paths) {
  if (grepl("^/", path)) path else file.path(paths$project$root, path)
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

read_required_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path)
  readr::read_tsv(path, show_col_types = FALSE)
}

summarise_leading_edge <- function(leading_edge_table, deg_table, reactome_only = FALSE) {
  if (reactome_only) {
    leading_edge_table <- leading_edge_table %>%
      filter(collection == "reactome")
  }

  leading_edge_table %>%
    filter(!is.na(leadingEdge), leadingEdge != "") %>%
    select(collection, contrast, pathway, Direction, FDR, NES, padj, leadingEdge) %>%
    separate_rows(leadingEdge, sep = ";") %>%
    rename(gene_symbol = leadingEdge) %>%
    mutate(gene_symbol = str_trim(gene_symbol)) %>%
    filter(gene_symbol != "") %>%
    group_by(gene_symbol) %>%
    summarise(
      recurrence = n_distinct(pathway),
      n_collections = n_distinct(collection),
      collections = paste(sort(unique(collection)), collapse = ";"),
      directions = paste(sort(unique(Direction)), collapse = ";"),
      min_pathway_fdr = min(FDR, na.rm = TRUE),
      max_abs_nes = max(abs(NES), na.rm = TRUE),
      pathways = paste(sort(unique(pathway)), collapse = ";"),
      .groups = "drop"
    ) %>%
    left_join(
      deg_table %>%
        select(gene_id, gene_symbol, gene_name, logFC, AveExpr, t, P.Value, adj.P.Val, B),
      by = "gene_symbol"
    ) %>%
    arrange(desc(recurrence), min_pathway_fdr, adj.P.Val)
}

paths <- load_paths()

leading_dir <- file.path(as_project_path(paths$outputs$enrichment, paths), "leading_edge")
deg_dir <- file.path(as_project_path(paths$outputs$limma_voom, paths), "deg")
out_dir <- ensure_dir(file.path(as_project_path(paths$outputs$leading_edge, paths), "driver_gene_tables"))

for (contrast in CONTRASTS) {
  leading_file <- file.path(leading_dir, paste0("driver_genes_leadingedge_", contrast, ".tsv"))
  deg_file <- file.path(deg_dir, paste0("deg_", contrast, "_annotated.tsv"))

  leading_edge_table <- read_required_tsv(leading_file)
  deg_table <- read_required_tsv(deg_file)

  required_leading_cols <- c("collection", "contrast", "pathway", "Direction", "FDR", "NES", "padj", "leadingEdge")
  required_deg_cols <- c("gene_id", "gene_symbol", "gene_name", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")

  if (!all(required_leading_cols %in% colnames(leading_edge_table))) {
    stop("Leading-edge table has unexpected columns: ", leading_file)
  }
  if (!all(required_deg_cols %in% colnames(deg_table))) {
    stop("DEG table has unexpected columns: ", deg_file)
  }

  recurrent_all <- summarise_leading_edge(
    leading_edge_table = leading_edge_table,
    deg_table = deg_table,
    reactome_only = FALSE
  )

  recurrent_reactome <- summarise_leading_edge(
    leading_edge_table = leading_edge_table,
    deg_table = deg_table,
    reactome_only = TRUE
  )

  readr::write_tsv(
    recurrent_all,
    file.path(out_dir, paste0("recurrent_leading_edge_", contrast, ".tsv"))
  )

  readr::write_tsv(
    recurrent_reactome,
    file.path(out_dir, paste0("recurrent_reactome_leading_edge_", contrast, ".tsv"))
  )

  message("Wrote leading-edge summaries for: ", contrast)
}

message("Leading-edge driver-gene summary complete.")
message("Output directory: ", out_dir)
