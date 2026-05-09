#!/usr/bin/env Rscript

# ============================================================
# Script: 06_reactome_communities_ko.R
#
# Purpose:
#   Builds Reactome pathway communities for the WT-vs-Gpnmb KO
#   injury interaction contrast using leading-edge gene overlap.
#
# Inputs:
#   - results/enrichment/leading_edge/driver_genes_leadingedge_Injury_interaction_WTvsKO.tsv
#
# Outputs:
#   - results/reactome_communities/ko/pathway_nodes.tsv
#   - results/reactome_communities/ko/pathway_edges.tsv
#   - results/reactome_communities/ko/community_summary.tsv
#   - results/reactome_communities/ko/community_recurrent_genes.tsv
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(yaml)
})

source("R/reactome_community_helpers.R")

load_paths <- function(config_file = "config/paths.yml") {
  yaml::read_yaml(config_file)
}

as_project_path <- function(path, paths) {
  if (grepl("^/", path)) path else file.path(paths$project$root, path)
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

paths <- load_paths()
contrast <- "Injury_interaction_WTvsKO"

input_file <- file.path(
  as_project_path(paths$outputs$enrichment, paths),
  "leading_edge",
  paste0("driver_genes_leadingedge_", contrast, ".tsv")
)

out_dir <- ensure_dir(file.path(
  as_project_path(paths$outputs$reactome_communities, paths),
  "ko"
))

if (!file.exists(input_file)) stop("Missing input file: ", input_file)

leading_edge <- readr::read_tsv(input_file, show_col_types = FALSE)

reactome_table <- prepare_reactome_table(
  leading_edge,
  fdr_cutoff = 0.05,
  min_ngenes = 10,
  max_ngenes = 300
)

result <- build_communities_for_table(
  reactome_table,
  jaccard_threshold = 0.20
)

readr::write_tsv(result$node_table, file.path(out_dir, "pathway_nodes.tsv"))
readr::write_tsv(result$edge_table, file.path(out_dir, "pathway_edges.tsv"))
readr::write_tsv(result$community_summary, file.path(out_dir, "community_summary.tsv"))
readr::write_tsv(result$community_genes, file.path(out_dir, "community_recurrent_genes.tsv"))

message("KO Reactome community analysis complete.")
message("Reactome pathways kept: ", nrow(reactome_table))
message("Output directory: ", normalizePath(out_dir))
