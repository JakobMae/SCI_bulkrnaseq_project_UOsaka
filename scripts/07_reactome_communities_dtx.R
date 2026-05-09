#!/usr/bin/env Rscript

# ============================================================
# Script: 07_reactome_communities_dtx.R
#
# Purpose:
#   Builds Reactome pathway communities for the DTX-vs-PBS
#   Gpnmb-DTR contrast using leading-edge gene overlap.
#
#   Upregulated and downregulated Reactome pathways are clustered
#   separately, matching the original notebook logic.
#
# Inputs:
#   - results/enrichment/leading_edge/driver_genes_leadingedge_DTR_DTX_vs_PBS_dpi14.tsv
#
# Outputs:
#   - results/reactome_communities/dtx/pathway_nodes.tsv
#   - results/reactome_communities/dtx/pathway_edges.tsv
#   - results/reactome_communities/dtx/community_summary.tsv
#   - results/reactome_communities/dtx/community_recurrent_genes.tsv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
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
contrast <- "DTR_DTX_vs_PBS_dpi14"

input_file <- file.path(
  as_project_path(paths$outputs$enrichment, paths),
  "leading_edge",
  paste0("driver_genes_leadingedge_", contrast, ".tsv")
)

out_dir <- ensure_dir(file.path(
  as_project_path(paths$outputs$reactome_communities, paths),
  "dtx"
))

if (!file.exists(input_file)) stop("Missing input file: ", input_file)

leading_edge <- readr::read_tsv(input_file, show_col_types = FALSE)

reactome_table <- prepare_reactome_table(
  leading_edge,
  fdr_cutoff = 0.05,
  min_ngenes = 10,
  max_ngenes = 300
)

results <- list()

for (direction in c("Up", "Down")) {
  direction_table <- reactome_table %>%
    filter(Direction == direction)

  if (nrow(direction_table) < 2) next

  direction_result <- build_communities_for_table(
    direction_table,
    jaccard_threshold = 0.20
  )

  direction_result$node_table <- direction_result$node_table %>%
    mutate(direction_clustered = direction)

  direction_result$edge_table <- direction_result$edge_table %>%
    mutate(direction_clustered = direction)

  direction_result$community_summary <- direction_result$community_summary %>%
    mutate(direction_clustered = direction)

  direction_result$community_genes <- direction_result$community_genes %>%
    mutate(direction_clustered = direction)

  results[[direction]] <- direction_result
}

node_table <- bind_rows(lapply(results, `[[`, "node_table"))
edge_table <- bind_rows(lapply(results, `[[`, "edge_table"))
community_summary <- bind_rows(lapply(results, `[[`, "community_summary"))
community_genes <- bind_rows(lapply(results, `[[`, "community_genes"))

readr::write_tsv(node_table, file.path(out_dir, "pathway_nodes.tsv"))
readr::write_tsv(edge_table, file.path(out_dir, "pathway_edges.tsv"))
readr::write_tsv(community_summary, file.path(out_dir, "community_summary.tsv"))
readr::write_tsv(community_genes, file.path(out_dir, "community_recurrent_genes.tsv"))

message("DTX Reactome community analysis complete.")
message("Reactome pathways kept: ", nrow(reactome_table))
message("Output directory: ", normalizePath(out_dir))
