#!/usr/bin/env Rscript

# ============================================================
# Script: 04_pathway_enrichment_camera_fgsea.R
#
# Purpose:
#   Runs pathway-level enrichment analysis on the limma-voom
#   bulk RNA-seq model.
#
#   camera() is used as the primary enrichment method because it
#   performs competitive gene set testing while accounting for
#   inter-gene correlation. fgsea is used supportively to provide
#   NES values and leading-edge genes.
#
# Inputs:
#   - results/limma_voom/objects/voom_object.rds
#   - results/limma_voom/objects/fit2_ebayes.rds
#   - results/limma_voom/objects/contrast_matrix.csv
#   - results/limma_voom/objects/logcpm_annotated.tsv
#   - data/references/msigdb/*.gmt
#
# Outputs:
#   - results/enrichment/camera/<collection>/camera_<contrast>.tsv
#   - results/enrichment/fgsea/<collection>/fgsea_<contrast>.tsv
#   - results/enrichment/summary/camera_allcollections_<contrast>.tsv
#   - results/enrichment/leading_edge/driver_genes_leadingedge_<contrast>.tsv
#   - results/figures/enrichment/<collection>/camera_top15_<contrast>.pdf/.png
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(limma)
  library(fgsea)
  library(yaml)
})

set.seed(1)

# -----------------------------
# Parameters
# -----------------------------

FDR_CUTOFF <- 0.05
MIN_GS_SIZE <- 15
MAX_GS_SIZE <- 500
TOP_N_PLOT <- 15

CONTRASTS_TO_RUN <- c(
  "WT_SCI_vs_SHAM_dpi7",
  "KO_SCI_vs_SHAM_dpi7",
  "DTR_DTX_vs_PBS_dpi14",
  "Injury_interaction_WTvsKO"
)

# -----------------------------
# Paths
# -----------------------------

load_paths <- function(config_file = "config/paths.yml") {
  if (!file.exists(config_file)) {
    stop("Missing config/paths.yml. Copy config/paths_template.yml and edit local paths.")
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

paths <- load_paths()

obj_dir <- file.path(as_project_path(paths$outputs$limma_voom, paths), "objects")
msigdb_dir <- as_project_path(paths$references$msigdb_dir, paths)

out_camera <- ensure_dir(file.path(as_project_path(paths$outputs$enrichment, paths), "camera"))
out_fgsea <- ensure_dir(file.path(as_project_path(paths$outputs$enrichment, paths), "fgsea"))
out_summary <- ensure_dir(file.path(as_project_path(paths$outputs$enrichment, paths), "summary"))
out_leading <- ensure_dir(file.path(as_project_path(paths$outputs$enrichment, paths), "leading_edge"))
out_fig <- ensure_dir(file.path(as_project_path(paths$outputs$figures, paths), "enrichment"))

f_voom <- file.path(obj_dir, "voom_object.rds")
f_fit2 <- file.path(obj_dir, "fit2_ebayes.rds")
f_contrasts <- file.path(obj_dir, "contrast_matrix.csv")
f_logcpm <- file.path(obj_dir, "logcpm_annotated.tsv")

stopifnot(file.exists(f_voom))
stopifnot(file.exists(f_fit2))
stopifnot(file.exists(f_contrasts))
stopifnot(file.exists(f_logcpm))
stopifnot(dir.exists(msigdb_dir))

gmt_files <- list(
  reactome = Sys.glob(file.path(msigdb_dir, "m2.cp.reactome.v*.Mm.symbols.gmt"))[1],
  hallmark = Sys.glob(file.path(msigdb_dir, "mh.all.v*.Mm.symbols.gmt"))[1],
  gobp = Sys.glob(file.path(msigdb_dir, "m5.go.bp.v*.Mm.symbols.gmt"))[1]
)

stopifnot(file.exists(gmt_files$reactome))
stopifnot(file.exists(gmt_files$hallmark))
stopifnot(file.exists(gmt_files$gobp))

# -----------------------------
# Load model objects
# -----------------------------

v <- readRDS(f_voom)
fit2 <- readRDS(f_fit2)
contrast_matrix <- as.matrix(read.csv(f_contrasts, row.names = 1, check.names = FALSE))
logcpm_annotated <- data.table::fread(f_logcpm) %>% as.data.frame()

stopifnot(!is.null(v$E))
stopifnot(!is.null(v$design))
stopifnot(!is.null(fit2$t))
stopifnot(!is.null(fit2$coefficients))
stopifnot(all(CONTRASTS_TO_RUN %in% colnames(fit2$coefficients)))
stopifnot(all(CONTRASTS_TO_RUN %in% colnames(contrast_matrix)))
stopifnot(all(c("gene_id", "gene_symbol") %in% colnames(logcpm_annotated)))

contrast_matrix <- contrast_matrix[colnames(v$design), , drop = FALSE]

gene_match <- match(rownames(v$E), logcpm_annotated$gene_id)
if (anyNA(gene_match)) {
  stop("Some voom rownames are missing from logcpm_annotated.tsv gene_id.")
}

symbol_for_row <- logcpm_annotated$gene_symbol[gene_match]
message("Symbol coverage among modelled genes: ",
        round(mean(!is.na(symbol_for_row)) * 100, 2), "%")

# -----------------------------
# Load and prepare gene sets
# -----------------------------

read_gmt <- function(path) {
  lines <- readLines(path)
  split_lines <- strsplit(lines, "\t", fixed = TRUE)
  gene_sets <- lapply(split_lines, function(x) unique(x[-c(1, 2)]))
  names(gene_sets) <- vapply(split_lines, `[`, character(1), 1)
  gene_sets
}

filter_gs_size <- function(gene_sets, min_size = MIN_GS_SIZE, max_size = MAX_GS_SIZE) {
  keep <- vapply(
    gene_sets,
    function(x) length(x) >= min_size && length(x) <= max_size,
    logical(1)
  )
  gene_sets[keep]
}

map_symbols_to_gene_ids <- function(symbols, gene_ids, row_symbols) {
  unique(gene_ids[!is.na(row_symbols) & row_symbols %in% symbols])
}

gene_sets_symbol <- list(
  reactome = read_gmt(gmt_files$reactome),
  hallmark = read_gmt(gmt_files$hallmark),
  gobp = read_gmt(gmt_files$gobp)
)

if (any(grepl("^R-HSA-", readLines(gmt_files$reactome, n = 5)))) {
  stop("Reactome GMT appears to contain human Reactome IDs. Expected mouse symbol GMT.")
}

gene_sets_symbol <- lapply(gene_sets_symbol, filter_gs_size)

gene_sets_ids <- lapply(gene_sets_symbol, function(collection) {
  mapped <- lapply(
    collection,
    map_symbols_to_gene_ids,
    gene_ids = rownames(v$E),
    row_symbols = symbol_for_row
  )
  mapped[vapply(mapped, length, integer(1)) > 0]
})

# -----------------------------
# Enrichment functions
# -----------------------------

ids_to_index <- function(ids, rownames_E) {
  idx <- match(ids, rownames_E)
  unique(idx[!is.na(idx)])
}

run_camera_one <- function(gene_sets_one_collection, contrast_name, collection_name) {
  index <- lapply(gene_sets_one_collection, ids_to_index, rownames_E = rownames(v$E))
  index <- index[vapply(index, length, integer(1)) > 0]

  result <- limma::camera(
    y = v$E,
    index = index,
    design = v$design,
    contrast = contrast_matrix[, contrast_name]
  )

  result %>%
    as.data.frame() %>%
    tibble::rownames_to_column("pathway") %>%
    mutate(collection = collection_name, contrast = contrast_name) %>%
    relocate(collection, contrast, pathway)
}

make_fgsea_stats <- function(contrast_name) {
  stats_table <- tibble(
    gene_id = rownames(v$E),
    gene_symbol = symbol_for_row,
    t = as.numeric(fit2$t[, contrast_name])
  ) %>%
    filter(!is.na(gene_symbol)) %>%
    group_by(gene_symbol) %>%
    slice_max(order_by = abs(t), n = 1, with_ties = FALSE) %>%
    ungroup()

  stats <- stats_table$t
  names(stats) <- stats_table$gene_symbol
  sort(stats, decreasing = TRUE)
}

run_fgsea_one <- function(gene_sets_one_collection, contrast_name, collection_name) {
  stats <- make_fgsea_stats(contrast_name)

  pathways <- lapply(gene_sets_one_collection, function(x) intersect(unique(x), names(stats)))
  pathways <- pathways[vapply(pathways, length, integer(1)) >= MIN_GS_SIZE]

  result <- fgsea::fgseaMultilevel(
    pathways = pathways,
    stats = stats,
    minSize = MIN_GS_SIZE,
    maxSize = MAX_GS_SIZE
  )

  result %>%
    as.data.frame() %>%
    mutate(collection = collection_name, contrast = contrast_name) %>%
    relocate(collection, contrast, pathway)
}

# -----------------------------
# Run enrichment
# -----------------------------

camera_results <- list()
fgsea_results <- list()

for (collection in names(gene_sets_ids)) {
  ensure_dir(file.path(out_camera, collection))
  ensure_dir(file.path(out_fgsea, collection))
  ensure_dir(file.path(out_fig, collection))

  for (contrast_name in CONTRASTS_TO_RUN) {
    message("camera: ", collection, " | ", contrast_name)
    camera_result <- run_camera_one(gene_sets_ids[[collection]], contrast_name, collection)
    readr::write_tsv(
      camera_result,
      file.path(out_camera, collection, paste0("camera_", contrast_name, ".tsv"))
    )
    camera_results[[paste(collection, contrast_name, sep = "__")]] <- camera_result

    message("fgsea: ", collection, " | ", contrast_name)
    fgsea_result <- run_fgsea_one(gene_sets_symbol[[collection]], contrast_name, collection)
    readr::write_tsv(
      fgsea_result,
      file.path(out_fgsea, collection, paste0("fgsea_", contrast_name, ".tsv"))
    )
    fgsea_results[[paste(collection, contrast_name, sep = "__")]] <- fgsea_result
  }
}

# -----------------------------
# Summary and leading-edge tables
# -----------------------------

for (contrast_name in CONTRASTS_TO_RUN) {
  camera_all <- bind_rows(lapply(names(gene_sets_ids), function(collection) {
    camera_results[[paste(collection, contrast_name, sep = "__")]]
  })) %>%
    arrange(FDR)

  readr::write_tsv(
    camera_all,
    file.path(out_summary, paste0("camera_allcollections_", contrast_name, ".tsv"))
  )

  driver_table <- bind_rows(lapply(names(gene_sets_ids), function(collection) {
    camera_result <- camera_results[[paste(collection, contrast_name, sep = "__")]]
    fgsea_result <- fgsea_results[[paste(collection, contrast_name, sep = "__")]]

    fgsea_small <- fgsea_result %>%
      select(pathway, NES, padj, leadingEdge) %>%
      mutate(
        leadingEdge = vapply(
          leadingEdge,
          function(x) paste(x, collapse = ";"),
          character(1)
        )
      )

    camera_result %>%
      filter(FDR < FDR_CUTOFF) %>%
      left_join(fgsea_small, by = "pathway") %>%
      mutate(
        NES = ifelse(is.na(NES), NA_real_, NES),
        padj = ifelse(is.na(padj), NA_real_, padj),
        leadingEdge = ifelse(is.na(leadingEdge), "", leadingEdge)
      )
  })) %>%
    arrange(FDR)

  readr::write_tsv(
    driver_table,
    file.path(out_leading, paste0("driver_genes_leadingedge_", contrast_name, ".tsv"))
  )
}

# -----------------------------
# Camera top-term plots
# -----------------------------

clean_pathway_label <- function(x) {
  x %>%
    str_replace("^HALLMARK_", "") %>%
    str_replace("^REACTOME_", "") %>%
    str_replace("^GOBP_", "") %>%
    str_replace("^GO_BP_", "") %>%
    str_replace("^GO_BIOLOGICAL_PROCESS_", "") %>%
    str_replace_all("_", " ")
}

plot_camera_top <- function(camera_result, contrast_name, collection_name) {
  plot_table <- camera_result %>%
    arrange(FDR) %>%
    slice_head(n = TOP_N_PLOT) %>%
    mutate(
      pathway_label = clean_pathway_label(pathway),
      pathway_label = factor(pathway_label, levels = rev(pathway_label)),
      signed_score = ifelse(Direction == "Up", -log10(FDR), log10(FDR))
    )

  ggplot(plot_table, aes(x = signed_score, y = pathway_label, fill = Direction)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = c("Up" = "#D73027", "Down" = "#4575B4")) +
    labs(
      title = paste0(collection_name, ": ", contrast_name),
      x = "Signed -log10 FDR",
      y = NULL
    ) +
    theme_classic(base_size = 9) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 9),
      axis.text.y = element_text(size = 7)
    )
}

for (collection in names(gene_sets_ids)) {
  for (contrast_name in CONTRASTS_TO_RUN) {
    plot_obj <- plot_camera_top(
      camera_results[[paste(collection, contrast_name, sep = "__")]],
      contrast_name,
      collection
    )

    out_base <- file.path(
      out_fig,
      collection,
      paste0("camera_top", TOP_N_PLOT, "_", contrast_name)
    )

    ggsave(paste0(out_base, ".pdf"), plot_obj, width = 7, height = 4.8)
    ggsave(paste0(out_base, ".png"), plot_obj, width = 7, height = 4.8, dpi = 300)
  }
}

message("Pathway enrichment complete.")
message("camera results: ", out_camera)
message("fgsea results: ", out_fgsea)
message("leading-edge tables: ", out_leading)
message("figures: ", out_fig)
