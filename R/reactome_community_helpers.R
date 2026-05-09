suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(igraph)
  library(readr)
})

parse_leadingedge <- function(x) {
  if (length(x) == 0 || is.na(x) || x == "") return(character(0))
  out <- unlist(strsplit(as.character(x), split = ";|,|\\t|\\r\\n|\\n|\\r"))
  out <- stringr::str_trim(out)
  out <- out[out != ""]
  unique(out)
}

jaccard_vec <- function(a, b) {
  a <- unique(as.character(a)); a <- a[a != ""]
  b <- unique(as.character(b)); b <- b[b != ""]
  if (length(a) == 0 || length(b) == 0) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

build_jaccard_matrix <- function(le_list) {
  n <- length(le_list)
  nm <- names(le_list)

  if (anyDuplicated(nm)) {
    stop("Pathway names must be unique before building the Jaccard matrix.")
  }

  M <- matrix(0, n, n, dimnames = list(nm, nm))

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      M[i, j] <- jaccard_vec(le_list[[i]], le_list[[j]])
    }
  }

  M
}

graph_from_similarity <- function(sim_mat, threshold = 0.20) {
  diag(sim_mat) <- 0

  idx <- which(sim_mat > threshold, arr.ind = TRUE)
  idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]

  if (nrow(idx) == 0) {
    g <- igraph::make_empty_graph(n = nrow(sim_mat), directed = FALSE)
    g <- igraph::set_vertex_attr(g, "name", value = rownames(sim_mat))
    return(g)
  }

  edges_df <- data.frame(
    from = rownames(sim_mat)[idx[, 1]],
    to = colnames(sim_mat)[idx[, 2]],
    weight = sim_mat[idx],
    stringsAsFactors = FALSE
  )

  igraph::graph_from_data_frame(
    edges_df,
    directed = FALSE,
    vertices = data.frame(name = rownames(sim_mat))
  )
}

cluster_similarity_graph <- function(g) {
  if (igraph::gorder(g) == 0) return(integer(0))

  if (igraph::gsize(g) == 0) {
    memb <- seq_len(igraph::gorder(g))
    names(memb) <- igraph::V(g)$name
    return(memb)
  }

  clu <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
  memb <- clu$membership
  names(memb) <- igraph::V(g)$name
  memb
}

prepare_reactome_table <- function(leading_edge_table,
                                   fdr_cutoff = 0.05,
                                   min_ngenes = 10,
                                   max_ngenes = 300) {
  required_cols <- c(
    "collection", "contrast", "pathway", "NGenes", "Direction",
    "PValue", "FDR", "NES", "padj", "leadingEdge"
  )

  if (!all(required_cols %in% colnames(leading_edge_table))) {
    stop("Leading-edge table has unexpected columns.")
  }

  leading_edge_table %>%
    mutate(
      collection = tolower(as.character(collection)),
      LE_set = lapply(leadingEdge, parse_leadingedge)
    ) %>%
    filter(collection == "reactome") %>%
    filter(!is.na(FDR), FDR < fdr_cutoff) %>%
    filter(!is.na(NES)) %>%
    filter(!is.na(NGenes), NGenes >= min_ngenes, NGenes <= max_ngenes) %>%
    filter(Direction %in% c("Up", "Down")) %>%
    filter(lengths(LE_set) > 0) %>%
    arrange(FDR, desc(abs(NES)))
}

build_communities_for_table <- function(reactome_table, jaccard_threshold = 0.20) {
  if (nrow(reactome_table) < 2) {
    stop("At least two Reactome pathways are required for community analysis.")
  }

  if (anyDuplicated(reactome_table$pathway)) {
    duplicated_pathways <- unique(reactome_table$pathway[duplicated(reactome_table$pathway)])
    stop(
      "Duplicate Reactome pathway names found after filtering: ",
      paste(head(duplicated_pathways, 10), collapse = "; ")
    )
  }

  le_list <- reactome_table$LE_set
  names(le_list) <- reactome_table$pathway

  sim_mat <- build_jaccard_matrix(le_list)
  graph <- graph_from_similarity(sim_mat, threshold = jaccard_threshold)
  membership <- cluster_similarity_graph(graph)

  node_table <- reactome_table %>%
    mutate(
      community = unname(membership[pathway]),
      degree = igraph::degree(graph)[pathway]
    ) %>%
    select(
      pathway, collection, contrast, NGenes, Direction,
      PValue, FDR, NES, padj, leadingEdge, community, degree
    ) %>%
    arrange(community, FDR)

  edge_table <- igraph::as_data_frame(graph, what = "edges") %>%
    rename(from = from, to = to, jaccard = weight)

  community_summary <- node_table %>%
    group_by(community) %>%
    summarise(
      directions = paste(sort(unique(Direction)), collapse = ";"),
      n_pathways = n(),
      n_up = sum(Direction == "Up"),
      n_down = sum(Direction == "Down"),
      median_NES = median(NES, na.rm = TRUE),
      min_FDR = min(FDR, na.rm = TRUE),
      representative_pathways = paste(head(pathway[order(FDR, -abs(NES))], 5), collapse = " | "),
      .groups = "drop"
    ) %>%
    arrange(min_FDR, desc(n_pathways))

  community_genes <- node_table %>%
    select(community, pathway, leadingEdge) %>%
    tidyr::separate_rows(leadingEdge, sep = ";") %>%
    mutate(gene_symbol = stringr::str_trim(leadingEdge)) %>%
    filter(gene_symbol != "") %>%
    count(community, gene_symbol, name = "recurrence") %>%
    arrange(community, desc(recurrence), gene_symbol)

  list(
    node_table = node_table,
    edge_table = edge_table,
    community_summary = community_summary,
    community_genes = community_genes
  )
}
