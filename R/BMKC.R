#' @importFrom WGCNA bicor
#' @importFrom igraph graph_from_adjacency_matrix max_cliques
NULL

#' Run the BMKC analysis to find differentially coexpressed gene modules.
#'
#' Implements the BMKC (Biweight Midcorrelation and k-Clique) method.
#' This method identifies modules of genes that are highly coexpressed in one
#' condition but not the other.
#'
#' @param expr_matrix A numeric matrix where rows are genes and columns are samples.
#' @param condition A factor or character vector with two levels, matching the columns of expr_matrix.
#' @param T1 The absolute correlation threshold for the first condition.
#' @param T2 The absolute correlation threshold for the second condition.
#' @param min_clique_size The minimum number of genes for a clique to be reported.
#' @return A list of character vectors, where each vector is a gene module (clique).
#' @export
run_bmkc <- function(expr_matrix, condition, T1 = 0.7, T2 = 0.3, min_clique_size = 4) {

  conditions <- unique(condition)
  if (length(conditions) != 2) stop("Condition vector must have exactly two levels.")

  cond1 <- conditions[1]; cond2 <- conditions[2]

  # Calculate bicor matrices for both conditions
  cor_1 <- WGCNA::bicor(t(expr_matrix[, condition == cond1]))
  cor_2 <- WGCNA::bicor(t(expr_matrix[, condition == cond2]))

  # Apply "Differential Coexpression Threshold" strategy
  # We look for high correlation in condition 1 and low in condition 2
  adj_matrix <- (abs(cor_1) >= T1) & (abs(cor_2) <= T2)

  # Convert to a numeric 0/1 matrix
  adj_matrix <- apply(adj_matrix, 2, as.numeric)
  rownames(adj_matrix) <- colnames(adj_matrix)

  # Use igraph to find maximal cliques
  graph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)
  cliques <- igraph::max_cliques(graph, min = min_clique_size)

  # Convert igraph vertex objects to gene names
  gene_modules <- lapply(cliques, function(clique) {
    names(clique)
  })

  return(gene_modules)
}
