#' Compute a squared matrix reporting the Jaccard similarity indices among couples of drug targets neighborhoods.
#' 
#' @param Graph An igraph object.
#' @param target_nodes A vector of gene symbols likely to be drug targets.
#' @param neighborhood The order of neighborhood of the provided target nodes.
#' @return A squared matrix reporting the Jaccard similarity indices among couples of drug targets neighborhoods.
#' @export

get_jacc_drug_targets_neighbors <- function(Graph, target_nodes, neighborhood) {
  
  drug_target_mat_L1000 <- matrix("NA", ncol = length(target_nodes), nrow = length(target_nodes), dimnames = list(target_nodes, target_nodes))
  jacc_index <- c()
  
  
  for (i in 1:length(as.character(unique(target_nodes)))) {
    degree_i <- igraph::ego(conGraph, order = neighborhood, nodes = as.character(unique(target_nodes)[i]))
    for (k in 1:length(as.character(unique(target_nodes)))) {
      degree_k <- igraph::ego(conGraph, order = neighborhood, nodes = as.character(unique(target_nodes)[k]))
      jacc_index <- length(intersect(as.character(names(degree_i[[1]])), as.character(names(degree_k[[1]]))))/length(union(as.character(names(degree_i[[1]])), as.character(names(degree_k[[1]]))))
      drug_target_mat_L1000[i, k] <- jacc_index
    }
    
  }
  
  return(drug_target_mat_L1000)
  
}

