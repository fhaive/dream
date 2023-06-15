#' Get a named list from an igraph object using ranked vectors
#'
#' This function takes an igraph object and a list of ranked vectors as input, and generates a named list based on the provided ranked vectors. 
#' The function assigns the ranked vectors as vertex attributes in the graph and then calculates a ranked gene list using specified network attributes. 
#' The resulting ranked gene list is then returned as a named list
#'
#' @param Graph An igraph object representing the graph.
#' @param ranked_vectors A list of one or more ranked vectors.
#' @return A named list based on the provided ranked vectors.
#' @export


get_graph_named_ranklist <- function(Graph, ranked_vectors){
  
  if(!class(Graph)=="igraph") {
    stop("Error: please provide an object of class igraph!")
  }
  
  if(is.null(Graph)) {
    stop("Error: please provide an igraph object!")
  }
  
  if(is.null(ranked_vectors)) {
    stop("Error: please provide a list of one or more ranked vectors!")
  }
  
  if(!class(ranked_vectors)=="list") {
    stop("Error: please provide a list of ranked vectors!")
  }
  
  
  if(length(ranked_vectors)==1) {
    Graph=igraph::set.vertex.attribute(Graph, "ranked_vectors_1", value = ranked_vectors[[1]][igraph::V(Graph)$name])
    
  }else if(length(ranked_vectors)>1) {
    for (rank in 1:length(ranked_vectors)) {
      Graph=igraph::set.vertex.attribute(Graph, paste0("ranked_vectors_", rank), value = ranked_vectors[[rank]][igraph::V(Graph)$name])
      
    }
    
  }
  
  net_attr <- c("betweenness", "cc", "degree", "closeness", "eigenvector", paste0("ranked_vectors_", 1:length(ranked_vectors)))
  
  ranked.graph=get_ranked_gene_list(iGraph = Graph, rank_list_attr = net_attr, debug_output = TRUE)
  named_ranklist=setNames(c(1:length(ranked.graph)), ranked.graph)
  
  return(named_ranklist)
  
}