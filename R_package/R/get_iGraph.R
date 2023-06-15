#' Create iGraph object from a symmetrical adjacency matrix, annotate it with centrality attributes and return the annotated iGraph.
#'
#' @importFrom igraph graph.adjacency vertex_attr betweenness closeness degree eigen_centrality transitivity edge_attr list.vertex.attributes
#'
#' @param adj_mat Binary consensus matrix computed by using get_ranked_consensus_binary_matrix() function.
#' @return Annotated igraph object representing the binary consensus matrix computed by get_ranked_consensus_binary_matrix() function.
#' @export
get_iGraph <- function(adj_mat){
  if(is.null(adj_mat)){
    stop("Input Adjacency Matrix is NULL!")
  }
  
  if(!isSymmetric(adj_mat)){
    print("Matrix is not symmetric. Getting symmetric matrix!")
    #adj_mat <- get_symmetric_matrix(adj_mat)
    adj_mat <- pmax(adj_mat, t(adj_mat))
  }
  
  iG <- igraph::graph.adjacency(adj_mat, mode="undirected", weighted=NULL)
  
  igraph::vertex_attr(iG, name="betweenness") <- as.vector(igraph::betweenness(iG))
  igraph::vertex_attr(iG, name="closeness") <- as.vector(igraph::closeness(iG, normalized=TRUE))
  igraph::vertex_attr(iG, name="degree") <- as.vector(igraph::degree(iG))
  igraph::vertex_attr(iG, name="eigenvector") <- as.vector(igraph::eigen_centrality(iG)$vector)
  igraph::vertex_attr(iG, name="cc") <- igraph::transitivity(iG, type="local", isolates="zero")
  
  igraph::vertex_attr(iG, name="color") <- "lightgray"
  igraph::vertex_attr(iG, name="highlightcolor") <- "darkgray"
  igraph::edge_attr(iG, name="color") <- "lightgray"
  igraph::edge_attr(iG, name="highlightcolor") <- "darkgray"
  
  #print(igraph::list.vertex.attributes(iG))
  iG
}
