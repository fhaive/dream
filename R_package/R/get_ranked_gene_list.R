#' Get top ranked candidates computed by Borda on the provided list of attributes and the cutoff.
#'
#' get_ranked_gene_list uses the annotation associated with each vertex to rank them and generated a list of vertices names ordered by rank. 
#' By default the annotations calculated by INfORM are "betweenness", "cc", "degree", "eccentricity", "closeness" & "eigenvector". The vertices are 
#' ranked by each annotation separately nad then the ranks are unified by means of Borda(). User must provide the annotations to use in ranking scheme, 
#' if the user has custom annotations such as "score" for differential gene expression then it can also be used for ranking.
#'
#' @importFrom igraph get.vertex.attribute
#' @importFrom TopKLists Borda
#' @importFrom utils head
#'
#' @param iGraph igraph object created from the binary consensus matrix computed by using get_ranked_consensus_binary_matrix() function.
#' @param rank_list_attr Vector of network attributes/scores to use for generating a combined gene rank
#' @param debug_output Print help and status messages to help debug the running of the function default:FALSE.
#' @return vector of genes ordered by rank, based on the selected attributes associated to each gene.
#' @export
get_ranked_gene_list <- function(iGraph, rank_list_attr, debug_output=FALSE){
  if(is.null(rank_list_attr)){
    stop("List of attributes for ranking is NULL!")
  }
  
  attrOrdMat <- list()
  for(a in 1:length(rank_list_attr)){
    val_ord=TRUE
    
    if(debug_output==TRUE)
      print(paste("Using attribute: ", rank_list_attr[a], sep=""))
    
    attrValueVector <- igraph::get.vertex.attribute(iGraph, rank_list_attr[a])
    
    if(rank_list_attr[a] == "score")
      attrValueVector <- abs(attrValueVector)
    
    attrOrdList <- cbind(igraph::V(iGraph)$name, attrValueVector)[order(attrValueVector, decreasing=val_ord)]
    
    if(debug_output==TRUE)
      print(utils::head(attrOrdList))
    
    attrOrdMat[[a]] <- attrOrdList
  }
  
  attrBorda <- TopKLists::Borda(attrOrdMat)
  
  if(debug_output==TRUE)
    print(utils::head(attrBorda$TopK, n=10))
  
  return(attrBorda$TopK$median)
}
