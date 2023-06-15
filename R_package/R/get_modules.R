#' Get modules from the main igraph.
#'
#' Extract modules from the main graph by using a specified community detection algorithm from igraph.
#'
#' @importFrom igraph cluster_walktrap cluster_spinglass cluster_louvain cluster_fast_greedy
#'
#' @param iGraph igraph object representing the main graph from which the subgraph should be extracted.
#' @param method community detection method from igraph.
#' @return communities object cotaining communities identified from igraph object by using a specific method
#' @export
get_modules <- function(iGraph, method="walktrap")
{
  switch(method,
         "walktrap" = {
           optimalStep <- c(2:10)[which.max(sapply(c(2:10), function(x){igraph::modularity(igraph::cluster_walktrap(iGraph, step=x))}))]
           igraph::cluster_walktrap(iGraph, step=optimalStep)
         },
         "spinglass" = igraph::cluster_spinglass(iGraph),
         "louvain" = igraph::cluster_louvain(iGraph),
         "greedy" = igraph::cluster_fast_greedy(iGraph)
  ) 
}
