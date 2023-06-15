#' Map OpenTargets taget drugs on a Gene Network
#'
#' This function maps OpenTargets target drugs onto the nodes of a given gene network.
#' 
#' @param Graph An igraph object.
#' @param parsed_opentargets A dataframe representing the OpenTargets database (e.g. obtained from the parse_opentargets_json).
#' @return A dataframe reporting the drugs and diseases mapped on the provided graph.
#' @export
map_opentargets_on_graph <- function(Graph, parsed_opentargets) {
  
  if(is.null(Graph)) {
    stop("Error: please provide a valid igraph object!")
  }
  
  if(!class(Graph)=="igraph") {
    stop("Error: the Graph object is not of class igraph!")
  }
  
  if(is.null(parsed_opentargets)) {
    stop("Error: please provide the parsed OpenTargets database (e.g. obtained from the parse_opentargets_json function)")
  }
  
  
  if(!class(parsed_opentargets)=="data.frame") {
    stop("Error: please provide an object of class data.frame!")
  }
  
  map_opentargets <- parsed_opentargets[which(parsed_opentargets$dat.target.gene_info.symbol %in% names(V(Graph))),]
  map_opentargets$dat.drug.molecule_name <- tolower(map_opentargets$dat.drug.molecule_name)
  map_opentargets$dat.drug.molecule_name <- gsub(map_opentargets$dat.drug.molecule_name, pattern = " ", replacement = ".")
  map_opentargets$dat.disease.efo_info.label <- gsub(map_opentargets$dat.disease.efo_info.label, pattern = " ", replacement = ".")
  
  return(map_opentargets)
  
}

