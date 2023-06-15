#' Creates a squared matrix of average shortes paths between couples of drugs. Note that each drug may have multiple targets, so the average shortest paths of all the targets of a certain drug is computed.
#' 
#' @param Graph An object of class "igraph".
#' @param drug_target_df Dataframe of at least two columns containing drug names and corresponding targets.
#' @param drugs_col Column index of the dataframe drug_target_df where the drugs are stored.
#' @param targets_col Column index of the dataframe drug_target_df where the drugs are stored.
#' @return A squared matrix of average shortest paths between couples of drug targets.
#' @export

get_avg_shortest_path <- function(Graph, drug_target_df, drugs_col, targets_col) {
  
  if(is.null(Graph) | !class(Graph)=="igraph"){
    stop("Error: please provide an igraph object!")
  }
  
  if(is.null(drug_target_df)) {
    stop("Error: please provide a dataframe with at least two columns: one for the drugs and another for the corresponding targets!")
  }
  
  if(is.null(drugs_col)) {
    stop("Error: please provide the column containing the drugs!")
  }
  
  if(is.null(targets_col)) {
    stop("Error: please provide the column containing the targets!")
  }
  
  if(!class(drugs_col)=="numeric") {
    stop("Error: please provide the index of the drug column!")
  }
  
  if(!class(targets_col)=="numeric") {
    stop("Error: please provide the index of the targets column!")
  }
  
  shortpath <- matrix(NA, nrow = length(unique(tolower(drug_target_df[, drugs_col]))), ncol = length(unique(tolower(drug_target_df[, drugs_col]))), dimnames = list(unique(tolower(drug_target_df[, drugs_col])), unique(tolower(drug_target_df[, drugs_col]))))
  
  
  for (drug1 in 1:length(unique(drug_target_df[, drugs_col]))) {
    print(paste("Computing avg shortest path for drug ", drug1 , ": ", unique(drug_target_df[, drugs_col])[drug1]))
    gene1 <- unique(as.vector(drug_target_df[, targets_col][which(drug_target_df[, drugs_col] %in% unique(tolower(drug_target_df[, drugs_col]))[drug1])]))
    
    for (drug2 in 1:length(unique(drug_target_df[, drugs_col]))){
      spathtotal <- c()
      gene2 <- unique(as.vector(drug_target_df[, targets_col][which(drug_target_df[, drugs_col] %in% unique(tolower(drug_target_df[, drugs_col]))[drug2])]))
      
      for(genes in 1:length(gene1)){
        spathgene <- igraph::get.shortest.paths(Graph, from = as.character(gene1)[genes], to = as.character(gene2))
        spathtotal <- c(spathtotal, mean(sapply(spathgene$vpath, function(x) length(x)-1)))
      }
      spathmean <- mean(spathtotal)
      shortpath[drug1,drug2] <- spathmean
    }
    
  }
  
  return(shortpath)
}

