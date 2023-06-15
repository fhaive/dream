#' Creates a squared matrix of average jaccard indices between subsets of neighbors of a certain order of couples of drugs. Note that each drug may have multiple targets, so the average shortest paths of all the targets of a certain drug is computed.
#' 
#' @param drug_target_df Dataframe of at least two columns containing drug names and corresponding targets.
#' @param drugs_col Column index of the dataframe drug_target_df where the drugs are stored.
#' @param targets_col Column index of the dataframe drug_target_df where the drugs are stored.
#' @param graph_rank An ordered named vector of a graph's nodes. Typically obtained from the get_graph_named_ranklist function.
#' @return A rank of drugs based on the average ranks of their targets.
#' @export
get_drug_avg_rank_from_graph <-  function(drug_target_df, drugs_col, targets_col, graph_rank) {
  
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
  
  if(is.null(graph_rank)) {
    stop("Error: please provide a named ranked vector of nodes!")
  }
  
  
  drug_inform_vec <- c()
  
  for(drugs in 1:length(unique(drug_target_df[, drugs_col]))) {
    targets <- as.vector(unique(drug_target_df[, targets_col][which(drug_target_df[, drugs_col]==unique(drug_target_df[, drugs_col])[drugs])]))
    inform_rank <- mean(which(names(graph_rank) %in% targets))
    drug_inform_vec <- c(drug_inform_vec, inform_rank)
    
  }
  
  names(drug_inform_vec) <- unique(drug_target_df[, drugs_col])
  
  return(drug_inform_vec)
  
}

