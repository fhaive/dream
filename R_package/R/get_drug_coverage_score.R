#' Creates a squared matrix of average jaccard indices between subsets of neighbors of a certain order of couples of drugs. Note that each drug may have multiple targets, so the average shortest paths of all the targets of a certain drug is computed.
#' 
#' @param Graph An object of class igraph.
#' @param drug_target_df Dataframe of at least two columns containing drug names and corresponding targets.
#' @param drugs_col Column index of the dataframe drug_target_df where the drugs are stored.
#' @param targets_col Column index of the dataframe drug_target_df where the drugs are stored.
#' @param graph_rank An ordered named vector of a graph's nodes. Typically obtained from the get_graph_named_ranklist function.
#' @param order Order of the neighborhood to be considered.
#' @return A rank of drugs based on the average ranks of their targets.
#' @export

get_drug_coverage_score <- function(Graph, drug_target_df, drugs_col, targets_col, graph_rank, order = 1) {
  darios_drug_score <- c()
  target_jacc_list <- list()
  ji <- c()
  
  for (drugs in 1:length(unique(drug_target_df[, drugs_col]))) {
    target_jacc_list <- list()
    targets <- as.vector(unique(drug_target_df[, targets_col][which(drug_target_df[, drugs_col]==unique(drug_target_df[, drugs_col])[drugs])]))
    print(paste("Computing graph coverage for", unique(drug_target_df[, drugs_col])[drugs], "based on", length(targets), "drug targets", sep = " "))
    inform_rank <- median(which(names(graph_rank) %in% targets))
    if(length(targets)==1) {
      ji <- 1
    }else{
      for (target in 1:length(targets)) {
        
        target_jacc_list[[target]] <- igraph::ego(Graph, order = order, nodes = targets[target])
        
      }
      
      for(jacc1 in 1:length(target_jacc_list)) {
        for(jacc2 in 1:length(target_jacc_list)) {
          ji <- c(ji, length(intersect(target_jacc_list[[jacc1]], target_jacc_list[[jacc2]]))/length(union(target_jacc_list[[jacc1]], target_jacc_list[[jacc2]])))
          
        }
        
      }
      
      
    }
    
    darios_drug_score <- c(darios_drug_score, length(targets)/inform_rank*median(ji))
    
  }
  
  names(darios_drug_score) <- unique(drug_target_df[, drugs_col])
  return(darios_drug_score)
  
}
