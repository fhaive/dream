#' Creates a squared matrix of average jaccard indices between subsets of neighbors of a certain order of couples of drugs. Note that each drug may have multiple targets, so the average shortest paths of all the targets of a certain drug is computed.
#' 
#' @param Graph An object of class "igraph".
#' @param drug_target_df Dataframe of at least two columns containing drug names and corresponding targets.
#' @param drugs_col Column index of the dataframe drug_target_df where the drugs are stored.
#' @param targets_col Column index of the dataframe drug_target_df where the drugs are stored.
#' @param order Order of the neighborhood to be considered.
#' @return A squared matrix of average shortest paths between couples of drug targets.
#' @export


get_avg_jacc_indices <- function(Graph, drug_target_df, drugs_col, targets_col, order = 1) {
  
  if(is.null(Graph) | !class(Graph)=="igraph") {
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
  
  if(!class(order)=="numeric") {
    stop("Error: please provide the neighborhood order as a numeric value!")
  }
  
  jacc_drugs <- matrix(NA, nrow = length(unique(tolower(drug_target_df[, drugs_col]))), ncol = length(unique(tolower(drug_target_df[, drugs_col]))), dimnames = list(unique(tolower(drug_target_df[, drugs_col])), unique(tolower(drug_target_df[, drugs_col]))))
  
  for (drug1 in 1:length(unique(drug_target_df[, drugs_col]))) {
    
    #target drug1
    print(paste("Computing avg Jaccard index for drug ", drug1 , ": ", unique(drug_target_df[, drugs_col])[drug1]))
    gene1 <- unique(as.vector(drug_target_df[, targets_col][which(tolower(drug_target_df[, drugs_col]) %in% unique(tolower(drug_target_df[, drugs_col]))[drug1])]))
    
    for (drug2 in 1:length(unique(drug_target_df[, drugs_col]))){
      
      #targets drug2
      gene2 <- unique(as.vector(drug_target_df[, targets_col][which(tolower(drug_target_df[, drugs_col]) %in% unique(tolower(drug_target_df[, drugs_col]))[drug2])]))
      
      jacc_vec <- c()
      #for each target of drug1 
      for(genes1 in 1:length(gene1)){
        # find close genes
        gene1_first_degree <- igraph::ego(Graph, order = order, nodes = gene1[genes1])
        #for each target of drug2
        for (genes2 in 1:length(gene2)) {
          #find close genes
          gene2_first_degree <- igraph::ego(Graph, order = order, nodes = gene2[genes2])
          
          # jaccard index between the genes
          jacc_index <- length(intersect(as.character(names(gene1_first_degree[[1]])), as.character(names(gene2_first_degree[[1]]))))/length(union(as.character(names(gene1_first_degree[[1]])), as.character(names(gene2_first_degree[[1]]))))
          
          jacc_vec <- c(jacc_vec, jacc_index)
          
        }
      }
      jaccmean <- mean(jacc_vec)
      jacc_drugs[drug1, drug2] <- jaccmean
      
    }
    
  }
  return(jacc_drugs)
}

