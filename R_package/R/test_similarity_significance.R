#' Performs a permutation test to assess the significance of the jaccard scores among couples of drug targets (e.g. obtained from the get_jacc_drug_targets_neighbors function).
#' 
#' @param drug_target_mat A squared matrix containing similarity scores between drug targets pairs (e.g. obtained from the get_jacc_drug_targets_neighbors function).
#' @param Graph An igraph object.
#' @param nIter Number of permutations to be performed.
#' @return A squared matrix reporting the statistical significance of the similarity measures (e.g. Jaccard indices) among couples of drug targets neighborhoods.
#' @export

test_similarity_significance <- function(drug_target_mat, Graph, nIter=1000) {
  
  if(is.null(drug_target_mat)){
    stop("Error: please provide a squared drug target matrix!")
  }
  
  if(!class(drug_target_mat)=="matrix") {
    stop("Error: please provide an object of class matrix!")
  }
  
  if(is.null(Graph) | !class(Graph)=="igraph"){
    stop("Error: please provide an igraph object!")
  }
  
  if(!class(nIter)=="numeric"){
    stop("Error: please provide the number of iterations!")
  }
  
  dist <- c()
  jacc_original <- c()
  nIter = nIter
  PvalueMatL1000 = matrix(NA, nrow = length(rownames(drug_target_mat)), ncol = length(rownames(drug_target_mat)) )
  
  
  pb = txtProgressBar(min = 0, max = length(rownames(drug_target_mat)), style = 3)
  
  for (i in 1:length(rownames(drug_target_mat))) {
    for (k in 1:length(colnames(drug_target_mat))) {
      print(paste(i, k, sep=","))
      jacc_permutated <-  c()
      jacc_original <- as.numeric(drug_target_mat[i,k])
      
      # start_time <- Sys.time()
      
      pbin = txtProgressBar(min = 1, max=nIter, style=2)
      for(nperm in 1:nIter){
        dist <- igraph::distances(graph = Graph, v = colnames(drug_target_mat)[i], to = colnames(drug_target_mat)[k], mode = "all")
        
        pick_first_node <- base::sample(names(V(Graph))[-c(which(names(V(Graph))==colnames(drug_target_mat)[i]), which(names(V(Graph))==colnames(drug_target_mat)[k]))], size = 1)
        first_node <- igraph::ego(Graph, order = dist, nodes = pick_first_node)
        pick_second_node <- base::sample(as.character(names(first_node[[1]])), size = 1)
        second_node <- igraph::ego(Graph, order = dist, nodes = pick_second_node)
        
        new_jacc <- length(intersect(as.character(names(first_node[[1]])), as.character(names(second_node[[1]]))))/length(union(as.character(names(first_node[[1]])), as.character(names(second_node[[1]]))))
        
        jacc_permutated <-  c(jacc_permutated, new_jacc)
        setTxtProgressBar(pbin,nperm)
      }
      close(pbin)
      # end_time = Sys.time()
      # end_time - start_time
      
      pvalue <- 1-sum(jacc_original > jacc_permutated)/nIter
      PvalueMatL1000[i,k] = pvalue
    }
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  return(PvalueMatL1000)
  #saveRDS(PvalueMatL1000, file = "PvalueMatL1000.rds")
  
}

