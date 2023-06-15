#' This function parses an edge rank matrix and performs edge selection based on different strategies. It returns a list containing a binary matrix representing selected edges and a vector of edge ranks.
#'
#' @importFrom TopKLists Borda
#'
#' @param edge_rank_matrix A symmetrix matrix with edge ranks as weight.
#' @param edge_selection_strategy Strategy to select edges for ensemble, "default" selects ranked edges untill all nodes have degree>=1, "top" swtiches to top N percentage ranked genes specfied by topN parameter.
#' @param mat_weights Type of scores in the ranked matrix; default:rank.
#' @param topN Top N percentage ranked edges to create ensemble if the edge_selection_strategy is "top"; default:10
#' @param debug_output Print help and status messages to help debug the running of the function; default:FALSE.
#' @param updateProgress Shiny application can request for update of progress from this function; default:NULL.
#' @return A list containing a vector of consensus edge ranks and a binary symmetrix matrix representing the edge rank matrix.
#' @export
parse_edge_rank_matrix <- function(edge_rank_matrix, edge_selection_strategy="default", mat_weights="rank", topN=10, debug_output=FALSE, updateProgress=NULL){
  bin_mat <- edge_rank_matrix
  #idx <- which(bin_mat>0)
  #if(length(idx)>0){
  #        print("Creating binary matrix...")
  #        bin_mat[which(bin_mat>0)] <- 1
  #}
  bin_mat[,] <- 0
  
  print("Getting edge list ordered by rank...")
  rank_matrix <- edge_rank_matrix
  rank_matrix[lower.tri(rank_matrix, diag=TRUE)] <- NA
  edge_df <- as.data.frame(as.table(rank_matrix))
  edge_df <- edge_df[-which(is.na(edge_df$Freq)),]
  edge_df <- data.frame(edge=paste0(edge_df$Var1,";",edge_df$Var2), weight=edge_df$Freq, stringsAsFactors=FALSE)
  edge_df <- edge_df[which(edge_df$weight>0),]
  ordOption <- ifelse(mat_weights=="rank", FALSE, TRUE)
  #print(paste0("ordOption : ", ordOption))
  #print(class(edge_df$weight))
  edge_rank <- edge_df$edge[order(edge_df$weight, decreasing=ordOption)]
  #print(paste0("edge rank before:", length(edge_rank)))
  #print(head(edge_rank))
  
  if(edge_selection_strategy=="default"){
    #input_genes <- nrow(rank_matrix)
    input_genes <- length(unique(unlist(strsplit(edge_rank, ";"))))
    #print(paste0("input genes : ", input_genes))
    genes <- NULL
    total_genes <- 0
    cutoffIdx <- NULL
    for(i in c(1:length(edge_rank))){
      if(total_genes<input_genes){
        local_genes <- strsplit(edge_rank[i], ";")[[1]]
        bin_mat[local_genes[1],local_genes[2]] <- 1
        bin_mat[local_genes[2],local_genes[1]] <- 1
        genes[local_genes[1]] <- 1
        genes[local_genes[2]] <- 1
        total_genes <- length(genes)
        cutoffIdx <- i
      }else{
        #print(paste0("Cutoff before break:", cutoffIdx))
        break
      }
    }
  }else if(edge_selection_strategy=="top"){
    cutoffIdx <- round((as.numeric(topN)*length(edge_rank))/100)
    for(i in c(1:length(edge_rank))){
      if(i<=cutoffIdx){
        local_genes <- strsplit(edge_rank[i], ";")[[1]]
        bin_mat[local_genes[1],local_genes[2]] <- 1
        bin_mat[local_genes[2],local_genes[1]] <- 1
      }else{
        break
      }
    }
  }
  #print(paste0("Cutoff:", cutoffIdx))
  edge_rank <- edge_rank[c(1:cutoffIdx)]
  #print(paste0("edge rank after:", length(edge_rank)))
  
  res_ll <- list(bin_mat=bin_mat, edge_rank=edge_rank)
  
  print("Returning!")
  return(res_ll)
}
