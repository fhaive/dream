#' Create edge ranked correlation matrix by combining information from different inference algorithms with the help of internal function calculate_correlation_matrix().
#'
#' get_ranked_consensus_matrix uses the internal function calculate_correlation_matrix() to get a single consensus matrix per inference algorithm. User can specify
#' the inference algorithms, correlation calculation methods and discretization methods, a combination of parameters will be created per inference algorithm to run calculate_correlation_matrix(), 
#' this will generate a consensus matrix per inference algorithm. The consesus matrices from different inference algorithms are used to create a single binary matrix by rank based selection of edges.
#'
#' @importFrom TopKLists Borda
#'
#' @param gx_table Gene expression table as a data frame.
#' @param iMethods Vector of valid inference algorithms for MINET package.
#' @param iEst Vector of valid correlation methods for MINET package.
#' @param iDisc Vector of valid discretization methods for MINET package.
#' @param ncores Number of cores for running instances of MINET in parallel default:2.
#' @param edge_selection_strategy How to select top ranked edges default:default. By default selectis top edges until all the nodes have at least one edge.
#' @param topN Top n percentage edges to select if edge_selection_strategy is 'top' default:10.
#' @param debug_output Print help and status messages to help debug the running of the function default:FALSE.
#' @param updateProgress Shiny application can request for update of progress from this function default:NULL.
#' @return A symmetrix matrix with median edge ranks representing the edge rank based consensus from different inference algorithms.
#' @export
get_ranked_consensus_matrix <- function(gx_table=NULL, iMethods=NULL, iEst=NULL, iDisc=NULL, ncores=2, matList=NULL, mat_weights="rank", ensemble_strategy="minet", debug_output=FALSE, updateProgress=NULL){
  mat_ll <- list()
  ranked_edges_ll <- list()
  
  if(length(grep("minet" ,ensemble_strategy))>0){
    mthdCount <- 1
    totalMthds <- length(iMethods)+1
    for(mthd in iMethods){
      if (is.function(updateProgress)) {
        text <- paste0("'", mthd, "' Median")
        value <- mthdCount / totalMthds
        updateProgress(detail = text, value = value)
      }
      mthdCount <- mthdCount + 1
      
      #print(paste0("Calculate correlation matrix for method : ", mthd))
      mat_ll[[mthd]] <- calculate_correlation_matrix(gx_table=gx_table, iMethods=mthd, iEst=iEst, iDisc=iDisc, ncores=ncores)
      
      #print(paste0("Get ranked edges for method : ", mthd))
      mat_ll[[mthd]][lower.tri(mat_ll[[mthd]], diag=TRUE)] <- NA
      
      edge_df <- as.data.frame(as.table(mat_ll[[mthd]]))
      #print("Minet edge_df before:")
      #print(str(edge_df))
      edge_df <- edge_df[-which(is.na(edge_df$Freq)),]
      edge_df <- data.frame(edge=paste0(edge_df$Var1,";",edge_df$Var2), weight=edge_df$Freq, stringsAsFactors=FALSE)
      #print("Minet edge_df after:")
      #print(str(edge_df))
      
      ranked_edges_ll[[mthd]] <- edge_df[order(edge_df$weight, decreasing=TRUE), "edge"]
    }
  }
  
  if(length(grep("user", ensemble_strategy))>0){
    matList <- as.list(matList)
    for(i in c(1:length(matList))){
      itmName <- paste0("user_", i)
      matList[[i]][lower.tri(matList[[i]], diag=TRUE)] <- NA
      
      edge_df <- as.data.frame(as.table(matList[[i]]))
      #print("User edge_df before:")
      #print(str(edge_df))
      edge_df <- edge_df[-which(is.na(edge_df$Freq)),]
      edge_df <- data.frame(edge=paste0(edge_df$Var1,";",edge_df$Var2), weight=edge_df$Freq, stringsAsFactors=FALSE)
      hasNulls <- FALSE
      if(mat_weights=="rank"){
        nullIdx <- which(edge_df$weight==0)
        if(length(nullIdx)>0){
          edge_df_null <- edge_df[nullIdx,]
          edge_df <- edge_df[-nullIdx,]
          hasNulls <- TRUE
        }
      }
      #print("User edge_df after:")
      #print(str(edge_df))
      
      ordOption <- ifelse(mat_weights=="rank", FALSE, TRUE)
      #print(paste0("ordOption : ", ordOption))
      edge_df <- edge_df[order(edge_df$weight,decreasing=ordOption),]
      if(hasNulls){
        edge_df <- as.data.frame(rbind(edge_df, edge_df_null), stringsAsFactors=FALSE)
      }
      #ranked_edges_ll[[itmName]] <- edge_df[order(edge_df$weight, decreasing=ordOption), "edge"]
      ranked_edges_ll[[itmName]] <- edge_df$edge
    }
  }
  #print("ranked_edges_ll:")
  #print(length(ranked_edges_ll))
  #print(names(ranked_edges_ll))
  #print(lapply(ranked_edges_ll, length))
  
  if (is.function(updateProgress)) {
    updateProgress(detail = "Consensus Binary", value = 1)
  }
  
  print("Perform Borda on list of list of ranked edges.")
  borda_res <- TopKLists::Borda(ranked_edges_ll)
  #Borda.plot(borda_res) #Is this the cause!!!?????
  
  print("Get a consensus binary matrix by selecting the most significant ranked edges from median rank of Borda result.")
  if(length(mat_ll)>0){
    rank_mat <- mat_ll[[1]]
  }else{
    rank_mat <- matList[[1]]
  }
  rank_mat[,] <- 0
  #print("rank_mat")
  #print(dim(rank_mat))
  
  median_list <- borda_res$TopK$median
  #Kendall.plot(ranked_edges_ll, median_list)
  
  #if(edge_selection_strategy=="default"){
  #	#input_genes <- dim(gx_table)[1]
  #	input_genes <- nrow(rank_mat)
  #	genes <- NULL
  #	total_genes <- 0
  #	cutoffIdx <- NULL
  #	for(i in c(1:length(median_list))){
  #		if(total_genes<input_genes){
  #			local_genes <- strsplit(median_list[i], ";")[[1]]
  #			rank_mat[local_genes[1],local_genes[2]] <- i
  #			rank_mat[local_genes[2],local_genes[1]] <- i
  #			genes[local_genes[1]] <- 1
  #			genes[local_genes[2]] <- 1
  #			total_genes <- length(genes)
  #		}else{
  #			cutoffIdx <- i-1
  #			break
  #		}
  #	}
  #}else if(edge_selection_strategy=="top"){
  #	cutOff <- round((as.numeric(topN)*length(median_list))/100)
  #	for(i in c(1:length(median_list))){
  #		if(i<=cutOff){
  #			local_genes <- strsplit(median_list[i], ";")[[1]]
  #			rank_mat[local_genes[1],local_genes[2]] <- i
  #			rank_mat[local_genes[2],local_genes[1]] <- i
  #		}else{
  #			break
  #		}
  #	}
  #}
  
  #print("median_list:")
  #print(str(median_list))
  #print(length(median_list))
  #print(head(median_list))
  for(i in c(1:length(median_list))){
    local_genes <- strsplit(median_list[i], ";")[[1]]
    rank_mat[local_genes[1],local_genes[2]] <- i
    rank_mat[local_genes[2],local_genes[1]] <- i
  }
  
  #print("Rank matrix computed, returning!")
  return(rank_mat)
}