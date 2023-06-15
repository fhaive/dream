#' Infer coexpression networks excluding the drugs for which less than minSamples samples are available
#' 
#' @param path The path where the expression matrices are stored.
#' @param pattern An optional regular expression. Only file names which match the regular expression will be returned.
#' @param minSamples Only the matrices with minSamples number of samples (columns) will be considered.
#' @param ncores Number of cores to infer the coexpression networks.
#' @return A list of adjacency matrices.
#' @export

get_multiple_coexpression_networks <- function(path, pattern=NULL, minSamples=5, ncores = 12) {
  
  if (is.null(path)) {
    path <- "."
  }
  
  myfiles <- list.files(path = path, pattern = pattern)
  netlist <- rep(list(NULL), length(myfiles))
  inputinform <- data.frame()
  
  for (exprmatrix in 1:length(myfiles)) {
    
    #for (exprmatrix in 1:5) {  
    print(paste("Processing: ", myfiles[exprmatrix]))
    print(paste("Drug", exprmatrix, "of", length(myfiles)))
    
    inputinform  <- read.delim(paste0(path, myfiles[exprmatrix]), header = TRUE, quote = "")
    
    if (is.null(minSamples)) {
      minSamples <- length(colnames(inputinform))
    }
    
    if (length(colnames(inputinform))>=minSamples) {    # Filter out the matrices having less than minSamples samples
      
      #controls <- inputinform[, grep(pattern = "CTRL", x = colnames(inputinform))]
      
      generatematrices=get_ranked_consensus_matrix(gx_table = inputinform, iMethods = c("clr","aracne","mrnet"),
                                                   iEst = c("pearson","spearman","kendall","mi.empirical","mi.mm","mi.shrink","mi.sg"),
                                                   iDisc=c("none","equalfreq","equalwidth","globalequalwidth"), ncores = ncores, debug_output = TRUE, updateProgress = TRUE)
      
      
      rankMat.parsed=parse_edge_rank_matrix(edge_rank_matrix = generatematrices, edge_selection_strategy = "default",
                                            mat_weights = "rank", topN = 10, debug_output = TRUE, updateProgress = TRUE)
      
      #conGraph <- get_iGraph(rankMat.parsed$bin_mat)
      
      
      netlist[[exprmatrix]] <- rankMat.parsed$bin_mat
      
      
    }else{print(paste("Number of samples is too low for ", myfiles[exprmatrix], ": skipping!"))}
    
    
  }
  
  return(netlist)
}

