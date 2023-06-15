#' Infer correlation matrix from gene expression table using mutual information.
#'
#' This function, calculate_correlation_matrix, utilizes the MINET package to create a correlation matrix based on the mutual information method.
#' It allows the user to specify the inference algorithms, correlation calculation methods, and discretization methods, enabling multiple combinations of parameters for multiple runs of minet().
#' The resulting correlation matrices are combined using the median to create a consensus matrix.
#'
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom minet minet
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach getDoParWorkers getDoParName %dopar%
#' @importFrom plyr aaply laply
#' @importFrom stats quantile median
#' @importFrom utils capture.output
#'
#' @param gx_table A data frame representing the gene expression table.
#' @param iMethods A vector of valid inference algorithms for the MINET package.
#' @param iEst A vector of valid correlation methods for the MINET package.
#' @param iDisc A vector of valid discretization methods for the MINET package.
#' @param ncores The number of cores to use for running instances of MINET in parallel. Default: 2.
#' @param debug_output If TRUE, it prints help and status messages to aid in debugging the function. Default: FALSE.
#' @return A binary symmetric matrix representing the median of mutual information correlations computed across various MINET combinations
#' @export
calculate_correlation_matrix <- function(gx_table, iMethods, iEst, iDisc, ncores=2, debug_output=FALSE){
  parList <- list()
  out.GX.list<- list()
  out.tmp.list<- list()
  out.list<- list()
  tGX <- t(gx_table)
  stdGX <- scale(tGX)
  
  if(is.null(iMethods)){
    stop("Please Select atLeast One Inference Algorithm!")
  }
  
  if(is.null(iEst)){
    stop("Please Select At Least One Correlation!")
  }
  
  if(is.null(iDisc)){
    stop("Please Select At Least One Discretization Method!")
  }
  
  cntr <- 0
  for(i in 1:length(iMethods)){
    #print(methods2[i])
    for(j in 1:length(iEst)){
      #print(paste("-",est.opt2[j]))
      for(k in 1:length(iDisc)){
        cntr <- cntr + 1
        parList[[cntr]] <- list()
        parList[[cntr]][["mt"]] <- iMethods[i]
        parList[[cntr]][["est"]] <- iEst[j]
        parList[[cntr]][["disc"]] <- iDisc[k]
      }
    }
  }
  #print(paste0("List of MINET Combinations: ", parList))
  
  if(ncores > parallel::detectCores()){
    ncores <- parallel::detectCores()
  }
  
  cl <- parallel::makeCluster(ncores, outfile="")
  doParallel::registerDoParallel(cl)
  #print(paste("Total Cores: ", parallel::detectCores(), sep=""))
  #print(paste("Total Workers: ", foreach::getDoParWorkers(), sep=""))
  #print(paste("DoPar Name: ",  foreach::getDoParName(), sep=""))
  
  #print(paste0("Before For Each, ", "Number of Combinations:", length(parList)))
  utils::capture.output(print(paste0("Starting Parallel Cluster For ", length(parList), " Combinations ", timestamp())), file="minet-log.txt", append=TRUE)
  utils::capture.output(print(paste0("Starting Parallel Cluster For ", length(parList), " Combinations ", timestamp())), file="minet-to-run.txt", append=TRUE)
  utils::capture.output(print(paste0("Starting Parallel Cluster For ", length(parList), " Combinations ", timestamp())), file="minet-completed.txt", append=TRUE)
  #out.tmp.list <- foreach(i=1:length(parList)) %do% {
  out.tmp.list <- foreach::foreach(i=1:length(parList)) %dopar% {
    #capture.output(print("Start For Each"), file="minet-log.txt", append=T)
    mt <- parList[[i]]$mt
    est <- parList[[i]]$est
    disc <- parList[[i]]$disc
    #print(parList[[i]]$mt)
    
    if((est == "mi.empirical" | est == "mi.mm" | est == "mi.shrink" | est == "mi.sg") & disc == "none") {
      miMat <- -1
      miMatName <- "np"
    }
    else{
      if(debug_output==TRUE)
        utils::capture.output(print(paste("----",mt,est,disc, sep="__")), file="minet-log.txt", append=TRUE)
      
      utils::capture.output(print(paste0("Iteration-", i, ": ", mt, "-", est, "-", disc)), file="minet-to-run.txt", append=TRUE)
      
      ptm <- proc.time()
      miMat <- minet::minet(stdGX, method=mt, estimator=est, disc=disc)
      utils::capture.output(print(paste0("Iteration-", i, ", ", mt, "-", est, "-", disc, ": ", "MINET Execution Time - ", round(proc.time() - ptm)[3], " sec")), file="minet-completed.txt", append=TRUE)
      #capture.output(print(proc.time() - ptm), file="minet-log.txt", append=T)
      miMatName <- paste(mt,est,disc,sep="__")
    }
    out.list <- list("mat"=miMat, "name"=miMatName)
  }
  utils::capture.output(print("For Each Finished, Stopping Cluster..."), file="minet-log.txt", append=TRUE)
  parallel::stopCluster(cl)
  
  for(i in 1:length(out.tmp.list)){
    out.GX.list[[i]] <- out.tmp.list[[i]]$mat
    names(out.GX.list)[[i]] <- out.tmp.list[[i]]$name
  }
  
  #if(debug_output==TRUE)
  #print(paste("Got minet list of lists --- ", str(out.GX.list)))
  
  #removing ones having no values
  #llGX <- out.GX.list[-(which(names(out.GX.list)=="np"))]
  idx <- which(names(out.GX.list)=="np")
  if(length(idx)>0){
    llGX <- out.GX.list[-idx]
  }else{
    llGX <- out.GX.list
  }
  
  if(length(llGX)>1){
    #print("Making Median Matrix...")
    arr1<-abind::abind(llGX,along=3)
    matGX <- apply(arr1,c(1,2),median)
  }else{
    #print("Only one matrix, not computing median!")
    matGX <- llGX[[1]]
  }
  
  #print("Return Matrix")
  return(matGX)
}
