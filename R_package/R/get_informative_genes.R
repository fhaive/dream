#' Identify informative genes for co-expression network construction.
#' 
#' The `get_informative_genes` function identifies informative genes, such as differentially expressed genes, that will be used to build a co-expressed network.
#' 
#' @param expr_mat An expression matrix where rows are genes and columns are samples.
#' @param group A factor indicating the biological conditions of the samples. For example, `group = as.factor(c(rep("T", length(tumor)), rep("N", length(normal))))` specifies two groups, "T" and "N", where "T" represents tumor samples and "N" represents normal samples.
#' @param test The statistical test to be performed for the identification of the informative genes. If test="wilcoxon" a Wilcoxon test is performed. If test="var.test" an F-test on the variance is performed. The default is test="wilcoxon + var.test" and both Wilcoxon and var.test are performed.
#' @param percentile An integer indicating the percentile of the ranked distribution of informative genes to be returned. The default is percentile=10.
#' @return An expression matrix for the informative genes to be used as input of the get_coexpression_network function.
#' @export

get_informative_genes <- function(expr_mat, group, test="wilcoxon + var.test", percentile=10) {
  
  if(is.null(expr_mat)){
    
    stop("Error: please provide an expression matrix!")
    
  }
  
  if(is.null(group)){
    stop("Error: please provide groups!")
    
  }
  
  if(length(group)!=ncol(expr_mat)){
    stop("Error: the group length is different from the sample size!")
  }
  
  if(!isClass(percentile, Class = "numeric")){
    stop("Error:please provide an integer")
  }
  
  if(percentile < 1  | percentile >100) {
    
    stop("Error: please provide a value between 1 and 100!")
    
  }
  
  
  if(!test %in% c("wilcoxon","var.test","wilcoxon + var.test")){
    stop("Error: please provide one of the following tests: wilcoxon, var.test, wilcoxon + var.test")
  }
  
  
  if (test == "wilcoxon") {
    
    wilres <- apply(expr_mat, 1, function(x){wilcox.test(x ~ group, data = expr_mat)})
    pvalwil <- sapply(wilres, function(x){x$p.value})
    pvalwil[pvalwil==0] <- 0.000000000000001
    
    logFC <- log(rowMeans(expr_mat[, which(group==names(table(group)[2]))])/rowMeans(expr_mat[, which(group==names(table(group)[1]))])) ### Change this according to factors levels order
    
    wilcox_score <- abs(logFC*(-log(pvalwil)))
    wilcox_score <- stats::setNames(object = wilcox_score, nm = rownames(expr_mat))
    
    uqua <- length(expr_mat[,1])*percentile/100
    
    inputinform <- list(expr_mat[names(sort(wilcox_score, decreasing=TRUE))[1:uqua],], wilcox_score)
    names(inputinform) <- c("inform_mat", "wilcox_score")
    
  }else if(test == "var.test") {
    
    varres <- apply(expr_mat, 1, function(x){var.test(x ~ group, data = expr_mat)})
    pvalvar <- sapply(varres, function(x){x$p.value})
    pvalvar[pvalvar == 0] <- 0.000000000000001
    
    logFC <- log(rowMeans(expr_mat[, which(group==names(table(group)[2]))])/rowMeans(expr_mat[, which(group==names(table(group)[1]))])) ### Change this according to factors levels order
    
    var_score <- abs(logFC*(-log(pvalvar)))
    
    var_score <- stats::setNames(object = var_score, nm = rownames(expr_mat))
    
    uqua <- length(expr_mat[,1])*percentile/100
    
    
    inputinform <- list(expr_mat[names(sort(var_score, decreasing=TRUE))[1:uqua],], var_score)
    names(inputinform) <- c("inform_mat", "var_score")
    
    
  }else if(test == "wilcoxon + var.test"){
    
    wilres <- apply(expr_mat, 1, function(x){wilcox.test(x ~ group, data = expr_mat)})
    pvalwil <- sapply(wilres, function(x){x$p.value})
    pvalwil[pvalwil==0] <- 0.000000000000001
    
    
    varres <- apply(expr_mat, 1, function(x){var.test(x ~ group, data = expr_mat)})
    pvalvar <- sapply(varres, function(x){x$p.value})
    pvalvar[pvalvar == 0] <- 0.000000000000001
    
    
    logFC <- log(rowMeans(expr_mat[, which(group==names(table(group)[2]))])/rowMeans(expr_mat[, which(group==names(table(group)[1]))])) ### Change this according to factors levels order
    
    wilcox_score <- abs(logFC*(-log(pvalwil)))
    
    wilcox_score <- stats::setNames(object = wilcox_score, nm = rownames(expr_mat))
    
    var_score <- abs(logFC*(-log(pvalvar)))
    var_score <- stats::setNames(object = var_score, nm = rownames(expr_mat))
    
    uqua <- length(expr_mat[,1])*percentile/100
    
    #all_samples_filtered_wilcoxon <- expr_mat[names(sort(wilcox_score, decreasing=TRUE))[1:uqua],]
    #all_samples_filtered_var <- expr_mat[names(sort(var_score, decreasing=TRUE))[1:uqua],]
    
    inputinform <- list(expr_mat[union(names(sort(wilcox_score, decreasing=TRUE))[1:uqua], names(sort(var_score, decreasing=TRUE))[1:uqua]),], var_score, wilcox_score)
    names(inputinform) <- c("inform_mat", "wilcox_score", "var_score")
    
  }
  return(inputinform)
}
