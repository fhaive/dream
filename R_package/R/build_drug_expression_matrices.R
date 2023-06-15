#' Converts expression matrices from L1000 data in tab delimited format.
#' 
#' This function converts expression matrices from the LINCS1000 dataset. 
#' It provides flexibility in including control samples and allows selection of specific genes for the resulting expression matrix.
#' 
#' @param L1000_data An object of class gct or gctx from the LINCS1000 study (e.g. parsed by the parse.gct function of the cmapR package).
#' @param L1000_metadata Object containing the metadata relative to the L1000 dataset (e.g file GSE70138_Broad_LINCS_sig_info_2017-03-06.txt).
#' @param include_controls A logical value indicating whether to include the plate-specific control samples in the expression matrix.
#' @param L1000_data_CTRL An object of class gct or gctx from the LINCS1000 study for control samples (e.g. DMSO).
#' @param L1000_metadata_CTRL Object containing the metadata relative to the L1000 dataset for control samples.
#' @param dest_path Folder to be created as a destination folder for the expression matrices. Default to current directory.
#' @param inputgenes A character vector of gene symbols to be included in the expression matrix (e.g. L1000 landmark genes). Default all the genes will be considered.
#' @return A squared matrix reporting the statistical significance of the similarity measures (e.g. Jaccard indices) among couples of drug targets neighborhoods.
#' @export

build_drug_expression_matrices <- function(L1000_data, L1000_metadata, include_controls = FALSE, L1000_data_CTRL=NULL, L1000_metadata_CTRL=NULL, dest_path = ".", inputgenes) {
  
  if(is.null(L1000_data)) {
    stop("Error: please provide a L1000 expression matrix!")
  }
  
  if(is.null(L1000_metadata)) {
    stop("Error: please provide a L1000 metadata file!")
  }
  
  if(!class(include_controls)=="logical") {
    stop("Error: please provide a logical value to the include_controls argument!")
  }
  
  if(include_controls==TRUE & is.null(L1000_data_CTRL)) {
    stop("Error: please provide a L1000 expression matrix for control samples!")
  }
  
  if(include_controls==TRUE & is.null(L1000_metadata_CTRL)) {
    stop("Error: please provide a L1000 metadata file for control samples!")
  }
  
  if (!class(dest_path)=="character") {
    stop("Error: please provide a valid path!")
  }
  
  if (is.null(inputgenes)) {
    inputgenes <- rownames(L1000_data@mat)
  }
  
  # if (!class(inputgenes)=="character") {
  #   stop("Error: please provide a vector of genes symbols!")
  # }
  
  if (!file.exists(dest_path)){
    dir.create(file.path(dest_path))
  } 
  
  
  distil_ids <- c()
  distil_ids_CTRL <- c()
  plate <- c()
  CTRL_expr <- data.frame()
  finalmat <- data.frame()
  
  for (drug in 1:length(unique(L1000_metadata$pert_iname))) {
    distil_ids <- L1000_metadata[which(L1000_metadata$pert_iname==L1000_metadata$pert_iname[drug]), "distil_id"]
    if(length(distil_ids)>1){
      drugs_expr <- L1000_data@mat[which(rownames(L1000_data@mat) %in% inputgenes), distil_ids] # filtra landmark genes
      rownames(drugs_expr) <- inputgenes[which(inputgenes %in% rownames(L1000_data@mat))]
    }
    
    if(include_controls == TRUE) {
      plate <- sapply(strsplit(distil_ids, ":", fixed=TRUE), "[", 1)
      CTRL_expr_all_plates <- data.frame(row.names = 1:length(rownames(drugs_expr)))
      for (distil_id in 1:length(plate)) {
        distil_ids_CTRL <- L1000_metadata_CTRL[grep(plate[distil_id], L1000_metadata_CTRL$distil_id, fixed = TRUE), "distil_id"]
        CTRL_expr <- L1000_data_CTRL@mat[which(rownames(L1000_data_CTRL@mat) %in% inputgenes), distil_ids_CTRL] # filtra landmark genes
        
        CTRL_expr_all_plates <- cbind(CTRL_expr_all_plates, CTRL_expr)
      }
      
      rownames(CTRL_expr_all_plates) <- inputgenes[which(inputgenes %in% rownames(L1000_data_CTRL))]
      colnames(CTRL_expr_all_plates) <- sapply(colnames(CTRL_expr_all_plates), function(x) paste0("CTRL_", x))
      
      finalmat <- cbind(drugs_expr, CTRL_expr_all_plates)
      write.table(finalmat, file = paste0(dest_path, "/", "Expression_matrix_", unique(L1000_metadata$pert_iname)[drug], ".txt"), quote = FALSE, sep = "\t", row.names=TRUE, col=NA)
      
    }else{
      write.table(drugs_expr, file = paste0(dest_path, "/", "Expression_matrix_", unique(L1000_metadata$pert_iname)[drug], ".txt"), quote = FALSE, sep = "\t", row.names=TRUE, col=NA)
    }
    
  }
  
}


