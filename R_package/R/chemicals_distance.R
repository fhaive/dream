#' Calculate pairwise chemical distances between SMILES strings
#'
#' Given a list of named SMILES strings, this function calculates pairwise
#' chemical distances between all possible pairs
#'
#' @param smiles A named list of SMILES strings.
#' @param method A character string specifying the method to use for calculating chemical distances. 
#'   "levenshtein", "fingerprint" and "mcs" are supported. Default is "fingerprint".
#' @param only_heavy_atoms A logical indicating whether to use only heavy atoms when calculating the MCS. 
#'   This parameter is only used when method = "mcs". Default is TRUE.
#' @param radius An integer specifying the radius for the ECFP fingerprint. This parameter is only used when 
#'   method = "fingerprint". Default is 4.
#' @param n_bits An integer specifying the length of the bit vector for the ECFP fingerprint. This parameter is 
#'   only used when method = "fingerprint". Default is 2048.
#' @param api_endpoint A character string specifying the URL of the API endpoint to use for calculating 
#'   chemical distances. Default is "http://127.0.0.1:5000/chemicals_distance".
#' @param python A character string specifying the path to the Python executable to use. Default is "python".
#' 
#' @return A data.frame containing three columns: drug1, drug2, and distance. The drug1 and drug2 columns contain 
#'   the names of the two drugs being compared, and the distance column contains the calculated 
#'   chemical distance between the two drugs.
#' @export
chemicals_distance <- function(smiles, method="fingerprint",
                               only_heavy_atoms=TRUE, radius=4, n_bits=2048,
                               api_endpoint="http://127.0.0.1:5000/chemicals_distance",
                               python="python") {
  if(method=="levenshtein"){
    return(smiles_edit_distance(smiles))
  }
  parameters <- list("method"=method)
  if(method=="fingerprint"){
    parameters[['radius']] <- radius
    parameters[['n_bits']] <- n_bits
  } else if (method=="mcs") {
    parameters[['only_heavy_atoms']] <- only_heavy_atoms
  }
  obj = list("smiles"=smiles)

  p <- start_server(python)
  Sys.sleep(10)
  resp = httr::POST(api_endpoint, body=obj, query=parameters, encode="json")
  stop_server(p)
  if(httr::status_code(resp) != 200) {
    httr::warn_for_status(resp)
    return(NA)
  }
  content_ <- httr::content(resp)
  dists <- content_$distances
  dists_df <- data.frame(do.call(rbind, dists))
  colnames(dists_df) <- c("drug1", "drug2", "distance")
  return(dists_df)
}

smiles_edit_distance <- function(smiles){
  smiles_int<-lapply(smiles, utf8ToInt)
  SMILES_dist <- stringdist::seq_distmatrix(a=smiles_int, method = "lv", nthread=4)
  SMILES_dist <- as.matrix(SMILES_dist)
  rownames(SMILES_dist) <- colnames(SMILES_dist) <- names(smiles)
  return(distance_to_dataframe(SMILES_dist))
}