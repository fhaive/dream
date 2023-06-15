#' Fetch SMILES from PubChem by compound name
#'
#' This function takes a vector of compound names and retrieves their SMILES
#' strings from PubChem using the PubChem REST API.
#'
#' @param compounds A character vector of compound names for which SMILES strings should be fetched.
#' @return A list of named SMILES strings corresponding to the input compounds.
#' @export
fetch_smiles <-function(compounds){
    my_smiles <- sapply(compounds, get_smiles_from_name, simplify=FALSE)
    return(my_smiles)
}

#' Retrieves the CanonicalSmiles property for a compound from its name using the PubChem API.
#'
#' This is an internal function that retrieves the CanonicalSmiles property for a given compound name using the PubChem API. It performs an HTTP GET request to the PubChem website and extracts the CanonicalSmiles property from the response.
#'
#' @param compound The name of the compound for which the CanonicalSmiles property should be retrieved.
#' @return The CanonicalSmiles property for the compound, or \code{NA} if the retrieval fails.
#'
#' @import httr
#'
get_smiles_from_name <- function(compound){
    url <- "https://pubchem.ncbi.nlm.nih.gov"
    endpoint <- paste0("/rest/pug/compound/name/", compound, "/property/CanonicalSmiles/txt")
    response <- httr::GET(url = url, path = endpoint)
    if (httr::status_code(response) != 200) {
        print(httr::status_code(response))
        return(NA)
    }
    smiles <- httr::content(response, encoding="UTF-8")
    return(gsub("\\n", "", smiles))
}