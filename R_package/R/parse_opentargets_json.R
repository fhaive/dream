#' Convenience function to parse OpenTargets JSON files
#' 
#' The parse_opentargets_json function is used to parse an OpenTargets JSON file and extract specific fields to create a final parsed data frame. The function takes the json_file as input and returns the parsed data frame.
#' 
#' @param json_file The OpenTargets JSON file.
#' @return A dataframe representing the OpenTargets database.
#' @export


parse_opentargets_json <- function(json_file) {
  
  if(is.null(json_file)) {
    stop("Error: please type in the OpenTargets JSON file to be processed!")
  }
  
  #json_file <- "19.02_evidence_data.json"
  out <- lapply(readLines(json_file), fromJSON) # load opentargets json file
  dat <- rbindlist(lapply(out, function(x) {as.list(unlist(x))}), fill=TRUE)
  final_opentargets_parsed <- data.frame(dat$target.gene_info.symbol, dat$target.gene_info.name, dat$target.gene_info.geneid, dat$disease.efo_info.label, dat$disease.name, dat$drug.molecule_name, dat$drug.id, dat$drug.molecule_type, dat$drug.max_phase_for_all_diseases.label, dat$unique_association_fields.pathway_id)
  
  return(final_opentargets_parsed)
}

