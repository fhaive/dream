#' Apply a genetic algorithm to optimize drug combinations.
#'
#' This function applies a genetic algorithm to find the best combinations of
#' drugs to treat a disease, based on multiple distance metrics.
#'
#' @param smiles a data frame with columns for drug 1, drug 2, and the distance between them based on their structure.
#' @param moa a data frame with columns for drug 1, drug 2, and the distance between them based on their mechanism of action (MOA)
#' @param graph a data frame with columns for drug 1, drug 2, and the distance between them based on the length of their shortest path on the disease-relevant gene network.
#' @param disease_network an interaction network of the disease-relevant genes.
#' @param rank a ranking of genes based on their degree in the disease-relevant gene network.
#' @param targets a data.frame in the format of OpenTargets specifying the target genes of each drug.
#' @param population_size an integer specifying the size of the population. Default is 5.
#' @param n_offsprings an integer specifying the number of offsprings produced in each generation. Default is 20.
#' @param attribute_init_prob a numeric value specifying the probability of choosing a drug during attribute initialization. Default is 0.3.
#' @param attribute_mutation_prob a numeric value specifying the probability of attribute mutation. Default is 0.1.
#' @param crossover_prob a numeric value specifying the probability of crossover. Default is 0.7.
#' @param mutation_prob a numeric value specifying the probability of mutation. Default is 0.3.
#' @param n_generations an integer specifying the number of generations to run. Default is 10.
#' @param api_endpoint the URL endpoint of the API. Default is "http://127.0.0.1:5000/genetic_algorithm".
#' @param python the path to the Python executable. Default is "python".
#'
#' @return a named list containing the following fields:
#'  drug_names: a list with the names of the drugs considered
#'  hall_of_fame: the list of best drug combinations identified during the execution of the algorithm, with their scores from the fitness function
#'  logbook: the average score of each fitness function for each generation
#'  population: a list of binary arrays representing the individuals in the population in the last generation of the algorithm.
#' @export
genetic_algorithm <- function(smiles, moa, graph, disease_network, rank, targets,
                              population_size=5,
                              n_offsprings=20,
                              attribute_init_prob=0.3,
                              attribute_mutation_prob=0.1,
                              crossover_prob=0.7,
                              mutation_prob=0.3,
                              n_generations=10,
                              api_endpoint="http://127.0.0.1:5000/genetic_algorithm",
                              python="python") {
  
  
  parameters <- list(
    population_size= population_size,
    n_offsprings=n_offsprings,
    attribute_init_prob=attribute_init_prob,
    attribute_mutation_prob=attribute_mutation_prob,
    crossover_prob=crossover_prob,
    mutation_prob=mutation_prob,
    n_generations=n_generations
  )
  
  obj = list(smiles_distances=smiles, moa_distances=moa, graph_distances=graph, ppi_network=disease_network, graph_rank=rank, drug_targets=targets)
  
  p <- start_server(python)
  resp = httr::POST(api_endpoint, body=obj, query=parameters, encode="json")
  stop_server(p)
  if(httr::status_code(resp) != 200) {
    httr::warn_for_status(resp)
    return(NA)
  }
  return(httr::content(resp))
}
