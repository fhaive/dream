import numpy as np
from deap import base
from deap import creator
from deap import tools
from deap import algorithms

from dream.genetic_algorithm import EvaluationFunction
from dream.genetic_algorithm import coverage_sum
from ..parmap import parmap
import functools

def genetic_algorithm(
    smiles_distances,
    moa_distances,
    graph_distances,
    ppi_network,
    graph_rank,
    drug_targets,
    population_size=100,
    n_offsprings=20,
    attribute_init_prob=0.3,
    attribute_mutation_prob=0.1,
    crossover_prob=0.7,
    mutation_prob=0.3,
    n_generations=2500,
    verbose=False,
    ):
    """
    API to a genetic algorithm to generate combinations of drugs that optimizes
    a multi-objective combination of distances of the drugs, including their
    SMILES representations, MOAs, and graph distances. The function also
    maximizes the coverage of disease targets, and minimizes the number of
    drugs used.

    Parameters
    ----------
    smiles_distances : dict
        A list of dictionaries of SMILES distances with keys: "drug1", "drug2",
        "distance".
    moa_distances : dict
        A list of dictionaries of MOA distances with keys: "drug1", "drug2",
        "distance".
    graph_distances : dict
        A list of dictionaries of graph distances with keys: "drug1", "drug2",
        "distance".
    ppi_network : dict
        A list of dictionaries representing the protein-protein interaction
        network with keys: "gene1", "gene2".
    graph_rank : dict
        A list of dictionaries representing the graph rank of the targets with
        keys: "gene", "rank", "_row".
    drug_targets : dict
        A list of dictionaries of drug targets with keys 
        "dat.drug.molecule_name" and "dat.target.gene_info.symbol".
    population_size : int, optional
        The size of the population of individuals, by default 100.
    n_offsprings : int, optional
        The number of offspring generated in each generation, by default 20.
    attribute_init_prob : float, optional
        The probability of adding a drug to the combination during the
        initialization of an individual, by default 0.3.
    attribute_mutation_prob : float, optional
        The probability of mutating each attribute of an individual, by default 0.1.
    crossover_prob : float, optional
        The probability of performing crossover between two individuals, by default 0.7.
    mutation_prob : float, optional
        The probability of performing mutation on an individual, by default 0.3.
    n_generations : int, optional
        The number of generations to run the algorithm for, by default 2500.
    verbose : bool, optional
        Whether to print verbose output during the algorithm, by default False.

    Returns
    -------
    drug_names : list
        A list of the names of the drugs in the dataset.
    population : list
        A list of the individuals in the final population.
    logbook : tools.Logbook
        A logbook containing statistics about the algorithm during each generation.
    solutions : list
        A list of tuples containing the non-dominated drug combinations and
        their fitness values encountered during the run.
    """
    # setup evaluation function
    evaluate = EvaluationFunction(
        smiles_distances=smiles_distances,
        moa_distances=moa_distances,
        graph_distances=graph_distances,
        ppi_network=ppi_network,
        graph_rank=graph_rank,
        drug_targets=drug_targets
    )

    toolbox, stats, hall_of_fame = setup_genetic_algorithm(
        evaluate,
        population_size,
        attribute_init_prob,
        attribute_mutation_prob
    )

    population = toolbox.population()
    population = toolbox.select(population, len(population))

    population, logbook = run_genetic_algorithm(
        population,
        toolbox,
        crossover_prob,
        mutation_prob,
        n_generations,
        n_offsprings,
        stats,
        hall_of_fame,
        verbose
    )
    solutions = parse_hall_of_fame(hall_of_fame, evaluate.drug_names)

    return evaluate.drug_names, population, logbook, solutions

def setup_genetic_algorithm(
    fitness_function,
    population_size,
    attribute_init_prob,
    attribute_mutation_prob
):
    """
    Set up the genetic algorithm for drug combination optimization using NSGA-II.

    Parameters
    ----------
    fitness_function : callable
        The fitness function that measures the quality of the drug combinations.
    population_size : int
        The number of individuals in each generation of the genetic algorithm.
    attribute_init_prob : float
        The probability of adding a drug to the combination during the
        initialization of an individual.
    attribute_mutation_prob : float
        The probability of mutating an attribute of an individual.

    Returns
    -------
    toolbox : deap.base.Toolbox
        The toolbox containing the operators for the genetic algorithm.
    stats : deap.tools.Statistics
        The statistics object to collect and report statistics on the population.
    hall_of_fame : deap.tools.ParetoFront
        The Pareto front object to store the best individuals seen during the optimization process.
    """
    ...

    creator.create(
        "Fitness", base.Fitness, weights=(1.0, 1.0, 1.0, -1.0, -1.0)
    )
    creator.create("Individual", list, fitness=creator.Fitness)

    IND_SIZE = len(fitness_function.drug_names)

    toolbox = base.Toolbox()
    toolbox.register("map", parmap, max_workers=16)
    toolbox.register(
        "init_attribute", lambda: np.random.random() < attribute_init_prob
    )
    toolbox.register(
        "individual",
        tools.initRepeat,
        creator.Individual,
        toolbox.init_attribute,
        n=IND_SIZE
    )
    toolbox.register(
        "population",
        tools.initRepeat,
        list,
        toolbox.individual,
        n=population_size
    )

    toolbox.register("evaluate", fitness_function)
    toolbox.register(
        "mutate", tools.mutShuffleIndexes, indpb=attribute_mutation_prob
    )
    toolbox.register("mate", tools.cxOnePoint)
    toolbox.register("select", tools.selNSGA2)

    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg_smile", lambda x: np.mean(x, axis=0)[0])
    stats.register("avg_moa", lambda x: np.mean(x, axis=0)[1])
    stats.register("avg_paths", lambda x: np.mean(x, axis=0)[2])
    stats.register("avg_coverage", lambda x: np.mean(x, axis=0)[3])
    stats.register("n_drugs", lambda x: np.mean(x, axis=0)[4])

    hall_of_fame = tools.ParetoFront()
    return toolbox, stats, hall_of_fame

def run_genetic_algorithm(
    population,
    toolbox,
    cxpb,
    mutpb,
    ngen,
    noffspring,
    stats=None,
    hall_of_fame=None,
    verbose=False
):
    """
    Runs a genetic algorithm to generate combinations of drugs that minimize a
    multi-objective combination of distances of the drugs; namely, the SMILES
    representations, the MOAs, and the graph distances of the drugs. It also
    maximizes the coverage of disease targets, and minimizes the number of
    drugs.

    Parameters
    ----------
    population : list
        A list of individuals comprising the initial population.
    toolbox : Toolbox
        A Toolbox object that contains the methods for performing various
        genetic operations (selection, crossover, mutation).
    cxpb : float
        The probability of mating two individuals.
    mutpb : float
        The probability of mutating an individual.
    ngen : int
        The number of generations to evolve.
    noffspring : int
        The number of offspring to produce at each generation.
    stats : deap.tools.Statistics, optional
        A Statistics object that accumulates statistics throughout the
        evolutionary process. Default is None.
    hall_of_fame : deap.tools.ParetoFront, optional
        A ParetoFront object that archives the best individuals generated during
        the evolutionary process. Default is None.
    verbose : bool, optional
        If True, information about the evolutionary process will be printed to
        the console. Default is False.

    Returns
    -------
    population : list
        A list of individuals comprising the final population.
    logbook : deap.tools.Logbook
        A Logbook object that contains the statistics of the evolutionary process.
    """
    ...
    logbook = tools.Logbook()
    logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])

    # Evaluate the individuals with an invalid fitness
    invalid_ind = [ind for ind in population if not ind.fitness.valid]
    evals = toolbox.map(toolbox.evaluate, invalid_ind)
    for ind, fit in zip(invalid_ind, evals):
        ind.fitness.values = fit

    if hall_of_fame is not None:
        hall_of_fame.update(population)

    record = stats.compile(population) if stats else {}
    logbook.record(gen=0, nevals=len(invalid_ind), **record)
    if verbose:
        print(logbook.stream)

    # Begin the generational process
    for gen in range(1, ngen + 1):
        # Generate new offspring
        offspring = algorithms.varOr(population, toolbox, noffspring, cxpb, mutpb)

        # evaluate the offspring
        evals = toolbox.map(toolbox.evaluate, offspring)
        for ind, fit in zip(offspring, evals):
            ind.fitness.values = fit

        # Select the next generation of fittest individuals
        population = toolbox.select(population+offspring, len(population))

        # Update the hall of fame with the generated individuals
        if hall_of_fame is not None:
            hall_of_fame.update(offspring)

        # Append the current generation statistics to the logbook
        record = stats.compile(population) if stats else {}
        logbook.record(gen=gen, nevals=len(invalid_ind), **record)
        if verbose:
            print(logbook.stream)

    return population, logbook

def parse_hall_of_fame(hall_of_fame, drug_names):
    """
    Extracts and returns the solutions from the Pareto front of a genetic
    algorithm optimization process.

    Parameters
    ----------
    hall_of_fame : deap.tools.ParetoFront
        A ParetoFront object that contains the individuals in the Pareto front.
    drug_names : list
        A list containing the names of the drugs being optimized.

    Returns
    -------
    list
        A list of tuples, each containing the names of the drugs and the
        fitness values of each solution.
    """
    solutions = []
    for individual in hall_of_fame:
        names = tuple(drug_names[i] for i, bit in enumerate(individual) if bit)
        fitness = individual.fitness
        solutions.append((names, parse_fitness(fitness.values)))
    return solutions

def parse_fitness(values):
    """
    Parse the fitness values of an individual into a dictionary with keys for
    different fitness components.

    Parameters
    ----------
    values : tuple or list
        Tuple or list of float values representing the fitness components of an
        individual.

    Returns
    -------
    dict
        Dictionary with keys for each of the five fitness components:
        "smiles_distance", "moa_distance", "shortest_path_distance",
        "coverage", and "n_drugs", and the corresponding values from the input
        `values` tuple or list.
    """
    keys = ("smiles_distance", "moa_distance", "shortest_path_distance", "coverage", "n_drugs")
    return dict(zip(keys, values))
