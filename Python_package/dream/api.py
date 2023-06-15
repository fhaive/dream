import io
import json
import pandas as pd
import numpy as np
from flask import Flask, jsonify, request
from dream.chemicals import pairwise_fingerprint_similarity
from dream.chemicals import pairwise_mcs_similarity


from dream.genetic_algorithm import genetic_algorithm

app = Flask(__name__)

@app.route("/chemicals_distance", methods=["GET", "POST"])
def chemicals_distance():
    parameters = request.args.to_dict(flat=True)
    method = parameters.pop("method", "fingerprint")
    data = request.get_json()
    smiles_dict = data["smiles"]
    rev_smiles_dict = {v:k for k,v in smiles_dict.items()}
    smiles = pd.Series(smiles_dict)

    try:
        if method == "mcs":
            only_heavy_atoms = (parameters.pop("only_heavy_atoms", "TRUE") == "TRUE")
            similarities = pairwise_mcs_similarity(
                smiles, only_heavy_atoms=only_heavy_atoms
            )
        elif method == "fingerprint":
            radius = int(parameters.pop("radius", 4))
            n_bits = int(parameters.pop("n_bits", 2048))
            similarities = pairwise_fingerprint_similarity(
                smiles, radius=radius, n_bits=n_bits
            )
        else:
            # http 400 error code: bad request: https://http.cat/400
            return jsonify({"error": "Unrecognized method: " + method}), 400
    except Exception as err:
        return jsonify({"error": str(err)}), 400
    # negative log transformation converts the similarity to a distance
    distances = [(rev_smiles_dict[mol1], rev_smiles_dict[mol2], -np.log(max(1e-5, sim))) for (mol1, mol2, sim) in similarities]
    return jsonify({"distances": distances})

@app.route("/genetic_algorithm", methods=["GET", "POST"])
def genetic_algorithm_api():
    parameters = request.args.to_dict(flat=True)

    population_size = int(parameters.pop("population_size", 100))
    n_offsprings = int(parameters.pop("n_offsprings", 20))
    attribute_init_prob = float(parameters.pop("attribute_init_prob", 0.3))
    attribute_mutation_prob = float(parameters.pop("attribute_mutation_prob", 0.1))
    crossover_prob = float(parameters.pop("crossover_prob", 0.7))
    mutation_prob = float(parameters.pop("mutation_prob", 0.3))
    n_generations = int(parameters.pop("n_generations", 2500))

    data = request.get_json()
    smiles_distances = data.pop("smiles_distances", None)
    moa_distances = data.pop("moa_distances", None)
    graph_distances = data.pop("graph_distances", None)
    ppi_network = data.pop("ppi_network", None)
    graph_rank = data.pop("graph_rank", None)
    drug_targets = data.pop("drug_targets", None)

    try:
        check_not_none(
            smiles_distances=smiles_distances,
            moa_distances=moa_distances,
            graph_distances=graph_distances,
            ppi_network=ppi_network,
            graph_rank=graph_rank,
            drug_targets=drug_targets
        )
    except ValueError as e:
        return jsonify(
            {"error": "{} is missing or not well formatted".format(str(e))}
        )

    drug_names, population, logbook, hall_of_fame = genetic_algorithm(
        smiles_distances=smiles_distances,
        moa_distances=moa_distances,
        graph_distances=graph_distances,
        ppi_network=ppi_network,
        graph_rank=graph_rank,
        drug_targets=drug_targets,
        population_size=population_size,
        n_offsprings=n_offsprings,
        attribute_init_prob=attribute_init_prob,
        attribute_mutation_prob=attribute_mutation_prob,
        crossover_prob=crossover_prob,
        mutation_prob=mutation_prob,
        n_generations=n_generations,
        verbose=False,
    )

    return jsonify({
        "drug_names": drug_names,
        "population": population,
        "logbook": logbook,
        "hall_of_fame":hall_of_fame,
    })

def check_not_none(**kwargs):
    for (k, v) in kwargs.items():
        if v is None:
            raise ValueError(k)

if __name__ == '__main__':
    app.debug = True
    app.run()
