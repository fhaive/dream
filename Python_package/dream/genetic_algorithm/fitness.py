import numpy as np
import pandas as pd
import igraph as ig

from dream.genetic_algorithm.coverage_sum import coverage_sum

class EvaluationFunction(object):
    """
    A class for evaluating drug combinations based on their similarities in
    SMILES, MOA, and graph distances, as well as their coverage of known drug
    targets.

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
        keys: "gene", "rank".
    drug_targets : dict
        A list of dictionaries of drug targets with keys 
        "dat.drug.molecule_name" and "dat.target.gene_info.symbol".

    Attributes
    ----------
    smiles_distances : pandas.DataFrame
        A DataFrame of SMILES distances between drugs.
    moa_distances : pandas.DataFrame
        A DataFrame of MOA distances between drugs.
    paths_distances : pandas.DataFrame
        A DataFrame of graph distances between drugs.
    ppi_network : pandas.DataFrame
        A DataFrame representing the protein-protein interaction network.
    graph_rank : pandas.DataFrame
        A DataFrame representing the graph rank of the targets.
    L1000_drug_targets : pandas.DataFrame
        A DataFrame of drug targets.
    drug_names : list
        A list of drug names.

    Methods
    -------
    __call__(individual)
        Evaluates the drug combination represented by the given binary array.

    Returns
    -------
    tuple
        A tuple of mean SMILES distance, mean MOA distance, mean graph distance,
        coverage p-value, and number of drugs in the combination.
    """
    def __init__(
        self,
        smiles_distances,
        moa_distances,
        graph_distances,
        ppi_network,
        graph_rank,
        drug_targets,
    ):
        """
        Initialize the EvaluationFunction object.
        """
        tempdf = pd.DataFrame.from_dict(smiles_distances)
        idx = sorted(set(tempdf['drug1']).union(tempdf['drug2']))
        tempdf = tempdf.pivot(index='drug1', columns='drug2', values='distance').reindex(index=idx, columns=idx).fillna(0, downcast='infer').pipe(lambda x: x+x.values.T)
        
        self.smiles_distances = tempdf
        
        tempdf = pd.DataFrame.from_dict(moa_distances)
        idx = sorted(set(tempdf['drug1']).union(tempdf['drug2']))
        tempdf = tempdf.pivot(index='drug1', columns='drug2', values='distance').reindex(index=idx, columns=idx).fillna(0, downcast='infer').pipe(lambda x: x+x.values.T)
        
        self.moa_distances = tempdf
        
        tempdf = pd.DataFrame.from_dict(graph_distances)
        idx = sorted(set(tempdf['drug1']).union(tempdf['drug2']))
        tempdf = tempdf.pivot(index='drug1', columns='drug2', values='distance').reindex(index=idx, columns=idx).fillna(0, downcast='infer').pipe(lambda x: x+x.values.T)
        
        self.paths_distances = tempdf
        
        self.ppi_network = pd.DataFrame.from_dict(ppi_network)
        
        
        tempdf = pd.DataFrame.from_dict(graph_rank) #set _row column as index (index from raw data)
        self.graph_rank = tempdf.set_index("gene")
        
        self.L1000_drug_targets = pd.DataFrame.from_dict(drug_targets)
        self.drug_names = self.smiles_distances.index.to_list() 
       
        
    def __call__(self, individual):
        """
        Evaluate the drug combination represented by the given binary array.

        Parameters
        ----------
        individual : array_like
            A binary array representing a drug combination.

        Returns
        -------
        tuple
            A tuple of mean SMILES distance, mean MOA distance, mean graph distance,
            coverage p-value, and number of drugs in the combination.
        """
        if np.sum(individual) <= 1:
            return (0, 0, 0, 1, len(self.drug_names))
        candidate_drugs = [
            self.drug_names[i] for i, bit in enumerate(individual) if bit == 1
        ]
       
        smiles_submat = self.smiles_distances.loc[candidate_drugs, candidate_drugs]
        moa_submat = self.moa_distances.loc[candidate_drugs, candidate_drugs]
        paths_submat = self.paths_distances.loc[candidate_drugs, candidate_drugs]
        _, coverage_pval = coverage_sum(
            candidate_drugs,
            self.L1000_drug_targets,
            self.ppi_network,
            self.graph_rank
        )
        n_drugs = np.sum(individual)
        linear_index = np.triu_indices(n_drugs, k=1)
        return (
            np.mean(smiles_submat.values[linear_index]),
            np.mean(moa_submat.values[linear_index]),
            np.mean(paths_submat.values[linear_index]),
            coverage_pval,
            n_drugs
        )
