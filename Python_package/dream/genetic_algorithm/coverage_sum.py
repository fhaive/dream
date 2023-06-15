import pandas as pd
import numpy as np
import itertools as it
from scipy.stats import norm
import igraph as ig 

def coverage_sum(candidate_drugs, drug_targets, ppi_network, graph_rank):
    """
    Compute coverage score for a set of candidate drugs based on the proportion
    of nodes in the protein-protein interaction network that are covered by the
    network neighborhood of the drugs' targets. Uses a permutation test to
    compute the significance of the coverage score.

    Parameters
    ----------
    candidate_drugs : list of str
        List of drug names to consider.
    drug_targets : pandas DataFrame
        DataFrame containing drug-target associations with columns
        "dat.drug.molecule_name" and "dat.target.gene_info.symbol".
    ppi_network : igraph.Graph
        The PPI network as an igraph.Graph object.
    graph_rank : pandas DataFrame
        DataFrame containing the graph rank of the nodes in the ppi network,
        with one column called "rank" and the row index corresponding to node
        names.

    Returns
    -------
    tuple
        A tuple containing the z-score and p-value for the coverage of
        candidate_drugs in the PPI network.

    """
    # take the targets of the drugs considered
    drug_target_df = drug_targets[
        drug_targets["dat.drug.molecule_name"].isin(candidate_drugs)
    ]
    # make a df with the nodes of the ppi network sorted by their degree
    #ppi needs to be a igraph object
    ppi_network = ig.Graph.TupleList(ppi_network.itertuples(index=False), directed=False, weights=False)
    deg_dist = pd.DataFrame(
        data=ppi_network.degree(),
        index=ppi_network.vs["name"],
        columns=["dist"]
    )
    deg_dist.sort_values(by="dist", inplace=True)
    deg_dist["bin"] = pd.qcut(deg_dist.dist, q=20, labels=range(20))
    return get_drug_input_coverage(
        candidate_drugs=candidate_drugs, 
        ppi_network=ppi_network, 
        drug_target_df=drug_target_df, 
        graph_rank=graph_rank, 
        deg_dist=deg_dist
    )

def get_drug_input_coverage(
    candidate_drugs,
    ppi_network,
    drug_target_df,
    graph_rank,
    deg_dist,
    order=1,
):
    """
    Permutation test to evaluate the significance of the coverage score for a
    set of candidate drugs based on the proportion of nodes in the PPI network
    that are covered by the network neighborhood of the drugs' targets. 

    Parameters
    ----------
    candidate_drugs : list of str
        List of drug names to consider.
    ppi_network : igraph.Graph
        The PPI network as an igraph.Graph object.
    drug_target_df : pandas DataFrame
        DataFrame containing drug-target associations with columns
        "dat.drug.molecule_name" and "dat.target.gene_info.symbol".
    graph_rank : pandas DataFrame
        DataFrame containing the graph rank of the nodes in the ppi network,
        with one column called "rank" and the row index corresponding to node
        names.
    deg_dist : pandas DataFrame
        DataFrame containing the degree distribution of nodes in the ppi
        network, with one column called "dist" and the row index corresponding
        to node names.
    order : int, optional
        The order of the neighborhood to consider. Default is 1.

    Returns
    -------
    tuple
        A tuple containing the standardized coverage score of candidate_drugs
        in the PPI network and the corresponding p-value computed with a
        permutation test.

    """
    # take the targets of selected drugs
    targets = set(
        drug_target_df
        .loc[
            drug_target_df["dat.drug.molecule_name"] .isin(candidate_drugs), 
            "dat.target.gene_info.symbol"
        ]
    )
    # consider only targets available in ppi network
    targets = list(set(ppi_network.vs["name"]) & targets)
    observed = compute_cov_score(targets, ppi_network, graph_rank, order=1)
    simulated = []
    for i in range(100):
        random_targets = generate_random_targets(targets, deg_dist)
        simulated.append(
            compute_cov_score(random_targets, ppi_network, graph_rank, order)
        )
    z_score = (observed - np.mean(simulated)) / np.std(simulated)
    return z_score, norm.cdf(-abs(z_score))

def compute_cov_score(targets, ppi_network, graph_rank, order=1):
    """
    Compute the coverage score for a set of drug targets based on the
    proportion of nodes in the protein-protein interaction network that are
    covered by the network neighborhood of the targets.

    Parameters
    ----------
    targets : list
        List of protein names corresponding to the drug targets to consider.
    ppi_network : igraph.Graph
        Protein-protein interaction network represented as an igraph object.
    graph_rank : pandas.DataFrame
        DataFrame containing the graph rank of the nodes in the ppi network,
        with one column called "rank" and the row index corresponding to node
        names.
    order : int, optional
        The order of the neighborhood of each target to consider. Default is 1.

    Returns
    -------
    float
        The coverage score normalized by the rank of the nodes covered by the
        drug targets.
    """
    target_neighbor_list = dict()
    for t in targets:
        target_neighbor_list[t] = [t] + [
            ppi_network.vs["name"][n] 
                for n in ppi_network.neighborhood(t, order=order)
        ]
    # each target neighbor is only considered once
    covered_nodes_overall = list(set(it.chain(*[v for k,v in target_neighbor_list.items()])))
    # take the median of all neighbors considered based on the graph_rank
    inform_rank = graph_rank.loc[covered_nodes_overall].median(numeric_only=True)[0]
    # evaluate ppi network coverage
    coverage_overall = float(len(covered_nodes_overall)) / len(ppi_network.vs)
    return coverage_overall/inform_rank

def generate_random_targets(targets, deg_dist):
    """
    Randomly generate a set of targets with a degree distribution similar to
    the input set.

    Parameters
    ----------
    targets : list
        A list of target nodes to be used as a reference for the degree
        distribution.
    deg_dist : pandas.DataFrame
        A dataframe containing the degree distribution of the ppi network.

    Returns
    -------
    list
        A list of nodes with a degree distribution similar to the input set.

    """
    drugNodes_bins = deg_dist.loc[targets].groupby(by="bin").count()
    random_targets = []
    for row in drugNodes_bins.itertuples():
        bin_id, freq = row
        if freq == 0:
            continue
        candidates = deg_dist[deg_dist["bin"] == bin_id]
        random_targets.extend(
            candidates.sample(n=freq, replace=False).index.to_list()
        )
    return random_targets


