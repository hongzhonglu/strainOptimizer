# -*- coding: utf-8 -*-
import pandas as pd
import networkx as nx
import pickle
import os

FILE_PATH = os.path.dirname(os.path.abspath(__file__))


def build_met_rxn_graph_from_GEM(model, exclusions=None):
    # 1. Build a directed bipartite graph
    G = nx.DiGraph()
    if exclusions is None:
        # exclude: cofactors, O2, H2O, H+, ATP
        exclusions = ['oxygen', 'H2O', 'H+', 'ATP', 'ADP', 'NADH', 'NADPH', 'NADP+', 'NAD+', 'FADH2', 'FAD',
                     'ferricytochrome c', 'ferrocytochrome c', 'ubiquinol-6', 'ubiquinone-6','diphosphate']

    # add metabolite nodes
    for met in model.metabolites:
        name = met.name.split(' [')[0]
        if name in exclusions:
            continue
        # Add metabolite nodes (with prefix to distinguish)
        G.add_node(f"{met.id}", bipartite=0)

    # add reaction nodes and directed edges
    for rxn in model.reactions:
        if rxn.bounds == (0, 0):
            # Skip reactions that are not active (bounds are zero)
            continue
        # Add reaction nodes
        G.add_node(f"{rxn.id}", bipartite=1)

        # if rxn is reversible, add edges in both directions
        if rxn.reversibility:
            for met, coeff in rxn.metabolites.items():
                name = met.name.split(' [')[0]
                if name in exclusions:
                    continue
                G.add_edge(f"{met.id}", f"{rxn.id}")
                G.add_edge(f"{rxn.id}", f"{met.id}")
        else:
            if rxn.lower_bound == 0 and rxn.upper_bound > 0:
                direction = 1
            elif rxn.lower_bound < 0 and rxn.upper_bound == 0:
                direction = -1
            # Connect reactions to their associated metabolites
            for met, coeff in rxn.metabolites.items():
                name = met.name.split(' [')[0]
                if name in exclusions:
                    continue
                coeff = coeff * direction  # Adjust coefficient based on reaction direction
                if coeff < 0:
                    # Substrate: edge from metabolite -> reaction
                    G.add_edge(f"{met.id}", f"{rxn.id}")
                elif coeff > 0:
                    # Product: edge from reaction -> metabolite
                    G.add_edge(f"{rxn.id}", f"{met.id}")
    return G


def calculate_shortest_reactions_distance(G, reaction1_id: str, reaction2_id: str):
    """
    Calculate the shortest path distance between two reactions in a metabolic network.

    Args:
        G (networkx.Graph): The bipartite graph representing the metabolic network.
        reaction1_id (str): The ID of the first reaction.
        reaction2_id (str): The ID of the second reaction.

    Returns:
        int: The distance between the two reactions (number of intermediate reactions).
             Returns float('inf') if no path exists.
        None: If the input is invalid.
    """


    # 1. Calculate the shortest path
    source_node = f"{reaction1_id}"
    target_node = f"{reaction2_id}"
    try:
        # networkx returns the path length as the number of edges
        path_length = nx.shortest_path_length(G, source=source_node, target=target_node)
        # record the shortest path
        path = nx.shortest_path(G, source=source_node, target=target_node)
        # print(f"Shortest path from {reaction1_id} to {reaction2_id}: {path}")

        # The path is in the form R-M-R-M...R.
        # Path length 2: R1-M-R2, distance is 1 reaction.
        # Path length 4: R1-M1-R_intermediate-M2-R2, distance is 2 reactions.
        # Therefore, metabolic distance = path_length / 2
        metabolic_distance = path_length / 2
        return int(metabolic_distance)

    except nx.NetworkXNoPath:
        # If the two nodes are not connected
        return float('inf')
    except nx.NodeNotFound:
        print(f"Error: One of the nodes {source_node} or {target_node} not in graph.")
        return None


def calculate_genetic_target_distance(model,genetic_target, substrate_ID, productID,networkfilepath=FILE_PATH+'/../../data/slim_Sce_met_rxn_diGraph.pkl'):
    '''Genetic target distance is defined as the minimum distance between the predicted genetic target and either the substrate or the product in the metabolic topological network.
    Args:
        genetic_target (str): gene ID
        substrate_ID (str): The reaction ID of the substrate uptake reaction.
        productID (str): The reaction ID of the product output reaction.
        model (cobra.Model): The GEM model object.
    Returns:
        float: The minimum distance between the genetic target and either the substrate or the product.
    '''
    try:
        rxnList=[rxn.id for rxn in model.genes.get_by_id(genetic_target).reactions]
    except KeyError:
        print(f"Error: Gene {genetic_target} not found in the model.")
        return None, None
    try:
        with open(networkfilepath, "rb") as f:
            G= pickle.load(f)
    except FileNotFoundError:
        # build the directed bipartite graph from the GEM model
        G=build_met_rxn_graph_from_GEM(model)

    # calculate the distance to the substrate
    substrate_distanceList = [calculate_shortest_reactions_distance(G,substrate_ID, rxn_id) for rxn_id in rxnList]
    # remove None values (no path found)
    substrate_distanceList = [d for d in substrate_distanceList if d is not None]
    # use average distance if multiple reactions are associated with the gene
    substrate_distance = sum(substrate_distanceList) / len(substrate_distanceList) if substrate_distanceList else float('inf')

    # calculate the distance to the product
    product_distanceList = [calculate_shortest_reactions_distance(G, rxn_id, productID) for rxn_id in rxnList]
    # remove None values (no path found)
    product_distanceList = [d for d in product_distanceList if d is not None]
    # use average distance if multiple reactions are associated with the gene
    product_distance = sum(product_distanceList) / len(product_distanceList) if product_distanceList else float('inf')

    return substrate_distance, product_distance
    # # return the minimum distance
    # target_distance = min(substrate_distance, product_distance)
    # if target_distance == float('inf'):
    #     warnings.warn(f"No path found for genetic target {genetic_target} to substrate {substrate_ID} or product {productID}.")
    #     return None
    # else:
    #     return target_distance


def extract_connected_pathway(model, check_rxnList, target_id=None):
    """
    Extracts a connected pathway from provided reaction List in a GEM model through bipartite graph.

    Args:
        model (cobra.Model): GEM model.
        check_rxnList (list): List of reaction IDs to check for connectivity.
        target_id (str): Optional; if provided, use it to extract the target pathway.

    Returns:
        rxns (list): List of reaction IDs in the connected pathway.
    """
    # 1. build genome-scale metabolic network as bipartite graph
    G = nx.Graph()
    for met in model.metabolites:
        # Add metabolite nodes (with prefix to distinguish)
        G.add_node(f"m_{met.id}", bipartite=0)
    for rxn in model.reactions:
        # Add reaction nodes
        G.add_node(f"r_{rxn.id}", bipartite=1)
        # Connect reactions to their associated metabolites
        for met in rxn.metabolites:
            G.add_edge(f"r_{rxn.id}", f"m_{met.id}")

    # 2. find connected components
    neighbor_mets = set()
    for rxn_id in check_rxnList:
        try:
            # Get the reaction node
            rxn_node = f"r_{rxn_id}"
            # Find neighbors (metabolites) of the reaction node
            neighbors = set(G.neighbors(rxn_node))
            neighbor_mets.update(neighbors)
        except nx.NodeNotFound:
            print(f"Warning: Reaction {rxn_id} not found in the graph.")
            continue

    # 3. Extract the connected pathway
    sub_nodes= neighbor_mets | {f"r_{rxn_id}" for rxn_id in check_rxnList}
    subgraph = G.subgraph(sub_nodes)

    # check connectivity: if the subgraph is connected, return the reactions
    if nx.is_connected(subgraph):
        return check_rxnList
    # If the subgraph is not connected, check if a target reaction is specified
    else:
        # If a target reaction is specified, find the connected component containing it
        if target_id is not None:
            target_node = f"r_{target_id}"
            for comp in nx.connected_components(subgraph):
                if target_node in comp:
                    # If the target reaction is in this component, return the reactions in this component
                    return [node[2:] for node in comp if node.startswith('r_')]
        # If no target reaction is specified, return the largest connected component
        else:
            largest_component = max(nx.connected_components(subgraph), key=len)
            return [node[2:] for node in largest_component if node.startswith('r_')]


def count_adjacent_reactions(G, reaction_id):
    """
    统计与指定反应相邻的反应数量（通过共享代谢物连接）。
    会保留方向性：即一个反应的产物连接到另一个反应的底物。
    """
    node_id = f"{reaction_id}"
    if node_id not in G:
        raise ValueError(f"Reaction {reaction_id} not found in graph.")

    adjacent_reactions = set()

    # 出边方向：reaction -> metabolite -> other reaction
    for met in G.successors(node_id):       # 反应的产物代谢物
        for rxn in G.successors(met):       # 以该代谢物为底物的反应
            if rxn.startswith("r_") and rxn != node_id:
                adjacent_reactions.add(rxn)

    # 入边方向：reaction <- metabolite <- other reaction
    for met in G.predecessors(node_id):     # 反应的底物代谢物
        for rxn in G.predecessors(met):     # 生成该代谢物的反应
            if rxn.startswith("r_") and rxn != node_id:
                adjacent_reactions.add(rxn)

    return len(adjacent_reactions), adjacent_reactions


def count_adjacent_reactions_from_gene(model,gene_id,networkfilepath=FILE_PATH+'/../../data/slim_Sce_met_rxn_diGraph.pkl'):
    """
    # could the sum of adjacent reactions of all reactions associated with a gene
    """
    # load the bipartite graph from file
    with open(networkfilepath, 'rb') as f:
        G = pickle.load(f)
    adjacent_reactions = set()
    for rxn in model.genes.get_by_id(gene_id).reactions:
        rxnID= rxn.id
        try:
            count_adjacent, adjacent_rxns = count_adjacent_reactions(G, rxnID)
            adjacent_reactions.update(adjacent_rxns)
        except ValueError as e:
            print(f"Error processing reaction {rxnID}: {e}")
            continue
    return len(adjacent_reactions), adjacent_reactions

