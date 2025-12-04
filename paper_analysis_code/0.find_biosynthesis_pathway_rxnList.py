# -*- coding: utf-8 -*-
'''Find linear biosynthesis pathway reactions for each product'''
from strainOptimizer.io import load_model
import networkx as nx
from strainOptimizer.manipulation.constraint.total_resource_allocation import constrain_enzymes
from strainOptimizer.strainDesign.ecFactory import run_ecFactory
import pandas as pd
import numpy as np


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

modelParams_dict={
    'ecGEM':{
    'model_type':'ecGEM',
    'c_source':"r_1714_REV",      # glucose exchange rxn
    'c_uptake':10,
    'growth_id':'r_2111',
    'total_enzymes':0.1,
        }
}

productParam_dict={
    '2-phenylethanol':{'productName':'2-phenylethanol',
                       'targetID':'r_1589',
                       'ecGEM_filepath':'examples/models/yeast/ecYeastGEM_batch.xml',
                       'etfl_filepath':'examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json'},
    # 'artemisinic_acid':{'productName':'artemisinic acid',
    #                'targetID':'SK_atemisinic_acid',
    #                'model_filepath':'examples/models/yeast/ecGEM_atemisinic.xml'},
    'heme_acid':{'productName':'heme',
                 'targetID':'EX_heme_a'
                 ,'ecGEM_filepath':'examples/models/yeast/heme_ecYeastGEM.xml',
                 'etfl_filepath':'examples/models/yeast/heme_cEFL.json'},
    'spermidine':{'productName':'spermidine',
                  'targetID':'r_2051',
                  'ecGEM_filepath':'examples/models/yeast/ecYeastGEM_batch.xml',
                'etfl_filepath':'examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json'
                  },
    'sclareol':{'productName':'sclareol',
                'targetID':'DM_sclareol_c',
                'ecGEM_filepath':'examples/models/yeast/sclareol_ecYeastGEM_batch.xml',
                'etfl_filepath':'examples/models/yeast/sclareol_cEFL.json'},
    # 'fatty_acid':{'productName':'free fatty acids',
    #                'targetID':'r_2189',
    #                'ecGEM_filepath':'examples/models/yeast/ecYeastGEM_batch.xml',
    #               'etfl_filepath':'examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json'},
}

products_pathway_dict={}
for product_key,productParam in productParam_dict.items():
    productName=productParam['productName']
    productID=productParam['targetID']
    ecGEM_filepath=productParam['ecGEM_filepath']

    modelParams= modelParams_dict['ecGEM']
    c_source = modelParams['c_source']
    c_uptake = modelParams['c_uptake']
    print('product:', product_key)

    model= load_model(ecGEM_filepath, model_type='ecGEM')
    with model:
        model.objective = productID
        sol = model.optimize()
        df_fluxes = sol.fluxes
        df_fluxes = df_fluxes.abs()
    potential_pathway_rxnList = df_fluxes[df_fluxes == df_fluxes[productID]].index.tolist()
    # keep only reactions with int value in df_multiple
    # Extract connected pathway reactions
    pathway_rxnList = extract_connected_pathway(model, potential_pathway_rxnList, target_id=productID)
    pathway_geneList = []
    for rxn_id in pathway_rxnList:
        rxn = model.reactions.get_by_id(rxn_id)
        pathway_geneList.extend([gene.id for gene in rxn.genes])
    pathway_geneList = list(set(pathway_geneList))  # remove duplicates
    products_pathway_dict[product_key] = pathway_geneList

# save it as json file
import json
with open('data/products_pathway_targets.json', 'w') as f:
    json.dump(products_pathway_dict, f)
