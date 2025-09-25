# -*- coding: utf-8 -*-
import pickle
from strainOptimizer.io import load_model
import cobra
import os
import networkx as nx
import warnings
from strainOptimizer.analysis.network import calculate_shortest_reactions_distance,build_met_rxn_graph_from_GEM,count_adjacent_reactions_from_gene

os.chdir('D:\code\github\strainOptimizer')



def calculate_genetic_target_distance(model,genetic_target, substrate_ID, productID,networkfilepath=r'data/Sce_met_rxn_diGraph.pkl'):
    '''Genetic target distance is defined as the minimum distance between the predicted genetic target and either the substrate or the product in the metabolic topological network.
    Args:
        genetic_target (str): gene ID
        substrate_ID (str): The reaction ID of the substrate uptake reaction.
        productID (str): The reaction ID of the product output reaction.
        model (cobra.Model): The GEM model object.
    Returns:
        float: The minimum distance between the genetic target and either the substrate or the product.
    '''
    rxnList=[rxn.id for rxn in model.genes.get_by_id(genetic_target).reactions]
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



model= load_model('examples/models/yeast/yeast-GEM.xml', model_type='gem')
with open('data/Sce_met_rxn_diGraph.pkl', "rb") as f:
    G = pickle.load(f)
# model=load_model('examples/models/yeast/sclareol_ecYeastGEM_batch.xml',model_type='ecGEM')
# target_id='r_1589'  # 2-phenylethanol
target_id='DM_sclareol_c'  # sclareol
c_source='r_1714'  # glucose

import pandas as pd
df=pd.read_csv('data/experiment_targets/2-phenylethanol_exp_targets.tsv',sep='\t')
geneList=df['geneID'].tolist()
gene_distance_dict={}
# for gene in model.genes:
for gene_id in geneList:
    distance=calculate_genetic_target_distance(genetic_target=gene_id, substrate_ID=c_source, productID=target_id, model=model)
    gene_distance_dict[gene_id]=distance
    # print(gene_id,distance)

    # calculate adjacent reactions

    count,adj_rxnList=count_adjacent_reactions_from_gene(G=G,model=model,gene_id=gene_id)
    print(gene_id,count,adj_rxnList)