import cobra
from cobra import Model, Reaction, Metabolite
import os
from os.path import join
import numpy as np
import pandas as pd
import copy
import re
from cobra.io import load_matlab_model, read_sbml_model

siwei = pd.read_excel('examples/siwei_heme_FSEOF/heme production enhance targets.xlsx')
ivan = pd.read_excel('examples/siwei_heme_FSEOF/genesResults_ecFSEOF.xlsx')
haoyu = pd.read_excel('examples/siwei_heme_FSEOF/haoyu.xlsx')

haoyu1 = haoyu[haoyu['k_score'] >1]['gene_IDs'].tolist()
ivan1 = ivan[ivan['K_score'] >1]['gene_IDs'].tolist()
siwei1 = siwei[siwei['k_scores'] >1]['gene'].tolist()

#compare difference

siwei_vs_ivan = list(set(siwei1)-set(ivan1))
ivan_vs_siwei = list(set(ivan1)-set(siwei1))
ivan_siwei_common = list(set(ivan1) & set(siwei1))


haoyu_vs_ivan = list(set(haoyu1)-set(ivan1))
ivan_haoyu_common = list(set(ivan1) & set(haoyu1))
ivan_vs_haoyu = list(set(ivan1)-set(haoyu1))

# then check which reaction connect with the genes which only exist in ivan procedures.
def getRxnByGene(model, gene0):
    """
    :param model: A metabolic model
    :param gene0: A gene
    :return:

    Example:
    getRxnByGene(model, gene0="YGR192C")

    """
    rxn_list2 = []
    for gene in model.genes:
        if gene.id == gene0:
            rxn_list = gene.reactions
            for x in rxn_list:
                rxn_list2.append(x.id)
    return rxn_list2



#model = cobra.io.load_matlab_model("examples/siwei_heme_FSEOF/y8.mat")

# heme model
ecYeast = read_sbml_model('examples/models/yeast/heme_ecYeastGEM.xml')
rxn_all = []
for i in ivan_vs_haoyu:
    #print(i)
    rxn0 = getRxnByGene(model=ecYeast, gene0=i)
    if len(rxn0) < 1:
        print(i)
    rxn_all.append(rxn0)
rxn_all = sum(rxn_all,[])
# input haoyu flux analysis
haoyu_flux = pd.read_excel('examples/siwei_heme_FSEOF/v_matrix.xlsx')
haoyu_k_score = pd.read_excel('examples/siwei_heme_FSEOF/k_matrix.xlsx')

rxn_all.append('EX_heme_a')
rxn_all.append('r_2111')
haoyu_flux_filter = haoyu_flux[haoyu_flux['rxnID'].isin(rxn_all)]
haoyu_flux_filter.to_excel("examples/siwei_heme_FSEOF/flux_check.xlsx")


haoyu_k_score_filter = haoyu_k_score[haoyu_k_score['rxnID'].isin(rxn_all)]
haoyu_k_score_filter.to_excel("examples/siwei_heme_FSEOF/haoyu_k_score_check.xlsx")






