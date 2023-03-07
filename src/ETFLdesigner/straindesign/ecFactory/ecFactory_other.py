# -*- coding: utf-8 -*-
# date : 2023/3/7 
# author : wangh
import numpy as np
import pandas as pd
from etfl.io.json import load_json_model


def find_leaks(candidates, targetID, model):
    '''function to find reactions that consume the target.
    :param candidates: a pandas dataframe with the following columns:
        1. geneID: gene ID
        2. k_score: k-score of the gene
        3. actions: action type for the gene
    :param targetID: target product exchange reaction ID
    :param model: ETFL model
    :return: a pandas dataframe with the following columns:

    '''
    # assume target objective rxn must be the exchange rxn
    met_ex_rxn = model.reactions.get_by_id(targetID)
    met_e_id=met_ex_rxn.reactants[0].id
    met_c_id=met_e_id.replace('e','c')
    met_c=model.metabolites.get_by_id(met_c_id)
    with model:
        # set the target exchange reaction as objective
        model.objective = targetID
        # optimize
        sol = model.optimize()

    # Find reactions that consume the product in the cytoplasm
    met_c_sum = met_c.summary(solution=sol)
    consuming_rxnIDs= met_c_sum.consuming_flux.index.tolist()

    target_genes=[]
    for rxnID in consuming_rxnIDs:
        rxn=model.reactions.get_by_id(rxnID)
        # remove rxn to the product target
        met_idList=[met.id for met in list(rxn.metabolites.keys())]
        if met_e_id in met_idList:
            consuming_rxnIDs.remove(rxnID)
            continue
        else:
            # find gene taget for the rxn
            geneList=list(rxn.genes)
            geneIDlist=[gene.id for gene in geneList]
            # remove gene that has been included in candidates
            geneIDlist=[geneID for geneID in geneIDlist if geneID not in candidates.index.tolist()]

            target_genes=target_genes+geneIDlist

    df_new=pd.DataFrame({'geneID':target_genes,'k_score':np.nan,'actions':np.nan})
    df_new.set_index('geneID',inplace=True)
    df_new['k_score']=[0]*len(target_genes)
    df_new['actions']=['KO']*len(target_genes)
    print('%s leak rxn has been found and added to candidates.\n'%len(target_genes))
    # add to candidates
    candidates=candidates.append(df_new)

    return candidates


def remove_essential_targets(candidates,essential_path=r'code_etfl/ETFLdesigner/data/essential_genes.txt'):
    '''function to remove essential genes from candidates.
    :param candidates: a pandas dataframe with the following columns:
        1. geneID: gene ID
        2. k_score: k-score of the gene
        3. actions: action type for the gene
    :param essential_path: path to essential genes data(default:path for S.cerevisiae essential genes table)
    :return: a updated candidates dataframe
    '''
    # remove essential genes
    essentials = pd.read_csv(essential_path, sep='\t').Ids.str.strip()
    essential_targets = candidates.loc[candidates.index.isin(essentials)].index.tolist()
    print(f"Removing {len(essential_targets)} essential targets\n")
    candidates = candidates.drop(essential_targets)
    return candidates

if __name__ == '__main__':
    test_model= load_json_model('models/ecoli_core.json')
    yefl=load_json_model('models/yeast8_cEFL_2584_enz_128_bins__20221115_120238.json')
    model=test_model
    geneTable=pd.DataFrame(columns=['geneID','k_score','actions'])
    # geneTable=pd.read_excel('code_etfl/ETFLdesigner/output/yefl_2PE_design_results.xlsx',sheet_name='geneTable')
    targetID='r_1589'
    candidates1=find_leaks(candidates=geneTable, targetID=targetID, model=model)

