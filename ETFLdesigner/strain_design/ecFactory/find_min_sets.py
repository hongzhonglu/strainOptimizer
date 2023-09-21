# -*- coding: utf-8 -*-
# date : 2023/3/14 
# author : wangh
# file : find_min_combined_sets.py
# project : etfl
'''used to find the minimal sets of genes that can be combined to achieve the optimal production performance
'''
import pandas as pd
import numpy as np
from ETFLdesigner.analysis import optimal_yield
from ETFLdesigner.manipulation.constraint import enzyme


def find_min_set(model,c_source,c_uptake,expYield,targetID,geneIDlist,gene_enz_fva_result,gene_enz_dict,model_type='etfl',tol_ratio=0.01):

    #step 1. construct optimal production mutant
    mutant_model=model.copy()
    target_genes_enzfva_result=gene_enz_fva_result.loc[geneIDlist]
    df_enz_bounds=pd.DataFrame(columns=['lb','ub'])
    for geneID in geneIDlist:
        enzID=gene_enz_dict[geneID][0]
        df_enz_bounds.loc[enzID,'lb']=target_genes_enzfva_result.loc[geneID,'prod_minprot']
        df_enz_bounds.loc[enzID,'ub']=target_genes_enzfva_result.loc[geneID,'prod_max']
    if model_type=='etfl':
        mutant_model=enzyme.ETFL_constrain_enz_conc(mutant_model,enzymes_bounds=df_enz_bounds,tol_ratio=tol_ratio)
    elif model_type=='ecGEM':
        mutant_model=enzyme.ecGEM_constrain_enz_conc(mutant_model,enzymes_bounds=df_enz_bounds,tol_ratio=tol_ratio)

    # calculate optimal production yield
    # fix carbon source uptake
    if model_type=='etfl':
        mutant_model.reactions.get_by_id(c_source).bounds=-c_uptake,-c_uptake
        growth_id = mutant_model.growth_reaction.id
    elif model_type=='ecGEM':
        mutant_model.reactions.get_by_id(c_source).bounds=c_uptake,c_uptake
        growth_id = 'r_2111'

    # fix experimental biomass yield
    c_source_MW=0.180156   # g/mmol
    exp_gr=expYield*c_uptake*c_source_MW
    mutant_model.reactions.get_by_id(growth_id).bounds=exp_gr,exp_gr
    # calculate optimal production yield,and optimal production rate
    opt_prod_yield, opt_prod_rate=optimal_yield.cal_max_yield(model=mutant_model,
                                                              targetID=targetID,
                                                              c_source=c_source,
                                                              c_uptake=c_uptake,
                                                              model_type=model_type)
    print('optimal production yield is: ',opt_prod_yield)
    print('optimal production rate is: ',opt_prod_rate)
    optimal_prod={'opt_prod_yield':opt_prod_yield,'opt_prod_rate':opt_prod_rate}

    #step 2. respectively convert each enzyme concentraction to WT-like constraints,and calculate the production yield ,and production rate
    df_min_set_result=pd.DataFrame(columns=['mod_prod_yield','mod_prod_rate','score'])
    for gene in geneIDlist:
        enzID = gene_enz_dict[gene][0]
        if enzID=='no enzyme':
            print('gene %s has no enzyme'%gene)
            continue
        # constrain the enzyme concentration into WT-like respectly
        if model_type=='etfl':
            enz_conc_constriant=mutant_model.constraints['MODC_enz_conc_'+enzID]
            prod_ub=enz_conc_constriant.ub
            prod_lb=enz_conc_constriant.lb
            # convert to WT-like constraint
            wt_ub=gene_enz_fva_result.loc[gene,'wt_max']
            wt_lb=gene_enz_fva_result.loc[gene,'wt_minprot']
            if wt_lb>wt_ub:
                wt_lb=wt_lb*(1-tol_ratio)
                wt_ub=wt_ub*(1+tol_ratio)
            enz_conc_constriant.bounds=wt_lb,wt_ub
        elif model_type=='ecGEM':
            prot_pseudo_rxn=mutant_model.reactions.get_by_id(enzID)
            prod_bounds=prot_pseudo_rxn.bounds
            wt_lb=gene_enz_fva_result.loc[gene,'wt_minprot']
            wt_ub=gene_enz_fva_result.loc[gene,'wt_max']
            if wt_lb>wt_ub:
                wt_ub=wt_ub*(1+tol_ratio)
                wt_lb=wt_lb*(1-tol_ratio)
            mutant_model.reactions.get_by_id(enzID).bounds=wt_lb,wt_ub
        # calculate mod_prod_yield and mod_production_rate
        mod_prod_yield, mod_prod_rate=optimal_yield.cal_max_yield(mutant_model,
                                                                  targetID,
                                                                  c_source,
                                                                  c_uptake,
                                                                  model_type=model_type)
        # calculate score
        score1=mod_prod_yield/opt_prod_yield
        score2=mod_prod_rate/opt_prod_rate
        score=np.mean([score1,score2])
        #round to four decimal
        score=round(score,4)
        df_min_set_result.loc[gene,'mod_prod_yield']=mod_prod_yield
        df_min_set_result.loc[gene,'mod_prod_rate']=mod_prod_rate
        df_min_set_result.loc[gene,'score']=score

        # convert back to original constraint
        if model_type=='etfl':
            enz_conc_constriant.ub=prod_ub
            enz_conc_constriant.lb=prod_lb
        elif model_type=='ecGEM':
            mutant_model.reactions.get_by_id(enzID).bounds=prod_bounds

    return df_min_set_result,optimal_prod




