# -*- coding: utf-8 -*-
# date : 2023/3/14 
# author : wangh
# file : find_min_combined_sets.py
# project : etfl
'''used to find the minimal sets of genes that can be combined to achieve the optimal production performance
'''
import pandas as pd
import numpy as np
from etfl.optim.constraints import ModelConstraint
from pytfa.optim.utils import symbol_sum
from etfl.optim.utils import safe_optim


def cal_max_yield(model,targetID,c_source,c_uptake,tol=1e-6):
    '''calculate the maximum yield of the target product
    parameters:
        model: ETFL model
        targetID: str, the objective reaction ID
        c_source: str, the carbon source uptake reaction ID
        c_uptake: float, the carbon source uptake rate
    return:
        max_yield: float, the maximum yield of the target product
    '''
    # 1. fix carbon source uptake
    model.reactions.get_by_id(c_source).bounds=-c_uptake,-c_uptake
    # 2. calculate optimal production rate
    model.objective=targetID
    model.objective_direction='max'
    sol=safe_optim(model)
    max_rate=sol.objective_value
    if model.solver.status == 'infeasible':
        return 0,0
    else:
        # 3. fix the max target product production and minimize the C source uptake
        model.reactions.get_by_id(targetID).bounds=max_rate-tol,max_rate
        model.objective=c_source
        model.objective_direction='max'
        opt_c_uptake=-model.slim_optimize()
        max_yield=max_rate/opt_c_uptake

        # 4. reset the original constraints
        model.reactions.get_by_id(c_source).bounds=-c_uptake,-c_uptake
        model.objective=model.growth_reaction.id
        model.objective_direction='max'

        return max_yield, max_rate


def constrain_enz_conc(model,enzymes_bounds):
    '''add constraint for enzyme concentration in ETFL model
    parameters:
        model: ETFL model
        enzymes_bounds: pd.DataFrame, columns=['lb','ub'], index=enzymeID.!! lb and ub should be the scaled values
    return:
        model: modified ETFL model
    '''
    enzymeIDlist=enzymes_bounds.index.tolist()
    for enzID in enzymeIDlist:
        enz_vars = model.get_variables_of_type('EnzymeVariable')
        enz_var=enz_vars.get_by_id(enzID)
        exp=symbol_sum([enz_var])
        model.add_constraint(
            kind=ModelConstraint,
            hook=model,
            expr=exp,
            id_='enz_conc_'+enzID,
            lb=enzymes_bounds.loc[enzID,'lb'],
            ub=enzymes_bounds.loc[enzID,'ub']
        )
        # enz_const=model.constraints.get_by_id('MODC_enz_conc_'+enzID)
    model.repair()
    return model


def find_min_set(model,c_source,c_uptake,expYield,targetID,geneIDlist,gene_enz_fva_result,gene_enz_dict,tol=1e-6):
    #step 1. construct optimal production mutant
    mutant_model=model.copy()
    target_genes_enzfva_result=gene_enz_fva_result.loc[geneIDlist]
    df_enz_bounds=pd.DataFrame(columns=['lb','ub'])
    for geneID in geneIDlist:
        # enzID=str(gene_enz_dict.loc[geneID,:].values).split('\'')[1]
        enzID=gene_enz_dict[geneID][0]
        df_enz_bounds.loc[enzID,'lb']=target_genes_enzfva_result.loc[geneID,'prod_minprot']
        df_enz_bounds.loc[enzID,'ub']=target_genes_enzfva_result.loc[geneID,'prod_max']

    mutant_model=constrain_enz_conc(mutant_model,enzymes_bounds=df_enz_bounds)

    # calculate optimal production yield
    # fix carbon source uptake
    mutant_model.reactions.get_by_id(c_source).bounds=-c_uptake,-c_uptake
    # fix experimental biomass yield
    c_source_MW=0.180156   # g/mmol
    exp_gr=expYield*c_uptake*c_source_MW
    mutant_model.growth_reaction.bounds=exp_gr,exp_gr
    # calculate optimal production yield,and optimal production rate
    opt_prod_yield, opt_prod_rate=cal_max_yield(mutant_model,targetID,c_source,c_uptake,tol=tol)
    print('optimal production yield is: ',opt_prod_yield)
    print('optimal production rate is: ',opt_prod_rate)
    optimal_prod={'opt_prod_yield':opt_prod_yield,'opt_prod_rate':opt_prod_rate}

    #step 2. respectively convert each enzyme concentraction to WT-like constraints,and calculate the production yield ,and production rate
    df_min_set_result=pd.DataFrame(columns=['mod_prod_yield','mod_prod_rate','score'])
    for gene in geneIDlist:
        # enzID = str(gene_enz_dict.loc[gene, :].values).split('\'')[1]
        enzID=gene_enz_dict[gene][0]
        enz_conc_constriant=mutant_model.constraints['MODC_enz_conc_'+enzID]
        prod_ub=enz_conc_constriant.ub
        prod_lb=enz_conc_constriant.lb
        # convert to WT-like constraint
        enz_conc_constriant.lb=gene_enz_fva_result.loc[gene,'wt_minprot']
        enz_conc_constriant.ub=gene_enz_fva_result.loc[gene,'wt_max']
        # calculate mod_prod_yield and mod_production_rate
        mod_prod_yield, mod_prod_rate=cal_max_yield(mutant_model,targetID,c_source,c_uptake,tol=tol)
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
        enz_conc_constriant.ub=prod_ub
        enz_conc_constriant.lb=prod_lb

    return df_min_set_result,optimal_prod





# # test
# geneIDlist=results['geneTable'][results['geneTable']['target_priority_leval']==2].index.tolist()
# gene_enz_fva_result=results['gene_enz_fva_result']
# gene_enz_dict=results['gene_enz_dict']