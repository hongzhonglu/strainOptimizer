# -*- coding: utf-8 -*-
# date : 2023/2/26 
# author : wangh
import sys
sys.path.append(r"D:\code\github\etfl\code_etfl\ETFLdesigner\ecFactory")

import pandas as pd
import fseof
from ecFactory_other import find_leaks,remove_essential_targets,getMetGeneMatrix,getGeneDepMatrix,getGenesGroups,enzymeFVA,genelist_to_enzymelist,minprotFBA_prot_conc
from etfl.optim.utils import safe_optim

def run_ecFSEOF_design(model, modelParam, expYield,action_thresholds=[0.05,0.5,1.05],remove_essential=False,model_type='etfl'):
    '''
    This function runs ecFSEOF method to identify gene targets for strain design
    :param model: ETFL model
    :param modelParam: a pandas series with the following parameters:
        targetID: target reaction ID
        productName: product name(used to find the leak reactions)
        c_source: carbon source exchange reaction ID
        c_uptake: carbon source uptake rate
    :param expYield: experimental yield of the target product
    :param action_thresholds: a list of three thresholds for gene targets:
        1. KO threshold
        2. KD threshold
        3. OE threshold
    :param remove_essential: a boolean value indicating whether to remove essential genes from the list of targets
    :param model_type: a string indicating the type of model ('etfl' or 'ecYeast')
    :return: a pandas dataframe with the following columns:
        1. geneID: gene ID
        2. k_score: k-score of the gene
        3. actions: action type for the gene
    '''
    # method parameters
    tol = 1E-13  # numeric tolerance for determining non-zero enzyme usages
    OEF = 2  # overexpression factor for enzyme targets
    KDF = 0.5  # down-regulation factor for enzyme targets
    step = 0

    # Parameters for FSEOF method
    Nsteps = 16  # number of FBA steps in ecFSEOF
    alphaLims = [0.5 * expYield, 2 * expYield]  # biomass yield limits for ecFSEOF
    # thresholds = [0.5, 1.05]  # K-score thresholds for valid gene targets
    # delLimit = 0.05  # K-score limit for considering a target as deletion

    # 1.- Run FSEOF to find gene candidates
    step = step + 1
    print(f'{step}.-  **** Running ecFSEOF method (ref: GECKO utilities) ****')
    results = fseof.run_FSEOF(model=model, targetID=modelParam['targetID'], c_source=modelParam['c_source'],c_uptake=modelParam['c_uptake'], alphaLims=alphaLims, Nsteps=Nsteps,model_type=model_type)
    # Format results table
    gene_result=results['geneTable']
    gene_result.loc[gene_result['k_score'] >= action_thresholds[2], 'actions'] = 'OE'
    gene_result.loc[gene_result['k_score'] <= action_thresholds[1], 'actions'] = 'KD'
    gene_result.loc[gene_result['k_score'] <= action_thresholds[0], 'actions'] = 'KO'
    # remove genes with no action
    gene_result = gene_result.loc[gene_result['actions'].notnull()]
    print(f'ecFSEOF returned {len(gene_result)} targets\n')
    results['geneTable'] = gene_result

    # 2.- Add flux leak targets (those genes not optimal for production that may consume the product of interest.
    # (probaly extend the approach to inmediate precurssors)
    step += 1
    print(f'{step}.-  **** Find flux leak targets to block ****')
    results['geneTable'] = find_leaks(candidates=results['geneTable'], targetID=modelParam['targetID'], model=model, product_name=modelParam['productName'])

    # 3.- discard essential genes from deletion targets
    if remove_essential:
        step += 1
        print(f'{step}.-  **** Removing essential targets ****')
        # print(type(results['geneTable']))
        results['geneTable'] = remove_essential_targets(candidates=results['geneTable'])

    # 4.- Construct Genes-metabolites network for classification of targets
    step += 1
    print(str(step) + '.-  **** Construct Genes-metabolites network for classification of targets ****')
    # Get Genes-metabolites network
    metGeneMatrix, metsConectivity, genesConectivity=getMetGeneMatrix(model=model,geneIDlist=results['geneTable'].index.tolist())
    # Get independent genes from GeneMetMatrix
    independant_genes,gene_equal_Matrix=getGeneDepMatrix(metGeneMatrix)
    # Get gene target groups (those connected to exactly the same metabolites)
    groups = getGenesGroups(gene_equal_Matrix)
    print('independent targets: ' + str(len(independant_genes[independant_genes==1])))
    print('target groups: ' + str(len(groups)))
    # print('\n')
    results['independent_genes'] = independant_genes
    results['groups'] = groups



    # 5.- enzyme usage variety analysis(EUVA)
    step += 1
    step = step + 1
    print(str(step) + '.-  **** Running EUVA analysis ****')

    # get target enzyme list
    gene_targetList=results['geneTable'].index.tolist()
    target_enz_list,gene_enz_dict=genelist_to_enzymelist(model=model,genelist=gene_targetList)
    results['gene_enz_dict']=gene_enz_dict
    # Fix unit C source uptake
    c_uptake = modelParam['c_uptake']
    c_source= modelParam['c_source']
    gluc_MW=0.180156  # g/mmol
    model.reactions.get_by_id(c_source).bounds=-c_uptake, -c_uptake
    print(f'  - Fixed unit glucose uptake rate at {c_uptake} mmol/gDW.h')

    print(str(step) + '.1-  **** Running EUVA for optimal production conditions ****')
    print('  - Fixed suboptimal biomass production, according to provided experimental yield')
    # Fix suboptimal experimental biomass yield conditions
    fix_gr = expYield * gluc_MW* c_uptake
    model.growth_reaction.bounds = fix_gr-tol, fix_gr+tol
    print(' Fix suboptimal experimental biomass = ' + str(fix_gr) + ' h-1')
    print('  - Maximize the production rate of the target product')
    # set maximize production rate as objective
    targetID= modelParam['targetID']
    model.objective = targetID
    model.objective_direction = 'max'
    sol=safe_optim(model)
    max_prod=sol.objective_value
    # run parsimonious protein usages FBA
    prod_minprotFBA_protconc=minprotFBA_prot_conc(model=model,target=model.growth_reaction.id,enzymeIDlist=target_enz_list,c_source=c_source,c_uptake=c_uptake)
    # run enzyme usage variety analysis
    prod_enz_fva_result=enzymeFVA(model=model,enzymeIDlist=target_enz_list)
    prod_enz_fva_result['minprotFBA']=prod_minprotFBA_protconc
    results['prod_enz_fva_result']=prod_enz_fva_result      # can be deleted

    # - EUVA for suboptimal biomasss production subject to a minimal (1%) production rate of the target product
    # and a unit CS uptake rate Get max biomass
    print(str(step) + '.2-  **** Running EUVA for optimal biomass growth condition ****')
    # fix production rate as 1% of the max production rate
    model.reactions.get_by_id(targetID).bounds = max_prod*0.01, max_prod*0.01
    print('  - Maximize biomass production')
    # set maximize biomass production as objective
    model.objective = model.growth_reaction
    model.objective_direction = 'max'
    sol=safe_optim(model)
    max_gr=sol.objective_value
    print(' Max biomass = ' + str(max_gr) + ' h-1')
    # run parsimonious protein usages FBA
    wt_minprotFBA_protconc=minprotFBA_prot_conc(model=model,target=targetID,enzymeIDlist=target_enz_list,c_source=c_source,c_uptake=c_uptake)
    # run enzyme usage variety analysis
    wt_enz_fva_result=enzymeFVA(model=model,enzymeIDlist=target_enz_list)
    wt_enz_fva_result['minprotFBA']=wt_minprotFBA_protconc
    results['wt_enz_fva_result']=wt_enz_fva_result     # can be deleted

    # discard some targets according to the EUVA results
    print(str(step) + '.3-  **** Discarding targets according to EUVA results ****')
    # build a dataframe with all candidate genes  EUVA results
    gene_enz_fva_result=pd.DataFrame(index=results['geneTable'].index.tolist(),columns=['prod_min','prod_max','prod_minprot','wt_min','wt_max','wt_minprot'])
    # gene_enz_dict=results['gene_enz_dict']
    # prod_enz_fva_result=results['prod_enz_fva_result']
    # wt_enz_fva_result=results['wt_enz_fva_result']
    for gene in results['geneTable'].index.tolist():
        if len(gene_enz_dict[gene])==1:
            enzID=gene_enz_dict[gene][0]
            enz_fva_result=list(prod_enz_fva_result.loc[enzID].values)+list(wt_enz_fva_result.loc[enzID].values)
            gene_enz_fva_result.loc[gene]=enz_fva_result
    results['gene_enz_fva_result']=gene_enz_fva_result

    # discard targets with no enzyme usage variety


    # 7.- combine candidate targets

    return results

