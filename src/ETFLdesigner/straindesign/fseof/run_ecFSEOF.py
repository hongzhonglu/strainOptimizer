# -*- coding: utf-8 -*-
# date : 2023/2/26 
# author : wangh
import sys
sys.path.append(r"D:\code\github\ETFLdesigner\src\ETFLdesigner\straindesign\fseof")

import pandas as pd
import fseof

def run_ecFSEOF_design(model, modelParam, expYield,action_thresholds=[0.05,0.5,1.05],remove_essential=False):
    '''
    This function runs ecFSEOF method to identify gene targets for strain design
    :param model: ETFL model
    :param modelParam: a pandas series with the following parameters:
        targetID: target reaction ID
        c_source: carbon source exchange reaction ID
        c_uptake: carbon source uptake rate
    :param expYield: experimental yield of the target product
    :param action_thresholds: a list of three thresholds for gene targets:
        1. KO threshold
        2. KD threshold
        3. OE threshold
    :param remove_essential: a boolean value indicating whether to remove essential genes from the list of targets
    :return: a pandas dataframe with the following columns:
        1. geneID: gene ID
        2. k_score: K-score of the gene
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

    # read file with essential genes list
    if remove_essential:
        essential = pd.read_csv('../data/essential_genes.txt', sep='\t').Ids.str.strip()

    # Get relevant rxn indexes
    # modelParam.targetIndx = model.rxns.index(modelParam.rxnTarget)
    # modelParam.CUR_indx = model.rxns.index(modelParam.CSrxn)
    # modelParam.prot_indx = model.rxns.index('prot_pool_exchange')
    # modelParam.growth_indx = model.rxns.index(modelParam.growthRxn)


    # 1.- Run FSEOF to find gene candidates
    step = step + 1
    print(f'{step}.-  **** Running ecFSEOF method (ref: GECKO utilities) ****')
    results = fseof.run_FSEOF(model=model, targetID=modelParam['targetID'], c_source=modelParam['c_source'],c_uptake=modelParam['c_uptake'], alphaLims=alphaLims, Nsteps=Nsteps)
    genes = results['geneTable'].index.tolist()
    print('\n')
    print(f'ecFSEOF returned {len(genes)} targets')
    print('\n')

    # Format results table
    gene_result=results['geneTable']
    gene_result.loc[gene_result['k_score'] >= action_thresholds[2], 'actions'] = 'OE'
    gene_result.loc[gene_result['k_score'] <= action_thresholds[1], 'actions'] = 'KD'
    gene_result.loc[gene_result['k_score'] <= action_thresholds[0], 'actions'] = 'KO'
    results['geneTable'] = gene_result

    # get candidate gene related enzyme information

    # 2.- Add flux leak targets (those genes not optimal for production that may consume the product of interest.
    # (probaly extend the approach to inmediate precurssors)

    # 3.- discard essential genes from deletion targets

    # 4.- Construct Genes-metabolites network for classification of targets

    # 5.- enzyme usage variety analysis(EUVA)

    # 6.- EUVA for suboptimal biomasss production subject to a minimal (1%) production rate of the target product
    # and a unit CS uptake rate Get max biomass

    # 7.- combine candidate targets

    return results


    # MWeigths = []
    # # Identify candidate genes in model enzymes
    # print(' Extracting enzymatic information for target genes')
    # iB = [model.enzGenes.index(gene) if gene in model.enzGenes else -1 for gene in genes]
    # candidates = pd.DataFrame({'genes': genes,
    #                            'enzymes': [model.enzymes[i] if i >= 0 else '' for i in iB],
    #                            'shortNames': geneShorts,
    #                            'MWs': [model.MWs[i] if i >= 0 else float('nan') for i in iB],
    #                            'pathways': [model.pathways[i] if i >= 0 else '' for i in iB],
    #                            'actions': actions,
    #                            'k_scores': k_scores})
    # # Keep results that comply with the specified K-score thresholds
    # print(f"Removing targets {thresholds[0]} < K_score < {thresholds[1]}")
    # toKeep = (candidates.k_scores >= thresholds[1]) | (candidates.k_scores <= thresholds[0])
    # candidates = candidates[toKeep]
    # print(f" * {len(candidates)} gene targets remain")
    # print('\n')
    #
    #
    # # 2.- Add flux leak targets (those genes not optimal for production that may consume the product of interest.
    # # (probaly extend the approach to inmediate precurssors)
    # step += 1
    # print(f'{step}.-  **** Find flux leak targets to block ****')
    # candidates = find_flux_leaks(candidates, modelParam.targetIndx, model)
    # print(f' * {len(candidates)} gene targets remain\n')
    #
    # # 3.- discard essential genes from deletion targets
    # step += 1
    # print(f'{step}.-  **** Removing essential genes from KD and KO targets list ****')
    # _, iB = np.isin(candidates['genes'], essential)
    # toRemove = iB & (candidates['k_scores'] <= delLimit)
    # candidates = candidates[~toRemove].reset_index(drop=True)
    # print(f' * {len(candidates)} gene targets remain\n')
    #
    # candidates.to_csv(f'{results_folder}/candidates_L1.txt', sep='\t', quoting=False, index=False)
    # proteins = 'draw_prot_' + candidates['enzymes']
    # _, enz_pos = np.isin(proteins, model.rxns)
    # candidates['enz_pos'] = enz_pos
    #
    #
    # # 4.- Construct Genes-metabolites network for classification of targets
    # step += 1
    # print(f'{step}.-  **** Construct Genes-metabolites network for classification of targets ****\n')
    #
    # # Get Genes-metabolites network
    # print('  Constructing genes-metabolites graph\n')
    # GeneMetMatrix, _, Gconect = getMetGeneMatrix(model, candidates['genes'])
    #
    # # Get independent genes from GeneMetMatrix
    # print('  Obtain redundant vectors in genes-metabolites graph (redundant targets)\n')
    # indGenes, G2Gmatrix, _ = getGeneDepMatrix(GeneMetMatrix)
    #
    # # Find unique targets (with no isoenzymes or not part of complexes)
    # candidates['unique'] = indGenes
    #
    # # Number of metabolites connected to each gene
    # candidates['conectivity'] = Gconect.mets_number
    #
    # # Get gene target groups (those connected to exactly the same metabolites)
    # _, groups = getGenesGroups(G2Gmatrix, candidates['genes'])
    # candidates['groups'] = groups