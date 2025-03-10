# -*- coding: utf-8 -*-
# date : 2023/2/26 
# author : wangh
import os
import numpy as np
import pandas as pd
from strainOptimizer.simulation import ppFBA


def k_matrix_filter(model, k_matrix, alpha, tol):
    # filter1.Take out rxns with no grRule:
    withGR = [rxn.id for rxn in model.reactions if rxn.gene_reaction_rule != '']
    k_matrix = k_matrix.loc[withGR,:]

    # filter2.remove out rxns that are always zero -> k=0/0=NaN:
    # all_nanList = [rxn for rxn in k_matrix.index if np.all(np.isnan(k_matrix.loc[rxn, :]))]
    all_nanList = k_matrix.apply(lambda x: np.all(pd.to_numeric(x, errors='coerce').isna()), axis=1)
    print('there are %d reactions with all NaNs' % np.sum(all_nanList))
    k_matrix = k_matrix.loc[~all_nanList]

    # 3.Replace remaining NaNs with 1s:
    # k_matrix[np.isnan(k_matrix)] = 1
    k_matrix = k_matrix.fillna(1)

    # 4.Replace any Inf value with 1000:
    k_matrix[np.abs(k_matrix) > 1000] = 1000

    # filter5.Filter out values that are inconsistent at different alphas:
    distinct_down = np.sum(k_matrix <= (1 - tol), axis=1) >= len(alpha)-2
    distinct_up = np.sum(k_matrix >= (1 + tol), axis=1) >=  len(alpha)-2
    # Identify those reactions with mixed patterns
    incons_rxns = distinct_down + distinct_up == 0
    print('there are %d reactions with mixed patterns' % np.sum(incons_rxns))
    k_matrix = k_matrix.loc[incons_rxns[~incons_rxns].index.tolist(), :]

    return k_matrix


def flux_scanning(model, targetID, c_source,c_uptake, alpha, filterG=False,model_type='etfl',tol=0.001):
    """
    Args:
        model : ecModel/ETFL model.
        targetID (str): Rxn ID for the production target reaction, a exchange reaction is recommended.
        c_source (str): Rxn ID for the main carbon source uptake reaction
        c_uptake(float): carbon source uptake rate
        alpha (np.array): scalling factor for production yield for enforced objective limits
        tol (float): numerical tolerance for fixing bounds
        filterG (bool, optional): TRUE if genes K_scores results should be filtered according to the alpha vector distribution
        model_type (str, optional): 'etfl' or 'ecGEM' model type

    Returns:
        dict: FC
            flux_WT (pd.Series): simulated flux for the wild-type (100% growth) condition
            alpha (np.array): scaling factors for production yield
            rxns (list): reaction information for each reaction in the model
            rxnNames (pd.Series): short name for each reaction in the model
            v_matrix (pd.DataFrame): matrix with simulated flux for each reaction and alpha value
            k_matrix (pd.DataFrame): matrix with normalized flux (k-score) for each reaction and alpha value
            grMatrix (pd.DataFrame): matrix with gene-reaction rules for each reaction
            k_rxns (pd.Series): median k-score across alpha values for each reaction
            genes (pd.Series): remaining genes after filtering inconsistent scores
            geneNames (pd.Series): short name for remaining genes
            k_genes (pd.DataFrame): k-score for each remaining gene
    """

    # check model_type
    if model_type not in ['etfl', 'ecGEM']:
        raise ValueError('model_type should be either "etfl" or "ecGEM"')
    if filterG is None:
        filterG = False

    if c_uptake is None:
        c_uptake = 1

    gluc_MW=0.180156 #g/mmol
    tol_ratio=0.01

    #step 1: build reactions k_matrix
    # Simulate WT (100% max growth):
    FC = {}
    if model_type == 'etfl':
        gr_rxnID = model.growth_reaction.id
    elif model_type == 'ecGEM':
        gr_rxnID = 'r_2111'
    # check if model has transcriptome attribute
    if hasattr(model, 'transcriptome'):
        print('integrating omic data to model')
        from strainOptimizer.manipulation.integration import integrate_omic_data_to_ecmodel
        # import numpy as np
        expression_threshold = np.percentile(model.transcriptome, 25)
        params = {'objective_reaction_id': gr_rxnID, 
                  'obj_frac': 0.4,
                  'expression_threshold':expression_threshold}
        with model:
            model = integrate_omic_data_to_ecmodel(model=model,
                                                    omic_data=model.transcriptome,
                                                    method='GIMME',
                                                    parameters=params)
        wt_sol = ppFBA(model=model,
                                   targetID=gr_rxnID,
                                   c_source=c_source,
                                   c_uptake=c_uptake,
                                   model_type=model_type,
                                   tol_ratio=tol_ratio)

    else:
        wt_sol = ppFBA(model=model,
                                   targetID=gr_rxnID,
                                   c_source=c_source,
                                   c_uptake=c_uptake,
                                   model_type=model_type,
                                   tol_ratio=tol_ratio)
    FC['flux_WT'] = wt_sol.fluxes
    # max_growth = FC['flux_WT'][gr_rxnID]
    # print('simulate WT-like:',max_growth)

    # simulate production in different suboptimal growth rate conditions
    FC['alpha'] = alpha
    # initialize fluxes and K_scores matrices
    rxnIDlist= [rxn.id for rxn in model.reactions]
    v_matrix = pd.DataFrame(index=rxnIDlist, columns=alpha)
    k_matrix = pd.DataFrame(index=rxnIDlist, columns=alpha)
    # simulate fluxes distribution in different growth rate conditions
    for i in range(len(alpha)):
        biomass_yield=alpha[i]
        growth=biomass_yield*c_uptake*gluc_MW
        model.reactions.get_by_id(gr_rxnID).bounds= growth*(1-tol_ratio), growth
        product_sol = ppFBA(model=model,
                                        targetID=targetID,
                                        c_source=c_source,
                                        c_uptake= c_uptake,
                                        model_type=model_type,
                                        tol_ratio=tol_ratio)
        FC['flux_MAX'] = product_sol.fluxes
        print('when growth rate is %s,production rate is %s'%(growth,product_sol.fluxes[targetID]))
        v_matrix.iloc[:, i] = FC['flux_MAX']
        k_matrix.iloc[:, i] = FC['flux_MAX'] / FC['flux_WT']

    # step 2: calculate reactions k scores
    # filter rxns according to k_matrix
    k_matrix_filtered = k_matrix_filter(model, k_matrix, alpha, tol)
    v_matrix_filtered = v_matrix.loc[k_matrix_filtered.index.tolist(), :]
    FC["v_matrix"] = v_matrix_filtered
    FC["k_matrix"] = k_matrix_filtered
    FC["rxns"] = k_matrix_filtered.index.tolist()
    # get median k_score for target rxn
    k_rxns=pd.Series(np.mean(k_matrix_filtered, axis=1),index=k_matrix_filtered.index.tolist())
    FC["k_rxns"] = k_rxns
    # Order from highest to lowest median k_score (across alphas)
    order = np.argsort(FC["k_rxns"])[::-1]
    FC["k_rxns"] = FC["k_rxns"][order]
    FC["rxns"] = FC['k_rxns'].index.tolist()
    FC['rxnNames'] = pd.Series([model.reactions.get_by_id(rxn).name for rxn in FC['rxns']], index=FC['rxns'])
    FC["v_matrix"] = FC["v_matrix"].loc[FC["rxns"], :]
    FC["k_matrix"] = FC["k_matrix"].loc[FC["rxns"], :]

    # step 3: calculate gene k scores
    # get gene scores according k_rxns
    # build gene-rxn matrix accoreding target rxns
    grMatrix = pd.DataFrame(index=[gene.id for gene in model.genes], columns=FC['rxns'])
    for rxn in FC['rxns']:
        related_genes = [gene.id for gene in model.reactions.get_by_id(rxn).genes]
        grMatrix[rxn] = np.array([gene.id in related_genes for gene in model.genes])
    genes_rxns_dict = {}
    k_genes=pd.Series()
    for gene in grMatrix.index.tolist():
        related_rxns= grMatrix.loc[gene, grMatrix.loc[gene, :] == True].index.tolist()
        if len(related_rxns) > 0:
            related_rxns_k=FC['k_rxns'].loc[related_rxns]
            # filter genes with inconsistent k_scores
            # check if all k_scores are below 1 or all above 1
            if np.sum(related_rxns_k <= (1 - tol)) == len(related_rxns) or np.sum(related_rxns_k >= (1 + tol)) == len(related_rxns):
                genes_rxns_dict[gene] = related_rxns
                k_genes[gene] = np.mean(k_matrix_filtered.loc[related_rxns, :].values.flatten())
            else:
                print('inconsistent k_scores for gene %s'%gene)
    # Order from highest to lowest median k_score (across alphas)
    order = np.argsort(k_genes)[::-1]
    k_genes = k_genes[order]
    FC["k_genes"] = k_genes
    FC["genes"]=FC["k_genes"].index.tolist()
    FC["geneNames"] = pd.Series([model.genes.get_by_id(geneID).name for geneID in FC["genes"]],index=FC["genes"])
    grMatrix= grMatrix.loc[FC["genes"], FC["rxns"]]
    FC["grMatrix"] = grMatrix

    # step 4(optional): # Filter any value between mean(alpha) and 1
    if filterG:
        unchanged = FC['k_genes'][(FC['k_genes'] >= np.mean(alpha) - tol) & (FC['k_genes'] <= 1 + tol)].index.tolist()
        # remove unchanged genes
        FC['genes'] = [g for g in FC['genes'] if g not in unchanged]
        FC['geneNames'] = FC['geneNames'][FC['genes']]
        FC['k_genes'] = FC['k_genes'][~unchanged]
        FC['grMatrix'] = FC['grMatrix'].loc[FC['genes'], :]
        # Update results for gene-related reactions (remove remaining reactions without any associated gene in rxnGeneM)
        FC['rxns'] = FC['grMatrix'].columns[FC['grMatrix'].any(axis=0)].tolist()
        FC['v_matrix'] = FC['v_matrix'].loc[FC['rxns'], :]
        FC['k_matrix'] = FC['k_matrix'].loc[FC['rxns'], :]
        FC['k_rxns'] = FC['k_rxns'][FC['rxns']]
        FC['rxnNames'] = [model.reactions.get_by_id(rxnID) for rxnID in FC['rxns']]

    return FC


def run_ecFSEOF(model, targetID, c_source,c_uptake, alphaLims, Nsteps,model_type='etfl'):
    """
    Run Flux-scanning with Enforced Objective Function for a specified production target.

    Arguments:
    - model:
        ETFL model
    - targetID: str
        Reaction ID for the production target reaction, a exchange reaction is recommended
    - c_source: str
        Reaction ID for the main carbon source uptake reaction (make sure that the correct directionality is indicated)
    - alphaLims: tuple of float
        Minimum and maximum biomass yield [gDw/mmol Csource] for enforced objective limits
    - Nsteps: int
        Number of steps for suboptimal objective in FSEOF
    - model_type: str
        Type of model (etfl or ecGEM)
    Returns:
    - dict
        Dictionary with the results of the FSEOF analysis
    """

    # Define alpha vector for suboptimal enforced objective values
    alphaV = np.linspace(alphaLims[0], alphaLims[1], Nsteps)
    # Run FSEOF analysis
    results = flux_scanning(model=model, targetID=targetID, c_source=c_source,c_uptake=c_uptake, alpha=alphaV,model_type=model_type)

    # Create gene table
    geneTable = pd.DataFrame(index=results['genes'], columns=["gene_names", "k_score"], dtype=float)
    geneTable.index.name = "gene_IDs"
    geneTable["gene_names"] = results["geneNames"]
    geneTable["k_score"] = results['k_genes']
    results["geneTable"] = geneTable

    # Create reaction table
    rxnTable = pd.DataFrame(index=results['rxns'], columns=["rxn_names", "k_score"], dtype=float)
    rxnTable.index.name = "rxn_IDs"
    rxnTable["rxn_names"] = results["rxnNames"]
    rxnTable["k_score"] = results['k_rxns']
    results["rxnTable"] = rxnTable

    # Remove redundant output fields
    results.pop('k_rxns', None)
    results.pop('rxns', None)
    results.pop('rxnNames', None)
    results.pop('genes', None)
    results.pop('geneNames', None)
    results.pop('k_genes', None)
    results.pop('flux_MAX', None)

    return results


