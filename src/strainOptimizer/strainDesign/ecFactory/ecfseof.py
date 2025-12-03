# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
from strainOptimizer.simulation import ppFBA,moma,mopa,pFBA
from strainOptimizer.strainDesign.ecFactory.ecFactory_other import default_scanning_range


def k_matrix_filter(model, k_matrix, tol):
    """
    Filter k_matrix to remove reactions with problematic patterns.
    
    Args:
        model: Metabolic model
        k_matrix (pd.DataFrame): Matrix with k-scores for reactions across alpha values
        tol (float): Tolerance threshold for filtering
        
    Returns:
        pd.DataFrame: Filtered k_matrix
    """
    # Filter 1: Remove reactions without gene-reaction rules
    reactions_with_gr = [rxn.id for rxn in model.reactions if rxn.gene_reaction_rule]
    k_matrix = k_matrix.loc[reactions_with_gr, :]
    
    # Filter 2: Remove reactions with all NaN values (zero fluxes)
    nan_mask = k_matrix.apply(lambda x: np.all(pd.to_numeric(x, errors='coerce').isna()), axis=1)
    nan_count = np.sum(nan_mask)
    print(f'Found {nan_count} reactions with all NaN values')
    k_matrix = k_matrix.loc[~nan_mask]
    
    # Filter 3: Replace remaining NaN values with 1 (neutral effect)
    k_matrix = k_matrix.fillna(1)
    
    # Filter 4: Cap extreme values to prevent numerical issues
    k_matrix = k_matrix.clip(-1000, 1000)
    
    # Filter 5: Remove reactions with inconsistent patterns across alpha values
    # A reaction is consistent if it's mostly up-regulated OR mostly down-regulated
    down_regulated = np.sum(k_matrix <= (1 - tol), axis=1) >= k_matrix.shape[1] - 2
    up_regulated = np.sum(k_matrix >= (1 + tol), axis=1) >= k_matrix.shape[1] - 2
    
    # Keep only reactions that show consistent regulation pattern
    consistent_mask = down_regulated | up_regulated
    inconsistent_count = np.sum(~consistent_mask)
    print(f'Found {inconsistent_count} reactions with inconsistent patterns')
    
    k_matrix = k_matrix.loc[consistent_mask]
    
    return k_matrix


def flux_scanning(model, target_id, c_source,c_uptake, alpha, substrate_MW,growth_id,filterG=False,model_type='etfl',tol=0.001,method='ppfba'):
    """
    Args:
        model : ecModel/ETFL model.
        target_id (str): Rxn ID for the production target reaction, a exchange reaction is recommended.
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
    if model_type not in ['etfl', 'ecGEM','GAN_ec']:
        raise ValueError('model_type should be either "etfl" or "ecGEM"')
    if filterG is None:
        filterG = False

    if c_uptake is None:
        c_uptake = 1

    # gluc_MW=0.180156 #g/mmol
    tol_ratio=0.01

    #step 1: build reactions k_matrix
    # Simulate WT (100% max growth):
    FC = {}
    # check if model has transcriptome attribute
    if hasattr(model, 'transcriptome'):
        print('integrating omic data to model')
        from strainOptimizer.manipulation.integration import integrate_omic_data_to_ecmodel
        # import numpy as np
        expression_threshold = np.percentile(model.transcriptome, 25)
        params = {'objective_reaction_id': growth_id, 
                  'obj_frac': 0.4,
                  'expression_threshold':expression_threshold}
        with model:
            model = integrate_omic_data_to_ecmodel(model=model,
                                                    omic_data=model.transcriptome,
                                                    method='GIMME',
                                                    parameters=params)
        wt_sol = ppFBA(model=model,
                                   target_id=growth_id,
                                   c_source=c_source,
                                   c_uptake=c_uptake,
                                   model_type=model_type,
                                   tol_ratio=tol_ratio)

    else:
        wt_sol = ppFBA(model=model,
                                   target_id=growth_id,
                                   c_source=c_source,
                                   c_uptake=c_uptake,
                                   model_type=model_type,
                                   tol_ratio=tol_ratio)
    FC['flux_WT'] = wt_sol.fluxes
    max_growth = FC['flux_WT'][growth_id]
    print('simulate WT-like:',max_growth)

    # simulate production in different suboptimal growth rate conditions
    FC['alpha'] = alpha
    # initialize fluxes and K_scores matrices
    rxnIDlist= [rxn.id for rxn in model.reactions]
    v_matrix = pd.DataFrame(index=rxnIDlist, columns=alpha)
    k_matrix = pd.DataFrame(index=rxnIDlist, columns=alpha)
    # simulate fluxes distribution in different growth rate conditions
    for i in range(len(alpha)):
        biomass_yield=alpha[i]
        growth=biomass_yield*c_uptake*substrate_MW
        growth=round(growth, 4)
        model.reactions.get_by_id(growth_id).bounds= growth*(1-tol_ratio), growth
        # model.reactions.get_by_id(growth_id).bounds = growth, growth
        if method == 'ppfba':
            with model:
                product_sol = ppFBA(model=model,
                                                target_id=target_id,
                                                c_source=c_source,
                                                c_uptake= c_uptake,
                                                model_type=model_type,
                                                tol_ratio=tol_ratio)
        elif method == 'pfba':
            with model:
                product_sol = pFBA(model=model,
                                                target_id=target_id,
                                                c_source=c_source,
                                                c_uptake=c_uptake,
                                                model_type=model_type,
                                                direction='max')
        elif method == 'moma':
            with model:
                if model_type == 'etfl':
                    model.reactions.get_by_id(c_source).bounds = -c_uptake, 0
                elif model_type == 'ecGEM':
                    model.reactions.get_by_id(c_source).bounds = 0, c_uptake
                elif model_type == 'GAN_ec':
                    model.reactions.get_by_id(c_source).bounds = -c_uptake, 0
                # 1.set the target objective
                model.reactions.get_by_id(target_id).bounds = 0, 1000
                model.objective = target_id
                model.objective_direction = 'max'
                product= model.slim_optimize()
                model.reactions.get_by_id(target_id).bounds= product,1000
                product_sol=moma(model=model,
                                 reference_solution=wt_sol,
                                 linear=True,
                                 model_type=model_type)
        elif method == 'mopa':
            with model:
                # open c uptake
                if model_type == 'etfl':
                    model.reactions.get_by_id(c_source).bounds = -c_uptake,0
                elif model_type == 'ecGEM':
                    model.reactions.get_by_id(c_source).bounds = 0, c_uptake
                elif model_type == 'GAN_ec':
                    model.reactions.get_by_id(c_source).bounds = -c_uptake, 0

                # 1.set the target objective
                model.reactions.get_by_id(target_id).bounds = 0, 1000
                model.objective = target_id
                model.objective_direction = 'max'
                product = model.slim_optimize()
                model.reactions.get_by_id(target_id).bounds = product, 1000
                # release growth rate constraint
                # model.reactions.get_by_id(growth_id)
                product_sol=mopa(model=model,
                                 reference_solution=wt_sol,
                                 linear=True,
                                 model_type=model_type,
                                 show=False)

        FC['flux_MUT'] = product_sol.fluxes
        try:
            print('when growth rate is %s,production rate is %s'%(growth,product_sol.fluxes[target_id]))
            # v_matrix.iloc[:, i] = FC['flux_MUT']
            v_matrix[alpha[i]] = FC['flux_MUT']
            k_matrix[alpha[i]] = FC['flux_MUT'] / FC['flux_WT']
        except:
            print('infeasible solution at alpha=%s'%alpha[i])
            # drop this alpha column
            v_matrix = v_matrix.drop(columns=[alpha[i]])
            k_matrix = k_matrix.drop(columns=[alpha[i]])

    # Step 2: Calculate reaction k-scores and apply filtering
    # Filter reactions based on k_matrix patterns
    k_matrix_filtered = k_matrix_filter(model, k_matrix, tol)
    v_matrix_filtered = v_matrix.loc[k_matrix_filtered.index, :]
    
    # Store filtered matrices in FC dictionary
    FC["v_matrix"] = v_matrix_filtered
    FC["k_matrix"] = k_matrix_filtered
    FC["rxns"] = k_matrix_filtered.index.tolist()
    
    # Calculate mean k-score for each reaction across alpha values
    k_rxns = pd.Series(np.mean(k_matrix_filtered.abs(), axis=1), 
                      index=k_matrix_filtered.index)
    
    # Filter out reactions with very low mean flux (< 1e-8)
    flux_threshold = 1e-8
    mean_flux_rxns = pd.Series(np.mean(v_matrix_filtered.abs(), axis=1), 
                              index=v_matrix_filtered.index)
    k_rxns = k_rxns[mean_flux_rxns > flux_threshold]

    FC["k_rxns"] = k_rxns
    # Order from highest to lowest median k_scor
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
                # print('inconsistent k_scores for gene %s'%gene)
                continue
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


def run_ecFSEOF(model, parameters=None,
                target_id=None,
                c_source=None,
                c_uptake=None,
                scanning_range=None,
                Nsteps=None,
                substrate_MW=None,
                growth_id=None,
                model_type='etfl',simulation_method='ppfba',action_thresholds=[0.05,0.3,1.05]):
    """
    Run Flux-scanning with Enforced Objective Function for a specified production target.

    Arguments:
    - model:
        ETFL model
    - target_id: str
        Reaction ID for the production target reaction, a exchange reaction is recommended
    - c_source: str
        Reaction ID for the main carbon source uptake reaction (make sure that the correct directionality is indicated)
    - scanning_range: tuple of float
        Minimum and maximum biomass yield [gDw/mmol Csource] for enforced objective limits
    - Nsteps: int
        Number of steps for suboptimal objective in FSEOF
    - model_type: str
        Type of model (etfl or ecGEM)
    Returns:
    - dict
        Dictionary with the results of the FSEOF analysis
    """
    if parameters is not None:
        target_id = parameters.strain['target_id']
        c_source = parameters.strain['c_source']
        c_uptake = parameters.strain['c_uptake']
        scanning_range = parameters.algorithm['scanning_range']
        Nsteps = parameters.algorithm.get('Nsteps',10)
        substrate_MW = parameters.strain['substrate_MW']
        model_type = parameters.model['model_type']
        simulation_method = parameters.algorithm['simulation_method']
        growth_id = parameters.model['growth_id']
        action_thresholds = parameters.algorithm['action_thresholds']
    else:
        if target_id is None:
            raise ValueError('target_id is required')
        if c_source is None:
            raise ValueError('c_source is required')
        if c_uptake is None:
            raise ValueError('c_uptake is required')
        if scanning_range is None:
            raise ValueError('scanning_range is required')
        if Nsteps is None:
            raise ValueError('Nsteps is required')
        if substrate_MW is None:
            raise ValueError('substrate_MW is required')
        if model_type is None:
            raise ValueError('model_type is required')
        if simulation_method is None:
            raise ValueError('simulation_method is required')

    if scanning_range is None:
        scanning_range = default_scanning_range(model=model,parameters=parameters)
    # Define alpha vector for suboptimal enforced objective values
    scanning_values = np.linspace(scanning_range[0], scanning_range[1], Nsteps)
    # Run FSEOF analysis
    results = flux_scanning(model=model,
                            target_id=target_id,
                            c_source=c_source,
                            c_uptake=c_uptake,
                            alpha=scanning_values,
                            growth_id=growth_id,
                            model_type=model_type,
                            substrate_MW=substrate_MW,
                            method=simulation_method)

    # Create gene table
    geneTable = pd.DataFrame(index=results['genes'], columns=["gene_name", "k_score"], dtype=float)
    geneTable.index.name = "gene_IDs"
    geneTable["gene_name"] = results["geneNames"]
    geneTable["k_score"] = results['k_genes']
    geneTable.loc[geneTable['k_score'] >= action_thresholds[2], 'action'] = 'OE'
    geneTable.loc[geneTable['k_score'] <= action_thresholds[1], 'action'] = 'KD'
    geneTable.loc[geneTable['k_score'] <= action_thresholds[0], 'action'] = 'KO'
    geneTable = geneTable.loc[geneTable['action'].isin(['OE','KD','KO'])]
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


