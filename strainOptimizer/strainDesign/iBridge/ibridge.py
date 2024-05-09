# -*- coding: utf-8 -*-
# date : 2024/4/24 
# author : wangh
from strainOptimizer.simulation import ppFBA2,moma,fba,mopa,pFBA
import numpy as np
from strainOptimizer.io import load_model
import pandas as pd


def generate_multiple_flux_profiles(model,targetID, c_source, c_uptake,model_type='etfl',method='moma',tol=0.01,linear=False,growthID=None):
    '''Simulate flux distribution at different production levels of targetID
    * Args:
        model: ecModel/ETFL model
        targetID (str): Rxn ID for the production target reaction, a exchange reaction is recommended.
        c_source (str): Rxn ID for the main carbon source uptake reaction
        c_uptake(float): carbon source uptake rate
        tol (float): numerical tolerance for fixing bounds
        model_type (str): 'etfl' or 'ecGEM' or 'GEM' model type
        method (str): Method for model simulation. 'moma' or 'ppfba' or 'mopa'
    * Returns:
        flux_profiles (dict): dictionary of flux distributions for each production level
    '''

    # check model_type
    if model_type not in ['etfl', 'ecGEM','GEM']:
        raise ValueError('model_type should be either "etfl" or "ecGEM"')

    if c_uptake is None:
        c_uptake = 1

    gluc_MW = 0.180156  # g/mmol
    tol_ratio = 0.01

    if growthID is None:
        if model_type == 'etfl':
            growthID = model.growth_reaction.id
        elif model_type == 'ecGEM':
            growthID = 'r_2111'
        elif model_type == 'GEM':
            growthID = 'r_2111'

    # 1.get wild-type flux distribution by minimizing total protein pool
    if model_type == 'GEM':
        if method != 'moma':
            print('Only moma method is supported for GEM model')
            method='moma'
        # get wild-type flux for GEM model by pFBA
        wt_sol=pFBA(model=model,
                targetID=growthID,
                c_source=c_source,
                c_uptake=c_uptake,
                model_type=model_type,
                direction='max')

    else:
        wt_sol=ppFBA2(model=model,
                    targetID=growthID,
                    c_source=c_source,
                    c_uptake=c_uptake,
                    model_type=model_type,
                    tol_ratio=tol_ratio)

    # 2. Get the range of production of targetID
    min_sol=fba(model=model,
                targetID=targetID,
                c_source=c_source,
                c_uptake=c_uptake,
                model_type=model_type,
                direction='min')
    min_prod=min_sol.fluxes[targetID]
    max_sol=fba(model=model,
                targetID=targetID,
                c_source=c_source,
                c_uptake=c_uptake,
                model_type=model_type,
                direction='max')
    max_prod=max_sol.fluxes[targetID]

    # 3. simulate flux distribution for each production level
    flux_profiles={}
    count=0
    for prod in np.linspace(min_prod,max_prod*0.95,10):
        count+=1
        print('Count: %d\tProduction: %0.2f' % (count, prod))
        # constrain the production of targetID
        model.reactions.get_by_id(targetID).bounds = round(prod*(1-tol),2), round(prod*(1+tol),2)

        # run simulation
        if method == 'moma':
            # if model_type == 'GEM':
            #     from cobra.flux_analysis import moma
            #     sol=moma(model=model,
            #             solution=wt_sol,
            #             linear=linear)
            # else:
            #     from strainOptimizer.simulation import moma
            #     sol=moma(model=model,
            #             reference_solution=wt_sol,
            #             linear=linear,
            #             model_type=model_type)
            sol=moma(model=model,
                    reference_solution=wt_sol,
                    linear=linear,
                    model_type=model_type)
        elif method == 'ppfba':
            sol=ppFBA2(model=model,
                    targetID=growthID,
                    c_source=c_source,
                    c_uptake=c_uptake,
                    model_type=model_type,
                    tol_ratio=tol_ratio)
        elif method == 'mopa':
            sol=mopa(model=model,
                    reference_solution=wt_sol,
                    linear=linear,
                    model_type=model_type)
        else:
            raise ValueError('method should be "moma" or "ppfba" or "mopa"')

        flux_profiles[prod]=sol.fluxes

    df_flux=pd.DataFrame.from_dict(flux_profiles)

    return df_flux


def calculate_MetScore_sum(model, covariance_data):
    """
    Calculate the SoC(sum of covariance) for each metabolite based on the provided covariance_data for each reactions.

    Parameters:
    - model: The cobra model containing metabolite and reaction information.
    - covariance_data: A dictionary containing covariance data for each reaction.

    Returns:
    - branch_metabolite_data: A dictionary mapping each metabolite ID to its calculated MetScore.
    """

    branch_metabolite_data = {}

    for each_metabolite in model.metabolites:
        tmp_met_id = each_metabolite.id
        branch_metabolite_data[tmp_met_id] = 0.0
        for each_reaction in each_metabolite.reactions:
            reactants = [met.id for met in each_reaction.reactants]
            products = [met.id for met in each_reaction.products]
            if each_reaction.id in covariance_data:
                if (tmp_met_id in reactants and covariance_data[each_reaction.id] >= 0):
                    branch_metabolite_data[tmp_met_id] += covariance_data[each_reaction.id]
                elif tmp_met_id in reactants and covariance_data[each_reaction.id] < 0:
                    branch_metabolite_data[tmp_met_id] += covariance_data[each_reaction.id]
                elif tmp_met_id in products and covariance_data[each_reaction.id] >= 0:
                    branch_metabolite_data[tmp_met_id] -= covariance_data[each_reaction.id]
                else:
                    branch_metabolite_data[tmp_met_id] -= covariance_data[each_reaction.id]

    return branch_metabolite_data


def calculate_rxn_correlation(df_flux,filter=False):
    # ignore reaction include both positive and negative flux
    if filter:
        df_flux = df_flux[(df_flux >= 0).all(1) | (df_flux <= 0).all(1)]

    df_flux_corr = df_flux.abs().T.corr()
    df_flux_cov = df_flux.abs().T.cov()

    return df_flux_corr,df_flux_cov


def select_candidates(model, targetID,metscore_df,
                      corr_df, cov_df,corr_threshold=0, cov_threshold=0.1):
    """
    Selects candidate reactions based on correlation and covariance thresholds, and calculate the positive and negative
    scores for each metabolite.

    Args:
        model: The model used for selection.
        target_reaction: The target reaction to be considered.
        output_file: The file to write the output to.
        metscore_df: The dataframe containing metabolite scores.
        corr_df: The dataframe containing correlation scores.
        cov_df: The dataframe containing covariance scores.
        corr_threshold: The correlation threshold for selection (default 0).
        cov_threshold: The civariance threshold for selection (default 0.1).

    Returns:
        df_candidates: The dataframe containing the selected candidate reactions.
    """
    # Filter the positive candidate reactions
    pcorr_df = corr_df[corr_df > corr_threshold].dropna()
    pcov_df = cov_df[cov_df > cov_threshold].dropna()
    positive_candidate_reactions = list(set(pcorr_df.index) & set(pcov_df.index))

    # Filter the negative candidate reactions
    ncorr_df = corr_df[corr_df < -corr_threshold].dropna()
    ncov_df = cov_df[cov_df < -cov_threshold].dropna()
    negative_candidate_reactions = list(set(ncorr_df.index) & set(ncov_df.index))

    print('Number of positive candidate reactions: %d' % (len(positive_candidate_reactions)))
    print('Number of negative candidate reactions: %d' % (len(negative_candidate_reactions)))

    # fp = open(output_file, 'w')
    col = ['Metabolite', 'Score',
              'No. of reactions', 'No. of positive reactions',
              'No. of negative reactions', 'candidate reactions',
              'positive reactions', 'negative reactions',
              'Positive score', 'Negative score']
    # fp.write('%s\n' % ('\t'.join(header)))
    df_candidates = pd.DataFrame(columns=col)

    for each_row, each_df in metscore_df.iterrows():
        cobra_metabolite = model.metabolites.get_by_id(each_row)
        candidate_reactions = []

        for each_reaction in cobra_metabolite.reactions:
            for each_reactant in each_reaction.reactants:
                if each_reactant.id == each_row:
                    candidate_reactions.append(each_reaction.id)

        candidate_reactions = list(set(candidate_reactions))
        pos_candidate_reactions = list(set(positive_candidate_reactions) & set(candidate_reactions))
        neg_candidate_reactions = list(set(negative_candidate_reactions) & set(candidate_reactions))

        positive_score = 0.0
        negative_score = 0.0
        for rxn in pos_candidate_reactions:
            positive_score += float(pcov_df.loc[rxn])
        for rxn in neg_candidate_reactions:
            negative_score += float(ncov_df.loc[rxn])

        df_candidates.loc[each_row] = [each_row, each_df[targetID],
                    len(candidate_reactions), len(pos_candidate_reactions),
                    len(neg_candidate_reactions), ';'.join(candidate_reactions),
                    ';'.join(pos_candidate_reactions), ';'.join(neg_candidate_reactions),
                    positive_score, negative_score]

    return df_candidates


def make_candidate_reaction_sets(df_candidate):
    '''Define all possible bridge reactions(Consume negative metabolite and produce positive metabolite)'''
    score_info = {}
    for each_met, each_df in df_candidate.groupby('Metabolite'):

        # if each_met[-2] == 'c':
            pos_reaction_num = each_df['No. of positive reactions'].values[0]
            neg_reaction_num = each_df['No. of negative reactions'].values[0]

            pos_score = each_df['Positive score'].values[0]
            neg_scroe = each_df['Negative score'].values[0]

            # normalize the score by the number of reactions to avoid reaction number bias
            final_pos_scroe = pos_score / np.sqrt(1 + pos_reaction_num)
            final_neg_scroe = neg_scroe / np.sqrt(1 + neg_reaction_num)

            # if both are 0, then ignore it
            if final_pos_scroe == 0 and final_neg_scroe == 0:
                continue

            score_info[each_met] = [final_pos_scroe, final_neg_scroe]

    negative_score_info = {}
    positive_score_info = {}

    for met in score_info:
        if abs(score_info[met][0]) > abs(score_info[met][1]):
            positive_score_info[met] = score_info[met][0]
        else:
            negative_score_info[met] = score_info[met][1]

    df_bridge = pd.DataFrame(columns=['Negative metabolite', 'Positive metabolite', 'Negative score', 'Positive score'])
    for negative_met in negative_score_info:
        negative_score = negative_score_info[negative_met]

        for positive_met in positive_score_info:
            positive_score = positive_score_info[positive_met]
            # add a new row
            df_bridge = df_bridge.append({'Negative metabolite': negative_met,
                                          'Positive metabolite': positive_met,
                                          'Negative score': negative_score,
                                          'Positive score': positive_score}, ignore_index=True)

    return df_bridge


def metabolite_set(cobra_model):
    met_info = []

    for each_reaction in cobra_model.reactions:
        reactants = [met.id for met in each_reaction.reactants]
        products = [met.id for met in each_reaction.products]
        met_info.append([reactants, products, each_reaction.id, each_reaction.reaction])
    return met_info


def get_endogenous_reaction_result(df_bridge_candidate, met_info, flux_corr_df, flux_cov_df,filter=False):
    rxn_result = dict()

    for each_row, each_df in df_bridge_candidate.iterrows():
        negative_met = each_df['Negative metabolite']
        positive_met = each_df['Positive metabolite']

        negative_score = each_df['Negative score']
        positive_score = each_df['Positive score']

        for each_met_set in met_info:
            # ignore arm reaction which contains 'arm' in rxnID
            if 'arm' in each_met_set[2]:
                continue
            # if negative_met in reactants and positive_met in products，then it is a OE targets
            if negative_met in each_met_set[0] and positive_met in each_met_set[1]:
                # up_reaction_list.append(each_met_set[2])
                rxn_score=positive_score-negative_score
                target_reaction = each_met_set[2]

                if target_reaction not in flux_corr_df.index:
                    corr_val = 0
                    cov_val = 0
                else:
                    corr_val = float(flux_corr_df.loc[target_reaction])
                    cov_val = float(flux_cov_df.loc[target_reaction])

                rxn_result[target_reaction] = [rxn_score,corr_val,cov_val]

            elif negative_met in each_met_set[1] and positive_met in each_met_set[0]:
                rxn_score=negative_score-positive_score
                target_reaction = each_met_set[2]
                if target_reaction not in flux_corr_df.index:
                    corr_val = 0
                    cov_val = 0
                else:
                    corr_val = float(flux_corr_df.loc[target_reaction])
                    cov_val = float(flux_cov_df.loc[target_reaction])

                rxn_result[target_reaction] = [rxn_score,corr_val,cov_val]

    col=['rxn_score','covariance','correlation']
    df_rxn_result=pd.DataFrame.from_dict(rxn_result,orient='index',columns=col)
    df_rxn_result['action']=['OE' if i>0 else 'KD' for i in df_rxn_result['rxn_score']]

    if filter:
        # fill nan with 0
        df_rxn_result['correlation']=df_rxn_result['correlation'].fillna(0)
        df_rxn_result['covariance']=df_rxn_result['covariance'].fillna(0)
        # remove OE reactions with negative correlation or negative covariance
        df_rxn_result=df_rxn_result[(df_rxn_result['action']=='KD') | (df_rxn_result['correlation']>0) | (df_rxn_result['covariance']>0)]
        # remove KD reactions with positive correlation or positive covariance
        df_rxn_result=df_rxn_result[(df_rxn_result['action']=='OE') | (df_rxn_result['correlation']<0) | (df_rxn_result['covariance']<0)]

    return df_rxn_result


def get_gene_result_from_rxn(rxn_result,model):

    gene_result=dict()
    # remove conflict genes that both have positive rxn and negative rxn
    conflict_genes=list()
    for rxnID in rxn_result.index:
        rxn=model.reactions.get_by_id(rxnID)
        gene_score=rxn_result.loc[rxnID,'rxn_score']
        gene_covariance=rxn_result.loc[rxnID,'covariance']
        gene_correlation=rxn_result.loc[rxnID,'correlation']
        gene_action=rxn_result.loc[rxnID,'action']
        for gene in rxn.genes:
            if gene.id not in gene_result:
                gene_result[gene.id]=[gene_score,gene_covariance,gene_correlation,gene_action]
            else:
                if gene_action!=gene_result[gene.id][3]:
                    conflict_genes.append(gene.id)
                else:
                    gene_score=gene_score+gene_result[gene.id][0]
                    if gene_action=='OE':
                        gene_covariance=max(gene_covariance,gene_result[gene.id][1])
                        gene_correlation=max(gene_correlation,gene_result[gene.id][2])
                    elif gene_action=='KD':
                        gene_covariance=min(gene_covariance,gene_result[gene.id][1])
                        gene_correlation=min(gene_correlation,gene_result[gene.id][2])

                    gene_result[gene.id]=[gene_score,gene_covariance,gene_correlation,gene_action]

    conflict_genes=list(set(conflict_genes))
    print(f'{len(conflict_genes)} conflict genes are removed')
    for gene in conflict_genes:
        gene_result.pop(gene)

    df_gene_result=pd.DataFrame.from_dict(gene_result,orient='index',columns=['gene_score','gene_covariance','gene_correlation','gene_action'])

    return df_gene_result


def run_iBridge_design(model,targetID,c_source,c_uptake,model_type,method,tol,linear,corr_threshold=0,cov_threshold=0.1):
    '''
    This function run iBridge method to identify gene/reaction targets for strain design.
    * Args:
        model: ETFL/ecGEM/GEM model
        targetID: the target product ID
        c_source: the carbon source uptake reaction ID
        c_uptake: the carbon source uptake amount
        model_type: the model type, 'etfl' or 'ecGEM' or 'GEM'
        method: the method for simulation, 'moma' or 'ppfba' or 'mopa'
        tol: the tolerance for simulation
        linear: whether the model is linear
        corr_threshold: the correlation threshold
        cov_threshold: the covariance threshold
    * Return:
        result (dict): the result of iBridge with the following keys:
            'bridge_candidates': a dataframe of bridge candidates
            'endogenous_reaction_result': a dataframe of endogenous reaction results
            'endogenous_gene_result': a dataframe of endogenous gene results

    ref:https://github.com/kaistsystemsbiology/iBridge

    '''
    print('Running iBridge for product:',targetID)
    # 1. generate multiple flux profiles
    print('step 1: generate multiple flux profiles...')
    df_flux = generate_multiple_flux_profiles(model=model,
                                             targetID=targetID,
                                             c_source=c_source,
                                             c_uptake=c_uptake,
                                             model_type=model_type,
                                             method=method,
                                             tol=tol,
                                             linear=linear)

    # 2. calculate the correlation and covariance of each reaction to the target reaction
    print('step 2: calculate the correlation and covariance of each reaction to the target reaction...')
    # calculate the correlation and covariance of each reaction to the target reaction
    df_flux_corr,df_flux_cov=calculate_rxn_correlation(df_flux,filter=True)
    df_target_corr=df_flux_corr[targetID]
    df_target_cov=df_flux_cov[targetID]
    # fill nan as 0
    df_target_corr=df_target_corr.fillna(0)
    df_target_cov=df_target_cov.fillna(0)

    # 3. calculate the MetScore for each metabolite
    print('step 3: calculate the MetScore for each metabolite...')
    fluxsum_dic = calculate_MetScore_sum(model,  dict(df_target_cov))
    df_metscore= pd.DataFrame.from_dict({targetID:fluxsum_dic})

    # 4. calculate all candidate reactions
    print('step 4: calculate all candidate bridge reactions...')
    df_candidates = select_candidates(model=model,
                                      targetID=targetID,
                                      metscore_df=df_metscore,
                                      corr_df=df_target_corr,
                                      cov_df=df_target_cov,
                                      corr_threshold=corr_threshold,
                                      cov_threshold=cov_threshold)

    df_bridge_candidates = make_candidate_reaction_sets(df_candidates)

    # 5. get final endogenous bridge reaction/gene result
    print('step 5: get final endogenous bridge reaction/gene result...')
    met_info=metabolite_set(model)
    df_endogenous_reaction_result=get_endogenous_reaction_result(df_bridge_candidate=df_bridge_candidates,
                                                                 met_info=met_info,
                                                                 flux_corr_df=df_target_corr,
                                                                 flux_cov_df=df_target_cov,
                                                                 filter=True)

    # get gene result from reaction result
    df_endogenous_gene_result=get_gene_result_from_rxn(rxn_result=df_endogenous_reaction_result,
                                                       model=model)


    if linear:
        method=method+'_linear'
    else:
        method=method+'_nonlinear'
    results={
        'bridge_candidates':df_bridge_candidates,
        'endogenous_reaction_result':df_endogenous_reaction_result,
        'endogenous_gene_result':df_endogenous_gene_result,
        'model_type':model_type,
        'method':method,
        'product':targetID,
        'c_source':c_source,
        'c_uptake':c_uptake,
        'flux_profiles':df_flux,
        'target_corr':df_target_corr,
        'target_cov':df_target_cov,

    }

    return results
