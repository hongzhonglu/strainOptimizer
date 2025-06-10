# -*- coding: utf-8 -*-
import sys
# sys.path.append(r"D:\code\github\etfl\code_etfl\strainOptimizer\ecFactory")
import pandas as pd
from strainOptimizer.strainDesign.ecFactory import ecfseof
from strainOptimizer.strainDesign.ecFactory.ecFactory_other import find_leaks,remove_essential_targets,getMetGeneMatrix,getGeneDepMatrix,getGenesGroups,genelist_to_enzymelist,compare_EUVR
from strainOptimizer.analysis.enzyme_variety_analysis import enzymeVA
from strainOptimizer.strainDesign.ecFactory import find_min_sets
from strainOptimizer.manipulation.constraint.total_resource_allocation import constrain_enzymes
from strainOptimizer.simulation.pprotFBA import pprotFBA_prot_conc

def run_ecFactory_design(model, modelParam, expYield,alphaLims,action_thresholds=[0.05,0.5,1.05],remove_essential=False,steps=123):
    '''
    This function runs ecFactory method to identify gene targets for strain design
    * Args:
    model: ETFL/ecGEM model
    modelParam: a pandas series with the following parameters:
        targetID: target reaction ID
        productName: product name(used to find the leak reactions)
        c_source: carbon source exchange reaction ID
        c_uptake: carbon source uptake rate
    expYield: experimental yield of the biomass growth
    action_thresholds: a list of three thresholds for gene targets:
        1. KO threshold
        2. KD threshold
        3. OE threshold
    remove_essential: a boolean value indicating whether to remove essential genes from the list of targets
    model_type: a string indicating the type of model ('etfl' or 'ecGEM')

    * Return:
    results: a dictionary with the following keys:
        geneTable: a pandas dataframe with the following columns:
            geneID: gene ID
            k_score: k-score value
            action: gene action
            minimal candidates set: a boolean value indicating whether the gene is in the minimal candidates set
    '''
    # model parameters
    targetID = modelParam['targetID']
    productName = modelParam['productName']
    c_source = modelParam['c_source']
    c_uptake = modelParam['c_uptake']
    model_type = modelParam['model_type']
    try:
        simulation_method = modelParam['simulation_method']
    except:
        simulation_method = 'ppfba'  # default simulation method is ppfba

    if model_type=='etfl':
        growth_rxnID = model.growth_reaction.id
    elif model_type=='ecGEM':
        growth_rxnID = 'r_2111'

    try:
        substrate_MW=modelParam['substrate_MW']  # g/mmol
    except:
        # if no substrate molecular weight, use glucose as default
        substrate_MW=0.180156

    step = 0

    # 1.- Run FSEOF to find gene candidates
    # Parameters for FSEOF method
    Nsteps = 10  # number of FBA steps in ecFSEOF
    step = step + 1
    print(f'{step}.-  **** Running ecFSEOF method (ref: GECKO utilities) ****')
    results = ecfseof.run_ecFSEOF(model=model,
                              targetID=targetID,
                              c_source=c_source,
                              c_uptake=c_uptake,
                              alphaLims=alphaLims,
                              Nsteps=Nsteps,
                              model_type=model_type,
                                  substrate_MW=substrate_MW,
                                  simulation_method=simulation_method)
    # Format results table
    gene_result=results['geneTable']
    gene_result.loc[gene_result['k_score'] >= action_thresholds[2], 'action'] = 'OE'
    gene_result.loc[gene_result['k_score'] <= action_thresholds[1], 'action'] = 'KD'
    gene_result.loc[gene_result['k_score'] <= action_thresholds[0], 'action'] = 'KO'
    # remove genes with no action
    # gene_result = gene_result.loc[gene_result['action'].notnull()]
    gene_result = gene_result.loc[gene_result['action'].isin(['OE','KD','KO'])]
    print(f'ecFSEOF returned {len(gene_result)} targets\n')
    results['geneTable'] = gene_result

    if steps==1:
        return results

    # 2.- Add flux leak targets (those genes not optimal for production that may consume the product of interest.
    # (probaly extend the approach to inmediate precurssors)
    step += 1
    print(f'{step}.-  **** Find flux leak targets to block ****')
    results['geneTable'] = find_leaks(candidates=results['geneTable'],
                                      targetID=targetID,
                                      model=model,
                                      product_name=productName)

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
    metGeneMatrix, metsConectivity, genesConectivity=getMetGeneMatrix(model=model,
                                                                      geneIDlist=results['geneTable'].index.tolist(),
                                                                      model_type=model_type)
    # Get independent genes from GeneMetMatrix
    independant_genes,gene_equal_Matrix=getGeneDepMatrix(metGeneMatrix)
    # Get gene target groups (those connected to exactly the same metabolites)
    groups = getGenesGroups(gene_equal_Matrix)
    print('independent targets: ' + str(len(independant_genes[independant_genes==1])))
    print('target groups: ' + str(len(groups)))
    results['independent_genes'] = independant_genes
    results['groups'] = groups
    results['metGeneMatrix'] = metGeneMatrix


    # 5.- enzyme usage variety analysis(EUVA)
    step += 1
    print(str(step) + '.-  **** Running EUVA analysis ****')
    # get target enzyme list
    gene_targetList=results['geneTable'].index.tolist()
    target_enz_list,gene_enz_dict=genelist_to_enzymelist(model=model,
                                                         genelist=gene_targetList,
                                                         model_type=model_type)
    df_gene_to_enz=pd.Series(gene_enz_dict)
    results['gene_enz_dict']=df_gene_to_enz

    # 5.1 - Run EUVA for optimal production conditions
    print(str(step) + '.1-  **** Running EUVA for optimal production conditions ****')
    # Fix suboptimal experimental biomass yield conditions
    fix_gr = expYield * substrate_MW* c_uptake
    model.reactions.get_by_id(growth_rxnID).bounds = fix_gr, fix_gr
    print(' Fix suboptimal experimental biomass = ' + str(fix_gr) + ' h-1')

    # calculate the protein abundance by parsimonious ptoteins FBA for optimal production condition
    prod_ppFBA_allprotconc = pprotFBA_prot_conc(model=model,
                                                  targetID=targetID,
                                                    c_source=c_source,
                                                  c_uptake=c_uptake,
                                                  model_type=model_type)

    prod_ppFBA_protconc=prod_ppFBA_allprotconc[target_enz_list]

    # fix total enzymes amount for ETFL model
    if model_type=='etfl':
        optm_fraction=1.01
        # calculate the total enzymes amount(exclude dummy enzyme)
        total_enzymes=prod_ppFBA_allprotconc.drop('dummy_enzyme',axis=0).sum()*optm_fraction
        print('  - Fix total enzymes amount to %s g/gDW'%total_enzymes)
        model=constrain_enzymes(model,total_enzymes,model_type=model_type)

    # calculate the max and min abundance of each candidate enzymes by enzyme usage variety analysis
    prod_enz_fva_result=enzymeVA(model=model,
                                 targetID=targetID,
                                 enzymeIDlist=target_enz_list,
                                 c_source=c_source,
                                 c_uptake=c_uptake,
                                 fraction_of_optimum=0.99,
                                 obj_direction='max',
                                 model_type=model_type)

    prod_enz_fva_result['minprotFBA']=prod_ppFBA_protconc
    results['prod_enz_fva_result']=prod_enz_fva_result

    # release biomass and production constraints
    model.reactions.get_by_id(growth_rxnID).bounds = 0, 1000

    # 5.2 - Run EUVA for optimal growth condition.
    print(str(step) + '.2-  **** Running EUVA for optimal biomass growth condition ****')
    # fix production rate as 1% of the max production rate
    # model.reactions.get_by_id(targetID).bounds = max_prod*0.01, max_prod*0.01
    print('  - Maximize biomass production')
    
    # set growth as objective
    with model:
        # check if model has transcriptome
        if hasattr(model, 'transcriptome'):
            from strainOptimizer.manipulation.integration import integrate_omic_data_to_ecmodel
            params = {'objective_reaction_id': growth_rxnID, 
                      'obj_frac': 0.4,
                      'expression_threshold':12}
            model=integrate_omic_data_to_ecmodel(model=model,
                                                  omic_data=model.transcriptome,
                                                  method='GIMME',
                                                  parameters=params)
            
        model.objective = growth_rxnID
        model.objective_direction = 'max'
        # run parsimonious protein usages FBA for max growth
        wt_minprotFBA_protconc=pprotFBA_prot_conc(model=model,
                                                targetID=growth_rxnID,
                                                enzymeIDlist=target_enz_list,
                                                c_source=c_source,
                                                c_uptake=c_uptake,
                                                model_type=model_type)
        # run enzyme usage variety analysis
        wt_enz_fva_result=enzymeVA(model=model,
                                targetID=growth_rxnID,
                                enzymeIDlist=target_enz_list,
                                c_source=c_source,
                                c_uptake=c_uptake,
                                fraction_of_optimum=0.99,
                                obj_direction='max',
                                model_type=model_type)

    wt_enz_fva_result['minprotFBA']=wt_minprotFBA_protconc
    results['wt_enz_fva_result']=wt_enz_fva_result     # can be deleted

    # release biomass constraints
    model.reactions.get_by_id(growth_rxnID).bounds = 0, 1000

    # 5.3 - Discard some targets according to the EUVA results
    print(str(step) + '.3-  **** Discarding targets according to EUVA results ****')
    # build a dataframe with all candidate genes  EUVA results
    gene_enz_fva_result=pd.DataFrame(index=results['geneTable'].index.tolist(),columns=['prod_min','prod_max','prod_minprot','wt_min','wt_max','wt_minprot'])
    for gene in results['geneTable'].index.tolist():
        if len(gene_enz_dict[gene])==1:
            enzID=gene_enz_dict[gene][0]
            enz_fva_result=list(prod_enz_fva_result.loc[enzID].values)+list(wt_enz_fva_result.loc[enzID].values)
            gene_enz_fva_result.loc[gene]=enz_fva_result
    results['gene_enz_fva_result']=gene_enz_fva_result

    # 5.3.1 Discard OE target candidates classified as futile for production(min=parsimonious=0)
    oe_candidates=results['geneTable'][results['geneTable']['action']=='OE'].index.tolist()
    futile_candidates=gene_enz_fva_result[(gene_enz_fva_result['prod_min']==0) & (gene_enz_fva_result['prod_minprot']==0)].index.tolist()
    futile_oe_toremove=list(set(oe_candidates).intersection(set(futile_candidates)))
    results['gene_enz_fva_result']=results['gene_enz_fva_result'].drop(futile_oe_toremove)
    print('  - Discard OE targets with min=parsimonious=0: ' + str(len(futile_oe_toremove)) + ' targets removed')

    # 5.3.2 Discard enzymes essential for production from KO target candidates(lb>0).
    ko_candidates=results['geneTable'][results['geneTable']['action']=='KO'].index.tolist()
    essential_candidates=gene_enz_fva_result[gene_enz_fva_result['prod_min']>0].index.tolist()
    essential_ko_toremove=list(set(ko_candidates).intersection(set(essential_candidates)))
    results['gene_enz_fva_result']=results['gene_enz_fva_result'].drop(essential_ko_toremove)
    print('  - Discard KO targets with min>0: ' + str(len(essential_ko_toremove)) + ' targets removed')

    # 5.3.3 Discard isoenzyme groups that contain an optimal isoform for biomass formation from KD and KO target candidates.
    kd_ko_candidates=results['geneTable'][results['geneTable']['action'].isin(['KD','KO'])].index.tolist()
    iso_group=[]
    for g_group in results['groups']:
        group_enz_fva_result=gene_enz_fva_result.loc[g_group]
        # check if there exist both >0 and =0 in prod_minprot
        if group_enz_fva_result['prod_minprot'].max()>0 and group_enz_fva_result['prod_minprot'].min()==0:
            iso_group=iso_group+g_group
    iso_group_toremove=list(set(kd_ko_candidates).intersection(set(iso_group)))
    # remove repeated genes in essential_ko_toremove
    iso_group_toremove=list(set(iso_group_toremove).difference(set(essential_ko_toremove)))
    results['gene_enz_fva_result']=results['gene_enz_fva_result'].drop(iso_group_toremove)
    print('  - Discard KD and KO targets that belong to isoenzyme groups with optimal isoform for biomass formation: ' + str(len(iso_group_toremove)) + ' targets removed')

    # 5.3.4  Discard enzymes with inconsistent result between EUVA and fseof
    # remove discarded genes according to EUVA results
    # results['gene_enz_fva_result']=results['gene_enz_fva_result'].loc[results['geneTable'].index.tolist(),:]
    gene_euvr_compare=compare_EUVR(gene_enz_fva_result=results['gene_enz_fva_result'])
    euva_up_List=gene_euvr_compare[gene_euvr_compare.str.contains('up_')].index.tolist()
    euva_down_List=gene_euvr_compare[gene_euvr_compare.str.contains('down_')].index.tolist()
    fseof_up_List=results['geneTable'][results['geneTable']['action']=='OE'].index.tolist()
    fseof_down_List=results['geneTable'][results['geneTable']['action'].isin(['KD','KO'])].index.tolist()
    # remove genes both in euva_up_List&fseof_down_List, euva_down_List&fseof_up_List
    euva_up_fseof_down=list(set(euva_up_List).intersection(set(fseof_down_List)))
    euva_down_fseof_up=list(set(euva_down_List).intersection(set(fseof_up_List)))
    euva_inconsistent_toremove=euva_up_fseof_down+euva_down_fseof_up
    results['gene_enz_fva_result']=results['gene_enz_fva_result'].drop(euva_inconsistent_toremove)
    gene_euvr_compare=gene_euvr_compare.drop(euva_inconsistent_toremove)
    print('  - Discard targets with inconsistent result between EUVA and fseof: ' + str(len(euva_inconsistent_toremove)) + ' targets removed')
    results['gene_euvr_compare']=gene_euvr_compare

    # 5.4 rank the remaining targets according to EUVA results
    print(str(step) + '.4-  ****  **** Rank targets by priority levels according to EUVA results ****')
    leval1_list=gene_euvr_compare[gene_euvr_compare.str.contains('distinct')].index.tolist()
    leval2_list=gene_euvr_compare[gene_euvr_compare.str.contains('overlaped')].index.tolist()
    leval3_list=gene_euvr_compare[gene_euvr_compare=='undistinguishable'].index.tolist()
    no_enzyme_list=results['gene_enz_dict'][results['gene_enz_dict']=='no enzyme'].index.tolist()
    leval3_list=list(set(leval3_list).difference(set(no_enzyme_list)))
    genetable=results['geneTable']
    # genetable['target_priority_leval']=0
    genetable.loc[leval1_list,'target_priority_leval']=1
    genetable.loc[leval2_list,'target_priority_leval']=2
    genetable.loc[leval3_list,'target_priority_leval']=3
    genetable.loc[no_enzyme_list,'target_priority_leval']='no enzyme'
    # mark the no enzyme gene target

    genetable['target_priority_leval'].fillna('removed by EUVA',inplace=True)
    results['geneTable']=genetable
    print('  - Rank targets by priority levels according to EUVA results: ' + str(len(leval1_list)) + ' targets in level 1, ' + str(len(leval2_list)) + ' targets in level 2, ' + str(len(leval3_list)) + ' targets in level 3')


    # return results
    if steps==12:
        return results

    # # 6.- combine candidate targets
    step=step+1
    tol_tatio=0.001
    print(str(step) + '.    **** Combine candidate targets ****')
    # # fix substrate uptake rate
    c_source=modelParam['c_source']
    c_uptake=modelParam['c_uptake']
    targetID=modelParam['targetID']

    # choose the targets range for combined set calculation
    l1_targets_numb=len(leval1_list)
    l2_targets_numb=len(leval2_list)
    l3_targets_numb=len(leval3_list)
    if l1_targets_numb>30:
        selection_range=[1]
    elif l1_targets_numb+l2_targets_numb>30:
        selection_range=[1,2]
    else:
        selection_range=[1,2,3]

    # get candidate
    candidatesID_list=results['geneTable'][results['geneTable']['target_priority_leval'].isin(selection_range)].index.tolist()
    print("Finding the minimal set of %s candidates to achieve the max target production yield" %len(candidatesID_list))
    # find minmal sets of targets
    min_set_analysis_result,optimal_prod_result=find_min_sets.find_min_set(model=model,
                                                                           c_source=c_source,
                                                                           c_uptake=c_uptake,
                                                                           expYield=expYield,
                                                                            targetID=targetID,
                                                                           geneIDlist=candidatesID_list,
                                                                           gene_enz_fva_result=results['gene_enz_fva_result'],
                                                                            gene_enz_dict=results['gene_enz_dict'],
                                                                           model_type=model_type,
                                                                           c_source_MW=substrate_MW)
    if not min_set_analysis_result.empty:
        print('  - Minimal set of targets: %s'%len(min_set_analysis_result[min_set_analysis_result['score']<(1-tol_tatio)]))
        results['min_set_analysis_result']=min_set_analysis_result
        results['optimal_prod_result']=optimal_prod_result
        min_set_IDlist=min_set_analysis_result[min_set_analysis_result['score']<(1-tol_tatio)].index.tolist()
        results['geneTable'].loc[min_set_IDlist,'minimal candidates set']=1
        results['geneTable']['minimal candidates set'].fillna(0,inplace=True)


    return results


