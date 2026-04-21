# -*- coding: utf-8 -*-
import pandas as pd
from strainOptimizer.strainDesign.ecFactory import ecfseof
from strainOptimizer.strainDesign.ecFactory.ecFactory_other import find_leaks,remove_essential_targets,getMetGeneMatrix,getGeneDepMatrix,getGenesGroups,compare_EUVR,default_scanning_range
from strainOptimizer.manipulation.variable import genelist_to_enzymelist
from strainOptimizer.analysis.enzyme_variety_analysis import enzymeVA
from strainOptimizer.strainDesign.ecFactory import find_min_sets
from strainOptimizer.manipulation.constraint.total_resource_allocation import constrain_enzymes
from strainOptimizer.simulation.pprotFBA import pprotFBA_prot_conc

# def run_ecFactory_design(model, modelParam, expYield,alphaLims,action_thresholds=[0.05,0.5,1.05],remove_essential=False,steps=123):
from strainOptimizer.analysis.FCC import calculate_FCC_by_abundance


def run_ecFactory_design(model, parameters):
    """
    Run ecFactory method to identify gene targets for strain design.
    
    The ecFactory algorithm analyzes enzyme usage variability to identify potential
    genetic engineering targets for metabolic strain optimization. It considers
    experimental yield constraints and applies thresholds to classify targets as
    knockout (KO), knockdown (KD), or overexpression (OE) candidates.
    
    Args:
        model: ETFL/ecGEM metabolic model object
        parameters: WorkflowParameters object containing three parameter dictionaries:
            - model: Model-specific parameters
                * model_type: Type of model ('etfl', 'ecGEM', etc.)
                * growth_id: Growth/biomass reaction ID

            - strain: Strain-specific parameters
                * target_id: Target reaction ID for product synthesis
                * product_name: Product name (used to identify leak reactions)
                * c_source: Carbon source exchange reaction ID
                * c_uptake: Carbon source uptake rate (mmol/gDW/h)
                * substrate_MW: Substrate molecular weight (g/mmol)
                * growth_id: Growth/biomass reaction ID
            - algorithm: Algorithm-specific parameters
                * experimental_yield: Experimental yield of biomass growth
                * simulation_method: Simulation method ('ppfba', 'pfba', 'moma', 'mopa')
                * action_thresholds: List of three thresholds [KO, KD, OE] for gene classification
                * remove_essential: Whether to exclude essential genes from targets (default: False)
                * steps: the ecFactory can be divided into 3 steps (default: 123)
    
    Returns:
        Dict[str, Any]: Dictionary containing design results with keys:
            - geneTable: pandas DataFrame with columns:
                * geneID: Gene identifier
                * k_score: K-score value indicating target priority
                * action: Recommended genetic action ('KO', 'KD', or 'OE')
                * minimal_candidates_set: Boolean indicating if gene is in minimal set
            - Additional analysis results and metadata
    """
    # prepare parameters
    model_type = parameters.model['model_type']
    
    target_id = parameters.strain['target_id']
    productName = parameters.strain['product_name']
    c_source = parameters.strain['c_source']
    substrate_MW = parameters.strain['substrate_MW']
    c_uptake = parameters.strain['c_uptake']

    
    simulation_method = parameters.algorithm['simulation_method']
    action_thresholds=parameters.algorithm['action_thresholds']
    steps=parameters.algorithm.get('steps',123) # default to 123 if not provided
    remove_essential=parameters.algorithm.get('remove_essential',False)
    calculate_fcc=parameters.algorithm.get('calculate_fcc',False)

    growth_id = parameters.model['growth_id']
    if growth_id is None:
        if model_type=='etfl':
            growth_id = model.growth_reaction.id
        elif model_type=='ecGEM':
            growth_id = 'r_2111'
        elif model_type=='GAN_ec':
            growth_id = 'r_2111'
        else:
            raise ValueError('Invalid model type!')

    expYield=parameters.algorithm['experimental_yield']
    scanning_range=parameters.algorithm['scanning_range']
    if scanning_range is None:
        scanning_range = default_scanning_range(model=model,parameters=parameters)
        if expYield is None:
            # if without experimental biomass yield, set as 1/2 of max theoretical yield
            expYield=(scanning_range[0]+scanning_range[1])/2
            expYield=round(expYield,4)

    print(f'Scanning range for FSEOF: {scanning_range}')
    print(f' Experimental biomass yield set to: {expYield} g/gSubstrate')
    step = 0
    # model_tmp=copy.deepcopy(model)

    # 1.- Run FSEOF to find gene candidates
    # Parameters for FSEOF method
    Nsteps = 16  # number of FBA steps in ecFSEOF
    step += 1
    print(f'{step}.-  **** Running ecFSEOF ****')
    results = ecfseof.run_ecFSEOF(model=model,
                              target_id=target_id,
                              c_source=c_source,
                              c_uptake=c_uptake,
                              scanning_range=scanning_range,
                              Nsteps=Nsteps,
                              model_type=model_type,
                                  substrate_MW=substrate_MW,
                                  simulation_method=simulation_method,
                                  growth_id=growth_id,
                                  action_thresholds=action_thresholds)
    # Format results table
    print(f'ecFSEOF returned {len(results["geneTable"])} targets')

    # 2.- Add flux leak targets (those genes not optimal for production that may consume the product of interest.
    # (probaly extend the approach to inmediate precurssors)
    step += 1
    print(f'{step}.-  **** Find flux leak targets to block ****')
    results['geneTable'] = find_leaks(candidates=results['geneTable'],
                                      target_id=target_id,
                                      model=model,
                                      product_name=productName)

    # 3.- discard essential genes from deletion targets
    if remove_essential:
        step += 1
        print(f'{step}.-  **** Removing essential targets ****')
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
    results['level 1 result'] = results['geneTable']

    if steps==1:
        return results

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
    fix_gr=round(fix_gr,2)
    model.reactions.get_by_id(growth_id).bounds = fix_gr, fix_gr
    print(' Fix suboptimal experimental biomass = ' + str(fix_gr) + ' h-1')

    # calculate the protein abundance by parsimonious ptoteins FBA for optimal production condition
    prod_ppFBA_allprotconc = pprotFBA_prot_conc(model=model,
                                                  target_id=target_id,
                                                    c_source=c_source,
                                                  c_uptake=c_uptake,
                                                  model_type=model_type)

    prod_ppFBA_protconc=prod_ppFBA_allprotconc[target_enz_list]

    # fix total enzymes amount for ETFL model
    if model_type=='etfl':
        optm_fraction=1.01
        # calculate the total enzymes amount(exclude dummy enzyme)
        total_enzymes=prod_ppFBA_allprotconc.drop('dummy_enzyme',axis=0).sum()*optm_fraction
        # 保留小数
        total_enzymes=round(total_enzymes, 5)
        print('  - Fix total enzymes amount to %s g/gDW'%total_enzymes)
        model=constrain_enzymes(model,total_enzymes,model_type=model_type)

    # calculate the max and min abundance of each candidate enzymes by enzyme usage variety analysis
    prod_enz_fva_result=enzymeVA(model=model,
                                 target_id=target_id,
                                 enzymeIDlist=target_enz_list,
                                 c_source=c_source,
                                 c_uptake=c_uptake,
                                 fraction_of_optimum=0.99,
                                 obj_direction='max',
                                 model_type=model_type)

    prod_enz_fva_result['minprotFBA']=prod_ppFBA_protconc
    results['prod_enz_fva_result']=prod_enz_fva_result

    # release biomass and production constraints
    model.reactions.get_by_id(growth_id).bounds = 0, 1000

    # 5.2 - Run EUVA for optimal growth condition.
    print(str(step) + '.2-  **** Running EUVA for optimal biomass growth condition ****')
    # fix production rate as 1% of the max production rate
    # model.reactions.get_by_id(target_id).bounds = max_prod*0.01, max_prod*0.01
    print('  - Maximize biomass production')
    
    # set growth as objective
    with model:
        # check if model has transcriptome
        if hasattr(model, 'transcriptome'):
            from strainOptimizer.manipulation.integration import integrate_omic_data_to_ecmodel
            params = {'objective_reaction_id': growth_id, 
                      'obj_frac': 0.4,
                      'expression_threshold':12}
            model=integrate_omic_data_to_ecmodel(model=model,
                                                  omic_data=model.transcriptome,
                                                  method='GIMME',
                                                  parameters=params)
            
        model.objective = growth_id
        model.objective_direction = 'max'
        # run parsimonious protein usages FBA for max growth
        wt_minprotFBA_protconc=pprotFBA_prot_conc(model=model,
                                                target_id=growth_id,
                                                enzymeIDlist=target_enz_list,
                                                c_source=c_source,
                                                c_uptake=c_uptake,
                                                model_type=model_type)
        # run enzyme usage variety analysis
        wt_enz_fva_result=enzymeVA(model=model,
                                target_id=growth_id,
                                enzymeIDlist=target_enz_list,
                                c_source=c_source,
                                c_uptake=c_uptake,
                                fraction_of_optimum=0.99,
                                obj_direction='max',
                                model_type=model_type)

    wt_enz_fva_result['minprotFBA']=wt_minprotFBA_protconc
    results['wt_enz_fva_result']=wt_enz_fva_result     # can be deleted

    # release biomass constraints
    model.reactions.get_by_id(growth_id).bounds = 0, 1000

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
    genetable.loc[leval1_list,'EUVA_priority_level']=1
    genetable.loc[leval2_list,'EUVA_priority_level']=2
    genetable.loc[leval3_list,'EUVA_priority_level']=3
    genetable.loc[no_enzyme_list,'EUVA_priority_level']='no enzyme'
    # mark the no enzyme gene target

    genetable['EUVA_priority_level'] = genetable['EUVA_priority_level'].fillna('discarded by EUVA')

    results['level 2 result'] = genetable[genetable['EUVA_priority_level']!='discarded by EUVA']
    results['geneTable']=genetable
    print('  - Rank targets by priority levels according to EUVA results: ' + str(len(leval1_list)) + ' targets in level 1, ' + str(len(leval2_list)) + ' targets in level 2, ' + str(len(leval3_list)) + ' targets in level 3')


    # # 6.- combine candidate targets
    step=step+1
    tol_tatio=0.001
    print(str(step) + '.    **** Combine candidate targets ****')

    # choose the targets range for combined set calculation
    l1_targets_numb=len(leval1_list)
    l2_targets_numb=len(leval2_list)
    l3_targets_numb=len(leval3_list)
    if l1_targets_numb>30:
        selection_range=[1]
    elif l1_targets_numb+l2_targets_numb>30:
        selection_range=[1,2]
    elif l1_targets_numb+l2_targets_numb+l3_targets_numb>30:
        selection_range=[1,2]
    else:
        selection_range=[1,2,3]

    # get candidate
    candidatesID_list=results['geneTable'][results['geneTable']['EUVA_priority_level'].isin(selection_range)].index.tolist()
    print("Finding the minimal set of %s candidates to achieve the max target production yield" %len(candidatesID_list))
    # find minmal sets of targets
    min_set_analysis_result,optimal_prod_result=find_min_sets.find_min_set(model=model,
                                                                           growth_id=growth_id,
                                                                           c_source=c_source,
                                                                           c_uptake=c_uptake,
                                                                           expYield=expYield,
                                                                            target_id=target_id,
                                                                           geneIDlist=candidatesID_list,
                                                                           gene_enz_fva_result=results['gene_enz_fva_result'],
                                                                            gene_enz_dict=results['gene_enz_dict'],
                                                                           model_type=model_type,
                                                                           c_source_MW=substrate_MW)
    if not min_set_analysis_result.empty and not optimal_prod_result is None:
        print('  - Minimal set of targets: %s'%len(min_set_analysis_result[min_set_analysis_result['score']<(1-tol_tatio)]))
        results['min_set_analysis_result']=min_set_analysis_result
        results['optimal_prod_result']=optimal_prod_result
        min_set_IDlist=min_set_analysis_result[min_set_analysis_result['score']<(1-tol_tatio)].index.tolist()
        results['geneTable'].loc[min_set_IDlist,'minimal_candidates_set']=1
        results['geneTable']['minimal_candidates_set'] = results['geneTable']['minimal_candidates_set'].fillna(0)
        results['level 3 result'] = results['geneTable'][results['geneTable']['minimal_candidates_set']==1]

    if calculate_fcc:
        step += 1
        print(str(step) + '.-  **** Calculating FCC scores for L1 genes ****')

        # resolve NGAM reaction id: try default first, fall back to auto-search
        from strainOptimizer.simulation.utils import find_ngam_reaction
        _DEFAULT_NGAM = 'r_4046'
        try:
            model.reactions.get_by_id(_DEFAULT_NGAM)
            ngam_id = _DEFAULT_NGAM
            print(f'  - NGAM reaction: {ngam_id} (default)')
        except KeyError:
            ngam_candidates = find_ngam_reaction(model)
            if not ngam_candidates:
                raise RuntimeError(
                    'Cannot find NGAM reaction in model. '
                    'Please set objective manually via calculate_FCC_by_abundance.'
                )
            ngam_id = ngam_candidates[0].id
            print(f'  - NGAM reaction not found by default id, auto-detected: {ngam_id}'
                  + (f' (candidates: {[r.id for r in ngam_candidates]})' if len(ngam_candidates) > 1 else ''))

        l1_genes = results['level 1 result'].index.tolist()
        fcc_records = {gene: {'FCCg': None, 'FCCp': None} for gene in l1_genes}
        calculated, skipped = 0, 0
        for gene in l1_genes:
            enz_list = gene_enz_dict.get(gene, [])
            if len(enz_list) != 1:
                skipped += 1
                continue
            protID = enz_list[0]
            try:
                FCCg, FCCp = calculate_FCC_by_abundance(
                    protID=protID,
                    model=model,
                    productID=target_id,
                    c_source=c_source,
                    c_uptake=c_uptake,
                    growthID=growth_id,
                    objective=ngam_id,
                    delta_conc=1
                )
                fcc_records[gene] = {'FCCg': FCCg, 'FCCp': FCCp}
                calculated += 1
            except Exception as e:
                print(f'  - FCC calculation failed for {gene} ({protID}): {e}')
        fcc_df = pd.DataFrame.from_dict(fcc_records, orient='index')
        results['fcc_result'] = fcc_df
        for key in ('geneTable', 'level 1 result', 'level 2 result', 'level 3 result'):
            if key in results and results[key] is not None:
                results[key] = results[key].join(fcc_df, how='left')
        print(f'  - FCC calculated: {calculated} genes, skipped (multi/no enzyme): {skipped}')

    if steps==123:
        return results


