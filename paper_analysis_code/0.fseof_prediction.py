# -*- coding: utf-8 -*-
import cobra
import sys
# add the path of the codebase
sys.path.append('reference/gengrong')
import codebase
from strainOptimizer.io import load_model

productParam_dict={
    '2-phenylethanol':{'productName':'2-phenylethanol',
                       'targetID':'r_1589',
                       'GEM_filepath':'examples/models/yeast/yeast-GEM.xml',
                       'etfl_filepath':'examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json'},
    # 'artemisinic_acid':{'productName':'artemisinic acid',
    #                'targetID':'SK_atemisinic_acid',
    #                'model_filepath':'examples/models/yeast/ecGEM_atemisinic.xml'},
    'heme_acid':{'productName':'heme',
                 'targetID':'EX_heme_a'
                 ,'GEM_filepath':'examples/models/yeast/heme_yeastGEM.xml',
                 'etfl_filepath':'examples/models/yeast/heme_cEFL.json'},
    'spermidine':{'productName':'spermidine',
                  'targetID':'r_2051',
                  'GEM_filepath':'examples/models/yeast/yeast-GEM.xml',
                'etfl_filepath':'examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json'
                  },
    'sclareol':{'productName':'sclareol',
                'targetID':'DM_sclareol_c',
                'GEM_filepath':'examples/models/yeast/sclareol_yeastGEM.xml',
                'etfl_filepath':'examples/models/yeast/sclareol_cEFL.json'},
    'fatty_acid':{'productName':'free fatty acids',
                   'targetID':'r_2189',
                   'GEM_filepath':'examples/models/yeast/yeast-GEM.xml',
                  'etfl_filepath':'examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json'},
}

modelParams_dict={'GEM':{
    'model_type':'GEM',
    'c_source':"r_1714",      # glucose exchange rxn
    'c_uptake':10,
    'growth_id':'r_2111',
        }}


for product_key,productParam in productParam_dict.items():
    productName=productParam['productName']
    productID=productParam['targetID']
    GEM_filepath=productParam['GEM_filepath']
    etfl_filepath=productParam['etfl_filepath']
    for model_key,modelParams in modelParams_dict.items():
        if modelParams['model_type']=='GEM':
            model_filepath=GEM_filepath
        elif modelParams['model_type']=='etfl':
            model_filepath=etfl_filepath

        c_source=modelParams['c_source']
        c_uptake=modelParams['c_uptake']
        growth_id=modelParams['growth_id']

        # load model
        model = load_model(model_filepath, model_type='gem')

        with model:
            model.objective = growth_id  # biomass rxn
            # model.objective=product
            model.reactions.get_by_id(c_source).bounds = -c_uptake, -c_uptake
            gluc_MW = 0.180156  # g/mmol
            max_yield = model.slim_optimize() / (c_uptake * gluc_MW)  # gDW / gGluc
            expYield = max_yield * 0.2
            alphaLims = (0.5 * expYield,
                         2 * expYield)  # The scanning range of the product yield can be adjusted according to the actual situation

        # set the action thresholds
        action_thresholds = [0.05, 0.2, 1]

        # Run the FSEOF algorithm to identify the potential targets
        Nsteps = 6 # number of FBA steps in ecFSEOF
        results = codebase.run_FSEOF(model=model,
                                     targetID=productID,
                                     c_source=c_source,
                                     c_uptake=c_uptake,
                                     alphaLims=alphaLims,
                                     Nsteps=Nsteps,
                                     model_type='GEM',
                                     biomass_id=growth_id)

        # Format results table
        final_result = results['geneTable']
        final_result.loc[final_result['k_score'] >= action_thresholds[2], 'action'] = 'OE'
        final_result.loc[final_result['k_score'] <= action_thresholds[1], 'action'] = 'KD'
        final_result.loc[final_result['k_score'] <= action_thresholds[0], 'action'] = 'KO'
        # remove genes with no action
        # final_result = final_result.loc[final_result['action'].notnull()]
        final_result = final_result.loc[final_result['action'].isin(['OE', 'KD', 'KO'])]

        # save result
        final_result.to_csv(f'analysis_code/results/moma_mopa_pfba/{model_key}_{productName}_gluc_{c_uptake}_fba_ecFSEOF_result.csv')
