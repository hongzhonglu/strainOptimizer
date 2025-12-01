# -*- coding: utf-8 -*-
from strainOptimizer.io import load_model
from strainOptimizer.manipulation.constraint.total_resource_allocation import constrain_enzymes
from strainOptimizer.strainDesign.ecFactory import run_ecFactory
import pandas as pd
import numpy as np


modelParams_dict={
    # 'ecGEM':{
    # 'model_type':'ecGEM',
    # 'c_source':"r_1714_REV",      # glucose exchange rxn
    # 'c_uptake':1,
    # 'growth_id':'r_2111',
    # 'total_enzymes':0.1,
    #     },
'etfl':{
    'model_type':'etfl',
    'c_source':"r_1714",
    'c_uptake':1,
    'growth_id':'r_2111',
    'total_enzymes':0.1,}
}

productParam_dict={
    '2-phenylethanol':{'productName':'2-phenylethanol',
                       'targetID':'r_1589',
                       'ecGEM_filepath':'examples/models/yeast/ecYeastGEM_batch.xml',
                       'etfl_filepath':'examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json'},
    # 'artemisinic_acid':{'productName':'artemisinic acid',
    #                'targetID':'SK_atemisinic_acid',
    #                'model_filepath':'examples/models/yeast/ecGEM_atemisinic.xml'},
    'heme_acid':{'productName':'heme',
                 'targetID':'EX_heme_a'
                 ,'ecGEM_filepath':'examples/models/yeast/heme_ecYeastGEM.xml',
                 'etfl_filepath':'examples/models/yeast/heme_cEFL.json'},
    'spermidine':{'productName':'spermidine',
                  'targetID':'r_2051',
                  'ecGEM_filepath':'examples/models/yeast/ecYeastGEM_batch.xml',
                'etfl_filepath':'examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json'
                  },
    'sclareol':{'productName':'sclareol',
                'targetID':'DM_sclareol_c',
                'ecGEM_filepath':'examples/models/yeast/sclareol_ecYeastGEM_batch.xml',
                'etfl_filepath':'examples/models/yeast/sclareol_cEFL.json'},
    'fatty_acid':{'productName':'free fatty acids',
                   'targetID':'r_2189',
                   'ecGEM_filepath':'examples/models/yeast/ecYeastGEM_batch.xml',
                  'etfl_filepath':'examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json'},
}

for product_key,productParam in productParam_dict.items():
    productName=productParam['productName']
    productID=productParam['targetID']
    ecGEM_filepath=productParam['ecGEM_filepath']
    etfl_filepath=productParam['etfl_filepath']
    for model_key,modelParams in modelParams_dict.items():
        if modelParams['model_type']=='ecGEM':
            model_filepath=ecGEM_filepath
        elif modelParams['model_type']=='etfl':
            model_filepath=etfl_filepath

        c_source=modelParams['c_source']
        c_uptake=modelParams['c_uptake']
        growth_id=modelParams['growth_id']
        total_enzymes=modelParams['total_enzymes']

        print('product:',product_key,'model:',model_key)
        output_file=f'{model_key}_{productName}_gluc_{c_uptake}_ecFactory_result.xlsx'
        # check if the output file already exists
        if pd.io.common.file_exists(f'analysis_code/results/ecGEM_vs_efl/{output_file}'):
            print(f'File {output_file} already exists, skipping...')
            continue

        # load model
        model=load_model(model_filepath,model_type=model_key)

        # # prepare model
        # constrain_enzymes(model,
        #                   total_prot=total_enzymes,
        #                   model_type=model_key)

        # calculate the max biomass yield
        model.objective = growth_id  # biomass rxn
        if model_key=='etfl':
            model.reactions.get_by_id(c_source).bounds=-c_uptake,0
        elif model_key=='ecGEM':
            model.reactions.get_by_id(c_source).bounds = 0, c_uptake
        gluc_MW = 0.180156  # g/mmol
        max_yield = model.slim_optimize() / (c_uptake * gluc_MW)  # gDW / gGluc
        expYield = max_yield * 0.45
        alphaLims = (0.5 * expYield, 2 * expYield)

        action_thresholds = [0.05, 0.5, 1.05]

        # run ecFactory
        modelParams.update(productParam)
        import time

        start_time = time.time()
        print('start time:', start_time)
        try:
            # run ecFSEOF
            results = run_ecFactory.run_ecFactory_design(model=model,
                                                         modelParam=modelParams,
                                                         expYield=expYield,
                                                         alphaLims=alphaLims,
                                                         action_thresholds=action_thresholds,
                                                         remove_essential=True,
                                                         steps=123
                                                         )
        except:
            continue
        end_time = time.time()
        print('end time:', end_time)
        print('time cost:', end_time - start_time)

        # save results
        for key in results.keys():
            if isinstance(results[key],np.ndarray):
                results[key]=pd.Series(results[key])
            # 检查是否为dict
            if isinstance(results[key],dict):
                results[key]=pd.Series(results[key])
        productName=modelParams['productName']
        with pd.ExcelWriter(f'analysis_code/results/ecGEM_vs_efl/{output_file}') as writer:
            for key in results.keys():
                results[key].to_excel(writer, sheet_name=key)

