# -*- coding: utf-8 -*-
from strainOptimizer.io import load_model
from strainOptimizer.manipulation.constraint.total_resource_allocation import constrain_enzymes
from strainOptimizer.strainDesign.ecFactory import run_ecFactory
import pandas as pd
import numpy as np
import os
os.chdir(r'D:\code\github\strainOptimizer')


modelParam_dict={
    'model_type':'GAN_ec',
    'c_source':"r_1714",
    'c_uptake':1,
    'growth_id':'r_2111',
    'total_enzymes':120,}

productParam_dict={
    '2-phenylethanol':{'productName':'2-phenylethanol',
                       'targetID':'r_1589',
                       'ecGEM_filepath':'examples/models/yeast/ecYeastGEM_batch.xml',
                       'etfl_filepath':'examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json',
                       'GAN_filepath':'examples/models/yeast/GAN_ecYeast/GAN_sen_v1.xml',
                       'GAN2_filepath':'examples/models/yeast/GAN_ecYeast/GAN_all_v2.xml'},
    # 'artemisinic_acid':{'productName':'artemisinic acid',
    #                'targetID':'SK_atemisinic_acid',
    #                'model_filepath':'examples/models/yeast/ecGEM_atemisinic.xml'},
    'heme_acid':{'productName':'heme',
                 'targetID':'EX_heme_a'
                 ,'ecGEM_filepath':'examples/models/yeast/heme_ecYeastGEM.xml',
                 'etfl_filepath':'examples/models/yeast/heme_cEFL.json',
                 'GAN_filepath':'examples/models/yeast/GAN_ecYeast/heme_GAN_sen_v1.xml',
                 'GAN2_filepath':'examples/models/yeast/GAN_ecYeast/heme_GAN_all_v2.xml'},
    'spermidine':{'productName':'spermidine',
                  'targetID':'r_2051',
                  'ecGEM_filepath':'examples/models/yeast/ecYeastGEM_batch.xml',
                'etfl_filepath':'examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json',
                'GAN_filepath':'examples/models/yeast/GAN_ecYeast/GAN_sen_v1.xml',
                'GAN2_filepath':'examples/models/yeast/GAN_ecYeast/GAN_all_v2.xml'
                  },
    'sclareol':{'productName':'sclareol',
                'targetID':'DM_sclareol_c',
                'ecGEM_filepath':'examples/models/yeast/sclareol_ecYeastGEM_batch.xml',
                'etfl_filepath':'examples/models/yeast/sclareol_cEFL.json',
                'GAN_filepath':'examples/models/yeast/GAN_ecYeast/sclareol_GAN_sen_v1.xml',
                'GAN2_filepath':'examples/models/yeast/GAN_ecYeast/sclareol_GAN_all_v2.xml'},
    'fatty_acid':{'productName':'free fatty acids',
                   'targetID':'r_2189',
                   'ecGEM_filepath':'examples/models/yeast/ecYeastGEM_batch.xml',
                  'etfl_filepath':'examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json',
                  'GAN_filepath':'examples/models/yeast/GAN_ecYeast/GAN_sen_v1.xml',
                  'GAN2_filepath':'examples/models/yeast/GAN_ecYeast/GAN_all_v2.xml'},
}

methodList=['mopa','ppfba']  # methods to be used for simulation

for product_key,productParam in productParam_dict.items():
    productName=productParam['productName']
    productID=productParam['targetID']
    model_filepath=productParam['GAN2_filepath']
    model_key='GAN_ec'

    c_source = modelParam_dict['c_source']
    c_uptake = modelParam_dict['c_uptake']
    growth_id = modelParam_dict['growth_id']
    total_enzymes = modelParam_dict['total_enzymes']

    for method in methodList:
        simulation_method = method
        modelParam_dict['simulation_method'] = simulation_method

        print('product:', product_key, 'model:', model_key, 'method:', simulation_method)
        output_file = f'GANall_{productName}_gluc_{c_uptake}_{simulation_method}_ecFSEOF_result.xlsx'
        # check if the output file already exists
        if pd.io.common.file_exists(f'analysis_code/results/moma_mopa_pfba/{output_file}'):
            print(f'File {output_file} already exists, skipping...')
            continue

        # load model
        model = load_model(model_filepath, model_type=modelParam_dict['model_type'])

        # # prepare model
        # constrain_enzymes(model,
        #                   total_prot=total_enzymes,
        #                   model_type='ecGEM')

        # calculate the max biomass yield
        model.objective = growth_id  # biomass rxn
        if model_key == 'etfl':
            model.reactions.get_by_id(c_source).bounds = -c_uptake, 0
        elif model_key == 'ecGEM':
            model.reactions.get_by_id(c_source).bounds = 0, c_uptake
        elif model_key == 'GAN_ec':
            model.reactions.get_by_id(c_source).bounds = -c_uptake, 0
        gluc_MW = 0.180156  # g/mmol
        max_yield = model.slim_optimize() / (c_uptake * gluc_MW)  # gDW / gGluc
        expYield = max_yield * 0.45
        alphaLims = (0.5 * expYield, 2 * expYield)

        action_thresholds = [0.05, 0.5, 1.05]

        # run ecFactory
        modelParam_dict.update(productParam)
        import time

        start_time = time.time()
        print('start time:', start_time)
        try:
            # run ecFSEOF
            results = run_ecFactory.run_ecFactory_design(model=model,
                                                         modelParam=modelParam_dict,
                                                         expYield=expYield,
                                                         alphaLims=alphaLims,
                                                         action_thresholds=action_thresholds,
                                                         remove_essential=True,
                                                         steps=1
                                                         )
        except:
            continue
        end_time = time.time()
        print('end time:', end_time)
        print('time cost:', end_time - start_time)

        # save results
        for key in results.keys():
            if isinstance(results[key], np.ndarray):
                results[key] = pd.Series(results[key])
            if isinstance(results[key], dict):
                results[key] = pd.Series(results[key])
        productName = modelParam_dict['productName']
        with pd.ExcelWriter(f'analysis_code/results/moma_mopa_pfba/{output_file}') as writer:
            for key in results.keys():
                results[key].to_excel(writer, sheet_name=key)