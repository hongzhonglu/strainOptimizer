# -*- coding: utf-8 -*-
from strainOptimizer.io import load_model
from strainOptimizer.manipulation.constraint.total_resource_allocation import constrain_enzymes
from strainOptimizer.strainDesign.ecFactory import run_ecFactory
import pandas as pd


modelParams_dict={'ecGEM':{
    'model_type':'ecGEM',
    'c_source':"r_1714_REV",      # glucose exchange rxn
    'c_uptake':10,
    'growth_id':'r_2111',
    'total_enzymes':0.1,
        }}

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

proteome_cost_dict = {}
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

        # load model
        model=load_model(model_filepath,model_type=model_key)

        # calculate theoretical maximal proteome cost
        with model:
            # set the tolerance
            model.tolerance = 1e-9

            # open substrate
            model.reactions.get_by_id(c_source).bounds = 0,1000

            # maximize the production
            model.objective = productID
            model.objective_direction = 'max'
            product=model.slim_optimize()

            # fix product and minimize the protein pool
            model.reactions.get_by_id(productID).bounds = product, product

            model.objective_direction = 'min'
            model.objective= 'prot_pool_exchange'

            min_prot_pool= model.slim_optimize()

            proteome_cost=min_prot_pool/product
            proteome_cost_dict[product_key]=proteome_cost

