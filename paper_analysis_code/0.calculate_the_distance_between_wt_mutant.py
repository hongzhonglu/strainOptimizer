# -*- coding: utf-8 -*-
from strainOptimizer.io import load_model
from strainOptimizer.manipulation.constraint.total_resource_allocation import constrain_enzymes
from strainOptimizer.simulation import ppFBA
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

wt_vs_mutant_dict={'protein':{},'flux':{}}
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
        model = load_model(model_filepath, model_type=model_key)

        # set the tolerance
        model.tolerance = 1e-9

        # set c uptake
        model.reactions.get_by_id(c_source).bounds=0,c_uptake

        # calculate the flux and proteome for the wild type
        with model:
            # set the objective to maximize growth
            model.objective = growth_id

            sol=ppFBA(model=model,
                      targetID=growth_id,
                      c_source=c_source,
                        c_uptake=c_uptake,
                      model_type=model_key
                      )
            max_growth=sol.fluxes[growth_id]
            df_all=sol.fluxes.to_frame(name='WT')
            # fill < 1e-9 with 0
            df_all[df_all.abs()<1e-9]=0
            df_flux = df_all[df_all.index.str.startswith('r_')]
            df_prot=df_all[df_all.index.str.startswith('draw_prot')]

        # calculate the flux and proteome for the production mutant
        with model:

            # fix growth as 50% of the max growth
            model.reactions.get_by_id(growth_id).bounds=0.2*max_growth,1000
            # set the objective to maximize production
            model.objective = productID

            sol=ppFBA(model=model,
                      targetID=productID,
                      c_source=c_source,
                        c_uptake=c_uptake,
                      model_type=model_key
                      )
            df_mutant=sol.fluxes
            # fill < 1e-9 with 0
            df_mutant[df_mutant.abs()<1e-9]=0
            df_prot['mutant']=df_mutant[df_mutant.index.str.startswith('draw_prot')]
            df_flux['mutant']=df_mutant[df_mutant.index.str.startswith('r_')]

        wt_vs_mutant_dict['protein'][f'{product_key}']=df_prot
        wt_vs_mutant_dict['flux'][f'{product_key}']=df_flux

# save as excel file
with pd.ExcelWriter('analysis_code/results/wt_vs_mutant_protein_flux_v2.xlsx') as writer:
    for key, df in wt_vs_mutant_dict['protein'].items():
        df.to_excel(writer, sheet_name=f'protein_{key}')
    for key, df in wt_vs_mutant_dict['flux'].items():
        df.to_excel(writer, sheet_name=f'flux_{key}')