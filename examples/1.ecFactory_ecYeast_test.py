# -*- coding: utf-8 -*-
# load packages
import pandas as pd
import numpy as np
from strainOptimizer.io import load_model
from strainOptimizer.strainDesign.ecFactory import run_ecFactory
from strainOptimizer.manipulation.constraint.total_resource_allocation import constrain_enzymes

# set tolerance
import cobra
cobra.Configuration().tolerance=1e-9

productParam_dict={
    # '2-phenylethanol':{'productName':'2-phenylethanol',
    #                    'targetID':'r_1589',
    #                    'model_filepath':'examples/models/yeast/ecYeastGEM_batch.xml'},
    # 'artemisinic_acid':{'productName':'artemisinic acid',
    #                'targetID':'SK_atemisinic_acid',
    #                'model_filepath':'examples/models/yeast/ecGEM_atemisinic.xml'},
    # 'heme_acid':{'productName':'heme',
    #              'targetID':'EX_heme_a'
    #              ,'model_filepath':'examples/models/yeast/heme_ecYeastGEM.xml'},
    # 'spermidine':{'productName':'spermidine',
    #               'targetID':'r_2051',
    #               'model_filepath':'examples/models/yeast/ecYeastGEM_batch.xml'}
    'sclareol':{'productName':'sclareol',
                'targetID':'DM_sclareol_c',
                'model_filepath':'examples/models/yeast/ecYeastGEM_sclareol_batch.xml'}
}

modelParams_dict={'ecGEM':{
    'model_type':'ecGEM',
    'c_source':"r_1714_REV",      # glucose exchange rxn
    'c_uptake':1,
    'growth_id':'r_2111',
    'total_enzymes':0.1,
        },
'etfl':{
    'model_type':'etfl',
    'c_source':"r_1714",
    'c_uptake':1,
    'growth_id':'r_2111',
    'total_enzymes':0.1,
}}

consistency_result=dict()
# load infomation
productParams=productParam_dict['sclareol']
model_key='ecGEM'
for model_key,modelParam in modelParams_dict.items():
    if model_key=='etfl':
        continue
    model_consistency_result=pd.DataFrame(index=['l1','l2','l3'])
    for product_key,productParam in productParam_dict.items():
        print('product:',product_key,'model:',model_key)
        productParams=productParam_dict[product_key]
        modelParams=modelParams_dict[model_key]

        modelParams.update(productParams)

        # load model
        model=load_model(filename=modelParams['model_filepath'],
                         model_type=modelParams['model_type'])

        # prepare model
        constrain_enzymes(model,
                        total_prot=modelParams['total_enzymes'],
                          model_type=modelParams['model_type'])

        # calculate the max biomass yield
        c_uptake=modelParams['c_uptake']
        growth_id=modelParams['growth_id']
        c_source=modelParams['c_source']
        model.objective=growth_id  # biomass rxn
        if modelParams['model_type']=='etfl':
            model.reactions.get_by_id(c_source).bounds=-c_uptake,-c_uptake
        elif modelParams['model_type']=='ecGEM':
            model.reactions.get_by_id(c_source).bounds=c_uptake,c_uptake
        gluc_MW=0.180156  # g/mmol
        max_yield=model.slim_optimize()/(c_uptake*gluc_MW) # gDW / gGluc
        expYield=max_yield*0.49
        alphaLims=(0.5*expYield,2*expYield)

        # maximize the production
        with model:
            model.objective=modelParams['targetID']
            # model.objective=growth_id
            constrain_enzymes(model,
                              total_prot=1,
                              model_type='ecGEM')
            for c in np.linspace(1,20,20):
                model.reactions.get_by_id(c_source).bounds=0, c
                print(model.slim_optimize())


        action_thresholds=[0.05,0.3,1.05]     # rules for overexpression, knockout and knockdown

        # set a timer
        import time
        start_time = time.time()
        print('start time:',start_time)
        # run ecFSEOF
        results=run_ecFactory.run_ecFactory_design(model=model,
                                                   modelParam=modelParams,
                                                   expYield=expYield,
                                                   alphaLims=alphaLims,
                                                   action_thresholds=action_thresholds,
                                                   remove_essential=True,
                                                   )
        end_time = time.time()
        print('end time:',end_time)
        print('time cost:',end_time-start_time)

        # # save result
        # save_result=False
        # if save_result:
        #     for key in results.keys():
        #         if isinstance(results[key],np.ndarray):
        #             results[key]=pd.Series(results[key])
        #         # 检查是否为dict
        #         if isinstance(results[key],dict):
        #             results[key]=pd.Series(results[key])
        #     productName=modelParams['productName']
        #     with pd.ExcelWriter(f'examples/result/ecYeast_{productName}_gluc_{c_uptake}_ecFactory_result.xlsx') as writer:
        #         for key in results.keys():
        #             results[key].to_excel(writer, sheet_name=key)
        #

        # Evaluate result

        # 1.check output target genes number
        # print FSEOF results
        print('FSEOF results:')
        genetable=results['geneTable']
        print('OE:',genetable[genetable['action']=='OE'].shape[0])
        print('KD:',genetable[genetable['action']=='KD'].shape[0])
        print('KO:',genetable[genetable['action']=='KO'].shape[0])

        # EUVA results
        print('EUVA results:')
        print('leval 1 candidates:',genetable[genetable['target_priority_leval']==1].shape[0])
        print('leval 2 candidates:',genetable[genetable['target_priority_leval']==2].shape[0])
        print('leval 3 candidates:',genetable[genetable['target_priority_leval']==3].shape[0])

        # minimal combined set
        min_result=results['min_set_analysis_result']
        min_set_list=min_result[min_result['score']<0.999].index.tolist()
        print('minimal combined set:',len(min_set_list))

        # 2. Compare with experiment data
        # load experiment data
        from strainOptimizer.analysis.dataset import load_experiment_targets,calculate_exp_consistency
        # load experiment data
        exp_data=load_experiment_targets(product=modelParams['productName'])

        level1_result=genetable
        level2_result=genetable[genetable['target_priority_leval'].isin([1,2,3])]
        level3_result=genetable[genetable['minimal candidates set']==1]

        print('level 1 result experiment consistency:')
        l1_eval=calculate_exp_consistency(predict_result=level1_result,exp_data=exp_data)
        print('level 2 result experiment consistency:')
        l2_eval=calculate_exp_consistency(predict_result=level2_result,exp_data=exp_data)
        print('level 3 result experiment consistency:')
        l3_eval=calculate_exp_consistency(predict_result=level3_result,exp_data=exp_data)

        l1_consistency=l1_eval['overall']['consistency']
        l2_consistency=l2_eval['overall']['consistency']
        l3_consistency=l3_eval['overall']['consistency']

        model_consistency_result[product_key]=[l1_consistency,l2_consistency,l3_consistency]

    consistency_result[model_key]=model_consistency_result

