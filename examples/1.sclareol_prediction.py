# -*- coding: utf-8 -*-
# date : 2024/12/23 
import pandas as pd
import numpy as np
from strainOptimizer.io import load_model
from strainOptimizer.strainDesign.ecFactory import run_ecFactory
from strainOptimizer.manipulation.constraint.total_resource_allocation import constrain_enzymes


# set patameters
productParams={'productName':'sclareol',
                'targetID':'DM_sclareol_c'}

modelParams_dict={'ecGEM':{
    'model_type':'ecGEM',
    'c_source':"r_1714_REV",      # glucose exchange rxn
    'c_uptake':1,
    'growth_id':'r_2111',
    'total_enzymes':0.4,
    'model_filepath':'examples/models/yeast/ecYeastGEM_sclareol_batch.xml',
        },
'etfl':{
    'model_type':'etfl',
    'c_source':"r_1714",
    'c_uptake':1,
    'growth_id':'r_2111',
    'total_enzymes':0.1,
    'model_filepath':'examples/models/yeast/cEFL_sclareol.json',
}}

# choose ecGEM or ETFL model
model_key='ecGEM'
# model_key='etfl'

modelParams=modelParams_dict[model_key]
modelParams.update(productParams)


# load model
model=load_model(filename=modelParams['model_filepath'],
                         model_type=modelParams['model_type'])

# prepare model
model=constrain_enzymes(model,
                total_prot=modelParams['total_enzymes'],
                  model_type=modelParams['model_type'])

# calculate the max biomass yield
c_uptake=modelParams['c_uptake']
growth_id=modelParams['growth_id']
c_source=modelParams['c_source']
productName=modelParams['productName']
model.objective=growth_id  # biomass rxn
if modelParams['model_type']=='etfl':
    model.reactions.get_by_id(c_source).bounds=-c_uptake,-c_uptake
elif modelParams['model_type']=='ecGEM':
    model.reactions.get_by_id(c_source).bounds=c_uptake,c_uptake
gluc_MW=0.180156  # g/mmol
max_yield=model.slim_optimize()/(c_uptake*gluc_MW) # gDW / gGluc
expYield=max_yield*0.49
alphaLims=(0.5*expYield,2*expYield)

action_thresholds=[0.05,0.5,1.05]     # rules for overexpression, knockout and knockdown
# set a timer
import time
start_time = time.time()
print('start time:',start_time)
# run ecFactory
results=run_ecFactory.run_ecFactory_design(model=model, modelParam=modelParams, expYield=expYield,alphaLims=alphaLims,action_thresholds=action_thresholds,remove_essential=True)
end_time = time.time()
print('end time:',end_time)
print('time cost:',end_time-start_time)

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


from strainOptimizer.analysis.dataset import load_experiment_targets,calculate_exp_consistency
# load experiment data
exp_data=load_experiment_targets(product=modelParams['productName'])

level1_result=genetable
level2_result=genetable[genetable['target_priority_leval'].isin([1,2])]
level3_result=genetable[genetable['minimal candidates set']==1]

print('level 1 result experiment consistency:')
l1_eval=calculate_exp_consistency(predict_result=level1_result,exp_data=exp_data)
print('level 2 result experiment consistency:')
l2_eval=calculate_exp_consistency(predict_result=level2_result,exp_data=exp_data)
print('level 3 result experiment consistency:')
l3_eval=calculate_exp_consistency(predict_result=level3_result,exp_data=exp_data)

# save result
# save results into excel file
for key in results.keys():
    if isinstance(results[key],np.ndarray):
        results[key]=pd.Series(results[key])
    # 检查是否为dict
    if isinstance(results[key],dict):
        results[key]=pd.Series(results[key])
with pd.ExcelWriter(f'examples/result/{model_key}_{productName}_gluc_{c_uptake}_ecFactory_result.xlsx') as writer:
    for key in results.keys():
        results[key].to_excel(writer, sheet_name=key)