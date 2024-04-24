# -*- coding: utf-8 -*-
# date : 2023/4/6 
# author : wangh
# file : ecYeast_ecFactory_heme_test.py
# project : etfl
import pandas as pd
import numpy as np
from strainOptimizer.io import load_model
from strainOptimizer.strainDesign.ecFactory import run_ecFactory

# load heme a product model
heme_ecYeast=load_model('examples/models/yeast/heme_ecYeastGEM.xml',model_type='ecGEM')
product_name='heme_a'
product_id='EX_heme_a'
model=heme_ecYeast

# calculate the max biomass yield
# set tolerance
model.tolerance=1e-9

# parameters for ecYeast
c_source="r_1714_REV"      # glucose exchange rxn
growth_id='r_2111'
model.objective=growth_id# biomass rxn
c_uptake=1
model.reactions.get_by_id('r_1714_REV').bounds=1,1
gluc_MW=0.180156  # g/mmol
max_yield=model.slim_optimize()/(c_uptake*gluc_MW) # gDW / gGluc
expYield=0.122
alphaLims=(0.5*expYield,2*expYield)

# prepare parameters for ecFactory algorithm
modelParam=pd.Series()
modelParam['targetID']=product_id
modelParam['c_source']=c_source
modelParam['c_uptake']=c_uptake
modelParam['productName']=product_name
action_thresholds=[0.05,0.3,1.05]     # rules for overexpression, knockout and knockdown

# set a timer
import time
start_time = time.time()
print('start time:',start_time)
# run ecFSEOF
results=run_ecFactory.run_ecFactory_design(model=model, modelParam=modelParam, expYield=expYield,alphaLims=alphaLims,action_thresholds=action_thresholds,remove_essential=True,model_type='ecGEM')
end_time = time.time()
print('end time:',end_time)
print('time cost:',end_time-start_time)

# print FSEOF results
print('FSEOF results:')
genetable=results['geneTable']
print('OE:',genetable[genetable['actions']=='OE'].shape[0])
print('KD:',genetable[genetable['actions']=='KD'].shape[0])
print('KO:',genetable[genetable['actions']=='KO'].shape[0])

# EUVA results
print('EUVA results:')
print('leval 1 candidates:',genetable[genetable['target_priority_leval']==1].shape[0])
print('leval 2 candidates:',genetable[genetable['target_priority_leval']==2].shape[0])
print('leval 3 candidates:',genetable[genetable['target_priority_leval']==3].shape[0])


# minimal combined set
min_result=results['min_set_analysis_result']
min_set_list=min_result[min_result['score']<0.999].index.tolist()
fseof_list=results['geneTable'].index.tolist()
euva_list=results['geneTable'][results['geneTable']['target_priority_leval'].isin([1,2,3])].index.tolist()


# save results into excel file
for key in results.keys():
    if isinstance(results[key],np.ndarray):
        results[key]=pd.Series(results[key])
    # 检查是否为dict
    if isinstance(results[key],dict):
        results[key]=pd.Series(results[key])
with pd.ExcelWriter(f'examples/result/ecYeast_{product_name}_gluc_{c_uptake}_ecFactory_result.xlsx') as writer:
# with pd.ExcelWriter(f'examples/result/ecYeast_{product_name}_gluc_{c_uptake}_ecFSEOF_result.xlsx') as writer:
    for key in results.keys():
        results[key].to_excel(writer, sheet_name=key)



