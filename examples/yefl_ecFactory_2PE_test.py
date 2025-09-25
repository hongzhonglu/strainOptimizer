# -*- coding: utf-8 -*-

# load packages
import pandas as pd
import numpy as np
from strainOptimizer.io import load_etfl_model
from strainOptimizer.strainDesign.ecFactory import run_ecFactory
from strainOptimizer.manipulation.constraint.total_resource_allocation import constrain_enzymes

# load model
yefl=load_etfl_model('examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json', solver='optlang-gurobi')
model=yefl

# # set total amount of enzymes
# total_enzymes=0.1 # g/gDW
# constrain_enzymes(model,total_enzymes,model_type='etfl')

# parameters for yefl siumlation
product_name='2-phenylethanol'
product_id='r_1589'    # 2-PE exchange rxn
c_source="r_1714"      # glucose exchange rxn
growth_id=model.growth_reaction.id     # biomass rxn
c_uptake=5
gluc_MW=0.180156  # g/mmol
model.reactions.get_by_id(c_source).bounds=-c_uptake,-c_uptake
model.objective=model.growth_reaction.id
sol=model.optimize()
max_yield=sol.objective_value/(c_uptake*gluc_MW) # gDW / gGluc
expYield=max_yield*0.49
alphaLims=(0.5*expYield,2*expYield)

# prepare parameters for ecFSEOF
modelParam=pd.Series()
modelParam['targetID']=product_id
modelParam['c_source']=c_source
modelParam['c_uptake']=c_uptake
modelParam['productName']=product_name
modelParam['model_type']='etfl'  # model type, 'ecGEM' or 'etfl'
action_thresholds=[0.05,0.2,1.1]     # rules for overexpression, knockout and knockdown
# set a timer
import time
start_time = time.time()
print('start time:',start_time)
# run ecFSEOF
results=run_ecFactory.run_ecFactory_design(model=model,
                                           modelParam=modelParam,
                                           expYield=expYield,
                                           alphaLims=alphaLims,
                                           action_thresholds=action_thresholds,
                                           remove_essential=True,
                                           steps=123)
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


# save results into excel file
for key in results.keys():
    if isinstance(results[key],np.ndarray):
        results[key]=pd.Series(results[key])
    # 检查是否为dict
    if isinstance(results[key],dict):
        results[key]=pd.Series(results[key])
with pd.ExcelWriter(f'examples/result/yefl_{product_name}_gluc_{c_uptake}_ecFactory_result.xlsx') as writer:
    for key in results.keys():
        results[key].to_excel(writer, sheet_name=key)

# test
