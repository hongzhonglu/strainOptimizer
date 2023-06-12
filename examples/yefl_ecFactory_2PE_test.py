# -*- coding: utf-8 -*-
# date : 2023/4/6 
# author : wangh
# file : yefl_ecFactory_2PE_test.py
# project : etfl
# load packages
import pandas as pd
import numpy as np
from ETFLdesigner.ETFLdesigner.io.ETFL import load_etfl_model
from ETFLdesigner.ETFLdesigner.strain_design.ecFactory import run_ecFactory

# load model
yefl=load_etfl_model('models/yeast8_cEFL_2584_enz_128_bins__20221228_090737.json', solver='optlang-gurobi')
model=yefl

# parameters for yefl siumlation
product_name='2-phenylethanol'
product_id='r_1589'    # 2-PE exchange rxn
c_source="r_1714"      # glucose exchange rxn
growth_id=model.growth_reaction.id     # biomass rxn
c_uptake=1
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
action_thresholds=[0.05,0.5,1.05]     # rules for overexpression, knockout and knockdown
# set a timer
import time
start_time = time.time()
print('start time:',start_time)
# run ecFSEOF
results=run_ecFactory.run_ecFactory_design(model=model, modelParam=modelParam, expYield=expYield,alphaLims=alphaLims,action_thresholds=action_thresholds,remove_essential=True,model_type='etfl')
end_time = time.time()
print('end time:',end_time)
print('time cost:',end_time-start_time)


# save results into excel file
for key in results.keys():
    if isinstance(results[key],np.ndarray):
        results[key]=pd.Series(results[key])
    # 检查是否为dict
    if isinstance(results[key],dict):
        results[key]=pd.Series(results[key])
with pd.ExcelWriter(f'code_etfl/ETFLdesigner/output/yefl_{product_name}_gluc_{c_uptake}_ecFactory_result.xlsx') as writer:
    for key in results.keys():
        results[key].to_excel(writer, sheet_name=key)
