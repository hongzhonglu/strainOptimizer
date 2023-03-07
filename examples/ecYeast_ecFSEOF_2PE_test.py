# -*- coding: utf-8 -*-
# date : 2023/3/5 
# author : wangh
# file : ecYeast_ecFSEOF_2PE_test.py
# project : etfl
from cobra.io import read_sbml_model
import pandas as pd
import numpy as np

import sys
sys.path.append('code_etfl/ETFLdesigner/ecFactory')
from code_etfl.ETFLdesigner.ecFactory import run_ecFSEOF

ecYeast=read_sbml_model('code_etfl/ETFLdesigner/models/ecYeastGEM_batch.xml')
# geckoYeast=GeckoModel('code_etfl/ETFLdesigner/models/ecYeastGEM_batch.xml')


# minimize total prot: total prot pseudoreaction: prot_pool_exchange
prot_pool_rxn=ecYeast.reactions.get_by_id('prot_pool_exchange')
met_prot_pool=ecYeast.metabolites.get_by_id('prot_pool[c]')

# calculate the max biomass yield
model=ecYeast
medium=model.medium
model.reactions.get_by_id('r_1714_REV').bounds=1,1
sol=model.optimize()

# parameters for 2-phenylethanol production design in yefl model
product_name='2-phenylethanol'
product_id='r_1589'    # 2-PE exchange rxn
c_source="r_1714_REV"      # glucose exchange rxn
growth_id='r_2111'    # biomass rxn
c_uptake=1
gluc_MW=0.180156  # g/mmol
max_yield=sol.objective_value/(c_uptake*gluc_MW) # gDW / gGluc
expYield=max_yield*0.5
alphaLims=(0.25*max_yield,0.75*max_yield)

# prepare parameters for ecFSEOF
modelParam=pd.Series()
modelParam['targetID']=product_id
modelParam['c_source']=c_source
modelParam['c_uptake']=c_uptake
action_thresholds=[0.05,0.3,1.05]     # rules for overexpression, knockout and knockdown
# set a timer
import time
start_time = time.time()
print('start time:',start_time)
# run ecFSEOF
results=run_ecFSEOF.run_ecFSEOF_design(model=model, modelParam=modelParam, expYield=expYield,action_thresholds=action_thresholds,remove_essential=False,model_type='ecGEM')
end_time = time.time()
print('end time:',end_time)
print('time cost:',end_time-start_time)

# save results into excel file
for key in results.keys():
    if isinstance(results[key],np.ndarray):
        results[key]=pd.Series(results[key])
with pd.ExcelWriter('code_etfl/ETFLdesigner/output/ecYeast_2PE_design_results.xlsx') as writer:
    for key in results.keys():
        results[key].to_excel(writer, sheet_name=key)