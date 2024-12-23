#!D:\softwares\programming\anoconda\envs\etfl\python.exe
# set workdir
import os
import numpy as np

os.chdir(r"D:/code/github/etfl")
# package path
import sys
sys.path.append(r"D:/code/github/etfl/code_etfl/strainOptimizer/ecFactory")
sys.path.append(r"D:/code/github/etfl")
# load packages
import pandas as pd
from etfl.io.json import load_json_model
import ecFactory.run_ecFSEOF as run_ecFSEOF

# load model
yefl=load_json_model('yetfl/models/yeast8_cEFL_2584_enz_128_bins__20221115_120238.json', solver='optlang-gurobi')

# calculate the max biomass yield
model=yefl
medium=model.medium
model.reactions.get_by_id('r_1714').bounds=-1,-1
model.objective=model.growth_reaction.id
sol=model.optimize()


# parameters for 2-phenylethanol production design in yefl model
product_name='2-phenylethanol'
product_id='r_1589'    # 2-PE exchange rxn
c_source="r_1714"      # glucose exchange rxn
growth_id=model.growth_reaction.id     # biomass rxn
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
action_thresholds=[0.05,0.5,1.05]     # rules for overexpression, knockout and knockdown
# set a timer
import time
start_time = time.time()
print('start time:',start_time)
# run ecFSEOF
results=run_ecFSEOF.run_ecFSEOF_design(model=model, modelParam=modelParam, expYield=expYield,action_thresholds=action_thresholds,remove_essential=False)
end_time = time.time()
print('end time:',end_time)
print('time cost:',end_time-start_time)

genetable=results['geneTable']
print('OE:',genetable[genetable['actions']=='OE'].shape[0])
print('KO:',genetable[genetable['actions']=='KO'].shape[0])
print('KD:',genetable[genetable['actions']=='KD'].shape[0])



# plot histogram of 0<k<1 k_scores distribution
import matplotlib.pyplot as plt
plt.hist(genetable[(genetable['k_score']<=1)&(genetable['k_score']>=0)]['k_score'],bins=20)
plt.xlabel('k_score')
plt.ylabel('targets count')
plt.title('yefl 2-PE design [0,1] gene k_score distribution')
plt.show()

# save results to excel file
# convert ndarray to Series
for key in results.keys():
    if isinstance(results[key],np.ndarray):
        results[key]=pd.Series(results[key])
with pd.ExcelWriter('code_etfl/strainOptimizer/output/yefl_2PE_design_results.xlsx') as writer:
    for key in results.keys():
        results[key].to_excel(writer, sheet_name=key)


#  save result to pickle file
# import pickle
# with open('code_etfl/strainOptimizer/output/yefl_2PE_design_results.pickle', 'wb') as handle:
#     pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)
# # load results dict
# with open('code_etfl/strainOptimizer/output/yefl_2PE_design_results.pickle', 'rb') as handle:
#     results = pickle.load(handle)
