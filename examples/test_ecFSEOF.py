# set workdir
import os
os.chdir(r"D:\code\github\etfl")
# load package path
import sys
sys.path.append(r"D:\code\github\ETFLdesigner\src\ETFLdesigner\straindesign\fseof")

# load packages
import pandas as pd
from etfl.io.json import load_json_model
import run_ecFSEOF

# load model
test_model=load_json_model('models/ecoli_core.json', solver='optlang-gurobi')

# parameters for succinate production design in ecoli core model
targetID="EX_succ_e"
c_source="EX_glc__D_e"
c_uptake=1
gluc_MW=0.180156  # g/mmol
max_yield=0.456   # gDW / gGluc
expYield=max_yield*0.5
alphaLims=(0.25*max_yield,0.75*max_yield)


# for ecoli core model
modelParam=pd.Series()
modelParam['targetID']="EX_succ_e"
modelParam['c_source']="EX_glc__D_e"
modelParam['c_uptake']=1.0
results=run_ecFSEOF.run_ecFSEOF_design(model=test_model, modelParam=modelParam, expYield=expYield,action_thresholds=[0.05,0.5,1.05],remove_essential=False)
