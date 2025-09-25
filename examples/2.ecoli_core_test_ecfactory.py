# -*- coding: utf-8 -*-
# date : 2023/2/26 
# author : wangh

# set workdir
# import os
# os.chdir(r"D:\code\github\etfl")
# # package path
# import sys
# sys.path.append(r"D:\code\github\etfl\code_etfl\strainOptimizer")

# load packages
import pandas as pd
from strainOptimizer.io import load_model
from strainOptimizer.strainDesign.ecFactory import run_ecFactory

# load model
test_model=load_model('examples/models/ecoli/ecoli_core_curated.json', solver='optlang-gurobi',model_type='etfl')

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
# modelParam['targetID']="EX_mal__L_e"
modelParam['c_source']="EX_glc__D_e"
modelParam['c_uptake']=1.0
modelParam['productName']='malate'
modelParam['model_type']='etfl'
results=run_ecFactory.run_ecFactory_design(model=test_model,
                                           modelParam=modelParam,
                                           expYield=expYield,
                                           alphaLims=alphaLims,
                                           action_thresholds=[0.05,0.3,1.05],
                                           remove_essential=False,)

print('OE candidate:',len(results['geneTable'][results['geneTable']['actions']=='OE']))
print('KD candidate:',len(results['geneTable'][results['geneTable']['actions']=='KD']))
print('KO candidate:',len(results['geneTable'][results['geneTable']['actions']=='KO']))

