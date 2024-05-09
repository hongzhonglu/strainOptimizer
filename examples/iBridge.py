# -*- coding: utf-8 -*-
# date : 2024/5/6 
# author : wangh
# file : iBridge2.py
# project : strainOptimizer
from strainOptimizer.strainDesign import run_iBridge_design
from strainOptimizer.io import load_model
import os
from cobra.io import read_sbml_model

os.chdir(r'D:\code\github\strainOptimizer')
solver = 'optlang-gurobi'

ecYeast = load_model('examples/models/yeast/ecYeastGEM_batch.xml',model_type='ecGEM',solver=solver)
efl = load_model('examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json', model_type='etfl',
                 solver=solver)
gem=read_sbml_model('examples/models/yeast/yeast-GEM.xml')
models_dict={'ecYeast':{'model':ecYeast,'model_type':'ecGEM'},
             'efl':{'model':efl,'model_type':'etfl'},
             'gem':{'model':gem,'model_type':'GEM'}}

for key,model_dict in models_dict.items():
    model = model_dict['model']
    model_type = model_dict['model_type']
    if model_type == 'etfl':
        c_source = 'r_1714'
        cov_threshold = 0.01
    elif model_type == 'ecGEM':
        c_source = 'r_1714_REV'
        cov_threshold = 0.1
    elif model_type == 'GEM':
        c_source = 'r_1714'
        cov_threshold = 0.1
    c_uptake = 10
    product='r_1589'
    results = run_iBridge_design(model=model,
                                  targetID='r_1589',
                                  c_source=c_source,
                                  c_uptake=c_uptake,
                                  model_type=model_type,
                                  method='moma',
                                  tol=0.01,
                                  linear=True,
                                  cov_threshold=0.01)

