# -*- coding: utf-8 -*-
# date : 2024/5/6 
# author : wangh
# file : iBridge.py
# project : strainOptimizer
from strainOptimizer.strainDesign import run_iBridge_design
from strainOptimizer.io import load_model
import os
from cobra.io import read_sbml_model

os.chdir(r'D:\code\github\strainOptimizer')
solver = 'optlang-gurobi'

# load models
ecYeast = load_model('examples/models/yeast/ecYeastGEM_batch.xml',model_type='ecGEM',solver=solver)
efl = load_model('examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json', model_type='etfl',
                 solver=solver)
gem=read_sbml_model('examples/models/yeast/yeast-GEM.xml')
models_dict={
    # 'ecYeast':{'model':ecYeast,'model_type':'ecGEM'},
    #          'efl':{'model':efl,'model_type':'etfl'},
             'gem':{'model':gem,'model_type':'GEM'}}

methods={'moma','mopa','ppfba'}
resuts={}
for key,model_dict in models_dict.items():
    resuts[key]={}
    model = model_dict['model']
    model_type = model_dict['model_type']
    if model_type == 'etfl':
        c_source = 'r_1714'
        cov_threshold = 0.01
        linear = True
    elif model_type == 'ecGEM':
        c_source = 'r_1714_REV'
        cov_threshold = 0.1
        linear = False
    elif model_type == 'GEM':
        c_source = 'r_1714'
        cov_threshold = 0.1
        linear = False
    c_uptake = 10
    # product='r_1589'   # 2-phenylethanol
    product='r_2051'   # spermidine
    for method in methods:
        with model:
            result = run_iBridge_design(model=model,
                                      targetID=product,
                                      c_source=c_source,
                                      c_uptake=c_uptake,
                                      model_type=model_type,
                                      method=method,
                                      tol=0.01,
                                      linear=linear,
                                      cov_threshold=cov_threshold)
        resuts[key][method]=result


# result evaluation
from strainOptimizer.analysis.dataset import load_experiment_targets,calculate_exp_consistency


if product=='r_1589':
    productName='2-phenylethanol'
    df_exp=load_experiment_targets(product=productName)
elif product=='r_2051':
    productName='spermidine'
    df_exp=load_experiment_targets(product=productName)

eval_results=dict()
for model_type,items in resuts.items():
    eval_results[model_type]={}
    for method,result in items.items():
        print(f'{model_type} {method}:')
        df_pred=result['endogenous_gene_result']
        # rename gene_action as action
        df_pred.rename(columns={'gene_action':'action'},inplace=True)
        exp_consistency=calculate_exp_consistency(df_pred,df_exp)
        eval_results[model_type][method]=exp_consistency


# load json data
import json

with open('data/universal_model.json', 'r') as f:
    universe = json.load(f)

for met in universe['metabolites']:
    annot_dict={}
    for annot in met['annotation']:
        key=annot[0]
        value=annot[1]
        annot_dict[key]=value
    met['annotation']=annot_dict

for rxn in universe['reactions']:
    annot_dict={}
    for annot in rxn['annotation']:
        key=annot[0]
        value=annot[1]
        annot_dict[key]=value
    rxn['annotation']=annot_dict

