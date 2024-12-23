# -*- coding: utf-8 -*-
# date : 2024/5/9 
# author : wangh
# file : plot_envelope.py
# project : strainOptimizer
from strainOptimizer.io import load_model
from matplotlib import pyplot as plt
import numpy as np

def extracti_from_list(l,indices):
    return [l[i] for i in indices]

if __name__ == '__main__':
    import os
    from strainOptimizer.analysis.dataset import load_experiment_targets
    from strainOptimizer.simulation import ppFBA
    from strainOptimizer.visualization import calculate_flux_space,plot_envelope
    os.chdir(r'D:\code\github\strainOptimizer')
    solver = 'optlang-gurobi'
    model= load_model('examples/models/yeast/ecYeastGEM_batch.xml',model_type='ecGEM',solver=solver)
    all_geneList=[gene.id for gene in model.genes]
    c_source='r_1714_REV'
    c_uptake=10
    oe_factor=2
    model.reactions.get_by_id(c_source).bounds = c_uptake, c_uptake

    # for 2-phenylethanol
    x='r_2111'
    y='r_1589'
    df_exp=load_experiment_targets('2-phenylethanol')
    df_exp=df_exp[df_exp.index.isin(all_geneList)]
    flux_spaceList=[]
    # calculate wild-tpe
    flux_spaceList.append(calculate_flux_space(model,x,y))
    # calculate mutant
    for geneID in df_exp.index:
        mutant=model.copy()
        if geneID not in all_geneList:
            continue
        action= df_exp.loc[geneID,'action']
        rxnidList=[rxn.id for rxn in model.genes.get_by_id(geneID).reactions if rxn.id.startswith('r_') ]
        with mutant:
            if action=='OE':
                wt_sol=ppFBA(model=mutant,
                              targetID=x,
                              c_source='r_1714_REV',
                              c_uptake=10,
                              model_type='ecGEM')
                for rxnid in rxnidList:
                    wt_flux=wt_sol[rxnid]
                    if wt_flux>0:
                        mutant.reactions.get_by_id(rxnid).bounds = wt_flux*oe_factor, 1000
                    elif wt_flux<0:
                        mutant.reactions.get_by_id(rxnid).bounds = -1000, wt_flux*oe_factor
            elif action in ['KD','KO']:
                mutant.genes.get_by_id(geneID).knock_out()
            try:
                with mutant:
                    flux_space=calculate_flux_space(mutant,x,y)
            except:
                flux_space=None
            flux_spaceList.append(flux_space)

    labels=['wild-type']+list(df_exp.index)
    fig=plot_envelope(flux_spaceList,labels,show=False)
    indices=[0,10,11]
    fig=plot_envelope(extracti_from_list(flux_spaceList,indices),
                  labels=extracti_from_list(labels,indices),
                  show=False)
    # set x,y axis label
    fig.axes[0].set_xlabel('growth(/h)',fontsize=16)
    fig.axes[0].set_ylabel('2-phenylethanol(mmol/gDW/h)',fontsize=16)
    plt.show()


    # for spermidine
    x='r_2111'
    y='r_2051'
    df_exp=load_experiment_targets('spermidine')
    df_exp=df_exp[df_exp.index.isin(all_geneList)]
    flux_spaceList=[]
    # calculate wild-tpe
    flux_spaceList.append(calculate_flux_space(model,x,y))
    # calculate mutant
    for geneID in df_exp.index:
        if geneID not in all_geneList:
            continue
        action= df_exp.loc[geneID,'action']
        rxnidList=[rxn.id for rxn in model.genes.get_by_id(geneID).reactions if rxn.id.startswith('r_') ]
        with model:
            if action=='OE':
                wt_sol=ppFBA(model=model,
                              targetID=x,
                              c_source='r_1714_REV',
                              c_uptake=10,
                              model_type='ecGEM')
                for rxnid in rxnidList:
                    wt_flux=wt_sol[rxnid]
                    if wt_flux>0:
                        model.reactions.get_by_id(rxnid).bounds = wt_flux*oe_factor, 1000
                    elif wt_flux<0:
                        model.reactions.get_by_id(rxnid).bounds = -1000, wt_flux*oe_factor
            elif action in ['KD','KO']:
                model.genes.get_by_id(geneID).knock_out()
            try:
                flux_space=calculate_flux_space(model,x,y)
            except:
                flux_space=None
            flux_spaceList.append(flux_space)

    labels=['wild-type']+list(df_exp.index)
    i=10
    fig=plot_envelope(flux_spaceList[:i],labels[:i],show=False)
    indices=[0,10,11]
    fig=plot_envelope(extracti_from_list(flux_spaceList,indices),
                  labels=extracti_from_list(labels,indices),
                  show=False)
    # set x,y axis label
    fig.axes[0].set_xlabel('growth(/h)',fontsize=16)
    fig.axes[0].set_ylabel('spermidine(mmol/gDW/h)',fontsize=16)
    plt.show()



