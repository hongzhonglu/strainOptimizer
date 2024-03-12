'''Compare the enzyme concentration range by EUVA analysis between ecModel and ETFL model'''
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from strainOptimizer.analysis.enzyme_variety_analysis import enzymeVA
from strainOptimizer.io import load_model
from strainOptimizer.manipulation.constraint.total_resource_allocation import constrain_enzymes
from strainOptimizer.manipulation.variable.enzyme import genelist_to_enzymelist
from strainOptimizer.simulation.pprotFBA import pprotFBA_prot_conc
from etfl.optim.constraints import ForwardCatalyticConstraint, BackwardCatalyticConstraint
import os
# set working directory
os.chdir(r'D:\github\strainOptimizer')


def saturate_enzymes(model,rxnList=None,sol=None,tol=1e-6):
    '''Saturate all catalytic enzymes in ETFL models'''
    if rxnList is None:
        rxnList = [reaction for reaction in model.reactions if type(reaction).__name__ == 'EnzymaticReaction']
    if sol is None:
        sol=model.optimize()
    fluxes=sol.fluxes
    # get all catalytic constraints
    for rxn in rxnList:
        rxnID = rxn.id
        ub=0
        lb=-tol
        flux=fluxes[rxnID]
        if flux>0:
            try:
                fwd_cons = model.forward_catalytic_constraint.get_by_id(rxnID)
            except:
                print('can not find forward_catalytic_constraint:%s'%rxnID)
                continue
            expr_fwd = fwd_cons.expr
            model.remove_constraint(fwd_cons)
            new_fwd_cons = ForwardCatalyticConstraint(reaction=rxn, expr=expr_fwd, ub=ub, lb=lb, queue=True)
            model._cons_dict[new_fwd_cons.name] = new_fwd_cons
        elif flux==0:
            continue
        else:
            try:
                bwd_cons = model.backward_catalytic_constraint.get_by_id(rxnID)
            except:
                print('can not find backward_catalytic_constraint:%s'%rxnID)
                continue
            expr_bwd = bwd_cons.expr
            model.remove_constraint(bwd_cons)
            new_bwd_cons = BackwardCatalyticConstraint(reaction=rxn, expr=expr_bwd, ub=ub, lb=lb, queue=True)
            model._cons_dict[new_bwd_cons.name] = new_bwd_cons

    model.repair()
    return model

def euva_analysis(model, c_source, c_uptake, targetID, geneIDlist, total_enzymes, model_type:['ecGEM','etfl']):

    # get enzymes from genes
    enzIDlist, gene_enz_dict = genelist_to_enzymelist(model, geneIDlist, model_type=model_type)

    # set total amount of enzymes
    model = constrain_enzymes(model, total_enzymes, model_type=model_type)

    # saturate enzymes for ETFL model
    if model_type=='etfl':
        model.reactions.get_by_id(c_source).bounds = -c_uptake, -c_uptake
        sol = model.optimize()
        # # option 1: Only saturate enzymes related to target geneList
        # rxnIDlist=[rxn.id for geneID in geneIDlist for rxn in model.genes.get_by_id(geneID).reactions]
        # rxnIDlist=list(set(rxnIDlist))
        # rxnList=[model.reactions.get_by_id(rxnID) for rxnID in rxnIDlist]

        # option 2: saturate all catalytic enzymes
        rxnList = [reaction for reaction in model.reactions if type(reaction).__name__ == 'EnzymaticReaction']
        rxnList = [reaction for reaction in rxnList if reaction.id.startswith('r_')]

        model=saturate_enzymes(model,rxnList=rxnList,sol=sol)
        gr=model.slim_optimize()
        if gr is None:
            print('Saturate enzymes failed!')
            return None
        else:
            print('Saturate enzymes successfully!')
            print('Growth rate:',gr)

    # run EUVA analysis
    euva_result = enzymeVA(model=model,
                              c_source=c_source,
                              c_uptake=c_uptake,
                              targetID=targetID,
                              model_type=model_type,
                              enzymeIDlist=enzIDlist)
    pprotFBA_result = pprotFBA_prot_conc(model=model,
                                            c_source=c_source,
                                            c_uptake=c_uptake,
                                            targetID=targetID,
                                            model_type=model_type,
                                            enzymeIDlist=enzIDlist)
    euva_result['pprotFBA'] = pprotFBA_result.values

    if model_type=='etfl':
        # extract molecular weight of enzymes
        enzMWs= [model.enzymes.get_by_id(enzID).molecular_weight for enzID in enzIDlist]
        euva_result['MW']=enzMWs  # 1 kDa=1000 g/mol= 1 g/mmol
        # convert all enzyme concentration from g/gDW to mmol/gDW
        euva_result['minimum'] = euva_result['minimum'] / euva_result['MW']
        euva_result['maximum'] = euva_result['maximum'] / euva_result['MW']
        euva_result['pprotFBA'] = euva_result['pprotFBA'] / euva_result['MW']

    index_genelist = []
    for enzID in euva_result.index:
        for g, enzlist in gene_enz_dict.items():
            if enzID in enzlist:
                index_genelist.append(g)
                break
    euva_result['geneID'] = index_genelist
    # set geneID as index
    euva_result.set_index('geneID', inplace=True)
    return euva_result


# set simulation parameters
c_uptake=1
total_enzymes=0.1 # g/gDW
growth_id='r_2111'


# EUVA analysis for ecModel
ecmodel=load_model('examples/models/yeast/ecYeastGEM_batch.xml',model_type='ecGEM')
model=ecmodel

# set model-specific parameters
c_source="r_1714_REV"      # glucose exchange rxn

# randomly sample 50 genes from all model genes
test_gene_numb=50
gene_list=[g.id for g in model.genes]
np.random.seed(1)
gene_list=np.random.choice(gene_list,test_gene_numb,replace=False)

# run EUVA analysis
ec_euva_result=euva_analysis(model=model,
                        c_source=c_source,
                        c_uptake=c_uptake,
                        targetID=growth_id,
                        geneIDlist=gene_list,
                        total_enzymes=total_enzymes,
                        model_type='ecGEM')


# EUVA analysis for ETFL model
efl=load_model('examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json',model_type='etfl',
               solver='optlang-gurobi')
model=efl

# set model-specific parameters
c_source='r_1714'

# saturate all enzymes

# model=saturate_enzymes(model)

# run euva analysis
etfl_euva_result=euva_analysis(model=model,
                        c_source=c_source,
                        c_uptake=c_uptake,
                        targetID=growth_id,
                        geneIDlist=gene_list,
                        total_enzymes=total_enzymes,
                        model_type='etfl')



# rename columns
ec_euva_result.rename(columns={'minimum':'ec_minimum','maximum':'ec_maximum','pprotFBA':'ec_pprotFBA'},inplace=True)
etfl_euva_result.rename(columns={'minimum':'efl_minimum','maximum':'efl_maximum','pprotFBA':'efl_pprotFBA'},inplace=True)
# remove duplicated index
etfl_euva_result=etfl_euva_result[~etfl_euva_result.index.duplicated(keep='first')]

# merge results
euva_result=pd.concat([ec_euva_result,etfl_euva_result],axis=1,join='inner')


# plot results
fig,ax=plt.subplots(1,3,figsize=(15,5))
# plot minimum
ax[0].scatter(euva_result['ec_minimum'],euva_result['efl_minimum'],s=15)
# ax[0].plot([0,0.1],[0,0.1],'k--')
ax[0].set_xlabel('ecModel minimum')
ax[0].set_ylabel('ETFL minimum')
# plot maximum
ax[1].scatter(euva_result['ec_maximum'],euva_result['efl_maximum'],s=15)
ax[1].plot([0,0.001],[0,0.001],'k--')
ax[1].set_xlabel('ecModel maximum')
ax[1].set_ylabel('ETFL maximum')
# set axis limits
ax[1].set_xlim([0,0.001])
ax[1].set_ylim([0,0.001])

# plot pprotFBA
ax[2].scatter(euva_result['ec_pprotFBA'],euva_result['efl_pprotFBA'],s=15)
ax[2].plot([0,0.00001],[0,0.00001],'k--')
ax[2].set_xlabel('ecModel pprotFBA')
ax[2].set_ylabel('ETFL pprotFBA')
# set axis limits
ax[2].set_xlim([0,0.00001])
ax[2].set_ylim([0,0.00001])
plt.tight_layout()
plt.show()

