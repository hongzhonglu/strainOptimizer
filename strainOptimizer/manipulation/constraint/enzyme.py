# -*- coding: utf-8 -*-
# date : 2023/3/18 
# author : wangh
# file : enzyme.py
# project : strainOptimizer
from etfl.optim.constraints import ModelConstraint
from pytfa.optim.utils import symbol_sum
from etfl.optim.constraints import ForwardCatalyticConstraint, BackwardCatalyticConstraint


def ETFL_constrain_enz_conc(model,enzymes_bounds,tol_ratio=0.01):
    '''add constraint for enzyme concentration in ETFL model
    parameters:
        model: ETFL model
        enzymes_bounds: pd.DataFrame, columns=['lb','ub'], index=enzymeID.!! lb and ub should be the scaled values
    return:
        model: modified ETFL model
    '''
    enzymeIDlist=enzymes_bounds.index.tolist()
    for enzID in enzymeIDlist:
        enz_vars = model.get_variables_of_type('EnzymeVariable')
        enz_var=enz_vars.get_by_id(enzID)
        exp=symbol_sum([enz_var])
        lb=enzymes_bounds.loc[enzID,'lb']
        ub=enzymes_bounds.loc[enzID,'ub']
        if lb>ub:
            print(enzID,"has a lower bound larger than upper bound")
            lb=lb*(1-tol_ratio)
            ub=ub*(1+tol_ratio)
        model.add_constraint(
            kind=ModelConstraint,
            hook=model,
            expr=exp,
            id_='enz_conc_'+enzID,
            lb=lb,
            ub=ub
        )
        # enz_const=model.constraints.get_by_id('MODC_enz_conc_'+enzID)
    model.repair()
    return model


def ecGEM_constrain_enz_conc(model,enzymes_bounds,tol_ratio=0.01):
    '''
    add constraint for enzyme concentration in ecGEM model
    parameters:
        model: ecGEM model
        enzymes_bounds: pd.DataFrame, columns=['lb','ub'], index=enzymeID.!! lb and ub should be the scaled values
    return:
        model: modified ecGEM model
    '''
    enzymeIDlist=enzymes_bounds.index.tolist()
    for enzID in enzymeIDlist:
        ub=enzymes_bounds.loc[enzID,'ub']
        lb=enzymes_bounds.loc[enzID,'lb']
        if ub>lb:
            model.reactions.get_by_id(enzID).bounds=lb,ub
        else:
            print(enzID,"has a lower bound larger than upper bound")
            lb=lb*(1-tol_ratio)
            ub=ub*(1+tol_ratio)
            model.reactions.get_by_id(enzID).bounds=lb,ub


    return model


def saturate_enzymes(model,rxnList=None,sol=None,tol=1e-6):
    '''Saturate all catalytic enzymes in ETFL models'''
    if rxnList is None:
        rxnList = [reaction for reaction in model.reactions if type(reaction).__name__ == 'EnzymaticReaction']
        rxnList = [reaction for reaction in rxnList if reaction.id.startswith('r_')]
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


def change_enz_conc_bounds(model, enzID, lb=None, ub=None):
    '''Change the enzyme concentration constraint bounds in ETFL model'''
    old_constraint = model.get_constraints_of_type('ModelConstraint').get_by_id('enz_conc_'+enzID)
    expression= old_constraint.expr
    old_lb=old_constraint.constraint.lb
    old_ub=old_constraint.constraint.ub
    if lb is None:
        lb=old_lb
    if ub is None:
        ub=old_ub
    # remove the original constraint
    model.remove_constraint(old_constraint)
    # add new constraint
    model.add_constraint(kind=ModelConstraint,
                         hook=model,
                         expr=expression,
                         id_='enz_conc_'+enzID,
                         lb=lb,
                         ub=ub)
    model.repair()
    return model