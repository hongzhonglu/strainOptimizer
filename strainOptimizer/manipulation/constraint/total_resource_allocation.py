# -*- coding: utf-8 -*-
from pytfa.optim.utils import symbol_sum
from etfl.optim.constraints import ConstantAllocation
from etfl.optim.variables import EnzymeVariable
import pandas as pd


def constrain_enzymes(model, total_prot, model_type: ['etfl', 'ecGEM'] = 'etfl'):
    '''Change total amount of enzymes resource for ecGEM or ETFL model'''
    if model_type == 'etfl':
        enz_vars = model.get_variables_of_type(EnzymeVariable)
        # we should first exclude dummy, ribosomes and rnaps
        exclusion = ['dummy_enzyme',  # 'rib', 'rib_mit', 'rnap', 'rnap_mit'
                     ]
        exclusion = ['EZ_{}'.format(x) for x in exclusion]
        enz_vars = [x for x in enz_vars if x.name not in exclusion]

        expr = symbol_sum([x for x in enz_vars])
        try:
            enz_fix = model.get_constraints_of_type(ConstantAllocation).get_by_id('enzyme_fix')
            model.remove_constraint(enz_fix)
        except:
            pass
        model.add_constraint(kind=ConstantAllocation,
                             hook=model,
                             expr=expr,
                             id_='enzyme_fix',
                             ub=total_prot,
                             lb=0)
        model.repair()  # update the new constraint
    elif model_type=='ecGEM':
        prot_pool_rxnID='prot_pool_exchange'
        try:
            model.reactions.get_by_id(prot_pool_rxnID).upper_bound=total_prot
        except:
            print('can not find prot_pool_exchange reaction:%s in the model'%prot_pool_rxnID)

    return model


def extract_geneList_from_compartment(compartment,annotaion_file=r'reference/proYeast9/data/compartment_annotation_refine.xlsx'):
    df_annotation=pd.read_excel(annotaion_file)
    # check if compartment in df_annotation
    if compartment not in df_annotation['compartment'].unique():
        print('compartment %s not in annotation file'%compartment)
        return []
    else:
        geneList=df_annotation[df_annotation['compartment']==compartment]['gene'].tolist()
        return geneList


def add_subproteome_constraint(model,geneList,id,ub,lb=None):
    prot_pool_rxnID='prot_pool_exchange'

    # extract enzyme list from gene list
    gene_to_prot_dict={}
    for g in model.genes:
        geneID = g.id
        for rxn in g.reactions:
            if 'draw_prot_' in rxn.id:
                draw_prot_rxnID = rxn.id
                gene_to_prot_dict[geneID] = draw_prot_rxnID
    enzymeList=[]
    for geneID in geneList:
        if geneID in gene_to_prot_dict.keys():
            enzymeList.append(gene_to_prot_dict[geneID])
    enzymeList=list(set(enzymeList))

    # build sub-proteome constraint
    prot_pool_met=model.metabolites.get_by_id('prot_pool[c]')
    coefficients = dict()
    for enzID in enzymeList:
        enzyme_var=model.reactions.get_by_id(enzID).forward_variable
        coefficient=float(model.reactions.get_by_id(enzID).metabolites[prot_pool_met])
        coefficients[enzyme_var] = coefficient
    # add sub-proteome constraint
    coefficients[model.reactions.get_by_id(prot_pool_rxnID).forward_variable] = float(ub)

    # set the constraint
    # whether set the lower bound of sub proteome constraint
    if lb is not None:
        range=float(ub-lb)
        constraint=model.problem.Constraint(range,lb=0,ub=range,name=id)
    else:
        # only set the upper bound
        constraint=model.problem.Constraint(0,lb=0,name=id)
    model.add_cons_vars([constraint])
    model.solver.update()
    constraint.set_linear_coefficients(coefficients=coefficients)

    # model.remove_cons_vars(constraint)
    return model


def add_organelleproteome_constraints(model, proteome_fractions, total_enzyme_pool=None, method='bounds', model_type='ecGEM'):
    # constrained_model=model.copy()
    if model_type=='ecGEM':
        constrained_model=model
        prot_pool_rxnID = 'prot_pool_exchange'
        if total_enzyme_pool is not None:
            constrained_model.reactions.get_by_id(prot_pool_rxnID).bounds=0,total_enzyme_pool
            print('total enzyme pool is set to %s' % total_enzyme_pool)

        for compartment,proteome_fraction in proteome_fractions.iterrows():
            if method=='bounds':
                ub=proteome_fraction['max']
                lb=proteome_fraction['min']
            elif method=='ub':
                ub=proteome_fraction['max']
                lb=None
            geneList=extract_geneList_from_compartment(compartment)
            # replase ' ' with '_'
            id=compartment.replace(' ','_')
            if len(geneList)>0:
                constrained_model=add_subproteome_constraint(model=constrained_model,
                                                             geneList=geneList,
                                                             id=id,
                                                             ub=ub,
                                                             lb=lb)
    else:
        raise ValueError('model_type should be ecGEM, but got %s' % model_type)
    return constrained_model