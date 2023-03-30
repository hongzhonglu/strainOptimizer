# -*- coding: utf-8 -*-
# date : 2023/3/18 
# author : wangh
from etfl.optim.constraints import ModelConstraint
from pytfa.optim.utils import symbol_sum

def ETFL_constrain_enz_conc(model, enzymes_bounds):
    '''add constraint for enzyme concentration in ETFL model
    parameters:
        model: ETFL model
        enzymes_bounds: pd.DataFrame, columns=['lb','ub'], index=enzymeID.!! lb and ub should be the scaled values
    return:
        model: modified ETFL model
    '''
    enzymeIDlist = enzymes_bounds.index.tolist()
    for enzID in enzymeIDlist:
        enz_vars = model.get_variables_of_type('EnzymeVariable')
        enz_var = enz_vars.get_by_id(enzID)
        exp = symbol_sum([enz_var])
        model.add_constraint(
            kind=ModelConstraint,
            hook=model,
            expr=exp,
            id_='enz_conc_' + enzID,
            lb=enzymes_bounds.loc[enzID, 'lb'],
            ub=enzymes_bounds.loc[enzID, 'ub']
        )
        # enz_const=model.constraints.get_by_id('MODC_enz_conc_'+enzID)
    model.repair()
    return model


# def ecGEM_constrain_enz_conc(model, enzymes_bounds):
#