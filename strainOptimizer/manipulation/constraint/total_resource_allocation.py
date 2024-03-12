# -*- coding: utf-8 -*-
from pytfa.optim.utils import symbol_sum
from etfl.optim.constraints import ConstantAllocation
from etfl.optim.variables import EnzymeVariable


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