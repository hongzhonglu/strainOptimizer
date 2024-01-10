'''change total protein constraint in ETFL model'''
from strainOptimizer.io import load_model
from pytfa.optim.utils import symbol_sum
from etfl.optim.constraints import ConstantAllocation
from etfl.optim.variables import EnzymeVariable

def constrain_enzymes(model, total_prot,model_type:['etfl','ecGEM']='etfl'):
    '''Change total amount of enzymes resource for ecGEM or ETFL model'''
    if model_type=='etfl':
        enz_vars = model.get_variables_of_type(EnzymeVariable)
        # we should first exclude dummy, ribosomes and rnaps
        exclusion = ['dummy_enzyme',  # 'rib', 'rib_mit', 'rnap', 'rnap_mit'
                        ]
        exclusion = ['EZ_{}'.format(x) for x in exclusion]
        enz_vars = [x for x in enz_vars if x.name not in exclusion]

        expr = symbol_sum([x for x in enz_vars])

        model.add_constraint(kind=ConstantAllocation,
                                 hook=model,
                                 expr=expr,
                                 id_='enzyme_fix',
                                 ub=total_prot)
        model.repair()    # update the new constraint
    #elif model_type=='ecGEM':

    return model


# load model
model=load_model(filename='examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json', model_type='etfl',solver='optlang-gurobi')

model.slim_optimize()
enzymes_ratio=0.1

# # get enzyme_fix constraint
# alloc_constraints=model.get_constraints_of_type(ConstantAllocation)
# enzymes_fix=alloc_constraints.get_by_id('enzyme_fix')
# # remove old constraint
# model.remove_constraint(enzymes_fix)
# add new constraint
constrain_enzymes(model,enzymes_ratio,model_type='etfl')

model.slim_optimize()
