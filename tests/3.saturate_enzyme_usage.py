'''Modify the enzyme catalytic constraint to saturate all enzyme usage by making v=kcat*[E].'''
from strainOptimizer.io import load_model
from etfl.optim.constraints import ForwardCatalyticConstraint, BackwardCatalyticConstraint
import os

# strategy 2: add a new saturate enzyme constraint for each reations
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


if __name__ == '__main__':
    # set working directory
    os.chdir(r'D:\github\strainOptimizer')
    # load model
    # test_model=load_model("examples/models/ecoli/ecoli_core_curated.json",model_type='etfl',solver='optlang-gurobi')
    # model=test_model

    model=load_model('examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json',model_type='etfl',solver='optlang-gurobi')

    # sol=model.optimize()
    # model=load_model("examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json",model_type='etfl',solver='optlang-gurobi')
    # all_rxns=[reaction for reaction in model.reactions if type(reaction).__name__ == 'EnzymaticReaction']
    # for rxn in all_rxns[:20]:
    #     with model:
    #         model=saturate_enzymes2(model,rxnList=[rxn],sol=sol)
    #         print(rxn.id,model.slim_optimize())
    all_rxns = [reaction for reaction in model.reactions if type(reaction).__name__ == 'EnzymaticReaction']
    all_rxns = [reaction for reaction in all_rxns if reaction.id.startswith('r_')]
    # split all rxns into 10 groups
    rxn_groups=[all_rxns[i:i + 100] for i in range(0, len(all_rxns), 100)]
    sol = model.optimize()
    i=0
    for rxnlist in rxn_groups:
        with model:
            # model.optimize()
            ## set glucose uptake rate as 1 mmol/gDW/h
            model.reactions.get_by_id('r_1714').bounds=-1,-1
            # print('Before saturate enzymes:',sol.objective_value)
            model=saturate_enzymes(model,rxnList=rxnlist,sol=sol)
            print('Saturate %s group enzyme'%i,model.slim_optimize())
            i+=1
            # sol=model.optimize()




