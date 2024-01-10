'''Modify the enzyme catalytic constraint to saturate all enzyme usage by making v=kcat*[E].'''
from strainOptimizer.io import load_model
from etfl.optim.constraints import ForwardCatalyticConstraint, BackwardCatalyticConstraint

def saturate_enzymes(model,rxnList=None,tol=1e-6):
    '''Saturate all catalytic enzymes in ETFL models'''
    if rxnList is None:
        rxnList = [reaction for reaction in model.reactions if type(reaction).__name__ == 'EnzymaticReaction']
    # get all catalytic constraints
    for rxn in rxnList:
        rxnID = rxn.id
        fwd_cons = model.forward_catalytic_constraint.get_by_id(rxnID)
        bwd_cons = model.backward_catalytic_constraint.get_by_id(rxnID)
        expr_fwd = fwd_cons.expr
        expr_bwd = bwd_cons.expr
        lb = -tol
        ub = 0
        model.remove_constraint(fwd_cons)
        new_fwd_cons = ForwardCatalyticConstraint(reaction=rxn, expr=expr_fwd, ub=ub, lb=lb, queue=True)
        model._cons_dict[new_fwd_cons.name] = new_fwd_cons

        # model.remove_constraint(bwd_cons)
        # new_bwd_cons = BackwardCatalyticConstraint(reaction=rxn, expr=expr_bwd, ub=ub, lb=lb, queue=True)
        # model._cons_dict[new_bwd_cons.name] = new_bwd_cons

    model.repair()
    return model


# strategy 2: add a new saturate enzyme constraint for each reations
from pytfa.optim import ReactionConstraint
class SaturateCatalyticConstraint(ReactionConstraint):
    """
    Class to represent a enzymatic constraint
    """

    prefix = 'SC_'

def saturate_enzymes2(model,rxnList=None,sol=None,tol=1e-6):
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
        if flux>=0:
            fwd_cons = model.forward_catalytic_constraint.get_by_id(rxnID)
            expr_fwd = fwd_cons.expr
            model.remove_constraint(fwd_cons)
            new_fwd_cons = ForwardCatalyticConstraint(reaction=rxn, expr=expr_fwd, ub=ub, lb=lb, queue=True)
            model._cons_dict[new_fwd_cons.name] = new_fwd_cons
        else:
            bwd_cons = model.backward_catalytic_constraint.get_by_id(rxnID)
            expr_bwd = bwd_cons.expr
            model.remove_constraint(bwd_cons)
            new_bwd_cons = BackwardCatalyticConstraint(reaction=rxn, expr=expr_bwd, ub=ub, lb=lb, queue=True)
            model._cons_dict[new_bwd_cons.name] = new_bwd_cons

    model.repair()
    return model


if __name__ == '__main__':
    # load model
    # test_model=load_model("examples/models/ecoli/ecoli_core_curated.json",model_type='etfl',solver='optlang-gurobi')
    # model=test_model

    model=load_model('examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json',model_type='etfl',solver='optlang-gurobi')

    sol=model.optimize()
    # model=load_model("examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json",model_type='etfl',solver='optlang-gurobi')
    all_rxns=[reaction for reaction in model.reactions if type(reaction).__name__ == 'EnzymaticReaction']
    for rxn in all_rxns[:20]:
        with model:
            model=saturate_enzymes2(model,rxnList=[rxn],sol=sol)
            print(rxn.id,model.slim_optimize())
    with model:
        # model.optimize()
        all_rxns = [reaction for reaction in model.reactions if type(reaction).__name__ == 'EnzymaticReaction']
        all_rxns= [reaction for reaction in all_rxns if reaction.id.startswith('r_')]
        model=saturate_enzymes2(model,rxnList=all_rxns[:50],sol=sol)
        print(model.slim_optimize())
        # sol=model.optimize()
    # model.reactions.get_by_id('r_1714').bounds=-1,-1
    model.slim_optimize()

    model.slim_optimize()



