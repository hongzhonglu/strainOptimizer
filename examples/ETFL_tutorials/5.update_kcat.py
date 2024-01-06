# -*- coding: utf-8 -*-
# date : 2023/2/8 
# author : wangh
# file : 8.update_kcat.py
# project : etfl
from etfl.io.json import load_json_model
import os
os.chdir(r"D:\code\github\etfl")

# 方案2:参考ETFL原函数：yeast.apply_enzyme_catalytic_constraint(new_constraint)
def update_single_kcat(model,enzymeID, new_kcat_bwd, new_kcat_fwd):
    '''update the kcat of a single target enzyme and update the related enzymatic constraints'''

    # update the kcat of the target enzyme in enzyme class variable
    enz=model.enzymes.get_by_id(enzymeID)
    old_kcat_bwd=enz.kcat_bwd
    old_kcat_fwd=enz.kcat_fwd
    enz.kcat_bwd=new_kcat_bwd
    enz.kcat_fwd=new_kcat_fwd
    print("old_kcat_bwd:",old_kcat_bwd)
    print("old_kcat_fwd:",old_kcat_fwd)
    print("new_kcat_bwd:",new_kcat_bwd)
    print("new_kcat_fwd:",new_kcat_fwd)

    # get all reactions that including the target enzyme
    geneID_List=[x for x in enz.composition.keys()]
    reationList=[]
    for geneID in geneID_List:
        gene=model.genes.get_by_id(geneID)
        # only update EnzymaticReaction constraint
        for reaction in gene.reactions:
            if reaction.__class__.__name__ == 'EnzymaticReaction':
                reationList.append(reaction)
        # print("%s rxns will update enzymatic constraint"%len(reationList))
        # update all reactions constraints involving the target enzyme
        # reationList=[x for x in enz.reactions]

    #update the related enzymatic constraint
    for reaction in reationList:
        integrate_enzyme_constraint(model, reaction)

    return model


def refoumulate_enzyme_expr(reaction):
    '''refoumulate the enzyme catalytic constraint expression of the target reaction.
    refer to: model.apply_enzyme_catalytic_constraint(reaction) in ETFL package.
    parameter:
    reaction: EnzymaticReaction class in ETFL model
    return: new enzyme catalytic constraint expressions:bwd_expr,fwd_expr'''
    v_max_fwd = dict()
    v_max_bwd = dict()

    # Write v_max constraint
    fwd_variable = reaction.forward_variable
    bwd_variable = reaction.reverse_variable

    for e, enz in enumerate(reaction.enzymes):
        # If the enzymes has the same kcat for both directions
        # v_fwd  <= kcat_fwd [E]
        # v_fwd - kcat_fwd [E] <= 0

        v_max_fwd[e] = enz.kcat_fwd * enz.concentration
        v_max_bwd[e] = enz.kcat_bwd * enz.concentration

    # Formulating the scaling factor on the max kcat
    k_f = max([x.kcat_fwd for x in reaction.enzymes])
    k_b = max([x.kcat_bwd for x in reaction.enzymes])
    # k_f = np.median([x.kcat_fwd for x in self.enzymes])
    # k_b = np.median([x.kcat_bwd for x in self.enzymes])
    E_m = max([x.scaling_factor for x in reaction.enzymes])

    # v_fwd <= sum(kcat_i*E_i)
    # for all i, E_i <= E_max (= 1g/gDW)
    # v_fwd / sum(kcat_i*E_i^max) <= sum(kcat_i*E_i) / sum(kcat_i*E_i^max) (<= 1)
    enz_constraint_expr_fwd = (fwd_variable - sum(v_max_fwd.values())) / (k_f * E_m)
    enz_constraint_expr_bwd = (bwd_variable - sum(v_max_bwd.values())) / (k_b * E_m)
    return enz_constraint_expr_fwd, enz_constraint_expr_bwd


def integrate_enzyme_constraint(model, reaction):
    '''integrate the enzyme catalytic constraint of the target reaction into the ETFL model.
    parameter:
    model: ETFL model
    reaction: EnzymaticReaction class in ETFL model'''
    from etfl.optim.constraints import ForwardCatalyticConstraint, BackwardCatalyticConstraint

    # formulate constraint expression
    new_expr_fwd, new_expr_bwd = refoumulate_enzyme_expr(reaction)

    # build the new enzyme catalytic constraint
    cons_fwd = ForwardCatalyticConstraint(reaction=reaction, expr=new_expr_fwd, ub=0, queue=True)
    cons_bwd = BackwardCatalyticConstraint(reaction=reaction, expr=new_expr_bwd, ub=0, queue=True)
    cons_fwd.change_expr(new_expr_fwd)
    cons_bwd.change_expr(new_expr_bwd)

    # add new constraints to model
    model._cons_dict[cons_fwd.name] = cons_fwd
    model._cons_dict[cons_bwd.name] = cons_bwd
    model.logger.info('Added constraint: {}'.format(cons_bwd.expr))
    model.logger.info('Added constraint: {}'.format(cons_fwd.expr))
    return model


# test
if __name__ == '__main__':
    # test_model=load_json_model("models/ME_ecoli_coreModel.json")
    yetfl = load_json_model("yetfl/models/yeast8_cEFL_2584_enz_128_bins__20221228_090737.json")
    model=yetfl
    sol=model.optimize()
    enzymeList = model.enzymes
    for enz in enzymeList:
        print(enz.id)
        print("old_enz_conc:", enz.X)
        with model:
            model = update_single_kcat(model=model, enzymeID=enz.id, new_kcat_bwd=1, new_kcat_fwd=1)
            sol2= model.optimize()
            print("sol2:", sol2.objective_value)
            # enzyme concentration
            print("new_enz_conc:", model.enzymes.get_by_id(enz.id).X)
