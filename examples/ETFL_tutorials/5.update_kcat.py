# -*- coding: utf-8 -*-
"""
Update kcat values in an ETFL model and re-apply enzymatic constraints.
"""
from pathlib import Path
import sys

def _resolve_project_root() -> Path:
    """Support both script mode and interactive mode."""
    start = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd().resolve()
    for candidate in [start, *start.parents]:
        if (candidate / "src" / "strainOptimizer").exists():
            return candidate
    return start


PROJECT_ROOT = _resolve_project_root()
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
MODEL_PATH = str(PROJECT_ROOT / 'examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json')

from strainOptimizer.etfl.io.json import load_json_model
from strainOptimizer.etfl.optim.constraints import ForwardCatalyticConstraint, BackwardCatalyticConstraint


def reformulate_enzyme_expr(reaction):
    """Reformulate forward and backward catalytic constraint expressions."""
    v_max_fwd = {}
    v_max_bwd = {}
    for e, enz in enumerate(reaction.enzymes):
        v_max_fwd[e] = enz.kcat_fwd * enz.concentration
        v_max_bwd[e] = enz.kcat_bwd * enz.concentration

    k_f = max(x.kcat_fwd for x in reaction.enzymes)
    k_b = max(x.kcat_bwd for x in reaction.enzymes)
    E_m = max(x.scaling_factor for x in reaction.enzymes)

    fwd_expr = (reaction.forward_variable - sum(v_max_fwd.values())) / (k_f * E_m)
    bwd_expr = (reaction.reverse_variable - sum(v_max_bwd.values())) / (k_b * E_m)
    return fwd_expr, bwd_expr


def integrate_enzyme_constraint(model, reaction):
    """Replace the catalytic constraints for a reaction with updated kcat expressions."""
    new_expr_fwd, new_expr_bwd = reformulate_enzyme_expr(reaction)
    cons_fwd = ForwardCatalyticConstraint(reaction=reaction, expr=new_expr_fwd, ub=0, queue=True)
    cons_bwd = BackwardCatalyticConstraint(reaction=reaction, expr=new_expr_bwd, ub=0, queue=True)
    cons_fwd.change_expr(new_expr_fwd)
    cons_bwd.change_expr(new_expr_bwd)
    model._cons_dict[cons_fwd.name] = cons_fwd
    model._cons_dict[cons_bwd.name] = cons_bwd
    return model


def update_single_kcat(model, enzymeID, new_kcat_fwd, new_kcat_bwd):
    """Update kcat for one enzyme and rebuild its catalytic constraints."""
    enz = model.enzymes.get_by_id(enzymeID)
    enz.kcat_fwd = new_kcat_fwd
    enz.kcat_bwd = new_kcat_bwd

    for geneID in enz.composition.keys():
        gene = model.genes.get_by_id(geneID)
        for rxn in gene.reactions:
            if rxn.__class__.__name__ == 'EnzymaticReaction':
                integrate_enzyme_constraint(model, rxn)
    return model


if __name__ == '__main__':
    model = load_json_model(MODEL_PATH)
    sol = model.optimize()
    print('Baseline growth:', sol.objective_value)

    # Test kcat update for first enzyme only
    enz = list(model.enzymes)[0]
    print(f'Updating kcat for {enz.id}  (fwd={enz.kcat_fwd}, bwd={enz.kcat_bwd})')
    with model:
        model = update_single_kcat(model, enz.id, new_kcat_fwd=1, new_kcat_bwd=1)
        sol2 = model.optimize()
        print(f'After kcat update — growth: {sol2.objective_value}')
        print(f'New enzyme concentration: {model.enzymes.get_by_id(enz.id).X}')
