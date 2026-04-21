# -*- coding: utf-8 -*-
"""
Add metabolic reactions to an ETFL model.
Demonstrates adding itaconate biosynthesis pathway to yeast cEFL model.
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
from cobra import Reaction, Metabolite

solver = 'optlang-gurobi'
growth_reaction_id = 'r_4041'

model = load_json_model(MODEL_PATH, solver=solver)
model.reactions.r_1714.lower_bound = -10
model.reactions.r_1714.upper_bound = 0
model.objective = growth_reaction_id
model.objective_direction = 'max'
sol = model.optimize()
print('Baseline growth:', sol.objective_value)

# Add itaconate pathway
model.add_metabolites(Metabolite('s_4349_c', compartment='c', formula='C5H4O4', name='itaconate'))
model.add_metabolites(Metabolite('s_4350_e', compartment='e', formula='C5H4O4', name='itaconate'))

# cis-aconitate decarboxylase: cis-aconitate[c] → itaconate[c] + CO2[c]
rxn_cad = Reaction('new_cis_aconitate_decarboxylase')
rxn_cad.name = 'cis-aconitate decarboxylase'
rxn_cad.lower_bound = 0.
rxn_cad.upper_bound = 1000.
rxn_cad.add_metabolites({
    model.metabolites.get_by_id('s_0516_c'): -1,
    model.metabolites.get_by_id('s_4349_c'): 1,
    model.metabolites.get_by_id('s_0456_c'): 1,
})
rxn_cad.gene_reaction_rule = 'AtCAD'

# itaconate transport: itaconate[c] → itaconate[e]
rxn_transport = Reaction('new_itaconate_transport')
rxn_transport.lower_bound = 0.
rxn_transport.upper_bound = 1000.
rxn_transport.add_metabolites({
    model.metabolites.get_by_id('s_4349_c'): -1,
    model.metabolites.get_by_id('s_4350_e'): 1,
})

# demand reaction: itaconate[e] →
rxn_exchange = Reaction('exchange_itaconate')
rxn_exchange.lower_bound = 0.
rxn_exchange.upper_bound = 1000.
rxn_exchange.add_metabolites({model.metabolites.get_by_id('s_4350_e'): -1})

model.add_reactions([rxn_cad, rxn_transport, rxn_exchange])

model.objective = 'exchange_itaconate'
model.objective_direction = 'max'
model.reactions.r_4041.lower_bound = 0.2
sol = model.optimize()
print('Itaconate production:', sol.objective_value)
print('Glucose uptake:', sol.fluxes['r_1714'])
