"""
Basic tutorial for ETFL models bundled in strainOptimizer.

Model types:
  cEFL  — expression + stoichiometry, constant biomass composition
  vEFL  — expression + stoichiometry, variable biomass composition
  cETFL — cEFL + thermodynamic constraints
  vETFL — vEFL + thermodynamic constraints
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

from strainOptimizer.etfl.io.json import save_json_model, load_json_model
from strainOptimizer.etfl.optim.config import standard_solver_config
from strainOptimizer.etfl.optim.utils import safe_optim
from strainOptimizer.etfl.optim.variables import EnzymeVariable
from strainOptimizer.simulation.pprotFBA import ppFBA
from pytfa.optim.utils import symbol_sum

MODEL_PATH = str(PROJECT_ROOT / 'examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json')

# 1. Load and save model
model = load_json_model(MODEL_PATH, solver='optlang-gurobi')
# save_json_model(model, 'your_output_path.json')

# Basic solver configuration
standard_solver_config(model)

# 2. Basic simulations
target_id = 'r_2111'
carbon_source = 'r_1714'

model.objective = target_id
solution = safe_optim(model)
print('FBA growth:', solution.objective_value)

# ppFBA — minimize total enzyme usage subject to optimal objective
solution = ppFBA(model, target_id=target_id, c_source=carbon_source, c_uptake=1, model_type='etfl')
print('ppFBA growth:', solution.objective_value)

# Customized objective: minimize sum of fluxes for a reaction list
object_rxnlist = []  # fill with reaction IDs of interest
if object_rxnlist:
    expr = symbol_sum([
        model.reactions.get_by_id(x).forward_variable + model.reactions.get_by_id(x).reverse_variable
        for x in object_rxnlist
    ])
    model.objective = expr
    model.objective_direction = 'min'
    model.optimize()

# 3. Extract macromolecule ratios from solution
model.objective = target_id
safe_optim(model)
prot_ratio         = model.interpolation_variable.prot_ggdw.variable.primal
mrna_ratio         = model.interpolation_variable.mrna_ggdw.variable.primal
dna_ratio          = model.interpolation_variable.dna_ggdw.variable.primal
lipid_ratio        = model.interpolation_variable.lipid_ggdw.variable.primal
carbohydrate_ratio = model.interpolation_variable.carbohydrate_ggdw.variable.primal
ion_ratio          = model.interpolation_variable.ion_ggdw.variable.primal
print(f'protein: {prot_ratio:.4f}  mRNA: {mrna_ratio:.4f}  lipid: {lipid_ratio:.4f}')

# 4. Access enzyme / mRNA variables
enz_vars = model.get_variables_of_type(EnzymeVariable)
