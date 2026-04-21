"""Build heme a production model in ETFL by adding a heme a exchange reaction."""
from pathlib import Path
import sys
import cobra
from cobra import Reaction

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

from strainOptimizer.io import load_model
from strainOptimizer.simulation.utils import set_carbon_source_bounds

model = load_model(
    filename=str(PROJECT_ROOT / 'examples/models/yeast/GAN_ecYeast/GAN_all_v2.xml'),
    model_type='ecGEM'
)

for heme_id in ('s_0811_c', 's_0811[m]'):
    try:
        heme_met = model.metabolites.get_by_id(heme_id)
        break
    except KeyError:
        continue

# add heme a demand reaction to simulate heme a production
reaction = Reaction('EX_heme_a')
reaction.name = 'heme a production'
reaction.lower_bound = 0
reaction.upper_bound = 1000
reaction.add_metabolites({heme_met: -1})
model.add_reactions([reaction])

model.objective = 'EX_heme_a'
model.objective_direction = 'max'

c_source = 'r_1714'
c_uptake = 1
set_carbon_source_bounds(model, c_source, c_uptake, model_type='ecGEM', fixed=True)

sol = model.optimize()
print('Max heme a production:', sol.objective_value)

cobra.io.write_sbml_model(
    model,
    str(PROJECT_ROOT / 'examples/models/yeast/GAN_ecYeast/heme_GAN_all_v2.xml')
)
