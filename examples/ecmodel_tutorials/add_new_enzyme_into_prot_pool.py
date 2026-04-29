"""Add a new enzyme into the protein pool constraint of an ecGEM model."""
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

from strainOptimizer.io import load_model
from cobra import Metabolite, Reaction

model = load_model(
    str(PROJECT_ROOT / 'examples/models/yeast/ecYeastGEM_batch.xml'),
    model_type='ecGEM',
)

met_prot_pool = model.metabolites.get_by_id('prot_pool[c]')
exchange_prot_pool = model.reactions.get_by_id('prot_pool_exchange')

# 1. Add new enzyme metabolite
new_enzyme = Metabolite('new_enzyme[c]', name='new enzyme', compartment='c')
model.add_metabolites(new_enzyme)
new_MW = 70  # kDa → g/mmol

# 2. Add draw_protein reaction to connect enzyme to the protein pool
new_draw_rxn = Reaction('draw_new_enzyme')
new_draw_rxn.add_metabolites({new_enzyme: 1, met_prot_pool: -new_MW})
model.add_reactions([new_draw_rxn])

print('Model reactions:', len(model.reactions))
print('draw_new_enzyme bounds:', new_draw_rxn.bounds)
