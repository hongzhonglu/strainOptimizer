"""Change the total protein pool constraint in an ETFL or ecGEM model."""
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
from strainOptimizer.manipulation.constraint.total_resource_allocation import constrain_enzymes

# Load model
model = load_model(
    filename=str(PROJECT_ROOT / 'examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json'),
    model_type='etfl',
    solver='optlang-gurobi',
)

print('Baseline growth:', model.slim_optimize())

enzymes_ratio = 0.1
constrain_enzymes(model, total_prot=enzymes_ratio, model_type='etfl')

print('Growth after enzyme constraint:', model.slim_optimize())
