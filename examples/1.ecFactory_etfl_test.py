# -*- coding: utf-8 -*-
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

from strainOptimizer import strainOptimizer_engine, WorkflowParameters

model_params = {
    'model_path': str(PROJECT_ROOT / 'examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json'),
    'model_type': 'etfl',
    'solver': 'optlang-gurobi',
    # 'solver': 'optlang-cplex',
    'growth_id': 'r_4041',
}

# Strain parameters - target product and growth conditions
strain_params = {
    'target_id': 'r_1589',
    'product_name': '2-phenylethanol',
    'c_source': 'r_1714',  # glucose exchange reaction
    'c_uptake': 1,  # glucose uptake rate (mmol/gDW/h)
}

# Algorithm control parameters - workflow and output settings
algorithm_params = {
    'design_algorithm': 'ecFactory',
    'remove_essential': True,
    'output_directory': str(PROJECT_ROOT / 'results'),
    'save_results': False,
    'steps': 123,
    'simulation_method': 'ppfba',
    # 'simulation_method': 'pfba',
    # 'scanning_range': [0.1, 0.4],
    # 'experimental_yield': 0.16,
}

# Create WorkflowParameters using the three-level structure
params = WorkflowParameters(
    model=model_params,
    strain=strain_params,
    algorithm=algorithm_params
)

engine = strainOptimizer_engine(params)

print(f"Engine created for {params.strain['product_name']} production")
print(f"Target reaction: {params.strain['target_id']}")
print(f"Carbon source: {params.strain['c_source']}")
print(f"Model type: {params.model['model_type']}")
print(f"Algorithm: {params.algorithm['design_algorithm']}")

# Load model
model = engine.load_model()
model.tolerance = 1e-7
model.solver.configuration.timeout = 1200

# Get model information
model_info = engine.get_model_info()
print(f"\nModel info: {model_info}")

# Run the design workflow
print("\nRunning strain design workflow...")
final_result = engine.run_design()

print("\n========== Design Summary ==========")
summary = engine.get_results_summary()
for k, v in summary.items():
    print(f"  {k}: {v}")

print("\n--- Final gene targets (geneTable) ---")
print(final_result.to_string())
