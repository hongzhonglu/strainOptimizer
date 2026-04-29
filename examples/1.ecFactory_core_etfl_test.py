# -*- coding: utf-8 -*-
"""
Test ecFactory algorithm for succinate production in E. coli using an ETFL core model.
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

from strainOptimizer import strainOptimizer_engine, WorkflowParameters

model_params = {
    'model_path': str(PROJECT_ROOT / 'examples/models/ecoli/ecoli_core_curated.json'),
    'model_type': 'etfl',
    'solver': 'optlang-gurobi',
    # 'solver': 'optlang-cplex',
    'growth_id': 'Biomass_Ecoli_core',
}

strain_params = {
    'target_id': 'EX_succ_e',
    'product_name': 'succinate',
    'c_source': 'EX_glc__D_e',
    'c_uptake': 1,
}

algorithm_params = {
    'design_algorithm': 'ecFactory',
    'remove_essential': False,
    'output_directory': str(PROJECT_ROOT / 'results'),
    'save_results': True,
    'steps': 123,
    'simulation_method': 'ppfba',
    # 'scanning_range': [0.1, 0.3],
    # 'experimental_yield': 0.1,
}

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

model = engine.load_model()
model.tolerance = 1e-9

model_info = engine.get_model_info()
print(f"\nModel info: {model_info}")

print("\nRunning strain design workflow...")
final_result = engine.run_design()

print("\n========== Design Summary ==========")
summary = engine.get_results_summary()
for k, v in summary.items():
    print(f"  {k}: {v}")

print("\n--- Final gene targets (geneTable) ---")
print(final_result.to_string())
