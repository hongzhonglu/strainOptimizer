# -*- coding: utf-8 -*-
"""
Predicting targets for tryptophan production in E. coli using the ecFactory algorithm in strainOptimizer.
"""
from pathlib import Path
import sys

def _resolve_project_root() -> Path:
    """Support both script mode and interactive mode."""
    if "__file__" in globals():
        return Path(__file__).resolve().parents[1]

    cwd = Path.cwd().resolve()
    for candidate in [cwd, *cwd.parents]:
        if (candidate / "src" / "strainOptimizer").exists():
            return candidate
    return cwd


PROJECT_ROOT = _resolve_project_root()
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from strainOptimizer import strainOptimizer_engine, WorkflowParameters

OUTPUT_DIR = PROJECT_ROOT / 'examples' / 'results'

model_params = {
    'model_path': str(PROJECT_ROOT / 'examples/models/ecoli/eciML1515_batch.xml'),
    'model_type': 'ecGEM',
    'solver': 'optlang-gurobi',
    # 'solver': 'optlang-cplex',
    'growth_id': 'BIOMASS_Ec_iML1515_core_75p37M',
}

strain_params = {
    'target_id': 'EX_trp__L_e',
    'product_name': 'Tryptophan',
    'c_source': 'EX_glc__D_e_REV',
    'c_uptake': 2,           # mmol/gDW/h
    'substrate_MW': 0.180156,  # g/mmol, glucose
}

algorithm_params = {
    'design_algorithm': 'ecFactory',
    'simulation_method': 'ppfba',
    'remove_essential': False,
    'steps': 123,
    'calculate_fcc': True,
    # 'scanning_range': [0.1, 0.3],
    # 'experimental_yield': 0.1,
    'output_directory': str(OUTPUT_DIR),
    'save_results': True,
}

params = WorkflowParameters(
    model=model_params,
    strain=strain_params,
    algorithm=algorithm_params
)

print(f"Product     : {params.strain['product_name']}")
print(f"Target rxn  : {params.strain['target_id']}")
print(f"Carbon source: {params.strain['c_source']}  ({params.strain['c_uptake']} mmol/gDW/h)")
print(f"Model type  : {params.model['model_type']}")
print(f"Algorithm   : {params.algorithm['design_algorithm']}")
print(f"Output dir  : {OUTPUT_DIR}")

engine = strainOptimizer_engine(params)

model = engine.load_model()
model.tolerance = 1e-9

model_info = engine.get_model_info()
print(f"\nModel info: {model_info}")

print("\nRunning strain design workflow...")
final_result = engine.run_design()

# Print summary
print("\n========== Design Summary ==========")
summary = engine.get_results_summary()
for k, v in summary.items():
    print(f"  {k}: {v}")

print("\n--- Final gene targets (geneTable) ---")
print(final_result.to_string())

all_results = engine.all_results
if 'level 3 result' in all_results:
    print("\n--- Level 3 result (minimal candidate set) ---")
    print(all_results['level 3 result'].to_string())
elif 'level 2 result' in all_results:
    print("\n--- Level 2 result ---")
    print(all_results['level 2 result'].to_string())
