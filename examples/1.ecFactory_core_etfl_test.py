# -*- coding: utf-8 -*-
import os
import sys
# sys.path.append(r'D:\code\github\strainOptimizer')
# sys.path.append(r'D:\code\github\strainOptimizer/strainOptimizer')
os.chdir(r'D:\code\github\strainOptimizer')
from strainOptimizer import strainOptimizer_engine,WorkflowParameters

# set tolerance
import cobra
# cobra.Configuration().tolerance=1e-9

model_params = {
    'model_path': 'examples/models/ecoli/ecoli_core_curated.json',
    'model_type': 'etfl',
    'solver': 'optlang-gurobi',
    # 'solver': 'optlang-cplex',
    'growth_id': 'Biomass_Ecoli_core',
    # 'total_enzymes': 0.1
}

# Strain parameters - target product and growth conditions
strain_params = {
    'target_id': 'EX_succ_e',
    'product_name': 'succinate',
    'c_source': 'EX_glc__D_e',  # glucose exchange reaction
    'c_uptake': 1,  # glucose uptake rate (mmol/gDW/h)
}

# Algorithm control parameters - workflow and output settings
algorithm_params = {
    'design_algorithm': 'ecFactory',
    'remove_essential': False,
    'output_directory': './results',
    'save_results': True,
    'steps': 123,
    'simulation_method': 'ppfba',
    # 'scanning_range':[0.1,0.3]
    # 'experimental_yield':0.1,
    # 'only_final_result': True,
    # Note: ecFactory-specific parameters like steps, action_thresholds, etc.
    # would need to be added to AlgorithmControl if they're used
}

# Create WorkflowParameters using the three-level structure
params = WorkflowParameters(
    model=model_params,
    strain=strain_params,
    algorithm=algorithm_params
)
parameters=params
engine = strainOptimizer_engine(params)

print(f"Engine created for {params.strain['product_name']} production")
print(f"Target reaction: {params.strain['target_id']}")
print(f"Carbon source: {params.strain['c_source']}")
print(f"Model type: {params.model['model_type']}")
print(f"Algorithm: {params.algorithm['design_algorithm']}")

# Load model
model = engine.load_model()
model.tolerance=1e-9
g = model.solver.problem

# Get model information
model_info = engine.get_model_info()
print(f"\nModel info: {model_info}")

# Run the design workflow
print("\nRunning strain design workflow...")
results = engine.run_design()