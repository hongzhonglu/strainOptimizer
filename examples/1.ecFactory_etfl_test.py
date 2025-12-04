# -*- coding: utf-8 -*-
import os
import sys
sys.path.append(r'D:\code\github\strainOptimizer')
os.chdir(r'D:\code\github\strainOptimizer')
from strainOptimizer import strainOptimizer_engine,WorkflowParameters

# set tolerance
import cobra
# cobra.Configuration().tolerance=1e-9

model_params = {
    'model_path': 'examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json',
    'model_type': 'etfl',
    'solver': 'optlang-gurobi',
    # 'solver': 'optlang-cplex',
    'growth_id': 'r_4041',
    # 'total_enzymes': 0.1
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
    'output_directory': './results',
    'save_results': False,
    'steps': 123,
    'simulation_method': 'ppfba',
    # 'simulation_method': 'pfba',
    # 'scanning_range':[0.1,0.4],
    # 'experimental_yield':0.16,
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

engine = strainOptimizer_engine(params)

print(f"Engine created for {params.strain['product_name']} production")
print(f"Target reaction: {params.strain['target_id']}")
print(f"Carbon source: {params.strain['c_source']}")
print(f"Model type: {params.model['model_type']}")
print(f"Algorithm: {params.algorithm['design_algorithm']}")

# Load model
model = engine.load_model()
model.tolerance=1e-7
model.solver.configuration.timeout = 1200
# model.solver.configuration.tolerances.optimality = 1e-4
# model.solver.problem.Params.MIPFocus = 1

# Get model information
model_info = engine.get_model_info()
print(f"\nModel info: {model_info}")

# Run the design workflow
print("\nRunning strain design workflow...")
results = engine.run_design()