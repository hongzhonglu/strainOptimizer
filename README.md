# strainOptimizer

Repo for strain design based on ETFL/ecGEM models.

## Installation
1. Clone the repo

```
git clone https://github.com/hongzhonglu/strainOptimizer.git
```

2. Set up the environment

```
conda create -n strainOpitimizer python=3.10
conda activate strainOpitimizer
conda env update --file environment.yml
pip install -e .
```

3. Solvers
- GUROBI (version 10.0)

We recommend obtaining a valid Gurobi license. Otherwise, adjust the code to use a different optimizer. To retrieve your license::
```
grbgetkey <YOUR-LICENSE-KEY>
```
For details, see the Gurobi licensing guide: https://www.gurobi.com/documentation/

- CPLEX
After installing CPLEX into the system, set the CPLEX Python API in your environment :
```
python [cplex_path]/python/setup.py install
```
## Usage guide
1. [ecFactory design example](examples/1.ecFactory_design_example.ipynb)
2. [ecFSEOF design example](examples/2.ecFSEOF_design_example.ipynb)

More examples can be found in the [examples](examples) folder.
## Example script
Here is an example script for 2-PE strain design:
```python
from strainOptimizer import strainOptimizer_engine,WorkflowParameters

# Model parameters - model path, type, solver, growth reaction
model_params = {
    'model_path': 'example/models\yeast\ecYeastGEM_batch.xml',
    'model_type': 'ecGEM',
    'solver': 'optlang-gurobi',
    'growth_id': 'r_2111',
}
# Strain parameters - target product and growth conditions
strain_params = {
    'target_id': 'r_1589',
    'product_name': '2-phenylethanol',
    'c_source': 'r_1714_REV',  # glucose exchange reaction
    'c_uptake': 5,  # glucose uptake rate (mmol/gDW/h)
}

# Algorithm control parameters - workflow and output settings
algorithm_params = {
    'design_algorithm': 'ecFactory',
    'simulation_method': 'ppfba',
    'experimental_yield': None, # if without experimental yield data, use the 1/2 
    'remove_essential': True,
    'output_directory': './results',
    'steps':123,
    'action_thresholds':[0.05,0.3,1.1]
    # 'save_results': False,
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

# Create workflow engine using the new framework
engine = strainOptimizer_engine(params)
# Load model
engine.load_model()
# Run the design workflow
engine.run_design()
```
## Contribution

* For contributors: Fork it to your Github account, and create a new branch from [`dev`](https://github.com/hongzhonglu/strainOptimizer/tree/dev).


