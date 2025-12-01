# strainOptimizer

Repo for strain design based on ETFL/ecGEM models.

## Installation
1. Clone the repo

```
git clone https://github.com/hongzhonglu/strainOptimizer.git
```

2. Dependencies installation from Pypi

```
conda create -n strainOpitimizer python=3.10
conda activate strainOpitimizer
pip install -r requirements.txt
```

3. Solvers
- GUROBI (version 10.0)
After installing GUROBI into the system, set the GUROBI Python API in your environment :
```
python [your_gurobi_path]/setup.py install
```

- CPLEX
After installing CPLEX into the system, set the CPLEX Python API in your environment :
```
python [cplex_path]/python/setup.py install
```
## Usage guide
1. [ecFactory design example](examples/1.ecFactory_design_example.ipynb)
2. [ecFSEOF design example](examples/2.ecFSEOF_design_example.ipynb)

More examples can be found in the [examples](examples) folder.

## Contribution

* For contributors: Fork it to your Github account, and create a new branch from [`dev`](https://github.com/hongzhonglu/strainOptimizer/tree/dev).


