# ETFLdesigner

Repo for strain design based on ETFL models.

## TODO📋

- [ ] Basic manipulation for ETFL model.

- [ ] Reproduce the ecFactory algorithm using Python.

- [ ] Build ecFactory algorithm for ETFL model.

- [ ] Build MOMA algorithm for ETFL model.

# Documentation

For documentation and API please check: 

# Installation
1. Clone the repo

```
git clone https://github.com/hongzhonglu/ETFLdesigner.git
```

2. Dependencies installation from Pypi

```
conda create -n env_name python=3.10
pip install -r requirements.txt
```

3. Solvers
- GUROBI
After installing GUROBI into the system, set the GUROBI Python API in your environment :
```
python [your_gurobi_path]/win64/setup.py install
```

- CPLEX
After installing CPLEX into the system, set the CPLEX Python API in your environment :
```
python [cplex_path]/python/setup.py install
```

# Contribution

* For contributors: Fork it to your Github account, and create a new branch from [`dev`](https://github.com/hongzhonglu/ETFLdesigner/tree/dev).


