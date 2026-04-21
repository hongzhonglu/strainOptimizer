# -*- coding: utf-8 -*-
"""iBridge strain design example for yeast 2-phenylethanol / spermidine production."""
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

from strainOptimizer.strainDesign import run_iBridge_design
from strainOptimizer.io import load_model
from cobra.io import read_sbml_model

solver = 'optlang-gurobi'

# load models — comment out models you don't need
models_dict = {
    # 'ecYeast': {
    #     'model': load_model(str(PROJECT_ROOT / 'examples/models/yeast/ecYeastGEM_batch.xml'),
    #                         model_type='ecGEM', solver=solver),
    #     'model_type': 'ecGEM',
    # },
    # 'efl': {
    #     'model': load_model(str(PROJECT_ROOT / 'examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json'),
    #                         model_type='etfl', solver=solver),
    #     'model_type': 'etfl',
    # },
    'gem': {
        'model': read_sbml_model(str(PROJECT_ROOT / 'examples/models/yeast/yeast-GEM.xml')),
        'model_type': 'GEM',
    },
}

methods = {'moma', 'mopa', 'ppfba'}
results = {}

for key, model_dict in models_dict.items():
    results[key] = {}
    model = model_dict['model']
    model_type = model_dict['model_type']

    if model_type == 'etfl':
        c_source = 'r_1714'
        cov_threshold = 0.01
        linear = True
    elif model_type == 'ecGEM':
        c_source = 'r_1714_REV'
        cov_threshold = 0.1
        linear = False
    elif model_type == 'GEM':
        c_source = 'r_1714'
        cov_threshold = 0.1
        linear = False

    c_uptake = 10
    # product = 'r_1589'  # 2-phenylethanol
    product = 'r_2051'    # spermidine

    for method in methods:
        with model:
            result = run_iBridge_design(
                model=model,
                targetID=product,
                c_source=c_source,
                c_uptake=c_uptake,
                model_type=model_type,
                method=method,
                tol=0.01,
                linear=linear,
                cov_threshold=cov_threshold,
            )
        results[key][method] = result

# result evaluation
from strainOptimizer.analysis.dataset import load_experiment_targets, calculate_exp_consistency

product_name = '2-phenylethanol' if product == 'r_1589' else 'spermidine'
df_exp = load_experiment_targets(product=product_name)

eval_results = {}
for model_key, items in results.items():
    eval_results[model_key] = {}
    for method, result in items.items():
        print(f'{model_key} {method}:')
        df_pred = result['endogenous_gene_result'].rename(columns={'gene_action': 'action'})
        exp_consistency = calculate_exp_consistency(df_pred, df_exp)
        eval_results[model_key][method] = exp_consistency
