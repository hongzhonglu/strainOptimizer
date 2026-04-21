'''
Expand the application of strainOptimizer to xylose-based 3-HP production in S. cerevisiae.

Case:
  3-Hydroxypropionate (3-HP) production on xylose in S. cerevisiae
  Model  : ecYeastGEM_batch.xml  (ecGEM)
  Pathway: Malonyl-CoA reductase (MCR) heterologous pathway
  Algorithm: ecFactory
  Experimental dataset source: 10.1038/s41467-025-59966-x
  Dataset: data/experiment_targets/yeast_xylose_3HP_exp_targets.tsv

Run prediction first:
  python paper_analysis_code/5.extend_usability_test_xylose.py
  (prediction and evaluation are combined in this script)
'''
from pathlib import Path
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from cobra import Reaction, Metabolite


def _resolve_project_root() -> Path:
    if '__file__' in globals():
        return Path(__file__).resolve().parents[1]
    cwd = Path.cwd().resolve()
    for candidate in [cwd, *cwd.parents]:
        if (candidate / 'src' / 'strainOptimizer').exists():
            return candidate
    return cwd


PROJECT_ROOT = _resolve_project_root()
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from strainOptimizer import strainOptimizer_engine, WorkflowParameters
from strainOptimizer.io import load_model
from strainOptimizer.analysis.dataset import load_experiment_targets, calculate_exp_consistency


PRODUCT = '3HP_xylose'
EXP_PRODUCT = 'yeast_xylose_3HP'
ALG = 'ecFactory'
RESULTS_DIR = PROJECT_ROOT / 'examples' / 'results'
OUTPUT_DIR = PROJECT_ROOT / 'paper_analysis_code' / 'results'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def add_mcr_3hp_pathway(model):
    """Add heterologous MCR pathway and 3-HP demand reaction to model."""
    try:
        hp3_c = model.metabolites.get_by_id('3hp_c')
    except KeyError:
        hp3_c = Metabolite('3hp_c', formula='C3H6O3', name='3-hydroxypropionate', compartment='c')

    if 'MCR_pathway' not in model.reactions:
        mcr_rxn = Reaction('MCR_pathway', name='Malonyl-CoA reductase pathway', lower_bound=0.0, upper_bound=1000.0)
        mcr_rxn.gene_reaction_rule = 'MCR'
        mcr_rxn.add_metabolites({
            model.metabolites.get_by_id('s_1101[c]'): -1.0,  # malonyl-CoA [cytoplasm]
            model.metabolites.get_by_id('s_1212[c]'): -2.0,  # NADPH [cytoplasm]
            model.metabolites.get_by_id('s_0794[c]'): -2.0,  # H+ [cytoplasm]
            hp3_c: 1.0,
            model.metabolites.get_by_id('s_0529[c]'): 1.0,   # CoA [cytoplasm]
            model.metabolites.get_by_id('s_1207[c]'): 2.0,   # NADP+ [cytoplasm]
        })
        model.add_reactions([mcr_rxn])

    if 'DM_3hp_c' not in model.reactions:
        model.add_boundary(hp3_c, type='demand', reaction_id='DM_3hp_c')

    return model


def fcc_consistent_subset(df):
    """Keep FCC-consistent targets: OE/FCCp>0 and KD|KO/FCCp<0."""
    if 'FCCp' not in df.columns:
        return pd.DataFrame(columns=df.columns)

    tmp = df[df['FCCp'].notna()].copy()
    if tmp.empty:
        return tmp

    tmp['FCCp'] = tmp['FCCp'].apply(lambda x: 0 if abs(x) < 1e-9 else x)
    oe_mask = (tmp['action'] == 'OE') & (tmp['FCCp'] > 0)
    kd_ko_mask = (tmp['action'].isin(['KD', 'KO'])) & (tmp['FCCp'] < 0)
    return tmp[oe_mask | kd_ko_mask]


# 1. Load model and add 3-HP heterologous pathway
print('Loading ecYeastGEM model...')
model = load_model(
    filename=str(PROJECT_ROOT / 'examples/models/yeast/ecYeastGEM_batch.xml'),
    model_type='ecGEM',
    solver='optlang-gurobi',
)
model.tolerance = 1e-9

print('Adding MCR (Malonyl-CoA reductase) 3-HP pathway...')
model = add_mcr_3hp_pathway(model)

# 2. Switch carbon source to xylose
model.reactions.get_by_id('r_1714_REV').bounds = (0.0, 0.0)   # glucose exchange (REV)
model.reactions.get_by_id('r_1718_REV').bounds = (0.0, 10.0)  # D-xylose exchange (REV)

with model:
    model.objective = 'DM_3hp_c'
    max_3hp = model.slim_optimize()
    print(f'Max 3-HP production on xylose: {max_3hp:.4f} mmol/gDW/h')
    model.objective = 'r_2111'
    max_gr = model.slim_optimize()
    print(f'Max growth on xylose: {max_gr:.4f} h-1')

# 3. Run ecFactory design
model_params = {
    'model_path': str(PROJECT_ROOT / 'examples/models/yeast/ecYeastGEM_batch.xml'),
    'model_type': 'ecGEM',
    'solver': 'optlang-gurobi',
    'growth_id': 'r_2111',
}

strain_params = {
    'target_id': 'DM_3hp_c',
    'product_name': PRODUCT,
    'c_source': 'r_1718_REV',
    'c_uptake': 5,
    'substrate_MW': 0.150130,
}

algorithm_params = {
    'design_algorithm': 'ecFactory',
    'simulation_method': 'ppfba',
    'remove_essential': True,
    'steps': 123,
    'action_thresholds': [0.05, 0.3, 1.1],
    'calculate_fcc': True,
    'output_directory': str(RESULTS_DIR),
    'save_results': True,
}

params = WorkflowParameters(
    model=model_params,
    strain=strain_params,
    algorithm=algorithm_params,
)

engine = strainOptimizer_engine(params)
engine.model = model

print(f'\nProduct     : {strain_params["product_name"]}')
print(f'Target rxn  : {strain_params["target_id"]}')
print(f'Carbon source: {strain_params["c_source"]}  ({strain_params["c_uptake"]} mmol/gDW/h)')
print(f'Algorithm   : {algorithm_params["design_algorithm"]}')
print(f'Output dir  : {RESULTS_DIR}')

print('\nRunning ecFactory workflow...')
final_result = engine.run_design()

print('\n--- Final gene targets (geneTable) ---')
print(final_result.to_string())

all_results = engine.all_results
if 'level 3 result' in all_results:
    print('\n--- Level 3 result (minimal candidate set) ---')
    print(all_results['level 3 result'].to_string())

# 4. Collect predicted results (xlsx first, fall back to in-memory)
prefix = f'{PRODUCT}_design_results_{ALG}'
_xlsx_keys = {
    'Level 1':      RESULTS_DIR / f'{prefix}_level_1_result.xlsx',
    'Level 2':      RESULTS_DIR / f'{prefix}_level_2_result.xlsx',
    'Level 3':      RESULTS_DIR / f'{prefix}_level_3_result.xlsx',
    'FCC-filtered': RESULTS_DIR / f'{prefix}_fcc_filtered_result.xlsx',
}
_memory_keys = {
    'Level 1':      'level 1 result',
    'Level 2':      'level 2 result',
    'Level 3':      'level 3 result',
    'FCC-filtered': 'fcc_filtered_result',
}

predict_levels = {}
for label, path in _xlsx_keys.items():
    if path.exists():
        predict_levels[label] = pd.read_excel(path, index_col=0)
        print(f'Loaded {label}: {len(predict_levels[label])} genes')
    elif _memory_keys[label] in all_results and all_results[_memory_keys[label]] is not None:
        predict_levels[label] = all_results[_memory_keys[label]]
        print(f'[in-memory] {label}: {len(predict_levels[label])} genes')
    else:
        print(f'[SKIP] {label}: not available')

# 5. Load experimental dataset and evaluate
print(f'\nLoading experimental dataset for: {EXP_PRODUCT}')
exp_data = load_experiment_targets(product=EXP_PRODUCT)
print(f'Experimental targets ({len(exp_data)}):')
print(exp_data.to_string())

print('\n' + '=' * 60)
print(' Prediction vs Experiment Consistency  (KO+KD merged)')
print('=' * 60)

MERGE_KO_KD = True
summary_rows = []
consistency_by_level = {}

for label, pred_df in predict_levels.items():
    print(f'\n--- {label} ({len(pred_df)} predicted genes) ---')
    consistency = calculate_exp_consistency(
        predict_result=pred_df,
        exp_data=exp_data,
        show=False,
        merge_ko_kd=MERGE_KO_KD,
    )
    if consistency is None:
        print('  No predicted targets match any action category.')
        continue
    consistency_by_level[label] = consistency

    overall = consistency['overall']
    print(f"  Predicted  : {overall['predict_num']}")
    print(f"  Exp targets: {overall['exp_num']}")
    print(f"  Hits       : {overall['hit_num']}")
    print(f"  Recall     : {overall['consistency']:.3f}  "
          f"({overall['hit_num']}/{overall['exp_num']})")
    print(f"  Precision  : {overall['precision']:.3f}  "
          f"({overall['hit_num']}/{overall['predict_num']})")

    action_labels = ('KD', 'OE') if MERGE_KO_KD else ('KO', 'KD', 'OE')
    for action in action_labels:
        if action in consistency:
            c = consistency[action]
            tag = 'KO+KD' if (MERGE_KO_KD and action == 'KD') else action
            print(f"  {tag:5s}  exp={c['exp_num']:2d}  "
                  f"hit={c['hit_num']:2d}  "
                  f"recall={c['consistency']:.2f}  "
                  f"hits: {c['hit']}")

    f1 = (2 * overall['consistency'] * overall['precision'] /
          (overall['consistency'] + overall['precision'])
          if (overall['consistency'] + overall['precision']) > 0 else 0)

    summary_rows.append({
        'Level': label,
        'Predicted': overall['predict_num'],
        'Exp': overall['exp_num'],
        'Hits': overall['hit_num'],
        'Recall': round(overall['consistency'], 4),
        'Precision': round(overall['precision'], 4),
        'F1': round(f1, 4),
    })

summary_df = pd.DataFrame(summary_rows)
if not summary_df.empty:
    summary_df = summary_df.set_index('Level')
else:
    summary_df = pd.DataFrame(columns=['Predicted', 'Exp', 'Hits', 'Recall', 'Precision', 'F1'])

print('\n' + '=' * 60)
print(' Summary Table')
print('=' * 60)
print(summary_df.to_string())

# 7. Hit gene details
print('\n' + '=' * 60)
print(' Hit gene details (all levels)')
print('=' * 60)

for label, consistency in consistency_by_level.items():
    hit_genes = []
    for action in ('KD', 'OE') if MERGE_KO_KD else ('KO', 'KD', 'OE'):
        if action in consistency:
            tag = 'KO+KD' if (MERGE_KO_KD and action == 'KD') else action
            for g in consistency[action]['hit']:
                hit_genes.append({'gene': g, 'action': tag})
    if hit_genes:
        hit_df = pd.DataFrame(hit_genes).set_index('gene')
        hit_df = hit_df.join(exp_data[['shortNames', 'enzymes']], how='left')
        print(f'\n{label} hits:')
        print(hit_df.to_string())

# 8. Save summary
out_xlsx = OUTPUT_DIR / f'{PRODUCT}_consistency_summary.xlsx'
with pd.ExcelWriter(out_xlsx) as writer:
    summary_df.to_excel(writer, sheet_name='Summary')
    exp_data.to_excel(writer, sheet_name='Exp_targets')
    for label, pred_df in predict_levels.items():
        pred_df.to_excel(writer, sheet_name=label.replace(' ', '_')[:31])
print(f'\nSummary saved to {out_xlsx}')

# 9. Bar chart
if not summary_df.empty:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    levels = summary_df.index.tolist()
    x = np.arange(len(levels))
    width = 0.35

    # colour: L1/L2/L3 in blue/orange, FCC-filtered in green
    bar_colors_r = ['#4C72B0' if 'FCC' not in lv else '#2ca02c' for lv in levels]
    bar_colors_p = ['#DD8452' if 'FCC' not in lv else '#98df8a' for lv in levels]

    ax = axes[0]
    ax.bar(x - width / 2, summary_df['Recall'],    width, color=bar_colors_r, label='Recall')
    ax.bar(x + width / 2, summary_df['Precision'], width, color=bar_colors_p, label='Precision', alpha=0.85)
    ax.set_xticks(x); ax.set_xticklabels(levels, rotation=15, ha='right')
    ax.set_ylim(0, 1.15); ax.set_ylabel('Score')
    ax.set_title(f'{PRODUCT} ({ALG}) – Recall & Precision')
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color='#4C72B0', label='Recall'),
                       Patch(color='#DD8452', label='Precision')], fontsize=8)
    for i, row in enumerate(summary_df.itertuples()):
        ax.text(i - width / 2, row.Recall + 0.02,    f'{row.Recall:.2f}',    ha='center', fontsize=8)
        ax.text(i + width / 2, row.Precision + 0.02, f'{row.Precision:.2f}', ha='center', fontsize=8)

    ax2 = axes[1]
    exp_total = summary_df['Exp'].iloc[0]
    hits   = summary_df['Hits'].values
    misses = exp_total - hits
    hit_colors = ['#55A868' if 'FCC' not in lv else '#2ca02c' for lv in levels]
    ax2.bar(levels, hits,   color=hit_colors, label='Hits')
    ax2.bar(levels, misses, bottom=hits, color='#C44E52', alpha=0.6, label='Missed')
    ax2.set_ylabel('Number of experimental targets')
    ax2.set_title(f'{PRODUCT} – Hits vs Missed (n={exp_total})')
    ax2.set_xticklabels(levels, rotation=15, ha='right')
    ax2.legend(fontsize=8)
    for i, h in enumerate(hits):
        if h > 0:
            ax2.text(i, h / 2, str(int(h)), ha='center', va='center',
                     color='white', fontweight='bold', fontsize=9)

    plt.tight_layout()
    out_fig = OUTPUT_DIR / f'{PRODUCT}_consistency_plot.png'
    plt.savefig(out_fig, dpi=150)
    print(f'Plot saved to {out_fig}')
    plt.show()
