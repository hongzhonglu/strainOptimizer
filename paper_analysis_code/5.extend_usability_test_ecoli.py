'''
Expand the application of strainOptimzier to E.coli cell design
example:
Tropotophan production in E. coli
model: eciML1515_batch.xml
algorithm: ecFactory
experiemtnal dataset source: 10.1101/2025.06.12.659423
dataset: data/experiment_targets/ecoli_tryptophan_exp_targets.tsv

run the prediction by: python examples/1.ecFactory_design_ecoli.py
predicted targets are saved in examples/results/Tryptophan_design_results_ecFactory*
'''

from pathlib import Path
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ── resolve project root ──────────────────────────────────────────────────────
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

from strainOptimizer.analysis.dataset import load_experiment_targets, calculate_exp_consistency

RESULTS_DIR  = PROJECT_ROOT / 'examples' / 'results'
PRODUCT      = 'Tryptophan'
EXP_PRODUCT  = 'ecoli_tryptophan'   # key used by load_experiment_targets
ALG          = 'ecFactory'
OUTPUT_DIR   = PROJECT_ROOT / 'paper_analysis_code' / 'results'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ── 1. Load predicted results ─────────────────────────────────────────────────
prefix = RESULTS_DIR / f'{PRODUCT}_design_results_{ALG}'

level_files = {
    # 'Level 1': prefix.parent / f'{PRODUCT}_design_results_{ALG}_level_1_result.xlsx',
    'Level 1': prefix.parent / f'{PRODUCT}_design_results_{ALG}_combined_result.xlsx',
    'Level 2': prefix.parent / f'{PRODUCT}_design_results_{ALG}_level_2_result.xlsx',
    'Level 3': prefix.parent / f'{PRODUCT}_design_results_{ALG}_level_3_result.xlsx',
    'FCC':     prefix.parent / f'{PRODUCT}_design_results_{ALG}_fcc_result.xlsx',
}

predict_levels = {}
for label, path in level_files.items():
    if path.exists():
        predict_levels[label] = pd.read_excel(path, index_col=0)
        print(f'Loaded {label}: {len(predict_levels[label])} genes')
    else:
        print(f'[WARN] {path.name} not found, skipping {label}')

if not predict_levels:
    raise FileNotFoundError(f'No result files found in {RESULTS_DIR}. '
                            'Run examples/1.ecFactory_design_ecoli.py first.')

# ── 2. Load experimental dataset ─────────────────────────────────────────────
print(f'\nLoading experimental dataset for: {EXP_PRODUCT}')
exp_data = load_experiment_targets(product=EXP_PRODUCT)
print(f'Experimental targets ({len(exp_data)}):')
print(exp_data.to_string())

# ── 3. Evaluate each level ────────────────────────────────────────────────────
print('\n' + '=' * 60)
print(' Prediction vs Experiment Consistency')
print('=' * 60)

summary_rows = []
consistency_by_level = {}

# KO and KD are treated as one down-regulation category for evaluation
MERGE_KO_KD = True

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

    # per-action breakdown (KD now covers merged KO+KD)
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
        'Level'    : label,
        'Predicted': overall['predict_num'],
        'Exp'      : overall['exp_num'],
        'Hits'     : overall['hit_num'],
        'Recall'   : round(overall['consistency'], 4),
        'Precision': round(overall['precision'], 4),
        'F1'       : round(f1, 4),
    })

summary_df = pd.DataFrame(summary_rows).set_index('Level')
print('\n' + '=' * 60)
print(' Summary Table')
print('=' * 60)
print(summary_df.to_string())

# ── 4. Hit gene details ───────────────────────────────────────────────────────
print('\n' + '=' * 60)
print(' Hit gene details (all levels)')
print('=' * 60)

for label, consistency in consistency_by_level.items():
    hit_genes = []
    for action in ('KO', 'KD', 'OE'):
        if action in consistency:
            for g in consistency[action]['hit']:
                hit_genes.append({'gene': g, 'action': action})
    if hit_genes:
        hit_df = pd.DataFrame(hit_genes).set_index('gene')
        # attach exp shortName
        hit_df = hit_df.join(exp_data[['shortNames', 'enzymes']], how='left')
        print(f'\n{label} hits:')
        print(hit_df.to_string())

# # ── 5. Save summary ───────────────────────────────────────────────────────────
# out_xlsx = OUTPUT_DIR / f'{PRODUCT}_ecoli_consistency_summary.xlsx'
# with pd.ExcelWriter(out_xlsx) as writer:
#     summary_df.to_excel(writer, sheet_name='Summary')
#     exp_data.to_excel(writer, sheet_name='Exp_targets')
#     for label, pred_df in predict_levels.items():
#         sheet = label.replace(' ', '_')
#         pred_df.to_excel(writer, sheet_name=sheet)
# print(f'\nSummary saved to {out_xlsx}')

# ── 6. Bar chart: Recall & Precision across levels ────────────────────────────
if not summary_df.empty:
    plt.rcParams['font.family'] = 'Arial'
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    fsize = 12
    levels = summary_df.index.tolist()
    x = np.arange(len(levels))
    width = 0.35

    # Recall / Precision grouped bar
    ax = axes[0]
    ax.bar(x - width/2, summary_df['Recall'],   width, label='Recall',    color='#4C72B0')
    ax.bar(x + width/2, summary_df['Precision'], width, label='Precision', color='#DD8452')
    ax.set_xticks(x)
    ax.set_xticklabels(levels,fontsize=fsize)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel('Score', fontsize=fsize)
    ax.set_title(f'{PRODUCT} – Recall & Precision')
    ax.legend(fontsize=fsize)
    for i, row in enumerate(summary_df.itertuples()):
        ax.text(i - width/2, row.Recall + 0.02, f'{row.Recall:.2f}', ha='center', fontsize=9)
        ax.text(i + width/2, row.Precision + 0.02, f'{row.Precision:.2f}', ha='center', fontsize=9)

    # Hit / Miss stacked bar (relative to exp targets)
    ax2 = axes[1]
    exp_total = summary_df['Exp'].iloc[0]
    hits   = summary_df['Hits'].values
    misses = exp_total - hits
    ax2.bar(levels, hits,   label='Hits',   color='#55A868')
    ax2.bar(levels, misses, bottom=hits, label='Missed', color='#C44E52', alpha=0.7)
    ax2.set_ylabel('Number of experimental targets', fontsize=fsize)
    ax2.set_xticks(x, labels=levels, fontsize=fsize)
    ax2.set_title(f'{PRODUCT} – Hits vs Missed (n={exp_total})', fontsize=fsize)
    ax2.legend(fontsize=fsize)
    for i, (h, m) in enumerate(zip(hits, misses)):
        ax2.text(i, h / 2, str(h), ha='center', va='center', color='white', fontweight='bold')

    plt.tight_layout()
    out_fig = OUTPUT_DIR / f'{PRODUCT}_ecoli_consistency_plot.png'
    plt.savefig(out_fig, dpi=500)
    # print(f'Plot saved to {out_fig}')
    plt.show()
