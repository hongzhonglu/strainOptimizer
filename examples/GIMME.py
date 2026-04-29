# -*- coding: utf-8 -*-
"""
GIMME transcriptome integration example.
Integrates RNA-seq data into a yeast ecGEM to generate a context-specific model.
Requires: troppo package  (pip install troppo)
"""
from pathlib import Path
import sys
import pandas as pd
import cobra

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

from strainOptimizer.io import load_model
from strainOptimizer.manipulation.integration import integrate_omic_data_to_ecmodel


if __name__ == "__main__":
    model = load_model(
        str(PROJECT_ROOT / 'examples/models/yeast/ecYeastGEM_batch.xml'),
        model_type='ecGEM',
    )

    df_transcriptomics = pd.read_csv(
        str(PROJECT_ROOT / 'sclareol_cell_factory_design/data/sce_sclareol_gene_fpkm.tsv'),
        index_col=0, sep='\t',
    )
    samples = ['SCX42-1', 'SCX42-2', 'SCX42-3']
    df_scx42 = df_transcriptomics[samples].mean(axis=1)

    parameters = {
        'objective_reaction_id': 'r_2111',
        'obj_frac': 0.8,
        'expression_threshold': 12,
    }

    specific_model = integrate_omic_data_to_ecmodel(
        model=model,
        omic_data=df_scx42,
        method='GIMME',
        parameters=parameters,
        copy_model=True,
    )
    print('Slim optimize:', specific_model.slim_optimize())
