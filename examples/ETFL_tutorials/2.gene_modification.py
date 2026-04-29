# -*- coding: utf-8 -*-
"""
Gene modification tutorial for ETFL models.
Demonstrates KO / knockdown / OE by constraining translation reaction bounds.
"""
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
MODEL_PATH = str(PROJECT_ROOT / 'examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json')

from strainOptimizer.etfl.io.json import load_json_model


def gene_knock_out_down_up(model, geneID, regulate_ratio):
    """Regulate gene expression by constraining the translation reaction.

    regulate_ratio > 1  → overexpression
    0 < regulate_ratio ≤ 1 → knockdown
    regulate_ratio = 0  → knockout
    """
    model.optimize()
    if model.solution.status != 'optimal':
        print("The model has no solution")
        return None
    trans_rxn = model.get_translation(geneID)
    initial_abundance = trans_rxn.flux
    target_abundance = initial_abundance * regulate_ratio
    if regulate_ratio > 1:
        model.logger.info("simulating %s overexpression in %s" % (geneID, model.id))
        trans_rxn.lower_bound = target_abundance
    else:
        model.logger.info("simulating %s knock-down/out in %s" % (geneID, model.id))
        trans_rxn.upper_bound = target_abundance
    return model


if __name__ == '__main__':
    model = load_json_model(MODEL_PATH)
    sol = model.optimize()
    print("original growth: %s" % sol.objective_value)

    coding_geneList = [mrna.id for mrna in model.mrnas]
    for geneID in coding_geneList[:5]:  # test first 5 genes
        with model:
            model_changed = gene_knock_out_down_up(model, geneID=geneID, regulate_ratio=2)
            if model_changed is not None:
                gr = model_changed.slim_optimize()
                print(f"{geneID} OE growth: {gr}")
