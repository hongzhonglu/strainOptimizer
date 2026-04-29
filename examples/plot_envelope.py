# -*- coding: utf-8 -*-
"""
Plot production envelopes for experimental gene targets.
Requires: strainOptimizer.visualization module.
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

from strainOptimizer.io import load_model
from strainOptimizer.analysis.dataset import load_experiment_targets
from strainOptimizer.simulation import ppFBA
from strainOptimizer.simulation.utils import set_carbon_source_bounds
from strainOptimizer.visualization import calculate_flux_space, plot_envelope
from matplotlib import pyplot as plt


def extracti_from_list(l, indices):
    return [l[i] for i in indices]


if __name__ == '__main__':
    solver = 'optlang-gurobi'
    model = load_model(
        str(PROJECT_ROOT / 'examples/models/yeast/ecYeastGEM_batch.xml'),
        model_type='ecGEM', solver=solver,
    )
    all_geneList = [gene.id for gene in model.genes]
    c_source = 'r_1714_REV'
    c_uptake = 10
    oe_factor = 2
    set_carbon_source_bounds(model, c_source, c_uptake, model_type='ecGEM', fixed=True)

    for product_id, product_label, exp_product_name in [
        ('r_1589', '2-phenylethanol(mmol/gDW/h)', '2-phenylethanol'),
        ('r_2051', 'spermidine(mmol/gDW/h)', 'spermidine'),
    ]:
        x = 'r_2111'
        y = product_id
        df_exp = load_experiment_targets(exp_product_name)
        df_exp = df_exp[df_exp.index.isin(all_geneList)]
        flux_spaceList = [calculate_flux_space(model, x, y)]

        for geneID in df_exp.index:
            action = df_exp.loc[geneID, 'action']
            rxnidList = [rxn.id for rxn in model.genes.get_by_id(geneID).reactions
                         if rxn.id.startswith('r_')]
            with model:
                if action == 'OE':
                    wt_sol = ppFBA(model=model, target_id=x, c_source=c_source,
                                   c_uptake=c_uptake, model_type='ecGEM')
                    for rxnid in rxnidList:
                        wt_flux = wt_sol.fluxes[rxnid]
                        if wt_flux > 0:
                            model.reactions.get_by_id(rxnid).bounds = wt_flux * oe_factor, 1000
                        elif wt_flux < 0:
                            model.reactions.get_by_id(rxnid).bounds = -1000, wt_flux * oe_factor
                elif action in ('KD', 'KO'):
                    model.genes.get_by_id(geneID).knock_out()
                try:
                    flux_space = calculate_flux_space(model, x, y)
                except Exception:
                    flux_space = None
                flux_spaceList.append(flux_space)

        labels = ['wild-type'] + list(df_exp.index)
        fig = plot_envelope(flux_spaceList, labels, show=False)
        fig.axes[0].set_xlabel('growth(/h)', fontsize=16)
        fig.axes[0].set_ylabel(product_label, fontsize=16)
        plt.show()
