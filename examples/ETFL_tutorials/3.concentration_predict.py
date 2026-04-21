# -*- coding: utf-8 -*-
"""
Concentration prediction tutorial for ETFL models.
Demonstrates extraction of enzyme / mRNA concentrations from ETFL solutions.
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
from strainOptimizer.etfl.analysis.utils import enzymes_to_peptides_conc


def extract_conc_results(model, mol_type, optimize=True):
    """Extract protein or mRNA concentrations after model optimization.

    Returns dict {mol_id: concentration (mmol/gDW)}.
    enzyme.scaled_X gives the relative fraction of total macromolecule.
    """
    if optimize:
        model.optimize()
    molecules = model.enzymes if mol_type == 'enzyme' else model.mrnas
    return {mol.id: mol.X for mol in molecules}


def gene_to_rxn(model, geneID):
    """Return reactions associated with a gene."""
    gene = model.genes.get_by_id(geneID)
    rxn_list = list(gene.reactions)
    if not rxn_list:
        print(f"Gene {geneID} is not involved in any reaction")
        return None
    return rxn_list[0] if len(rxn_list) == 1 else rxn_list


if __name__ == '__main__':
    model = load_json_model(MODEL_PATH, solver='optlang-gurobi')

    enzyme_conc = extract_conc_results(model, 'enzyme')
    mrna_conc   = extract_conc_results(model, 'mRNA', optimize=False)  # already optimal

    print(f"Extracted {len(enzyme_conc)} enzyme concentrations")
    print(f"Extracted {len(mrna_conc)} mRNA concentrations")

    # Peptide concentrations from enzyme concentrations
    peps_conc = enzymes_to_peptides_conc(model, enzyme_conc)
    print(f"Derived {len(peps_conc)} peptide concentrations")

    # Show gene → reaction mapping for first 3 coding genes
    for geneID in list(mrna_conc.keys())[:3]:
        rxns = gene_to_rxn(model, geneID)
        print(f"  {geneID}: {rxns}")
