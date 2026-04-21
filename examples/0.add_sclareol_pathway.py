# -*- coding: utf-8 -*-
"""Add sclareol biosynthesis pathway to yeast ecGEM and test production."""
from pathlib import Path
import sys
import cobra
from cobra import Reaction, Metabolite

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

model = load_model(
    filename=str(PROJECT_ROOT / 'examples/models/yeast/GAN_ecYeast/GAN_all_v2.xml'),
    model_type='ecGEM'
)

# build non-exist metabolites
ggdp_c = Metabolite('ggdp_c', formula='C20H33O7P2', name='Geranylgeranyl diphosphate', compartment='c')
lpp_c = Metabolite('lpp_c', formula='C20H35O8P2', name='(13E)-8a-hydroxylabden-15-yl diphosphate', compartment='c')
sclareol_c = Metabolite('sclareol_c', formula='C20H36O2', name='sclareol', compartment='c')

# load existing metabolites (try different ID conventions for ecGEM / ETFL)
for suffix in ('', '[c]', '_c'):
    try:
        frdp_c = model.metabolites.get_by_id(f's_0190{suffix}')
        ipdp_c = model.metabolites.get_by_id(f's_0943{suffix}')
        ppi_c  = model.metabolites.get_by_id(f's_0633{suffix}')
        h_c    = model.metabolites.get_by_id(f's_0794{suffix}' if suffix != '_c' else 's_0793_c')
        h2o_c  = model.metabolites.get_by_id(f's_0803{suffix}')
        break
    except KeyError:
        continue

print('frdp_c:', frdp_c.name)
print('ipdp_c:', ipdp_c.name)
print('ppi_c:',  ppi_c.name)
print('h_c:',    h_c.name)
print('h2o_c:',  h2o_c.name)

# build pathway reactions
frtt = Reaction('frtt', name='farnesyltranstransferase', lower_bound=0, upper_bound=1000)
frtt.add_metabolites({frdp_c: -1, ipdp_c: -1, ggdp_c: 1, ppi_c: 1})
frtt.gene_reaction_rule = 'FRTT'

lpps = Reaction('lpps', name='(13E)-8a-hydroxylabden-15-yl diphosphate synthase', lower_bound=0, upper_bound=1000)
lpps.add_metabolites({ggdp_c: -1, h_c: -1, lpp_c: 1, h2o_c: 1})
lpps.gene_reaction_rule = 'SsLPPS'

tps = Reaction('tps', name='sclareol synthase', lower_bound=0, upper_bound=1000)
tps.add_metabolites({lpp_c: -1, h2o_c: -1, sclareol_c: 1, h_c: 1, ppi_c: 1})
tps.gene_reaction_rule = 'SsTPS'

model.add_reactions([frtt, lpps, tps])
model.add_boundary(sclareol_c, type='demand')

with model:
    model.objective = 'DM_sclareol_c'
    print('Max sclareol production:', model.slim_optimize())

cobra.io.write_sbml_model(
    model,
    str(PROJECT_ROOT / 'examples/models/yeast/GAN_ecYeast/sclareol_GAN_all_v2.xml')
)
