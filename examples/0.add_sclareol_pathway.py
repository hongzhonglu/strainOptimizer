# -*- coding: utf-8 -*-
# date : 2024/12/23 
import cobra
from cobra import Model, Reaction, Metabolite
from strainOptimizer.io import load_model

# load model
# model=load_model(filename='examples/models/yeast/ecYeastGEM_batch.xml',model_type='ecGEM')
# model=load_model(filename='examples/models/yeast/yeast8_cEFL_2584_enz_64_bins__20231221_083715.json',model_type='etfl')
# model=load_model(filename='examples/models/yeast/yeast-GEM.xml',model_type='gem')
model=load_model(filename='examples/models/yeast/GAN_ecYeast/GAN_all_v2.xml',model_type='ecGEM')

# build non-exist metabolite
ggdp_c=cobra.Metabolite('ggdp_c',formula='C20H33O7P2',name='Geranylgeranyl diphosphate',compartment='c')
lpp_c=cobra.Metabolite('lpp_c',formula='C20H35O8P2',name='(13E)-8a-hydroxylabden-15-yl diphosphate',compartment='c')
sclareol_c=cobra.Metabolite('sclareol_c',formula='C20H36O2',name='sclareol',compartment='c')

# load exist metabolite
try:
    frdp_c=model.metabolites.get_by_id('s_0190')   # Farnesyl diphosphate
    ipdp_c=model.metabolites.get_by_id('s_0943')   # Isopentenyl diphosphate
    ppi_c=model.metabolites.get_by_id('s_0633')    # diphosphate
    h_c=model.metabolites.get_by_id('s_0794')      # H+
    h2o_c=model.metabolites.get_by_id('s_0803')    # H2O
except:
    try:
        # for ecModel
        # add [c]
        frdp_c = model.metabolites.get_by_id('s_0190[c]')  # Farnesyl diphosphate
        ipdp_c = model.metabolites.get_by_id('s_0943[c]')  # Isopentenyl diphosphate
        ppi_c = model.metabolites.get_by_id('s_0633[c]')  # diphosphate
        h_c = model.metabolites.get_by_id('s_0794[c]')  # H+
        h2o_c = model.metabolites.get_by_id('s_0803[c]')  # H2O
    except:
        # for ETFL model
        # add _c
        frdp_c = model.metabolites.get_by_id('s_0190_c')  # Farnesyl diphosphate
        ipdp_c = model.metabolites.get_by_id('s_0943_c')  # Isopentenyl diphosphate
        ppi_c = model.metabolites.get_by_id('s_0633_c')  # diphosphate
        h_c = model.metabolites.get_by_id('s_0793_c')  # H+
        h2o_c = model.metabolites.get_by_id('s_0803_c')  # H2O

print('frdp_c:',frdp_c.name)
print('ipdp_c:',ipdp_c.name)
print('ppi_c:',ppi_c.name)
print('h_c:',h_c.name)
print('h2o_c:',h2o_c.name)



# build non-exist reactions
frtt=cobra.Reaction('frtt',name='farnesyltranstransferase',lower_bound=0,upper_bound=1000)
frtt.add_metabolites({frdp_c:-1,ipdp_c:-1,ggdp_c:1,ppi_c:1})
frtt.gene_reaction_rule='FRTT'

lpps=cobra.Reaction('lpps',name='(13E)-8a-hydroxylabden-15-yl diphosphate synthase',lower_bound=0,upper_bound=1000)
lpps.add_metabolites({ggdp_c:-1,h_c:-1,lpp_c:1,h_c:1,h2o_c:1})
lpps.gene_reaction_rule='SsLPPS'

tps=cobra.Reaction('tps',name='sclareol synthase',lower_bound=0,upper_bound=1000)
tps.add_metabolites({lpp_c:-1,h2o_c:-1,sclareol_c:1,h_c:1,ppi_c:1})
tps.gene_reaction_rule='SsTPS'

# add reactions to model
model.add_reactions([frtt,lpps,tps])

# add demand reaction for sclareol: DM_sclareol_c
model.add_boundary(sclareol_c,type='demand')

with model:
    model.objective='DM_sclareol_c'
    print(model.slim_optimize())

# save model
cobra.io.write_sbml_model(model,'examples/models/yeast/GAN_ecYeast/sclareol_GAN_all_v2.xml')

# from etfl.io.json import save_json_model
# save_json_model(model,'examples/models/yeast/cEFL_sclareol.json')