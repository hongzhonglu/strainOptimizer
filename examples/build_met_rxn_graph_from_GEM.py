# -*- coding: utf-8 -*-
import networkx as nx
import cobra


def build_met_rxn_graph_from_GEM(model, exclusions=None):
    # 1. Build a directed bipartite graph
    G = nx.DiGraph()
    if exclusions is None:
        # exclude: cofactors, O2, H2O, H+, ATP
        exclusions = ['oxygen', 'H2O', 'H+', 'ATP', 'ADP', 'NADH', 'NADPH', 'NADP+', 'NAD+', 'FADH2', 'FAD',
                     'ferricytochrome c', 'ferrocytochrome c', 'ubiquinol-6', 'ubiquinone-6','diphosphate','phosphate']

    # add metabolite nodes
    for met in model.metabolites:
        name = met.name.split(' [')[0]
        if name in exclusions:
            continue
        # Add metabolite nodes (with prefix to distinguish)
        G.add_node(f"{met.id}", bipartite=0)

    # add reaction nodes and directed edges
    for rxn in model.reactions:
        if rxn.bounds == (0, 0):
            # Skip reactions that are not active (bounds are zero)
            continue
        # Add reaction nodes
        G.add_node(f"{rxn.id}", bipartite=1)

        # if rxn is reversible, add edges in both directions
        if rxn.reversibility:
            for met, coeff in rxn.metabolites.items():
                name = met.name.split(' [')[0]
                if name in exclusions:
                    continue
                G.add_edge(f"{met.id}", f"{rxn.id}")
                G.add_edge(f"{rxn.id}", f"{met.id}")
        else:
            if rxn.lower_bound == 0 and rxn.upper_bound > 0:
                direction = 1
            elif rxn.lower_bound < 0 and rxn.upper_bound == 0:
                direction = -1
            # Connect reactions to their associated metabolites
            for met, coeff in rxn.metabolites.items():
                name = met.name.split(' [')[0]
                if name in exclusions:
                    continue
                coeff = coeff * direction  # Adjust coefficient based on reaction direction
                if coeff < 0:
                    # Substrate: edge from metabolite -> reaction
                    G.add_edge(f"{met.id}", f"{rxn.id}")
                elif coeff > 0:
                    # Product: edge from reaction -> metabolite
                    G.add_edge(f"{rxn.id}", f"{met.id}")
    return G



model=cobra.io.read_sbml_model('examples/models/yeast/yeast-GEM.xml')

# add heterogenous pathway
### sclareol biosynthesis pathway
ggdp_c=cobra.Metabolite('ggdp_c',formula='C20H33O7P2',name='Geranylgeranyl diphosphate',compartment='c')
lpp_c=cobra.Metabolite('lpp_c',formula='C20H35O8P2',name='(13E)-8a-hydroxylabden-15-yl diphosphate',compartment='c')
sclareol_c=cobra.Metabolite('sclareol_c',formula='C20H36O2',name='sclareol',compartment='c')
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

### heme a exchange reaction
try:
    heme_metID='s_0811_c' # for ETFL model
    heme_met=model.metabolites.get_by_id(heme_metID)
except:
    heme_metID = 's_0811' # for yeast-GEM 9 model
    heme_met = model.metabolites.get_by_id(heme_metID)

# add heme a demand reaction to simulate heme a production
from cobra import Reaction
reaction=Reaction('EX_heme_a')
reaction.name='heme a production'
reaction.subsystem='heme a production'
reaction.lower_bound=0
reaction.upper_bound=1000
reaction.add_metabolites({heme_met:-1})
model.add_reactions([reaction])

# build slim model:remove cofactors , H20 ... from network
import pandas as pd
df_met_info=pd.DataFrame(columns=['met_id','met_name','C','rxn_count'])
for met in model.metabolites:
    met_id=met.id
    met_name=met.name
    try:
        c_numb=met.elements['C']
    except:
        c_numb=0
    rxn_count=len(met.reactions)
    df_met_info=df_met_info.append({'met_id':met_id,'met_name':met_name,'C':c_numb,'rxn_count':rxn_count},ignore_index=True)
# rebuild the df according to met_name, and sum rxn_count for same met_name
df_met_info=df_met_info.groupby(['met_name']).agg({'rxn_count':'sum','C':'first'}).reset_index()

noC_list=df_met_info[df_met_info['C']==0]['met_name'].tolist()
cofactor_list=df_met_info[df_met_info['rxn_count']>80]['met_name'].tolist()
exclusion_metList=list(set(noC_list+cofactor_list))


# build the directed bipartite graph
# 1. build complete network
G=build_met_rxn_graph_from_GEM(model,exclusions=[])
# save as pickle file
import pickle
with open("data/Sce_met_rxn_diGraph.pkl", "wb") as f:
    pickle.dump(G, f)

# 2. build slim network
slim_model=model.copy()
# remove reactions without gene association except for exchange reactions and transport reactions
original_model=cobra.io.read_sbml_model('examples/models/yeast/yeast-GEM.xml')
remove_rxnIDList=[]
for rxn in original_model.reactions:
    if 'exchange' in rxn.name or 'transport' in rxn.name:
        continue
    if len(rxn.genes)==0:
        if 'pseudoreaction' in rxn.name or 'SLIME' in rxn.name:
            print(rxn.name)
            remove_rxnIDList.append(rxn.id)
slim_model.remove_reactions(remove_rxnIDList,remove_orphans=True)

slim_G=build_met_rxn_graph_from_GEM(slim_model, exclusions=exclusion_metList)
with open("data/slim_Sce_met_rxn_diGraph.pkl", "wb") as f:
    pickle.dump(slim_G, f)

# # load file
# with open("data/Sce_met_rxn_diGraph.pkl", "rb") as f:
#     G_loaded = pickle.load(f)

