import cobra
from cobra import Model, Reaction, Metabolite
import os
'''from ETFLdesigner.ETFLdesigner.io.ecModel import load_ecmodel'''
model = cobra.io.read_sbml_model(r'D:\Modelling\Basic cobra\model\ecYeastGEM_batch860.xml')

'''这一步可以开始导入青蒿酸模型进入我们的现有模型'''
# 第一个反应
reaction = Reaction('Amorpha_diene_sys')
reaction.name = 'Amorpha_diene systhesis'
reaction.subsystem = 'Cell cytoplasm Biosynthesis'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.

Amorpha_diene = Metabolite(
    'Amorpha_diene',
    formula='C15H24',
    name='Amorpha_diene',
    compartment='c')

# 酶代谢物
'''prot_ADS = Metabolite(
    'prot_ADS',
    formula='',
    name='ADS enzyme',
    compartment='c')'''

reaction.add_metabolites({
    model.metabolites.get_by_id('s_0190[c]'): -1.0,  # FPP
    Amorpha_diene: 1.0,                              # 紫穗槐二烯
    model.metabolites.get_by_id('s_0633[c]'): 1.0,   # 焦磷酸
    # prot_ADS: -1.9305/3600                    # 这里我们姑且让他和后面的kcat一致，试一试
})
model.add_reactions([reaction])

# 第二个反应
reaction = Reaction('artemisinic_alcohol_sys')
reaction.name = 'artemisinic_alcohol systhesis'
reaction.subsystem = 'Cell cytoplasm Biosynthesis'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.

artemisinic_alcohol = Metabolite(
    'artemisinic_alcohol',
    formula='C15H24O',
    name='artemisinic_alcohol',
    compartment='c')

# 酶代谢物
'''prot_ALDH1 = Metabolite(
    'prot_ALDH1',
    formula='',
    name='prot_ALDH1',
    compartment='c')
'''
reaction.add_metabolites({
    Amorpha_diene: -1.0,                                 # 紫穗槐二烯
    model.metabolites.get_by_id('s_1275[c]'): -1.0,      # O2
    model.metabolites.get_by_id('s_1275[c]'): -1.0,      # NADPH
    model.metabolites.get_by_id('s_1207[c]'): 1.0,       # NADP+
    model.metabolites.get_by_id('s_0794[c]'): -1.0,      # H+
    model.metabolites.get_by_id('s_0803[c]'): 1.0,       # H2O
    artemisinic_alcohol: 1.0,                            # 青蒿醇
    # prot_ALDH1: -0.65359/3600/4                        # 第2步酶的1/kcat
})
model.add_reactions([reaction])

# 第三个反应
reaction = Reaction('atemisinic_aldehyde_sys')
reaction.name = 'atemisinic_aldehyde systhesis'
reaction.subsystem = 'Cell cytoplasm Biosynthesis'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.

atemisinic_aldehyde = Metabolite(
    'atemisinic_aldehyde',
    formula='C15H22O',
    name='atemisinic_aldehyde',
    compartment='c')

# 酶代谢物
'''prot_ALDH1 = Metabolite(
    'prot_ADS',
    formula='',
    name='ADS enzyme',
    compartment='c')'''

reaction.add_metabolites({
    artemisinic_alcohol: -1.0,                           # 青蒿醇
    model.metabolites.get_by_id('s_1275[c]'): -1.0,      # O2
    model.metabolites.get_by_id('s_1275[c]'): -1.0,      # NADPH
    model.metabolites.get_by_id('s_1207[c]'): 1.0,       # NADP+
    model.metabolites.get_by_id('s_0794[c]'): -1.0,      # H+
    model.metabolites.get_by_id('s_0803[c]'): 2.0,       # H2O
    atemisinic_aldehyde: 1.0,                            # 青蒿醛
    # prot_ALDH1: -0.65359/3600/4                                # 第3步酶的1/kcat
})
model.add_reactions([reaction])

# 第四个反应
reaction = Reaction('atemisinic_acid_sys')
reaction.name = 'atemisinic_acid systhesis'
reaction.subsystem = 'Cell cytoplasm Biosynthesis'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.

atemisinic_acid = Metabolite(
    'atemisinic_acid',
    formula='C15H21O2',
    name='atemisinic_acid',
    compartment='c')

# 酶代谢物
'''prot_ALDH1 = Metabolite(
    'prot_ADS',
    formula='',
    name='ADS enzyme',
    compartment='c')'''

reaction.add_metabolites({
    atemisinic_aldehyde: -1.0,  # 青蒿醛
    model.metabolites.get_by_id('s_1275[c]'): -1.0,      # O2
    model.metabolites.get_by_id('s_1275[c]'): -1.0,      # NADPH
    model.metabolites.get_by_id('s_1207[c]'): 1.0,       # NADP+
    model.metabolites.get_by_id('s_0794[c]'): -1.0,      # H+
    model.metabolites.get_by_id('s_0803[c]'): 1.0,       # H2O
    atemisinic_acid: 1.0,                                # 青蒿酸
    # prot_ALDH1: -0.65359/3600/4                                # 第4步酶的1/kcat
})
model.add_reactions([reaction])

model.add_boundary(model.metabolites.get_by_id("atemisinic_acid"), type="sink")  # SK_atemisinic_acid

'''# ADS加入蛋白池
reaction = Reaction('draw_ADS')
reaction.name = 'draw_prot_ADS'
reaction.subsystem = 'Cell cytoplasm Biosynthesis'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.

reaction.add_metabolites({ model.metabolites.get_by_id('prot_pool[c]'): -63.933, prot_ADS: 1})
model.add_reactions([reaction])

reaction = Reaction('draw_ALDH1')
reaction.name = 'draw_prot_ALDH1'
reaction.subsystem = 'Cell cytoplasm Biosynthesis'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.

reaction.add_metabolites({ model.metabolites.get_by_id('prot_pool[c]'): -53.8, prot_ALDH1: 1})
model.add_reactions([reaction])
'''

'''model.remove_reactions(model.reactions.get_by_id('draw_EGT2'))'''

output_path = "C:/Users/Administrator/Desktop/ecGEM with atemisinic acid.xml"
cobra.io.write_sbml_model(model, output_path)