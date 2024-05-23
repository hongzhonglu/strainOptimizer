# -*- coding: utf-8 -*-
# file : pyTFA.py
# project : strainOptimizer
'''Tutorial codes for pyTFA'''
from cobra.io import read_sbml_model
from strainOptimizer.simulation.TFA import tfa,build_tfa_model
import pandas as pd
import os

os.chdir(r'D:\code\github\strainOptimizer')


# load model
ecModel=read_sbml_model('examples/models/yeast/ecYeastGEM_batch.xml')
yeast9=read_sbml_model('examples/models/yeast/yeast-GEM.xml')

c_source="r_1714_REV"      # glucose exchange rxn
growth_id='r_2111'
c_uptake=5


# build TFA model
tfa_model=build_tfa_model(model=ecModel,model_type='ecGEM')
tfa_yeast9=build_tfa_model(model=yeast9,model_type='GEM')

# run TFA
sol=tfa(model=ecModel,
        c_source=c_source,
        targetID=growth_id,
        c_uptake=c_uptake,
        model_type='ecGEM',
        direction='max')



from pytfa.io.plotting import plot_fva_tva_comparison
from bokeh.plotting import show, output_file
from bokeh.layouts import column
import pytfa.io.viz




