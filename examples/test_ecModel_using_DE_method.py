# newly added
import sys
import os
from cobra.io import load_matlab_model, read_sbml_model
sys.path.append('/Users/xluhon/Documents/GitHub/yetfl/code')
sys.path.append('/Users/xluhon/Documents/GitHub/etfl')
os.chdir('/Users/xluhon/Documents/GitHub/ETFLdesigner')

# load the GEM model
ecGEM = read_sbml_model("examples/models/yeast/heme_ecYeastGEM.xml")
ecGEM.solver ='glpk' # if not setting as this, the model correction will produce error
model = ecGEM.copy()



