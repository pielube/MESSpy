"""
MESSpy - Run_test

Do not change this script
Create your own run_dev.py and input_dev/ 
"""

#%% ###########################################################################
"""
PRE PROCESSING
==============
"""

# Import modules
from rec import REC
import os
import json

# Selecting simulation names
study_case = 'REC_test' # str name for results file.pkl
reference_case = 'buiseness as usual' # str name for results file.pkl

# Selecting input files:
path = r'./input_test' # change the path with r'./input_dev' if you are working on your own run_dev 67
# path = r'./input_dev'

file_structure = 'structure.json'
file_general   = 'general.json'
file_eco       = 'economics.json'
file_refcase   = 'refcase.json'

# Opening input files:
with open(os.path.join(path,file_structure),'r') as f: structure = json.load(f)
with open(os.path.join(path,file_general),'r') as   f: general = json.load(f)
with open(os.path.join(path,file_eco),'r') as       f: economic_data = json.load(f)
with open(os.path.join(path,file_refcase),'r') as   f: structure0 = json.load(f)

# Edit input files:
import preprocess_test as pre
#import preprocess_dev as pre
# Instead of modifying the original input files.json we suggest to modify the dictionary variables 
# structure, structure0, general and economic_data using specific functions that you can define in preprocess_dev. 
# An example:
# structure = pre.change_peakP(structure, 'p1', 5) 

#%% ###########################################################################
"""
SOLVER
======
"""

rec = REC(structure,general,path) # create REC object
rec.REC_energy_simulation() # simulate REC enegy balances
rec.save(study_case) # save results in 'study_case.pkl'
  
#%% ###########################################################################
"""
POST PROCESS - ECONOMIC ANALYSIS
================================
"""

from economics import NPV
    
# Reference case simulation (run only if changed)
rec0 = REC(structure0,general,path) # create REC
rec0.REC_energy_simulation() # simulate REC 
rec0.save(reference_case) # save results in 'reference_case.pkl'

#%% Actual economic analysis
NPV(structure,structure0,study_case,reference_case,economic_data,general['simulation years'],path) 


#%% ###########################################################################
"""
POST PROCESS - PLOTTING
================================
"""

import postprocess_test as pp
# import postprocess_dev as pp

"""
some post-process are alredy avaiable as examples in postprocess_test
you should create your own postprocess_dev.py
"""

import postprocess_test as pp
# import postprocess_dev as pp

"""
some post-process are alredy avaiable as examples in postprocess_test
you should create your own postprocess_dev.py
"""

# pp.total_balances(study_case,'p1','hydrogen')
# pp.total_balances(study_case,'p1','heat')
# pp.total_balances(study_case,'p2','hydrogen')
# pp.total_balances(study_case,'c1')
# pp.REC_electricity_balance(study_case)

# pp.LOC_plot(study_case)
#pp.storage_control(study_case)

# pp.NPV_plot(study_case)

# pp.hourly_balances_electricity(study_case,'p1', 0, 5)
# pp.hourly_balances_heat(study_case,'p1',0, 5)

pp.ele_param(study_case,0,40)
pp.fc_param(study_case,0,40)
# pp.Flows(study_case)

# pp.load_profile_FC(study_case, 'p1', 0, 50)
# pp.load_profile_ele(study_case, 'p1', 0, 50 )

# pp.csc_allocation_sum(study_case)

# rec.locations['p1'].energy_balance['electricity']['PV'].sum()

#%% ##########################################################################
"Parametric analysis - workflow example"

# for parameter in list of values:
    # open input files
    # edit input files using parameter and function you can define in preprocess_test 
    # create and simulate rec
    # if you need it, run the economic analysis (remember that you need to create and simulate rec0 just ones)
    # save the results
# use results to make what you want 

# enjoy ;)






