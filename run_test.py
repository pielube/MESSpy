"""
MESSpy - Run_test

don't work on this script:
    you should create your own run_dev.py, input_dev/, postprocess_dev.py and preprocess_dev.py
"""

#%% ###########################################################################
"""
PRE PROCESSING
==============
"""

# Import modules
from rec import REC
from economics import NPV
import postprocess_test as pp
import preprocess_test as pre
#import postprocess_dev as pp
#import preprocess_dev as pp
import os
import json

# Selecting simulation names
study_case = 'REC_test' # str name for results file.pkl
reference_case = 'buiseness as usual' # str name for results file.pkl

# Selecting input files:
path = r'./input_test' # change the path with r'./input_dev' if you are working on your own run_dev
#path = r'./input_dev'

file_structure = 'structure.json'
file_general = 'general.json'
file_eco = 'economics.json'
file_refcase = 'refcase.json'

# Opening input files:
with open(os.path.join(path,file_structure),'r') as f: structure = json.load(f)
with open(os.path.join(path,file_general),'r') as f: general = json.load(f)
with open(os.path.join(path,file_eco),'r') as f: economic_data = json.load(f)
with open(os.path.join(path,file_refcase),'r') as f: structure0 = json.load(f)

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

some post-process are alredy avaiable as examples in postprocess_test
you should create your own postprocess_dev.py
"""

pp.total_balances(study_case,'prosumer_1','electricity')
pp.total_balances(study_case,'prosumer_2','electricity')
pp.total_balances(study_case,'prosumer_2','hydrogen')
pp.total_balances(study_case,'consumer_1','electricity')

#pp.total_balances(reference_case,'consumer_2','electricity')
#pp.total_balances(reference_case,'consumer_2','heating water')
#pp.total_balances(reference_case,'consumer_2','gas')
pp.total_balances(study_case,'consumer_2','electricity')
pp.total_balances(study_case,'consumer_2','heating water')


pp.REC_electricity_balance(study_case)

pp.LOC_plot(study_case)

pp.NPV_plot(study_case)

pp.hourly_balances_electricity(study_case,'prosumer_1', 2, 3)
pp.hourly_balances_electricity(study_case,'prosumer_2', 2, 3)
#pp.hourly_balances_electricity(study_case,'consumer_1', 2, 3)
#pp.hourly_balances_electricity(study_case,'consumer_2', 2, 3)

#pp.csc_allocation_sum(study_case)
#pp.storage_control(study_case)
#pp.ele_param(study_case, 2, 3)
#pp.fc_param(study_case, 2, 3)


#%% ##########################################################################
"Sensitivity analysis - practical example"
# click on the following code and press Ctrl+5 to discomment all together

# =============================================================================
# import numpy as np
# import pickle
# import matplotlib.pyplot as plt
# 
# pv_size = np.arange(1,11)
# npv = []
# sc = [] # self-consumption
# ss = [] # self-sufficiency
# 
# for pv in pv_size:
#     study_case = f"PV size = {pv}"
#     new_structure = pre.change_peakP(structure, 'prosumer_1', pv)
#     rec = REC(new_structure,general,path) # create REC object
#     rec.REC_energy_simulation() # simulate REC enegy balances
#     rec.save(study_case) # save results in 'study_case.pkl'
#     
#     with open('results/balances_'+study_case+'.pkl', 'rb') as f: balances = pickle.load(f)
#     demand = -balances['prosumer_1']['electricity']['demand'].sum() # read from saved results .pkl
#     production = balances['prosumer_1']['electricity']['PV'].sum()
#     into_grid = balances['prosumer_1']['electricity']['grid'].sum(where=balances['prosumer_1']['electricity']['grid']<0)
#     from_grid = balances['prosumer_1']['electricity']['grid'].sum(where=balances['prosumer_1']['electricity']['grid']>0)
#     
#     # you can also read values from the python object rec. (in this case you do not need rec.save)
#     # demand = -rec.locations['prosumer_1'].energy_balance['electricity']['demand'].sum() # read from python 
#     
#     sc.append((production+into_grid)/production*100)
#     ss.append((demand-from_grid)/demand*100)
#     
# plt.figure(dpi=1000)
# plt.plot(pv_size,sc,label='Self-consumption')
# plt.plot(pv_size,ss,label='Self-sufficiency')
# plt.xlabel("PV peak power [kWp]")
# plt.ylabel("[%]")
# plt.grid()
# plt.legend()
# plt.title('prosumer_1')
# plt.show()
# =============================================================================





    






