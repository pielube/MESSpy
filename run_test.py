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
from economics import NPV, LCOH
import postprocess_test as pp
import preprocess_test as pre
#import postprocess_dev as pp
#import preprocess_dev as pre
import os
import json
import pickle

# Selecting simulation names
name_studycase = 'Rec' # str name for energy_balances_results file.pkl
name_refcase = 'Rec0' # str name for energy_balances_results file.pkl
name_economic = 'From Rec0 to Rec' # str name for economic_assesment_results file.pkl

# Selecting input files:
path = r'./input_test' # change the path with r'./input_dev' if you are working on your own run_dev
#path = r'./input_dev'

# file_studycase      = 'studycase'
# file_refcase        = 'refcase'
file_studycase      = 'studycase_hydrogen'
file_refcase        = 'refcase_hydrogen'
file_general        = 'general'
file_tech_cost      = 'tech_cost'
file_energy_market  = 'energy_market'
file_emissions      = 'emissions'

# Opening input files:
with open(os.path.join(path,f"{file_studycase}.json"),'r') as f: studycase = json.load(f)
with open(os.path.join(path,f"{file_refcase}.json"),'r') as f: refcase = json.load(f)
with open(os.path.join(path,f"{file_emissions}.json"),'r') as f: emissions = json.load(f)
with open(os.path.join(path,f"{file_general}.json"),'r') as f: general = json.load(f)
with open(os.path.join(path,f"{file_tech_cost}.json"),'r') as f: tech_cost = json.load(f)
with open(os.path.join(path,f"{file_energy_market}.json"),'r') as f: energy_market = json.load(f)


#%% ###########################################################################
"""
SOLVER - studycase simulation
======
"""

rec = REC(studycase,general,file_studycase,file_general,path) # create REC object
rec.REC_energy_simulation() # simulate REC enegy balances
rec.tech_cost(tech_cost) # calculate the cost of all technologies 
rec.save(name_studycase) # save results in 'name_studycase.pkl'
    

#%% ###########################################################################
"""
POST PROCESS - refcase simulation and investment assessment
================================
"""
  
# Reference case simulation (run only if changed)
rec0 = REC(refcase,general,file_refcase,file_general,path) # create REC object
rec0.REC_energy_simulation() # simulate REC 
rec0.tech_cost(tech_cost) # calculate the cost of all technologies 
rec0.save(name_refcase) # save results in 'name_refcase.pkl'

    
#%% Net present value calculation to asses the investment comparing refcase and studycase
NPV(file_studycase,name_studycase,name_refcase,energy_market,general['simulation years'],path,name_economic) 


#%% ###########################################################################
"""
POST PROCESS - PLOTTING
================================

some post-process are alredy avaiable as examples in postprocess_test
you should create your own postprocess_dev.py
"""

# =============================================================================
# Here you can read the main results (balances and economic):
# =============================================================================

# with open('results/economic_assessment_'+name_economic+'.pkl', 'rb') as f: economic = pickle.load(f)
# with open('results/balances_'+name_studycase+'.pkl', 'rb') as f: balances = pickle.load(f)


# =============================================================================
# studycase postprocess useful functions
# =============================================================================

# pp.total_balances(name_studycase,'prosumer_1','electricity')
# pp.total_balances(name_studycase,'prosumer_2','electricity')
# pp.total_balances(name_studycase,'prosumer_2','hydrogen')
# pp.total_balances(name_studycase,'consumer_1','electricity')

# pp.total_balances(name_refcase,'consumer_2','electricity')
# pp.total_balances(name_refcase,'consumer_2','heating water')
# pp.total_balances(name_refcase,'consumer_2','gas')
# pp.total_balances(name_studycase,'consumer_2','electricity')
# pp.total_balances(name_studycase,'consumer_2','heating water')

# pp.REC_electricity_balance(name_studycase)
 
# pp.hourly_balances_electricity(name_studycase,'prosumer_1', 20, 21)
# pp.hourly_balances_electricity(name_studycase,'prosumer_2', 2, 3)

# pp.hourly_balances_electricity(name_studycase,'consumer_1', 2, 3)
# pp.hourly_balances_electricity(name_studycase,'consumer_2', 2, 3)

# pp.csc_allocation_sum(name_studycase)
# pp.storage_control(name_studycase)
# pp.ele_param(name_studycase, 2, 3)
# pp.fc_param(name_studycase, 2, 3)

# balances_pp = pp.renewables_results(studycase,name_studycase,general['simulation years'],'prosumer_2', emissions, '2025',235,237, print_ = True, plot = True, plot_cum = True, ghg = True,save=True)
# balances_pp = pp.renewables_results(studycase,name_studycase,general['simulation years'],'prosumer_2', emissions, '2025',120,125, print_ = False, plot = True, plot_cum = False, ghg = False,save=False)
# balances_pp = pp.renewables_results(studycase,name_studycase,general['simulation years'],'prosumer_2', emissions, '2025',220,240, print_ = False, plot = True, plot_cum = False, ghg = False,save=False)

# LCOH('prosumer_2',balances_pp,studycase,name_studycase,energy_market,general['simulation years'],path,name_economic, revenues = False, refund = True)

# pp.LOC_plot(name_studycase)

# pp.NPV_plot(name_economic)

# =============================================================================
# studycase_hydrogen postprocess useful functions
# =============================================================================

balances_pp = pp.renewables_results(studycase,name_studycase,general['simulation years'],'prosumer_1', emissions, '2025',235,237, print_ = True, plot = True, plot_cum = True, ghg = True,save=True)
balances_pp = pp.renewables_results(studycase,name_studycase,general['simulation years'],'prosumer_1', emissions, '2025',120,125, print_ = False, plot = True, plot_cum = False, ghg = False,save=False)
balances_pp = pp.renewables_results(studycase,name_studycase,general['simulation years'],'prosumer_1', emissions, '2025',220,240, print_ = False, plot = True, plot_cum = False, ghg = False,save=False)

LCOH('prosumer_1',balances_pp,studycase,name_studycase,energy_market,general['simulation years'],path,name_economic, revenues = False, refund = True)

pp.ele_param(name_studycase, 2, 3)
pp.fc_param(name_studycase, 2, 3)

pp.LOC_plot(name_studycase)

pp.NPV_plot(name_economic)


#%% ##########################################################################
"Sensitivity analysis - practical example"
# click on the following code and press Ctrl+5 to discomment all together

# =============================================================================
# import numpy as np
# import matplotlib.pyplot as plt
# 
# pv_size = np.arange(1,11)
# npv = []
# sc = [] # self-consumption
# ss = [] # self-sufficiency
# 
# for pv in pv_size:
#     print(pv)
#     name_studycase = f"PV size = {pv}"
#     new_studycase = pre.change_peakP(studycase, 'prosumer_1', pv) # change PV size 
#     rec = REC(new_studycase,general,file_studycase,file_general,path) # create REC object
#     rec.REC_energy_simulation() # simulate REC enegy balances
#     rec.tech_cost(tech_cost) # calculate the cost of all technologies 
#     rec.save(name_studycase) # save results in 'name_studycase.pkl'
#     
#     with open('results/balances_'+name_studycase+'.pkl', 'rb') as f: balances = pickle.load(f)
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










