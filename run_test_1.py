"""
MESSpy - Run_test_1 - Assessing a small Energy Community composed by 2 consumers and one prosumer with PV and battery.

We advice against working on this script:
    you should create your own run_xxx.py, input_xxx/, postprocess_xxx.py and preprocess_xx.py
"""

#%% ###########################################################################
"""
PRE PROCESSING
==============
"""

# Import modules: don't change it
from core import rec
from core import economics as eco
import os
import json
import pickle

# Import modules: change it if you want to run your your own simulation
import preprocess_test as pre
import postprocess_test as pp
#import postprocess_xxx as pp
#import preprocess_xxx as pre

# Selecting input files folder: change it if you want to run your your own simulation
path = r'./input_test_1' 
#path = r'./input_xxx'

# Selecting simulation names: change them to your liking
name_studycase = 'Post' # just a name for saving results of the studycase simulation
name_refcase = 'Pre' # just a name for saving results of the refcase simulation
name_economic = 'Post vs Pre' # just a name for saving results of economic evaluation

# Selecting the names of the .json files to be read as simulations input: change it if you want to run your your own simulation
file_tech_cost      = 'tech_cost'
file_energy_market  = 'energy_market'
file_general        = 'general'
file_studycase      = 'studycase'
file_refcase        = 'refcase'

# If interested in testing simulations with 15 minutes timestep resolution:
# =============================================================================
# # file_general        = 'general15' 
# # file_studycase      = 'studycase15' 
# # file_refcase        = 'refcase15'
# =============================================================================

### NOW you can Run file (F5) ###

# Opening input files:
with open(os.path.join(path,f"{file_studycase}.json"),'r')      as f: studycase     = json.load(f)
with open(os.path.join(path,f"{file_refcase}.json"),'r')        as f: refcase       = json.load(f)
with open(os.path.join(path,f"{file_general}.json"),'r')        as f: general       = json.load(f)
with open(os.path.join(path,f"{file_tech_cost}.json"),'r')      as f: tech_cost     = json.load(f)
with open(os.path.join(path,f"{file_energy_market}.json"),'r')  as f: energy_market = json.load(f)


#%% ###########################################################################
"""
SOLVER - studycase simulation
======
"""

sim = rec.REC(studycase,general,file_studycase,file_general,path) # create REC object
sim.REC_power_simulation() # simulate REC power balances
sim.tech_cost(tech_cost) # calculate the cost of all technologies 
sim.save(name_studycase,'pkl') # saving results in .pkl format (useful if you make postprocess using python)
sim.save(name_studycase,'csv',sep=';',dec=',') # saving results in .csv (useful if you make postprocess using other languages/programmes)
    
#%% ###########################################################################
"""
SOLVER - refcase simulation
================================
"""
  
sim0 = rec.REC(refcase,general,file_refcase,file_general,path) # create REC object
sim0.REC_power_simulation() # simulate REC power balances
sim0.tech_cost(tech_cost) # calculate the cost of all technologies 
sim0.save(name_refcase,'pkl') # saving results in .pkl format (useful if you make postprocess using python)
sim0.save(name_refcase,'csv',sep=';',dec=',') # saving results in .csv (useful if you make postprocess using other languages/programmes)

#%% ###########################################################################
"""
POST PROCESS - Investment assessment comparing refcase and studycase
================================
"""
# Net present value calculation to asses the investment comparing refcase and studycase (saves results in both .pkl and .csv)
eco.NPV(file_studycase,file_refcase,name_studycase,name_refcase,energy_market,path,name_economic,'pkl')
eco.NPV(file_studycase,file_refcase,name_studycase,name_refcase,energy_market,path,name_economic,'csv',sep=';',dec=',')

#%% ###########################################################################
"""
POST PROCESS - PLOTTING
================================
some post-process are alredy avaiable as examples in postprocess_test_1
you should create your own postprocess_dev.py and create your own graphs
"""

# Here the main simulation results are read: balances, balances0 and economic. Balances are also available in the Variable Explorer panel: sim and sim0. Results are also available in .csv format in results/csv folder.
with open('results/pkl/balances_'+name_studycase+'.pkl', 'rb')              as f: balances  = pickle.load(f)
with open('results/pkl/balances_'+name_refcase+'.pkl', 'rb')                as f: balances0 = pickle.load(f)
with open('results/pkl/economic_assessment_'+name_economic+'.pkl', 'rb')    as f: economic  = pickle.load(f)

# Here some examples of graphs
pp.location_balance(name_refcase,'prosumer','electricity')
pp.location_balance(name_studycase,'prosumer','electricity')
pp.location_balance(name_refcase,'consumer 1','electricity')
pp.location_balance(name_studycase,'consumer 1','electricity')
#pp.location_balance(name_refcase,'consumer 2','electricity')
#pp.location_balance(name_studycase,'consumer 2','electricity')

pp.hourly_balances_electricity(name_studycase,'prosumer', 2, 3)
pp.hourly_balances_electricity(name_studycase,'consumer 1', 2, 3)

b = pp.REC_electricity_balance(name_studycase)
pp.hist_12_balances_pc(name_studycase,1100)

pp.csc_allocation_sum(name_studycase)
pp.LOC_plot(name_studycase)
pp.NPV_plot(name_economic)
    

#%% ##########################################################################
"Sensitivity analysis - practical example"
# Run MESS several times as the installed PV power varies

import numpy as np
import matplotlib.pyplot as plt

pv_size = np.arange(1,11)

npv = []
sc = [] # self-consumption
ss = [] # self-sufficiency

print('\n Sensitivity analysis running:')
for pv in pv_size:
    print(pv)
    name_studycase = f"PV size = {pv}"
    new_studycase = pre.change_peakP(studycase, 'prosumer', pv) # change PV size 
    cs = rec.REC(new_studycase,general,file_studycase,file_general,path) # create REC object
    cs.REC_power_simulation() # simulate REC enegy balances
    cs.tech_cost(tech_cost) # calculate the cost of all technologies 
    cs.save(name_studycase,'pkl') # save results in 'name_studycase.pkl'
    
    with open('results/pkl/balances_'+name_studycase+'.pkl', 'rb') as f: balances = pickle.load(f)
    demand = -balances['prosumer']['electricity']['demand'].sum() # read from saved results .pkl
    production = balances['prosumer']['electricity']['PV'].sum()
    into_grid = balances['prosumer']['electricity']['grid'].sum(where=balances['prosumer']['electricity']['grid']<0)
    from_grid = balances['prosumer']['electricity']['grid'].sum(where=balances['prosumer']['electricity']['grid']>0)
    
    # you can also read values from the python object rec. (in this case you do not need rec.save)
    # demand = -rec.locations['prosumer_1'].energy_balance['electricity']['demand'].sum() # read from python 
    
    sc.append((production+into_grid)/production*100)
    ss.append((demand-from_grid)/demand*100)
    
plt.figure(dpi=1000)
plt.plot(pv_size,sc,label='Self-consumption')
plt.plot(pv_size,ss,label='Self-sufficiency')
plt.xlabel("PV peak power [kWp]")
plt.ylabel("[%]")
plt.grid()
plt.legend()
plt.title('prosumer_1')
plt.show()










