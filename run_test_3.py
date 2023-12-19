"""
MESSpy - Run_test_3 - Analysis of an illustrative case study involving hydrogen technologies

don't work on this script:
    you should create your own run_dev.py, input_dev/, postprocess_dev.py and preprocess_dev.py
"""

#%% ###########################################################################
"""
PRE PROCESSING
==============
"""
# Import modules: do not edit
from core import rec
from core import economics as eco
import os
import json
import pickle

# Import modules: change if willing to use different pre- and post-process modules
import preprocess_test_3 as pre
import postprocess_test_3 as pp
#import postprocess_xxx as pp
#import preprocess_xxx as pre

# Selecting input files folder: change it if you want to run your your own simulation
path = r'./input_test_3' 
#path = r'./input_xxx'

# Selecting simulation names: change them to your liking
name_studycase  = 'Post' # illustrative name for saving the results of the studycase simulation
name_refcase    = 'Pre' # illustrative name for saving the results of the refcase simulation
name_economic   = 'Post vs Pre' # illustrative name for saving the results of economic evaluation

# Selecting the names of the .json files to be read as simulations input: change it if you want to run your your own simulation
file_tech_cost      = 'tech_cost'
file_energy_market  = 'energy_market'
file_general        = 'general'
file_studycase      = 'studycase'
file_refcase        = 'refcase'

# =============================================================================
# # If you are interested in testing 15 minutes timestep simulations:
# file_general        = 'general15' 
# file_studycase      = 'studycase15' 
# file_refcase        = 'refcase15'
# =============================================================================

### NOW you can Run the simulation (F5) ###

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
sim0.save(name_refcase,'pkl') # saving results in .pkl format (usefull if you make postprocess using python)
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
some post-process are alredy avaiable as examples in postprocess_test
you should create your own postprocess_dev.py and create your own graphs
"""

# Here the main simulation results are read: balances, balances0 and economic. Balances are also available in the Variable Explorer panel: sim and sim0. Results are also available in .csv format in results/csv folder.
with open('results/pkl/balances_'+name_studycase+'.pkl', 'rb')              as f: balances  = pickle.load(f)
with open('results/pkl/balances_'+name_refcase+'.pkl', 'rb')                as f: balances0 = pickle.load(f)
with open('results/pkl/economic_assessment_'+name_economic+'.pkl', 'rb')    as f: economic  = pickle.load(f)

# Total balances figures and hydrogen-related ghg emissions calculation in post-process. 'balance_pp': dictionary containing total balances calculation useful for LCOH calculation, NPV calculation and post process plots
balances_pp = pp.energy_balance_results(studycase,name_studycase,'industrial_facility',print_=True,plot=True)
ghg         = pp.ghg_emissions(name_studycase,'industrial_facility',energy_market,print_= True)

# Levelised Cost of Hydrogen calculation
LCOH = eco.LCOH('industrial_facility',balances_pp,studycase,name_studycase,energy_market,path,name_economic,revenues=False,refund=True,plot=True,print_=True)

# Some plot examples
pp.hydrogen_production(name_studycase,'industrial_facility')   
pp.plot_post_process(balances_pp,studycase,'industrial_facility',20,24)
pp.LOC_plot(name_studycase)
pp.NPV_plot(name_economic)
    

#%% ##########################################################################
"Sensitivity analysis - practical example"
# Run MESS several times as the installed PV power varies

# =============================================================================
# import numpy as np
# import matplotlib.pyplot as plt
# 
# intervals   = 11
# wind_size   = np.linspace(100000,300000,intervals)
# 
# # Creating lists to store parameters of the sensitivity analisys
# lcoh        = [] # levelized cost of hydrogen
# ghg         = [] # emission intensity
# 
# print('\nSensitivity analysis running:')
# print('Wind Farm Size:')
# 
# for wind in wind_size: # varying wind power size
#     print('                 '+str(int(wind/1000))+ ' MW' )
#     
#     name_studycase = f"Wind size = {wind/1000} MW"
#     new_studycase = pre.change_peakW(studycase, 'industrial_facility', wind) # change wind size 
#     cs = rec.REC(new_studycase,general,file_studycase,file_general,path) # create REC object
#     cs.REC_power_simulation() # simulate REC enegy balances
#     cs.tech_cost(tech_cost) # calculate the cost of all technologies 
#     cs.save(name_studycase,'pkl') # save results in 'name_studycase.pkl'
#     
#     with open('results/pkl/balances_'+name_studycase+'.pkl', 'rb') as f: balances = pickle.load(f)
#     balances_pp = pp.energy_balance_results(studycase,name_studycase,'industrial_facility')
#     lcoh.append(eco.LCOH('industrial_facility',balances_pp,studycase,name_studycase,energy_market,path,name_economic)) 
#     ghg.append(pp.ghg_emissions(name_studycase,'industrial_facility',energy_market))
#     
# # Plotting the results
# fig, ax1 = plt.subplots(dpi=1000)   
# ax1.plot(wind_size/1000,lcoh,label='LCOH',color='tab:blue')
# ax1.set_xlabel("Wind Farm Size [MW]")
# ax1.set_ylabel("LCOH [€/kgH$_\mathregular{2}$]", color='tab:blue')
# ax1.grid()
# ax1.set_title('industrial_facility')
# 
# ax2 = ax1.twinx()
# ax2.plot(wind_size/1000,ghg,label='GHG',color='tomato')
# ax2.set_ylabel("GHG [kgH$_\mathregular{2}$/kgCO$_\mathregular{2}$]",color='tomato')
# 
# lines, labels = ax1.get_legend_handles_labels()
# lines2, labels2 = ax2.get_legend_handles_labels()
# ax2.legend(lines + lines2, labels + labels2, loc='best')  
# 
# plt.show()
# =============================================================================







