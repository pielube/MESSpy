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
from economics import NPV,LCOH
import postprocess_test as pp
import preprocess_test as pre
# import postprocess_dev as pp
# import preprocess_dev as pre
import os
import json
import pickle
import pandas as pd


# Selecting simulation names
name_studycase  = 'Industrial Facility Hydrogen' # str name for energy_balances_results file.pkl
name_refcase    = 'Industrial Facility Fossil Fuels' # str name for energy_balances_results file.pkl
name_economic   = 'From Fossil Fuels to H2' # str name for economic_assesment_results file.pkl

path = r'./input_test' # change the path with r'./input_dev' if you are working on your own run_dev
# path = r'./input_dev'

file_studycase      = 'studycase'
file_refcase        = 'refcase'
file_general        = 'general'
file_tech_cost      = 'tech_cost'
file_energy_market  = 'energy_market'
file_emissions      = 'emissions'

# Opening input files:
with open(os.path.join(path,f"{file_refcase}.json"),'r') as f: refcase = json.load(f)
with open(os.path.join(path,f"{file_general}.json"),'r') as f: general = json.load(f)
with open(os.path.join(path,f"{file_emissions}.json"),'r') as f: emissions = json.load(f)
with open(os.path.join(path,f"{file_studycase}.json"),'r') as f: studycase = json.load(f)
with open(os.path.join(path,f"{file_energy_market}.json"),'r') as f: energy_market = json.load(f)
with open(os.path.join(path,f"{file_tech_cost}.json"),'r') as f: tech_cost = json.load(f)


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
# """
# POST PROCESS - refcase simulation and investment assessment
# ================================
# """
  
# Reference case simulation (run only if changed)
# rec0 = REC(refcase,general,file_refcase,file_general,path) # create REC object
# rec0.REC_energy_simulation() # simulate REC 
# rec0.tech_cost(tech_cost) # calculate the cost of all technologies 
# rec0.save(name_refcase) # save results in 'name_refcase.pkl'


#%% Net present value calculation to asses the investment comparing refcase and studycase
# NPV(name_studycase,name_refcase,energy_market,general['simulation years'],path,name_economic) 


#%% ###########################################################################
"""
POST PROCESS - PLOTTING
================================

some post-process are alredy avaiable as examples in postprocess_test
you should create your own postprocess_dev.py
"""

pp.total_balances(name_studycase,'Industrial_Facility_1','electricity')
pp.total_balances(name_studycase,'Industrial_Facility_1','hydrogen')

#pp.total_balances(name_refcase,'consumer_2','electricity')
#pp.total_balances(name_refcase,'consumer_2','heating water')
#pp.total_balances(name_refcase,'consumer_2','gas')
#pp.total_balances(name_studycase,'consumer_2','electricity')
#pp.total_balances(name_studycase,'consumer_2','heating water')

# pp.REC_electricity_balance(name_studycase)

pp.LOC_plot(name_studycase)
pp.hydrogen_production(name_studycase, 'Industrial_Facility_1')

# pp.NPV_plot(name_economic)

pp.renewables(name_studycase,general['simulation years'],'Industrial_Facility_1',190,200)#,plot=True)
LCOH(studycase,name_studycase,energy_market,general['simulation years'],path,name_economic)#, revenues = ['oxygen','electricity']) 
pp.ghg_emissions(name_studycase,'Industrial_Facility_1', emissions, '2025')
pp.hourly_balances_electricity(name_studycase,'Industrial_Facility_1',40,85)

# pp.storage_control(name_studycase)
#pp.ele_param(name_studycase, 100, 150, plot=True)
   
# # here you can read the main results (balances and economic):
# #with open('results/economic_assessment_'+name_economic+'.pkl', 'rb') as f: economic = pickle.load(f)
# with open('results/balances_'+name_studycase+'.pkl', 'rb') as f: balances = pickle.load(f)

#%% ##########################################################################
"Sensitivity analysis - practical example"
# click on the following code and press Ctrl+5 to uncomment all together

# =============================================================================
# import numpy as np
# import matplotlib.pyplot as plt
# import pandas as pd
# 
# intervals = 5
# 
# electrolyzer_size = np.round(np.linspace(50,150,intervals)).astype(int)
# wind_size         = np.linspace(50000,150000,intervals)
# pv_size           = np.linspace(25000,75000,intervals)
# #battery_size     = np.linspace(25000,75000,intervals)
# #tank_size        = np.linspace(250000,600000,intervals)
# 
# # Creating dictionaries to store parameters of the sensitivity analisys
# 
# lcoh    = {} # levelized cost of hydrogen
# ghg     = {} # emission intensity
# cf      = {} # electorlyzers capacity factor
# 
# for wind in wind_size:
#     
#     print('Wind Farm Size '+str(wind/1000)+ ' MW' )
#     pre.change_peakW(studycase, 'Industrial_Facility_1', wind) # change Wind size 
#     lcoh[wind]  = []
#     ghg[wind]   = []
#     cf[wind]    = []
#     # for pv in pv_size:
#     for ele in electrolyzer_size:   
#         
#         print('Electrolyzer Size '+str(ele)+ ' MW' )
# 
#         # name_studycase = f"Wind size {wind} PV size = {pv}"
#         name_studycase = f"Wind size {wind} Ele size = {ele}"
#         new_studycase = pre.change_Elesize(studycase, 'Industrial_Facility_1', ele) # change PV size 
#         # new_studycase = pre.change_peakP(studycase, 'Industrial_Facility_1', pv) # change PV size 
#         rec = REC(new_studycase,general,file_studycase,file_general,path) # create REC object
#         rec.REC_energy_simulation() # simulate REC enegy balances
#         rec.tech_cost(tech_cost) # calculate the cost of all technologies 
#         rec.save(name_studycase) # save results in 'name_studycase.pkl'
#         
#         pp.renewables(name_studycase,general['simulation years'],'Industrial_Facility_1',198,200)
#         lcoh[wind].append(LCOH(studycase,name_studycase,energy_market,general['simulation years'],path,name_economic,plot=False)) 
#         ghg[wind].append(pp.ghg_emissions(name_studycase,'Industrial_Facility_1', emissions, '2025'))
#         cf[wind].append(pp.ele_param(name_studycase,2, 100))
#         
#         
#         with open('results/balances_'+name_studycase+'.pkl', 'rb') as f: balances = pickle.load(f)
#         # demand = -balances['Industrial_Facility_1']['electricity']['demand'].sum() # read from saved results .pkl
#         # production = balances['Industrial_Facility_1']['electricity']['PV'].sum()
#         # into_grid = balances['Industrial_Facility_1']['electricity']['grid'].sum(where=balances['Industrial_Facility_1']['electricity']['grid']<0)
#         # from_grid = balances['Industrial_Facility_1']['electricity']['grid'].sum(where=balances['Industrial_Facility_1']['electricity']['grid']>0)
#         
#         # you can also read values from the python object rec. (in this case you do not need rec.save)
#         # demand = -rec.locations['Industrial_Facility_1'].energy_balance['electricity']['demand'].sum() # read from python 
#         
#         # sc.append((production+into_grid)/production*100)
#         # ss.append((demand-from_grid)/demand*100)
#         print()
# 
# #%%
# 
# df = pd.DataFrame.from_dict(lcoh, orient='index', columns=electrolyzer_size)
# df.index.name = 'Wind size'
# # df4 = pd.DataFrame.from_dict(cf, orient='index', columns=pv_size)
# df1 = pd.DataFrame.from_dict(ghg, orient='index', columns=electrolyzer_size)
# df2 = pd.DataFrame.from_dict(cf, orient='index', columns=electrolyzer_size)
# 
#        
# 'LCOH'
# # # Creating 2-D grid of features
# [X, Y] = np.meshgrid(df.columns, df.index)
# fig, ax = plt.subplots(dpi=600, figsize=(7,6))
# Z = df.values
# # contourf_ = ax.contourf(X/1000, Y/1000, Z ,cmap='Blues')
# contourf_ = ax.contourf(X, Y/1000, Z ,cmap='Blues')
# # ax.set_title(str(List_CGS[k]) + ' - PI ()')
# # ax.set_xlabel('PV farm size [MW]')
# ax.set_xlabel('Electorlyser size [MW]')
# ax.set_ylabel('Wind farm size [MW]')
# cbar = fig.colorbar(contourf_) #, orientation='horizontal')
# cbar.set_label(r'$\mathbf{LCOH}$ [€/kgH$_\mathregular{2}$]')
# 
# 'GHG'
# [X, Y] = np.meshgrid(df1.columns, df1.index)
# fig, ax = plt.subplots(dpi=600, figsize=(7,6))
# Z = df1.values
# # contourf_ = ax.contourf(X/1000, Y/1000, Z ,cmap='Blues')
# contourf_ = ax.contourf(X, Y/1000, Z ,cmap='Greens')
# # ax.set_title(str(List_CGS[k]) + ' - PI ()')
# # ax.set_xlabel('PV farm size [MW]')
# ax.set_xlabel('Electorlyser size [MW]')
# ax.set_ylabel('Wind farm size [MW]')
# cbar = fig.colorbar(contourf_)#, orientation='horizontal')
# cbar.set_label(r'$\mathbf{GHG}$ [kgCO$_\mathregular{2}$/kgH$_\mathregular{2}$]')
# 
# 'Capacity factor'
# [X, Y] = np.meshgrid(df2.columns, df2.index)
# fig, ax = plt.subplots(dpi=600, figsize=(7,6))
# Z = df2.values
# # contourf_ = ax.contourf(X/1000, Y/1000, Z ,cmap='Blues')
# contourf_ = ax.contourf(X, Y/1000, Z ,cmap='Oranges')
# # ax.set_title(str(List_CGS[k]) + ' - PI ()')
# # ax.set_xlabel('PV farm size [MW]')
# ax.set_xlabel('Electorlyser size [MW]')
# ax.set_ylabel('Wind farm size [MW]')
# cbar = fig.colorbar(contourf_)#, orientation='horizontal')
# cbar.set_label(r'Capacity Factor $_\mathregular{ele}$  [%]')
# 
# 
# # [X, Y] = np.meshgrid(df4.columns, df4.index)
# # fig, ax = plt.subplots(dpi=600, figsize=(7,6))
# # Z = df4.values
# # # contourf_ = ax.contourf(X/1000, Y/1000, Z ,cmap='Blues')
# # contourf_ = ax.contourf(X/1000, Y/1000, Z ,cmap='Oranges')
# # # ax.set_title(str(List_CGS[k]) + ' - PI ()')
# # # ax.set_xlabel('PV farm size [MW]')
# # ax.set_xlabel('PV field size [MW]')
# # ax.set_ylabel('Wind farm size [MW]')
# # cbar = fig.colorbar(contourf_)#, orientation='horizontal')
# # cbar.set_label(r'Capacity Factor $_\mathregular{ele}$  [%]')
# 
# 
# 
# =============================================================================






