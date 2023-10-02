#%% ###########################################################################
"""
PRE PROCESSING
==============
"""
# Import modules
from rec import REC
from economics import NPV
import preprocess_CND as pre
import postprocess_CND as ppdev
import postprocess_test as pp
import os
import json
import pickle
import pandas as pd
import matplotlib.dates as mdates  
from matplotlib.dates import MonthLocator, DateFormatter
import matplotlib.pyplot as plt

# Selecting simulation names
name_studycase = 'Studicase_GEC_20kW_pv' # str name for energy_balances_results file.pkl
name_refcase = 'Refcase_GEC_8kW_pv' # str name for energy_balances_results file.pkl
name_economic = 'From Refcase to Studycase' # str name for economic_assesment_results file.pkl

# Selecting input files:
path = r'./input_CND' # change the path with r'./input_dev' if you are working on your own run_dev
file_studycase = 'studycase'
file_refcase = 'refcase'
file_general = 'general'
file_tech_cost = 'tech_cost'
file_energy_market  = 'energy_market'

# Opening input files:
with open(os.path.join(path,f"{file_studycase}.json"),'r') as f: studycase = json.load(f)
with open(os.path.join(path,f"{file_refcase}.json"),'r') as f: refcase = json.load(f)
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
SOLVER - refcase simulation
================================
"""
rec0 = REC(refcase,general,file_refcase,file_general,path) # create REC object
rec0.REC_energy_simulation() # simulate REC 
rec0.tech_cost(tech_cost) # calculate the cost of all technologies 
rec0.save(name_refcase) # save results in 'name_refcase.pkl'

#%% ###########################################################################
"""
POST PROCESS - Investment assessment
================================
"""
# Net present value calculation to asses the investment comparing refcase and studycase
NPV(file_studycase,file_refcase,name_studycase,name_refcase,energy_market,general['simulation years'],path,name_economic)

#%% ###########################################################################
"""
POST PROCESS - PLOTTING
================================
some post-process are alredy avaiable as examples in postprocess_test
you should create your own postprocess_dev.py
"""
pp.total_balances(name_studycase,'GEC','electricity')

pp.total_balances(name_refcase,'GEC','electricity')

pp.REC_electricity_balance(name_studycase)
 
pp.hourly_balances_electricity(name_studycase,'GEC', 20, 21)

pp.hourly_balances_electricity(name_studycase,'GEC', 2, 3)

# pp.csc_allocation_sum(name_studycase)
# pp.storage_control(name_studycase)

pp.NPV_plot(name_economic)

#%% ##########################################################################
"Sensitivity analysis"

scenari = []
for pv in [25,30,35]:
    name_studycase = f"{pv} kWp"
    scenari.append(name_studycase)
    with open(os.path.join(path,f"{file_studycase}.json"),'r') as f: studycase = json.load(f)
    
    studycase = pre.change_peakP(studycase, 'GEC', pv)
    
    rec = REC(studycase,general,file_studycase,file_general,path) # create REC object
    rec.REC_energy_simulation() # simulate REC enegy balances
    rec.tech_cost(tech_cost) # calculate the cost of all technologies 
    rec.save(name_studycase) # save results in 'name_studycase.pkl'
    
    name_economic = name_studycase # str name for economic_assesment_results file.pkl
    file_energy_market = 'energy_market'
    with open(os.path.join(path,f"{file_energy_market}.json"),'r') as f: energy_market = json.load(f)
    NPV(file_studycase,file_refcase,name_studycase,name_refcase,energy_market,general['simulation years'],path,name_economic) 

# grafici

name_studycase = scenari[0] 
ppdev.typedays(name_studycase,[19,164,197,287],['17 gennaio','17 aprile','15 luglio','15 ottobre'],30)  
ppdev.NPV_scenarios(scenari, ['GEC'], 'Valore attuale netto')
df = ppdev.bilanci_scenari([name_refcase]+scenari)
print(df)
#df.to_excel("bilanci energetici.xlsx")
df2 = ppdev.flussi_di_cassa_scenari(scenari)
print(df2) 
#df2.to_excel("flussi di cassa.xlsx")
       
#a = pp.REC_electricity_balance(name_studycase)
for name_studycase in [name_refcase]+scenari:
    ppdev.hist_12_balances_pc(name_studycase,7000)   
    
with open('results/economic_assessment_'+name_economic+'.pkl', 'rb') as f: economic = pickle.load(f)
with open('results/balances_'+name_studycase+'.pkl', 'rb') as f: balances = pickle.load(f)


#%% self consumption series (to match with EV series) 

# e-mobility staff load data reading 
path = r'./input_CND/loads'
ev_staff_load = r'e_mobility_load_staff'
csv_file_path = os.path.join(path, f"{ev_staff_load}.csv")
df = pd.read_csv(csv_file_path)
df.set_index('Time stamp', inplace=True)
df.rename(columns={'ee897394-71fd-4179-8056-4787a0b0c225':'Bay 8','8733ec0a-c87b-4e98-8202-5a2ce395332c':'Bay 13'}, inplace= True)
df['Total staff']= df.sum(axis=1)

# location self-consumption
demand = -balances['GEC']['electricity']['demand']
from_grid = balances['GEC']['electricity']['grid']
from_grid[from_grid<0] = 0
self_consumption = demand-from_grid

  
fig, ax = plt.subplots(figsize=(8, 6), dpi=600)
ax.plot(df.index,df['Total staff'], linewidth=0.6)
ax.plot(df.index,self_consumption, linewidth=0.6)
# ax.set_title('Power [kW]')
ax.set_xlabel('Time [h]')
ax.set_ylabel('Production [kW]')
ax.xaxis.set_major_locator(MonthLocator())
ax.xaxis.set_major_formatter(DateFormatter('%b'))
# plt.xticks(rotation=45) 
ax.set_xlim(min(df.index), max(df.index))
ax.set_ylim(0, None)
plt.grid(alpha=0.3)
plt.show()



