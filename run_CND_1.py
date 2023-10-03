#%% ###########################################################################
"""
PRE PROCESSING
==============
"""
# Import modules
from rec import REC
from economics import NPV
import preprocess_CND as pre
import postprocess_CND as pp
import numpy as np
import os
import json
import pickle
import pandas as pd
import matplotlib.dates as mdates  
from matplotlib.dates import MonthLocator, DateFormatter
import matplotlib.pyplot as plt

# Selecting input files:
path = r'./input_CND' # change the path with r'./input_dev' if you are working on your own run_dev
file_studycase = 'refcase'
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

#%% ##########################################################################
"RUN Sensitivity analysis"

scenarios = []
pv_sizes = [8,20,30]
#pv_sizes = np.arange(8,31,2)
bess_sizes = [0,30]

for pv in pv_sizes:
    for bess in bess_sizes:

        with open(os.path.join(path,f"{file_studycase}.json"),'r') as f: studycase = json.load(f)
        
        if pv > 8: # if pv >8, change studycase peakP
            studycase = pre.change_peakP(studycase, 'GEC', pv)
        name_studycase = f"{pv} kWp"
        
        if bess >0: # if there is a BESS, add it to the studycase 
           studycase = pre.add_battery(studycase,'GEC',bess)
           name_studycase += f" {bess} kWh"
            
        scenarios.append(name_studycase)
        
        rec = REC(studycase,general,file_studycase,file_general,path) # create REC object
        rec.REC_energy_simulation() # simulate REC enegy balances
        rec.tech_cost(tech_cost) # calculate the cost of all technologies 
        rec.save(name_studycase) # save results in 'name_studycase.pkl'
              
        
#%% graphs and post-processing

df = pp.bilanci_scenarios(scenarios)
print(df)
#df.to_excel("bilanci energetici.xlsx")
       
for name_studycase in scenarios:
    pp.hist_12_balances_pc(name_studycase,6200)   
for name_studycase in scenarios:  
    pp.typedays(name_studycase,[19,164,197,287],['17 January','17 April','15 July','15 October'],32) 


#%% Self-Sufficiency and Self-Consumption vs PV peakP. With and without battery.

# preparing data for plot: calculate ss e sc of each scenarios and divide scenarios with and withous battery
SS = {}
SC = {}
for s in scenarios:
    res = pp.REC_electricity_balance(s,noprint=True)
    SS[s] = res['Value / demand [%]']['SC']
    SC[s] = res['Value / production [%]']['SC']

ss = []
sc = []
ssb = []
scb = []
for s in scenarios:
    if 'kWh' in s:
        ssb.append(SS[s])
        scb.append(SC[s])
    else:
        ss.append(SS[s])
        sc.append(SC[s])
        
plt.figure(dpi=1000)
plt.plot(pv_sizes,sc,label='SC')
plt.plot(pv_sizes,scb,label='SC with battery')
plt.plot(pv_sizes,ss,label='SS')
plt.plot(pv_sizes,ssb,label='SS with battery')
plt.ylim(0,100)
plt.xlim(8,30)
plt.xticks([8,12,16,20,24,28],['8', '12', '16', '20', '24', '28'])
plt.grid()
plt.xlabel("PV peak power [kWp]")
plt.ylabel('SC e SS [%]')
plt.legend()
plt.title("Self-Consumption and Self-Sufficiency")
plt.show()


#%% Assessing EV self-sufficiency

# e-mobility staff load data reading 
path = r'./input_CND/loads'
ev_staff_load = r'e_mobility_load_staff'
csv_file_path = os.path.join(path, f"{ev_staff_load}.csv")
df = pd.read_csv(csv_file_path)
df.set_index('Time stamp', inplace=True)
df.rename(columns={'ee897394-71fd-4179-8056-4787a0b0c225':'Bay 8','8733ec0a-c87b-4e98-8202-5a2ce395332c':'Bay 13'}, inplace= True)
df['Total staff']= df.sum(axis=1)
ev = np.array(df['Total staff'])

# location self-consumption
ss = []
ssb = []
for s in scenarios:
    name_studycase = s
    with open('results/balances_'+name_studycase+'.pkl', 'rb') as f: balances = pickle.load(f)
    demand = -balances['GEC']['electricity']['demand']
    from_grid = balances['GEC']['electricity']['grid']
    from_grid[from_grid<0] = 0
    self_consumption = demand-from_grid
    ev_ss = np.minimum(ev,self_consumption)
    ev_ss_index = round(ev_ss.sum()/ev.sum() *100,2)
    if 'kWh' in s:
        ssb.append(ev_ss_index)
    else:
        ss.append(ev_ss_index)

plt.figure(dpi=1000)
plt.plot(pv_sizes,ss,label='SS')
plt.plot(pv_sizes,ssb,label='SS with battery')
plt.ylim(0,100)
plt.xlim(8,30)
plt.grid()
plt.xlabel("PV peak power [kWp]")
plt.ylabel('SS [%]')
plt.xticks([8,12,16,20,24,28],['8', '12', '16', '20', '24', '28'])
plt.legend()
plt.title("EV Self-Sufficiency")
plt.show()    
  

