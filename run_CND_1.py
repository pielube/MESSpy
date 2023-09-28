#%% ###########################################################################

# Import modules
from rec import REC
from economics import NPV
import preprocess_CND as pre
import os
import json
import pickle
import pandas as pd

# Selecting input files:
path = r'./input_CND' # change the path with r'./input_dev' if you are working on your own run_dev
file_studycase = 'refcase'
file_refcase = 'refcase'
file_general = 'general'
file_tech_cost = 'tech_cost'

# Opening input files:
with open(os.path.join(path,f"{file_studycase}.json"),'r') as f: studycase = json.load(f)
with open(os.path.join(path,f"{file_refcase}.json"),'r') as f: refcase = json.load(f)
with open(os.path.join(path,f"{file_general}.json"),'r') as f: general = json.load(f)
with open(os.path.join(path,f"{file_tech_cost}.json"),'r') as f: tech_cost = json.load(f)

# refcase simulation
name_refcase = 'Without PV' # str name for energy_balances_results file.pkl
rec0 = REC(refcase,general,file_refcase,file_general,path) # create REC object
rec0.REC_energy_simulation() # simulate REC 
rec0.tech_cost(tech_cost) # calculate the cost of all technologies 
rec0.save(name_refcase) # save results in 'name_refcase.pkl'

# studycase simulation:

scenari = []
for pv in [20,25,30]:
    name_studycase = f"{pv} kWp"
    scenari.append(name_studycase)
    with open(os.path.join(path,f"{file_studycase}.json"),'r') as f: studycase = json.load(f)
    
    studycase = pre.add_PV(studycase, 'GEC', pv, 10 , 0)
    
    rec = REC(studycase,general,file_studycase,file_general,path) # create REC object
    rec.REC_energy_simulation() # simulate REC enegy balances
    rec.tech_cost(tech_cost) # calculate the cost of all technologies 
    rec.save(name_studycase) # save results in 'name_studycase.pkl'
    
    name_economic = name_studycase # str name for economic_assesment_results file.pkl
    file_energy_market = 'energy_market'
    with open(os.path.join(path,f"{file_energy_market}.json"),'r') as f: energy_market = json.load(f)
    NPV(file_studycase,file_refcase,name_studycase,name_refcase,energy_market,general['simulation years'],path,name_economic) 

# grafici
import postprocess_CND as pp

name_studycase = scenari[0] 
pp.typedays(name_studycase,[19,164,197,287],['17 gennaio','17 aprile','15 luglio','15 ottobre'],30)  
df = pp.bilanci_scenari(scenari)
pp.NPV_scenarios(scenari, ['GEC'], 'Valore attuale netto')
print(df)
#df.to_excel("bilanci energetici.xlsx")
df2 = pp.flussi_di_cassa_scenari(scenari)
print(df2) 
#df2.to_excel("flussi di cassa.xlsx")
       
#a = pp.REC_electricity_balance(name_studycase)
for name_studycase in scenari:
    pp.hist_12_balances_pc(name_studycase,6000)   
    
with open('results/economic_assessment_'+name_economic+'.pkl', 'rb') as f: economic = pickle.load(f)
with open('results/balances_'+name_studycase+'.pkl', 'rb') as f: balances = pickle.load(f)