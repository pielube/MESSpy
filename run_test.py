from rec import REC
import time
import os
import json

#%%

"""
MESSpy - Run_test

doon't work on this script:
    you should create your own run_dev.py and input_dev/ 
"""

path = r'./input_test' # change the put with r'./input_dev' if you are working on your own run_dev

study_case = 'REC_test' # str name for results file.pkl
reference_case = 'buiseness as usual' # str name for results file.pkl


"""
Input files
"""

file = 'structure.json'
filepath = os.path.join(path,file)
with open(filepath,'r') as f:
    structure = json.load(f)

file = 'general.json'
filepath = os.path.join(path,file)
with open(filepath,'r') as f:
    general = json.load(f)

time1 = time.time()
print('Creating structure..')

# Creating initial structure
rec = REC(structure,general,path) # create REC structure

time2 = time.time()
print('Structure created in {:.2f} seconds'.format(time2-time1))

#%% ###########################################################################
print('Running the model..')
time2 = time.time()

# Running the model
#rec.reset() # reset REC energy balances
rec.REC_energy_simulation() # simulate REC structure
rec.save(study_case) # save results in 'study_case.pkl'

time3 = time.time()
print('Model runned in {:.2f} seconds'.format(time3-time2))
  
#%% ###########################################################################
from economics import NPV
print('Economic analysis..') 
time3 = time.time()

file = 'economics.json'
filepath = os.path.join(path,file)
with open(filepath,'r') as f:
    economic_data = json.load(f)

file = 'refcase.json'
filepath = os.path.join(path,file)
with open(filepath,'r') as f:
    structure0 = json.load(f)
    
# Reference case simulation (run only if changed)
rec0 = REC(structure0,general,path) # create REC
rec0.REC_energy_simulation() # simulate REC 
rec0.save(reference_case) # save results in 'reference_case.pkl'

# Actual economic analysis (It has no sense if simulation_years = 1)
NPV(structure,structure0,study_case,reference_case,economic_data,general['simulation years'],path) 

time4 = time.time()  
print('Eonomic analysis performend in {:.2f} seconds'.format(time4-time3))

#%% post process
import postprocess_test as pp
#import postprocess_dev as pp

"""
some post-process are alredy avaiable as examples in postprocess_test
you should create your own postprocess_dev.py
"""

print('Post processing..')
time4 = time.time()

#pp.total_balances(study_case, 'p1')
#pp.total_balances(study_case,'p2')
#pp.total_balances(study_case,'c1')

pp.LOC_plot(study_case)

pp.NPV_plot(study_case)

pp.hourly_balances(study_case,'p1', 2, 3)
pp.hourly_balances(study_case,'p2', 2, 3)
pp.hourly_balances(study_case,'c1', 2, 3)

pp.csc_allocation_sum(study_case)

pp.Flows(study_case) # if it doesn't work try to open the file.html directrly from the results/ folder
pp.Flows(study_case) # if it doesn't work try to open the file.html directrly from the results/ folder

time5 = time.time()  
print('Post process performend in {:.2f} seconds'.format(time5-time4))



