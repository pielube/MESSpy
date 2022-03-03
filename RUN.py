from REC import REC
from Economy import NPV
import Post_processing as pp
from Input import*
import time

" RECS Renewable Energy Communities Simulator "

study_case = 'study case' # str name for results file.pkl
reference_case = 'reference case' # str name for results file.pkl

#%% study case simulation 
time1 = time.time()  

rec = REC(structure,simulation_years) # create REC structure
#rec.reset() # reset REC energy balances

time2 = time.time() 
 
rec.REC_energy_simulation() # simulate REC structure
rec.save(study_case) # save results in 'study_case.pkl'
pp.total_balances(study_case)
# pp.SOC_plot(study_case)

time3 = time.time()  

#%% reference case simulation (run only if changed)
rec0 = REC(structure0,simulation_years) # create REC
rec0.REC_energy_simulation() # simulate REC 
rec0.save(reference_case) # save results in 'reference_case.pkl'

#%% economic assesment (It has no sense if simulation_years = 1)
NPV(structure,structure0,study_case,reference_case,economic_data,simulation_years) 
pp.NPV_plot()

time4 = time.time()  

print('Execution time {:.2f}'.format(time2-time1))
print('Execution time {:.2f}'.format(time3-time2))
print('Execution time {:.2f}'.format(time4-time3))

###################################################### work in progress...
### version 1.0
# togliere vettori inutile che si possono calcolare a posteriori (es. self.consumption)
# inserire bilanci input output singole tech
# economic assessment
    # add not float electricity price (series)
    # add redistribution of REC incentives
# add battery ageing 
# update / totally rewrite Electolyzer and Fuel cell modells 
# general production from .csv? Usefull f.i for wind production

### TO SPPED UP
# parametric anlysis and optimizations
    # don't recreate REC object every time but change only the parameter of interest and then resimulate the REC
        # use REC.reset() 
    # don't resimulate REC object if the parameter of interest is an economic one
# Casching
    # pvlib standardised serie (it changes only is a PV parameter different from Ppeak changes)
    # load
    # reference_case
# if every year is the same, f.i without battery aging
    # simulate only one year and duplicate it for the economic analysis
# parallelizations of locationss that can be solved simultaneously

### version 1.0S (for students)
# .exe
# input from yaml/jason/xml
# instruction
# idea: ask them to find load info, give them only balances without economic analysis. 

### version 2.0 and following
# HP Guglielmo
# enviromental analysis
# add new logic to simulate smart batteries and collective batteries (single location needs information of entire REC, simulation order?)
# graphic and web interface
# add gas and heat/cool components ...
# add warming.py
# add test.py
# different timestep? (collective self-consumption is defined by law on an hourly basis)







