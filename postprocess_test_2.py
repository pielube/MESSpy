"""
MESSpy - postprocessing
"""
import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from core import constants as c
from matplotlib.patches import Patch

        
def location_balance(simulation_name,loc,var=None):
    """
    Total balances figures
    
    simulationa_name : str 
    loc : str 
    var : str -> possibility of specifying a single energy carrier
    
    """
    
    with open('results/pkl/balances_'+simulation_name+'.pkl', 'rb') as f: balances = pickle.load(f)

    ###### load analysis
    
    ###### total energy balances
    
    carriers = ['electricity','heating water','gas','hydrogen','HP hydrogen']
    units    = {'electricity': 'kJ', 'hydrogen': 'kg', 'HP hydrogen': 'kg', 'gas':'kJ', 'heating water':'kJ'}
    
    if var: 
        
        carriers = [var]
        units    = {var : units[var]}

        if var == 'hydrogen' and 'mechanical compressor' in balances[loc][var]:
            balances[loc][var].pop('mechanical compressor')  # dict.values() to be removed as they alter the total hydrogen balance and ar enot supposed to do so. 
            
    for carrier in carriers:
        balance  = 0   # initializing the variable to visualize the balance at the end of the simulation period
        print('\n\n'+loc+' '+carrier+' balance: '+simulation_name+'\n') 

        for b in balances[loc][carrier]:
            
            positiv=balances[loc][carrier][b][balances[loc][carrier][b]>0].sum()*c.P2E
            negativ=balances[loc][carrier][b][balances[loc][carrier][b]<0].sum()*c.P2E
            
            if positiv != 0:
                print(b+' '+str(round(positiv,1))+' '+units[carrier])
                balance += positiv
                
            if negativ != 0:
                print(b+' '+str(round(negativ,1))+' '+units[carrier])
                balance += negativ
            

def NPV_plot(study_case):
    ##### economic
    with open('results/pkl/economic_assessment_'+study_case+'.pkl', 'rb') as f: economic = pickle.load(f)
    
    plt.figure(dpi=1000)
    
    for loc in economic:    
        y = economic[loc]['NPV']
        x = np.linspace(0,len(y)-1,len(y))
        plt.plot(x,y,label=loc)
    plt.plot(x,np.zeros(len(x)),color='k')
        
    plt.legend()
    plt.title("Investments")
    plt.grid()
    plt.ylabel('Net Present Value [€]')
    plt.xlabel('Time [years]')
    plt.xlim(0,len(y)-1)
    plt.show()

def LOC_plot(simulation_name):
      
    with open('results/pkl/LOC_'+simulation_name+'.pkl', 'rb') as f:
        LOC = pickle.load(f)
           
    unit = {'H tank': '[kg]', 'HPH tank': '[kg]', 'battery': '[kJ]', 'inertial TES': '[°C]'}
    
    for location_name in LOC:
        for tech in LOC[location_name]:
            
            plt.figure(dpi=600)                               
            y = LOC[location_name][tech]
            x = np.linspace(0,len(y)-1,len(y))        
            plt.plot(x,y,label=location_name)
            plt.grid()
            plt.ylabel('LOC '+unit[tech])
            plt.xlabel('Time [hours]')
            xticks = list(np.linspace(0, len(x) - 1, 13).astype(int))
            xticklabels = ['          Jan','         Feb','          Mar','         Apr','         May','          Jun','         Jul','          Aug','           Sep','          Oct','          Nov','           Dec','']
            plt.xticks(xticks,xticklabels,rotation=45)
            plt.title(location_name+' '+tech)
            plt.xlim(0,x[-1])
            plt.show()

def satisfaction_story(sati):
    mode = {0: 'no demand', 
            1: 'demand satisfied by iTES without switching on the HP ',
            2: 'demand satisfied by switching on the HP',
            3: 'demand satisfied by switching on the HP and iTES overheated',
            -1: 'unsatisfied demand because there is not enough power',
            -2: 'unsatisfied demand because i TES cant reach the minimum required temperature in one step',
            -3: 'unsatisfied demand becasue t amb is too cold to heat water to the desired temperature'       
            }
    
    print('\n history of heating demand satisfaction:')
    for m in mode:
        print(f"{mode[m]}: {(sati == m).sum()} steps")
  
def cop(cop):
    plt.figure(dpi=1000)
    plt.plot(np.arange(len(cop)),cop)
    plt.ylabel('COP [-]')
    xticks = list(np.linspace(0, len(cop) - 1, 13).astype(int))
    xticklabels = ['          Jan','         Feb','          Mar','         Apr','         May','          Jun','         Jul','          Aug','           Sep','          Oct','          Nov','           Dec','']
    plt.xticks(xticks,xticklabels,rotation=45)
    plt.grid()
    plt.title('HP Coefficient of Performance')
    plt.xlim(0,len(cop))
    plt.show()
    
def heating_demand(demand):
    plt.figure(dpi=1000)
    plt.plot(np.arange(len(demand)),demand,color='red')
    xticks = list(np.linspace(0, len(demand) - 1, 13).astype(int))
    xticklabels = ['          Jan','         Feb','          Mar','         Apr','         May','          Jun','         Jul','          Aug','           Sep','          Oct','          Nov','           Dec','']
    plt.xticks(xticks,xticklabels,rotation=45)
    plt.ylabel('Demand [kW]')
    plt.grid()
    plt.xlim(0,len(demand))
    plt.title('Heating water demand')
    plt.show()