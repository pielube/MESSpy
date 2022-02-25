"""
Post processing
"""
import pickle
import numpy as np
import matplotlib.pyplot as plt

def total_balances(simulation_name):

    with open('Results/balances_'+simulation_name+'.pkl', 'rb') as f:
        balances = pickle.load(f)
    
    ###### load analysys
    
    ##### total electricity balances
    carriers = ['electricity','hydrogen']
    units = {'electricity': 'kWh', 'hydrogen': 'kg'}
    
    for carrier in carriers:
        print('\nTotal '+carrier+' balances:\n')
        for loc in balances:
            print('\n'+loc)   
            for b in balances[loc][carrier]:
                print(b+' '+str(round(sum(balances[loc][carrier][b]),1))+' '+units[carrier])
    
def NPV_plot():
    ##### economic
    simulation_name = 'economic_assessment'
    with open('Results/'+simulation_name+'.pkl', 'rb') as f:
        economic = pickle.load(f)
    
    plt.figure(dpi=1000)
    
    for loc in economic:    
        y = economic[loc]['NPV']
        x = np.linspace(0,len(y)-1,len(y))
        plt.plot(x,y,label=loc)
        
    plt.legend()
    plt.grid()
    plt.ylabel('Net Present Value [€]')
    plt.xlabel('Time [years]')
    plt.show()
    
def SOC_plot(simulation_name):
      
    with open('Results/SOC_'+simulation_name+'.pkl', 'rb') as f:
        SOC = pickle.load(f)
           
    unit = {'H tank': '[m3]', 'battery': '[kWh]'}
    
    for loc in SOC:
        for tech in SOC[loc]:
            
            y = SOC[loc][tech]
            x = np.linspace(0,len(y)-1,len(y))        
            plt.plot(x,y,label=loc)
            plt.grid()
            plt.ylabel('SOC '+unit[tech])
            plt.xlabel('Time [hours]')
            plt.title(loc+' '+tech)
            plt.show()
        

           
    
    
        