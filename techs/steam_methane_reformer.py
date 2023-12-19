# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 17:27:00 2023

@author: Andrea
"""

import numpy as np
import os
import sys 
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),os.path.pardir)))   # temporarily adding constants module path 
from core import constants as c
import matplotlib.pyplot as plt

class SMR:    
    
    def __init__(self,parameters,timestep=False):
        """
        Create a SMR object 
    
        parameters : dictionary
            'Ppeak': float peak hydrogen output thermal power [kWp] 
            'efficiency': float overall SMR plant efficiency [-]

        outputs : SMR object able to:
            consume Natural gas and produce hydrogen.use(step,hyd)
        """
        
        self.Ppeak = parameters['Ppeak']
        self.efficiency = parameters['efficiency']              # Efficiency value taken from https://www.sciencedirect.com/science/article/pii/S0360319914014372
        self.cost = False # will be updated with tech_cost()    # Cost of 280 €/kW of output hydrogen thermal power taken from https://www.sciencedirect.com/science/article/pii/S2666790822001574

        try:
            self.timestep   = c.timestep            # [min]       simulation timestep if launched from main
        except AttributeError:
            self.timestep   = timestep              # [min]       simulation timestep if launched from Steam methane reformer.py
        
    def use(self,hyd):
        """
        Compute thermal power consumption (natural gas consumption)
        
        inputs :
            hyd float hydrogen demand in step [kg/s]
            
        outputs : 
            thermal power consumption float [kW]
        """
        
        hyd_kW = - hyd*c.LHVH2*1e3    # kW hydrogen thermal power output needed
        if hyd_kW > self.Ppeak:
            print(f"Warning: the Steam methane reformer peak power is too low to cover an hydrogen demand of {-hyd.round(2)} kg/s equivalent to {hyd_kW.round(0)} kW")
        NG_consumption_kW = min(hyd_kW/self.efficiency,self.Ppeak/self.efficiency)
        hyd_produced =  (min(hyd_kW,self.Ppeak))/(c.LHVH2*1e3)
        
        return( -NG_consumption_kW, hyd_produced)
    
    def tech_cost(self,tech_cost):
        """
        Parameters
        ----------
        tech_cost : dict
            'cost per unit': float [€/kW]
            'OeM': float, percentage on initial investment [%]
            'refud': dict
                'rate': float, percentage of initial investment which will be rimbursed [%]
                'years': int, years for reimbursment
            'replacement': dict
                'rate': float, replacement cost as a percentage of the initial investment [%]
                'years': int, after how many years it will be replaced

        Returns
        -------
        self.cost: dict
            'total cost': float [€]
            'OeM': float, percentage on initial investment [%]
            'refud': dict
                'rate': float, percentage of initial investment which will be rimbursed [%]
                'years': int, years for reimbursment
            'replacement': dict
                'rate': float, replacement cost as a percentage of the initial investment [%]
                'years': int, after how many years it will be replaced
        """
        
        tech_cost = {key: value for key, value in tech_cost.items()}        
        
        size = self.Ppeak # kW
        
        if tech_cost['cost per unit'] == 'default price correlation':
            C0 = 160 # €/kW
            scale_factor = 0.8 # 0:1
            C = C0 * size ** scale_factor
        else:
            C = size * tech_cost['cost per unit']

        tech_cost['total cost'] = tech_cost.pop('cost per unit')
        tech_cost['total cost'] = C
        tech_cost['OeM'] = tech_cost['OeM'] *C /100 # €


        self.cost = tech_cost
        
###########################################################################################################################################################

if __name__ == "__main__":
    
    """
    Test
    """
    
    inp_test_SMR = {'Ppeak': 10000, 'efficiency': 0.687}   # Efficiency value taken from https://www.sciencedirect.com/science/article/pii/S0360319914014372
                                                            # Cost of 280 €/kW of output hydrogen thermal power taken from https://www.sciencedirect.com/science/article/pii/S2666790822001574

    SMR = SMR(inp_test_SMR)
    
    hyd_dem = -np.arange(0.02,0.15,0.01) # [kg/s]
    
    Nsteps = len(hyd_dem)
    NG_consumption   = np.zeros(Nsteps)
    hyd_produced   = np.zeros(Nsteps)
    
    for i in range(len(hyd_dem)):
        
        NG_consumption[i],hyd_produced[i] = SMR.use(hyd_dem[i])
        
    print('/n total SMR NG consumption = {} kWh '.format(np.sum(-NG_consumption)))
    print('total SMR H2 production = {} kg \n'.format(np.sum(hyd_produced)))
    
    plt.figure(dpi=200)
    plt.plot(-NG_consumption,hyd_produced)
    plt.xlabel('NG consumption [kW]')
    plt.ylabel('H2 produced [kg/s]')
    plt.show()
    
    from collections import Counter
        # Count occurrences of each unique point
    point_counts = Counter(zip(NG_consumption, hyd_produced))
    
    # Create a scatter plot with transparency
    plt.figure(dpi=200)
    plt.scatter(-NG_consumption, hyd_produced, alpha=0.5, label='Data points')
    plt.xlabel('NG consumption [kW]')
    plt.ylabel('H2 produced [kg/s]')
    
    # Add count labels for each point
    for (x, y), count in point_counts.items():
        plt.text(-x, y, str(count), fontsize=8, ha='right', va='bottom', color='red')
    
    plt.legend()
    plt.show()

    