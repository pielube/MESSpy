# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 15:46:34 2022

@author: pietro
"""

import numpy as np
import os
import sys 
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),os.path.pardir)))   # temorarily adding constants module path 

from core import constants as c

class boiler: 
    def __init__(self, parameters,timestep=False):
        """
        Create a general boiler object, serving as the parent class for different specific boiler models.
    
        parameters : dictionary
            'Ppeak': float peak thermal power [kWp] 
            'efficiency': float boiler efficiency [-]
            

        outputs : boiler object able to:
            consume fuel or electricity and produce heat .use(demand,timestep)
        """
        
        self.Ppeak          = parameters['Ppeak']
        self.efficiency     = parameters['efficiency']
        self.cost           = False  # to be updated via tec_cost() function
        if timestep == False: 
            self.timestep   = c.timestep              # [min]       simulation timestep if launched from main
        else:
            self.timestep   = timestep                # [min]       simulation timestep if launched from boiler.py
        
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
            C0              = 160 # €/kW
            scale_factor    = 0.8 # 0:1
            C               = C0 * size ** scale_factor
        else:
            C               = size * tech_cost['cost per unit']

        tech_cost['total cost'] = tech_cost.pop('cost per unit')
        tech_cost['total cost'] = C
        tech_cost['OeM']        = tech_cost['OeM'] *C /100 # €


        self.cost = tech_cost
    
    
class boiler_el(boiler):    
    
    def __init__(self,parameters):
        """
        Create an electric boiler object 
    
        parameters : dictionary
            'Ppeak': float peak thermal power [kWp] 
            'efficiency': float boiler efficiency [-]
            
        outputs : boiler object able to:
            consume electricity and produce heat .use(timestep,demand)
        """
        super().__init__(parameters)
        
    def use(self,step,demand):
        """
        Computes electricity consumption and heat produced
        
        inputs :
            considered timestep []
            float energy demand in the considered timestep [kW]
            
        outputs : 
            consumption : float energy consumption [kW]
            heatprod    : float heat produced [kW] 
        """

        if -demand/self.efficiency > self.Ppeak:
            print(f'Warning: the boiler nominal power is too low to cover heat demand - simulation step: {step}')
        heatprod    = min(-demand,self.Ppeak*self.efficiency)
        consumption = - heatprod/self.efficiency 
        
        return(consumption,heatprod)
       
     
class boiler_ng(boiler):    
    
    def __init__(self,parameters):
        """
        Create a natual gas fuelled boiler object 
    
        parameters : dictionary
            'Ppeak': float peak thermal power [kWp] 
            'efficiency': float boiler efficiency [-]

        outputs : boiler object able to:
            consume natural gas and produce heat .use(timestep,demand)
        """
        super().__init__(parameters)
        self.LHVNG      = c.LHVNG                   # [MJ/kg]      Natural Gas Lower Heating Value
        self.LHVNGVOL   = c.LHVNGVOL*1000           # [kJ/Sm^3]    Natural Gas Lower Heating Value - Volumetric
        
    def use(self,step,demand):
        """
        Compute natural gas consumption and heat produced
        
        inputs :
            considered timestep []
            demand float energy demand in timestep [kW] (-)
            
        outputs : 
            consumption : float natural gas consumption [Sm3/s]
            heatprod    : float heat produced [kW] 
        """

        if demand < 0:  # [kW] heat required
            if -demand/self.efficiency > self.Ppeak:
                print(f'Warning: the boiler nominal power is too low to cover heat demand - simulation step: {step}')            
            heatprod        = min(-demand,self.Ppeak*self.efficiency)       # [kW] produced heat at the considered timestep
            ng_mflowrate    = (- heatprod/self.efficiency)/(self.LHVNGVOL)  # [Sm^3/s] natural gas consumption to produce the required heat  
            
            return(ng_mflowrate,heatprod)   

              
class boiler_h2(boiler):    

    def __init__(self,parameters):
        """
        Create an hydrogen-fuelled boiler object  
    
        parameters : dictionary
            'Ppeak'     : float peak thermal power [kWp] 
            'efficiency': float boiler efficiency [-]
            
        outputs : boiler object able to:
            consume hydrogen and produce heat .use(timestep,demand)
        """
        super().__init__(parameters)
        self.LHVH2      = c.LHVH2*1000          # [kJ/kg]      Hydrogen Lower Heating Value

    def use(self,step,demand,available_hyd):
        """
        Compute consumption and heat produced
        
        inputs :
            considered timestep []
            demand float energy demand in timestep [kW] (-)
            
        outputs : 
            consumption : float hydrogen consumption [kg/s]
            heatprod    : float heat produced [kW]    
        """

        max_available_hyd = available_hyd/(self.timestep*60) # [kg/s] available hydrogne mass flow 
        
        if demand < 0: # heat required [kW]
            if -demand/self.efficiency > self.Ppeak:
                print(f'Warning: the boiler nominal power is too low to cover heat demand - simulation step: {step}')
            heatprod        = min(-demand,self.Ppeak*self.efficiency)   # [kW] produced heat at the considered timestep
            input_heat      = - heatprod/self.efficiency                # [kW] boiler input gross heat needed to satisfy demand
            h2_mflowrate    = input_heat/self.LHVH2                     # [kg/s] hydrogen consumption to produce the required heat
        
            if -h2_mflowrate > max_available_hyd:   # partial load operation, if not enough hydrogen is available to meet demand (H tank is nearly empty)
                heatprod        = (max_available_hyd*self.LHVH2)*self.efficiency    # [kW] produced heat at the considered timestep
                h2_mflowrate    = max_available_hyd                                 # [kg/s] hydrogen consumption to produce the required heat
            
            return(h2_mflowrate,heatprod)
        
###########################################################################################################################################################

if __name__ == "__main__":
    
    """
    Test
    """
    
    inp_test_NG_noncond = {'Ppeak': 24., 'efficiency': 0.85}
    inp_test_el         = {'Ppeak': 24., 'efficiency': 0.92}  # P: 3-6 kW @ 230 V, 9, 12, 14 kW or >36 kW @ 400 V (triphase)
    inp_test_NG_cond    = {'Ppeak': 12., 'efficiency': 1.00} 
    inp_test_h2         = {'Ppeak': 24., 'efficiency': 0.92}
  
    # boiler_NG_noncond = boiler_ng(inp_test_NG_noncond)
    boiler_NG_cond    = boiler_ng(inp_test_NG_cond)
    # boiler_el         = boiler_el(inp_test_el)
    # boiler_h2         = boiler_h2(inp_test_h2)
    
    # timestep = 60  # [min]
    # Nsteps = 10
    
    # demand = -np.arange(0,6,1) # [kWh]
    # available_hyd=np.linspace(10,2,6)
    
    # Qprod_NG   = np.zeros(Nsteps)
    # Qcons_NG   = np.zeros(Nsteps)
    # Qprod_el   = np.zeros(Nsteps)
    # Qcons_el   = np.zeros(Nsteps)
    # Qprod_H2   = np.zeros(Nsteps)
    # Qcons_H2   = np.zeros(Nsteps)
    # Qconskg_H2   = np.zeros(Nsteps)
    
    # for i in range(len(demand)):
        
    #     Qcons_NG[i],Qprod_NG[i] = boiler_NG_noncond.use(demand[i],timestep)
    #     Qcons_el[i],Qprod_el[i] = boiler_el.use(demand[i],timestep)
    #     Qcons_H2[i],Qconskg_H2[i],Qprod_H2[i] = boiler_h2.use(demand[i],available_hyd[i],timestep)
        
        
    # print('Boiler NG generation = {} kWh '.format(np.sum(Qprod_NG)))
    # print('// // // consumption = {} kWh \n'.format(np.sum(-Qcons_NG)))
    # print('Boiler Ele generation = {} kWh '.format(np.sum(Qprod_el)))
    # print('// // // consumption = {} kWh \n'.format(np.sum(-Qcons_el)))
    # print('Boiler H2 generation = {} kWh '.format(np.sum(Qprod_H2)))
    # print('// // // consumption = {} kWh '.format(np.sum(-Qcons_H2)))
    # print('// // // consumption = {} kg '.format(np.sum(-Qconskg_H2)))

    