# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 15:46:34 2022

@author: pietro
"""

import numpy as np
import sys 
sys.path.append('..')
import constants as c

class boiler_el:    
    
    def __init__(self,parameters):
        """
        Create a boiler object 
    
        parameters : dictionary
            'Ppeak': float peak thermal power [kWp] 
            'efficiency': float boiler efficiency [-]
            

        outputs : boiler object able to:
            consume fuel or electricity and produce heat .use(demand,timestep)
        """
        
        self.Ppeak = parameters['Ppeak']
        self.efficiency = parameters['efficiency']
        
        
    def use(self,demand,timestep):
        """
        Compute consumption and heat produced
        
        inputs :
            demand float energy demand in timestep [kWh]
            timestep float timestep in hours [h]
            
        outputs : 
            consumption float energy consumption [kWh]
            heatprod float heat produced [kWh] 
        """
        
        heatprod = min(demand,self.Ppeak*timestep)
        consumption = - heatprod/self.efficiency #
        
        return(consumption,heatprod)
    
class boiler_ng:    
    
    def __init__(self,parameters):
        """
        Create a boiler object 
    
        parameters : dictionary
            'Ppeak': float peak thermal power [kWp] 
            'efficiency': float boiler efficiency [-]
            

        outputs : boiler object able to:
            consume fuel or electricity and produce heat .use(demand,timestep)
        """
        
        self.Ppeak = parameters['Ppeak']
        self.efficiency = parameters['efficiency']
        
        
    def use(self,demand,timestep):
        """
        Compute consumption and heat produced
        
        inputs :
            demand float energy demand in timestep [kWh] (-)
            timestep float timestep in hours [h] 
            
        outputs : 
            consumption float energy consumption [kWh]
            heatprod float heat produced [kWh] 
        """
        
        if demand < 0: # heat required
        
            heatprod = min(-demand,self.Ppeak*timestep)
            consumption = - heatprod/self.efficiency 
            
            return(consumption,heatprod)
        
        else:
            return(0,0)       
        
###########################################################################################################################################################

if __name__ == "__main__":
    
    """
    Test
    """
    
    inp_test_NG_noncond = {'Ppeak': 24., 'efficiency': 0.85}
    inp_test_NG_cond    = {'Ppeak': 24., 'efficiency': 0.92}
    inp_test_el         = {'Ppeak': 12., 'efficiency': 1.00} # P: 3-6 kW @ 230 V, 9, 12, 14 kW or >36 kW @ 400 V (triphase)
  
    boiler_NG_noncond = boiler_ng(inp_test_NG_noncond)
    boiler_NG_cond    = boiler_ng(inp_test_NG_cond)
    boiler_el         = boiler_el(inp_test_el)
    
    timestep = 1 # h
    Nsteps = 10
    
    demands = np.arange(0,6,1) # kWh
    
    Qth   = np.zeros(Nsteps)
    Qfuel = np.zeros(Nsteps)
    
    for i in range(len(demands)):
        
        Qth[i],Qfuel[i] = boiler_NG_noncond.use(demands[i],timestep)
        # Qth[i],Qfuel[i] = boiler_NG_cond.use(demands[i],timestep)
        # Qth[i],Qfuel[i] = boiler_el.use(demands[i],timestep)
        
    print(np.sum(Qth))
    print(np.sum(Qfuel))

    