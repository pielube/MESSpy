# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 17:11:30 2022

@author: pietro
"""

import numpy as np
import pvlib

class heatpump:    
    
    def __init__(self,parameters,general):
        """
        Create a heatpump object
        Priority given to heating, can do both heating and cooling
        
        inputs:
            parameters : dictionary
                'Ppeak_heat': float peak heating power [kWp]
                'Ppeak_cool': float peak cooling power [kWp] 
                'COP_heat': float coefficient of performance [-]
                'COP_cool': float coefficient of performance [-]
                'Tsupply_heat': float temp supply heat [°C]
                'Tsupply_cool': float temp supply cool [°C]
   
        outputs : heatpump object able to:
            consume electricity to heat or cool .use()
        """
        
        if 'Ppeak_heat' in parameters:
            self.Ppeak_heat = parameters['Ppeak_heat']
        else:
            self.Ppeak_heat = 0.
                       
        if 'Ppeak_cool' in parameters:
            self.Ppeak_cool = parameters['Ppeak_cool']
        else:
            self.Ppeak_cool = 0.
            
        if 'COP_heat' in parameters:
            self.COP_heat = parameters['COP_heat']
        else:
            self.COP_heat = 0.
        
        if 'COP_cool' in parameters:
            self.COP_cool = parameters['COP_cool']
        else:
            self.COP_cool = 0.
            
        if 'Tsupply_heat' in parameters:
            self.Tsupply_heat = parameters['Tsupply_heat']
        else:
            self.Tsupply_heat = 0.
        
        if 'Tsupply_cool' in parameters:
            self.Tsupply_cool = parameters['Tsupply_cool']
        else:
            self.Tsupply_cool = 0.
            
        
        latitude = general['latitude']
        longitude = general['longitude']
        tmybool = general['TMY']
        
        if tmybool:
            weather = pvlib.iotools.get_pvgis_tmy(latitude, longitude, map_variables=True)[0]
        else:
            print('WARNING: heat pump requires Tamb')
            
        self.Tamb = weather['temp_air'].to_numpy()

        
        
    def use(self,dem_heat,dem_cool,h,timestep):
        """
        Compute electricity consumption and heating or cooling produced
        
        inputs :
            dem_heat float heating demand in timestep [kWh] (neg when there is demand to be covered)
            dem_cool float cooling demand in timestep [kWh] (neg when there is demand to be covered)
            h index which timestep we are in [-]
            timestep timestep [hours]
            
        outputs : 
            consumption float electricity consumption [kWh]
            heatprod float heat produced [kWh]
            coolprod float cool produced [kWh]
        """
        
        dem_heat = - dem_heat
        dem_cool = - dem_cool
        
        consumption = 0
        heatprod = 0
        coolprod = 0
               
        COP_ref = 3.91592 # see below
        
        if dem_heat > 0: # heating
            heatprod = min(dem_heat,self.Ppeak_heat*timestep)
            deltaT = max(0, self.Tsupply_heat - self.Tamb[h]) # determine temp difference, if neg set to 0
            COP = (6.81 - 0.121 * deltaT + 0.000630 * deltaT**2)/COP_ref*self.COP_heat # Eq (4) in Staffell et al.
            consumption = - heatprod / COP
        
        elif dem_heat == 0 and dem_cool > 0: # cooling
            coolprod = min(dem_cool,self.Ppeak_cool*timestep)
            deltaT = max(0, self.Tamb[h] - self.Tsupply_cool) # determine temp difference, if neg set to 0
            COP = (6.81 - 0.121 * deltaT + 0.000630 * deltaT**2)/COP_ref*self.COP_cool # Eq (4) in Staffell et al.
            consumption = - coolprod / COP
        
        return(consumption,heatprod,coolprod)


###########################################################################################################################################################

if __name__ == "__main__":
    
    """
    Test
    """

    inp_test = {'Ppeak_heat': 10., 'Ppeak_cool': 10.,'COP_heat': 4., 'COP_cool': 3.5, 'Tsupply_heat': 35.,'Tsupply_cool': 10.}
    
    general = {
    "simulation years": 30,
    "latitude": 43.41,
    "longitude": 11.15,
    "reference year": 2015,
    "TMY": True
    }
  
    hp = heatpump(inp_test,general)
    
    timestep = 1 # h
    modes = np.array([0,1,1,2,2,0,0]) # array 0 - off 1 - heat 2 - cool
    Nsteps = len(modes)

    Tamb = np.zeros(Nsteps)
    Tamb[np.where(modes==1)[0]] = 5.
    Tamb[np.where(modes==2)[0]] = 30.
    
    demheat = - np.ones(Nsteps)*5.
    demcool = - np.ones(Nsteps)*10.
    
    heatprod = np.zeros(Nsteps)
    coolprod = np.zeros(Nsteps)
    consump  = np.zeros(Nsteps)
    
    for i in range(Nsteps):
        
        consump[i],heatprod[i],coolprod[i] = hp.use(demheat[i],demcool[i],3,timestep)
        
        
        
"""
Air-Water heat pump. Ambient temperature as reservoir temperature.
COP based off regression analysis of manufacturers data
Source: "A review of domestic heat pumps, Iain Staffell, Dan Brett, Nigel Brandonc and Adam Hawkes"
http://pubs.rsc.org/en/content/articlepdf/2012/ee/c2ee22653g
TODO: Validate this methodology

1) hp.
   15 < deltaT < 60
   only heating

   Here we are using it for all deltaT and for heating and cooling

2) It is common for heat pumps in Europe to be
   tested at 7°C external temperature (dry bulb) and 20°C indoor
   temperature, with 35 °C output from the pump
   
   7°C and 35°C temperatures used to calculate ref avg COP to adim the curve
   and use it for different ref COPs
   deltaTref = 35. -7.
   COPref = 6.81 - 0.121 * deltaTref + 0.000630 * deltaTref**2 = 3.91592

3) Ref values for Theatsupply
   direct air heating: 25–35°C
   underfloor heating: 30–45°C
   large-area radiators: 45–60°C
   conventional radiators: 60–75°C
"""


    