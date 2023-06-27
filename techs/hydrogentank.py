import numpy as np
import os
import sys 
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),os.path.pardir)))   # temorarily adding constants module path 
import constants as c
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt

class H_tank:    
    
    def __init__(self,parameters,simulation_hours):
        
        """
        Create a H_tank object
    
        parameters : dictionary
            'max capacity': float [kg]
            'pressure': float [bar]
            'self discharge': ?
                      
        output : H tank object able to:
            supply or abrosrb hydrogen .use(h,hyd)
            record the level of charge .LOC
            calculate its own volume (pressure) .volume(pressure)
        """
        
        self.cost = False # will be updated with tec_cost()

        self.pressure       = parameters['pressure']                # [bar] H tank storage pressure
        self.LOC            = np.zeros(simulation_hours+1)          # array H tank level of Charge 
        self.max_capacity   = parameters['max capacity']            # H tank max capacity [kg]
        self.used_capacity  = 0                                     # H tank used capacity <= max_capacity [kg]
        temperature         = 273.15 + 15                           # [K] temperature at which hydrogen is stored
        self.density        = PropsSI('D', 'P', self.pressure*100000, 'T', temperature, 'hydrogen')  # [kg/m^3] hydrogen density for selected density and temperature
        if self.max_capacity:
            self.tank_volume = round(self.max_capacity/self.density,2)   # [m^3] tank volume
        
    def use(self,h,hyd,constant_demand=False):
        """
        The H tank can supply or absorb hydrogen
     
        h: int hour to be simulated
        hyd: hydrogen requested (hyd<0) or provided (hyd>0) [kg]
      
        output : hydrogen supplied or absorbed that hour [kg]
        """
        if self.max_capacity:
            
            if hyd >= 0:                                         # charge H tank
                
                charge = min(hyd,self.max_capacity-self.LOC[h])  # how much hydrogen can H tank absorb? Minimum between h2 produced and available capacity in tank -> maximum capacity - Level od Charge (always a positive value, shifted above 0 every time it goes beyond)
                self.LOC[h+1] = self.LOC[h]+charge               # charge H tank
                
                
                # self.used_capacity -> parameter in which is stored the memory of the tank charging story.
                #                       It is not representative of the used capacity at time h.
                #                       It represents the maximum level reached inside the tank up to time h of the simulation.
                #                       Once self.used_capacity reaches the value of self.max.capacity it is no longer possible to allow 
                #                       the tank's LOC to move towards negative values depending on the hydrogen demand.
                #                       From this point onwards, the value of self.used_capacity remains the same until the end of the simulation. 
                                       
                if self.LOC[h+1] > self.used_capacity: # update used capacity
                    self.used_capacity = self.LOC[h+1]      
                return(-charge) # return hydrogen absorbed
                
            else: # discharge H tank (this logic allows to back-calculate LOC[0], it's useful for long term storage systems)
                
                if (self.used_capacity==self.max_capacity):  # the max_capacity has been reached, so LOC[h+1] can't become negative 
                       
                    discharge = min(-hyd,self.LOC[h]) # how much hydrogen can H tank supply?
                    self.LOC[h+1] = self.LOC[h]-discharge # discharge H tank
            
                else:   # the max_capacity has not yet been reached, so LOC[h+1] may become negative and then the past LOC can be shifted upwards to positive values. 
                        # The history of LOC during simulation is created by taking the self.used_capacity parameter into account. 
                        # Once the maximum capacity is reached for the first time, no more shifts to negative values are permitted and LOC
                        # remains the only parameter to represent the actual hydrogen amount in side the tank.
                                                      
                    discharge = min(-hyd,self.LOC[h]+self.max_capacity-self.used_capacity) # how much hydrogen can H tank supply?
                    self.LOC[h+1] = self.LOC[h]-discharge                                  # discharge H tank
                    if self.LOC[h+1] < 0:                                                  # if the level of charge has become negative
                        self.used_capacity += - self.LOC[h+1]                              # incrase the used capacity
                        self.LOC[:h+2] += - self.LOC[h+1]                                  # traslate the past LOC array
                return(discharge) # return hydrogen supplied
            
        else:   # Tjis option is activated when hydrogen tank is sized at the end of simulation as a result of 'supply-led' operation strategy.
                # A costant mass-flow rate demand is created based on the clumulative production of electrolyzers througout the year. 
                # Tank size in this case smoothes surplus or deficit of production during operation, allowing for a constant rate deliver. 
            charge = hyd - constant_demand        # how much hydrogen can H tank absorb?
            self.LOC[h+1] = self.LOC[h]+charge    # charge H tank
            if h == (len(self.LOC)-2):            # at the end of simulation. LOC array has self.simulation_hours + 1 values
                self.max_capacity   = max(self.LOC)+abs(min(self.LOC))  # [kg] max tank capacity
                self.shift          = abs(min(self.LOC))                # Hydrogen amount in storage at time 0
                self.LOC            = self.LOC + self.shift             # shifting the Level Of Charge curve to avoid negative minimum value (minimum is now at 0kg)
                                                                        # It is now possible to define how much H2 must be present in storage at the beginning of simulation. 
                self.tank_volume = round(self.max_capacity/self.density,2)   # [m^3] tank volume   
                                                     
            return(charge)
        
        
    def tech_cost(self,tech_cost):
        """
        Parameters
        ----------
        tech_cost : dict
            'cost per unit': float [€/kg]
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

        size = self.max_capacity
        # print('H2 Tank Max Capacity [tonnes]', self.max_capacity/1000)
        
        if tech_cost['cost per unit'] == 'default price correlation':
            C0 = 1500 # €/kg
            scale_factor = 0.4 # 0:1
            C = size * C0 **  scale_factor
        else:
            C = size * tech_cost['cost per unit']

        tech_cost['total cost'] = tech_cost.pop('cost per unit')
        tech_cost['total cost'] = C
        tech_cost['OeM'] = tech_cost['OeM'] *C /100 # €

        self.cost = tech_cost   
        
class HPH_tank:    
    
    def __init__(self,parameters,simulation_hours):
        
        """
        Create a H_tank object
    
        parameters : dictionary
            'max capacity': float [kg]
            'pressure': float [bar]
            'self discharge': ?
                      
        output : H tank object able to:
            supply or abrosrb hydrogen .use(h,hyd)
            record the level of charge .LOC
            calculate its own volume (pressure) .volume(pressure)
        """
        
        self.cost = False # will be updated with tec_cost()

        self.pressure = parameters['pressure']          # H tank storage pressure
        self.LOC = np.zeros(simulation_hours+1)         # array H tank level of Charge 
        self.max_capacity = parameters['max capacity']  # H tank max capacity [kg]
        self.used_capacity = 0                          # H tank used capacity <= max_capacity [kg]
        temperature         = 273.15 + 15                           # [K] temperature at which hydrogen is stored
        self.density        = PropsSI('D', 'P', self.pressure*100000, 'T', temperature, 'hydrogen')  # [kg/m^3] hydrogen density for selected density and temperature
        if self.max_capacity:
            self.tank_volume = round(self.max_capacity/self.density,2)   # [m^3] tank volume
        
    def use(self,h,hyd,constant_demand=False):
        """
        The H tank can supply or absorb hydrogen
     
        h: int hour to be simulated
        hyd: hydrogen requested (hyd<0) or provided (hyd>0) [kg]
      
        output : hydrogen supplied or absorbed that hour [kg]
        """
        if self.max_capacity:
            if hyd >= 0:                                         # charge H tank
                
                charge = min(hyd,self.max_capacity-self.LOC[h])  # how much hydrogen can H tank absorb?
                self.LOC[h+1] = self.LOC[h]+charge               # charge H tank
                
                if self.LOC[h+1] > self.used_capacity: # update used capacity
                    self.used_capacity = self.LOC[h+1]      
                
                return(-charge) # return hydrogen absorbed
                
            else: # discharge H tank (this logic allows to back-calculate the LOC[0], it's useful for long term storage systems)
                
                if (self.used_capacity==self.max_capacity):  # the max_capacity has been reached, so LOC[h+1] can't become negative 
                       
                    discharge = min(-hyd,self.LOC[h]) # how much hydrogen can H tank supply?
                    self.LOC[h+1] = self.LOC[h]-discharge # discharge H tank
            
                else: # the max_capacity has not yet been reached, so LOC[h+1] may become negative and then the past LOC may be shifted upwards  
                                                      
                    discharge = min(-hyd,self.LOC[h]+self.max_capacity-self.used_capacity) # how much hydrogen can H tank supply?
                    self.LOC[h+1] = self.LOC[h]-discharge                                  # discharge H tank
                    if self.LOC[h+1] < 0:                                                  # if the level of charge has become negative
                        self.used_capacity += - self.LOC[h+1]                              # incrase the used capacity
                        self.LOC[:h+2] += - self.LOC[h+1]                                  # traslate the past LOC array
                        
                return(discharge) # return hydrogen supplied
            
        else: 
            charge = hyd - constant_demand        # how much hydrogen can H tank absorb?
            self.LOC[h+1] = self.LOC[h]+charge    # charge H tank
            if h == (len(self.LOC)-2):            # at the end of simulation. LOC array has self.simulation_hours + 1 values
                self.max_capacity   = max(self.LOC)+abs(min(self.LOC))  # [kg] max tank capacity
                self.shift          = abs(min(self.LOC))                # Hydrogen amount in storage at time 0
                self.LOC            = self.LOC + self.shift             # shifting the Level Of Charge curve to avoid negative minimum value (minimum is now at 0kg)
                                                                        # It is now possible to define how much H2 must be present in storage at the beginning of simulation. 
                self.tank_volume = round(self.max_capacity/self.density,2)   # [m^3] tank volume
                
            return(charge)
        
        
    def tech_cost(self,tech_cost):
        """
        Parameters
        ----------
        tech_cost : dict
            'cost per unit': float [€/kg]
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

        size = self.max_capacity
        
        if tech_cost['cost per unit'] == 'default price correlation':
            C0 = 1500 # €/kg
            scale_factor = 0.4 # 0:1
            C = size * C0 **  scale_factor
        else:
            C = size * tech_cost['cost per unit']

        tech_cost['total cost'] = tech_cost.pop('cost per unit')
        tech_cost['total cost'] = C
        tech_cost['OeM'] = tech_cost['OeM'] *C /100 # €


        self.cost = tech_cost    