import numpy as np
import os
import sys 
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),os.path.pardir)))   # temorarily adding constants module path 
from core import constants as c
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt

class H_tank:    
    
    def __init__(self,parameters,timestep_number):
        
        """
        Create a H_tank object. Hydrogen storage system in tanks. 
    
        parameters : dictionary
            'max capacity': float tank capacity [kg]
            'pressure': float storage pressure [bar]
                      
        output : H tank object able to:
            supply or abrosrb hydrogen .use(step,hyd)
            keep track of the level of charge .LOC
        """
        
        self.cost = False # will be updated with tec_cost()
        self.timestep       = c.timestep                            # [min] selected timestep for simulation
        self.pressure       = parameters['pressure']                # [bar] H tank storage pressure
        self.LOC            = np.zeros(timestep_number+1)           # [kg] array keeping trak hydrogen tank level of charge 
        self.max_capacity   = parameters['max capacity']            # [kg] H tank max capacity 
        self.used_capacity  = 0                                     # [kg] H tank used capacity <= max_capacity 
        temperature         = 273.15 + 15                           # [K] temperature at which hydrogen is stored
        self.density        = PropsSI('D', 'P', self.pressure*100000, 'T', temperature, 'hydrogen')  # [kg/m^3] hydrogen density for selected density and temperature
        if self.max_capacity:
            self.tank_volume = round(self.max_capacity/self.density,2)   # [m^3] tank volume
        
    def use(self,step,hyd,constant_demand=False):
        """
        Hydrogen tank can supply or absorb hydrogen
     
        step: step to be simulated
        hyd: hydrogen requested (hyd<0) or provided (hyd>0) [kg/s]
      
        output : hydrogen supplied or absorbed in  the timestep [kg/s]
        """
        
        hyd = hyd*self.timestep*60      # [kg] conversion from kg/s to kg for the considered timestep
        
        if self.max_capacity:           # if tank size has been defined when setting up the case study
            if hyd >= 0:                # [kg] if hyd has positive value it means excess hydrogen is being produced by the system. Hydrogen tank is charged
                charge = min(hyd,self.max_capacity-self.LOC[step])  # [kg] how much hydrogen can H tank absorb? Minimum between hydrogen produced and available capacity in tank -> maximum capacity - Level of Charge (always a positive value, promptly shifted above 0 every time it goes beyond)
                self.LOC[step+1] = self.LOC[step]+charge            # [kg] charge H tank
                
                """
                EXPLANATION self.used_capacity parameter 
                self.used_capacity -> parameter in which is stored the memory of the tank charging story.
                                      It is not representative of the used capacity at the considered timestep.
                                      It represents the maximum level reached inside the tank up to time step of the simulation.
                                      Once self.used_capacity reaches the value of self.max.capacity it is no longer possible to allow 
                                      the tank's LOC to move towards negative values depending on the hydrogen demand.
                                      From this point onwards, the value of self.used_capacity remains the same until the end of the simulation.
                """
                                       
                if self.LOC[step+1] > self.used_capacity:           # update used capacity
                    self.used_capacity = self.LOC[step+1]           # [kg]
                    
                charge_flow_rate = charge/(self.timestep*60)        # [kg/s] converting the amount of hydrogen stored into a flow rate
                return(-charge_flow_rate)                           # [kg/s] return hydrogen absorbed 
                
            else: # if hyd value is negative. Discharge H tank (this logic allows to back-calculate LOC[0], it's useful for long term storage systems)
                if (self.used_capacity==self.max_capacity):         # the tank max_capacity has been reached in this step, so LOC[step+1] can't become negative anymore in case requested hydrogen is greater than available value
                    discharge           = min(-hyd,self.LOC[step])  # how much hydrogen can H tank supply?
                    self.LOC[step+1]    = self.LOC[step]-discharge  # discharge H tank
            
                else:   # the max_capacity has not yet been reached, so LOC[step+1] may become negative and then the past LOC can be shifted upwards to positive values. 
                        # The history of LOC during simulation is created by taking the self.used_capacity parameter into account. 
                        # Once the maximum capacity is reached for the first time, no more shifts to negative values are permitted and LOC
                        # remains the only parameter to represent the actual hydrogen amount in side the tank.
                                                      
                    discharge = min(-hyd,self.LOC[step]+self.max_capacity-self.used_capacity)   # [kg] hydrogen that can be supplied by H tank at the considered timestep
                    self.LOC[step+1] = self.LOC[step]-discharge                                 # [kg] discharge H tank
                    if self.LOC[step+1] < 0:                                                    # if the level of charge has become negative
                        self.used_capacity  += - self.LOC[step+1]                               # incrase the used capacity
                        self.LOC[:step+2]   += - self.LOC[step+1]                               # shift the past LOC array
                
                discharge_flow_rate = discharge/(self.timestep*60)                              # [kg/s] converting the amount of hydrogen to be dischrged into a flow rate
                return(discharge_flow_rate)                                                     # [kg/s] return hydrogen supplied 
            
        else:           # This option is activated when hydrogen tank is sized at the end of simulation as a result of 'supply-led' operation strategy.
                        # A costant mass-flow rate demand is created based on the clumulative production of electrolyzers througout the year. 
                        # Tank size in this case smoothes surplus or deficit of production during operation, allowing for a constant rate deliver. 
            
            constant_demand     = constant_demand*self.timestep*60              # [kg] transforming kg/s into kg for hydrogen  demand
            charge              = hyd - constant_demand                         # [kg] for the considered strategy hyd[step] is either positive or 0, representing the hydrogen produced by the elctrolyzer at each timestep. 
            self.LOC[step+1]    = self.LOC[step]+charge                         # [kg] charge H tank
            if step == (len(self.LOC)-2):                                       # at the end of simulation. LOC array has self.sim_timesteps+1 values
                self.max_capacity   = max(self.LOC)+abs(min(self.LOC))          # [kg] max tank capacity
                self.shift          = abs(min(self.LOC))                        # [kg] hydrogen amount in storage at time 0
                self.LOC            = self.LOC + self.shift                     # shifting the Level Of Charge curve to avoid negative minimum value (minimum is now at 0kg)
                                                                                # It is now possible to define how much hydrogen must be present in storage at the beginning of simulation. 
                self.tank_volume    = round(self.max_capacity/self.density,2)   # [m^3] tank volume   
             
            charge_flow_rate = charge/(self.timestep*60)                        # [kg/s] converting the amount of hydrogen stored into a flow rate                                     
            return(charge_flow_rate)                                            # [kg/s] return hydrogen absorbed 
        
        
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
    
    def __init__(self,parameters,timestep_number):
        
        """
        Create a H_tank object
    
        parameters : dictionary
            'max capacity': float [kg]
            'pressure': float [bar]
            'self discharge': ?
                      
        output : H tank object able to:
            supply or abrosrb hydrogen .use(step,hyd)
            record the level of charge .LOC
            calculate its own volume (pressure) .volume(pressure)
        """
        
        self.cost = False # will be updated with tec_cost()
        
        self.timestep = c.timestep
        
        self.pressure = parameters['pressure']          # H tank storage pressure
        self.LOC = np.zeros(timestep_number+1)         # array H tank level of Charge 
        self.max_capacity = parameters['max capacity']  # H tank max capacity [kg]
        self.used_capacity = 0                          # H tank used capacity <= max_capacity [kg]
        temperature         = 273.15 + 15                           # [K] temperature at which hydrogen is stored
        self.density        = PropsSI('D', 'P', self.pressure*100000, 'T', temperature, 'hydrogen')  # [kg/m^3] hydrogen density for selected density and temperature
        if self.max_capacity:
            self.tank_volume = round(self.max_capacity/self.density,2)   # [m^3] tank volume
        
    def use(self,step,hyd,constant_demand=False):
        """
        The H tank can supply or absorb hydrogen
     
        step: int hour to be simulated
        hyd: hydrogen requested (hyd<0) or provided (hyd>0) [kg/s]
      
        output : hydrogen supplied or absorbed that hour [kg/s]
        """
        
        hyd = hyd*self.timestep*60          # Conversion from kg/s to kg for the considered timestep
        
        if self.max_capacity:
            if hyd >= 0:                                         # charge H tank
                
                charge = min(hyd,self.max_capacity-self.LOC[step])  # how much hydrogen can H tank absorb?
                self.LOC[step+1] = self.LOC[step]+charge               # charge H tank
                
                if self.LOC[step+1] > self.used_capacity: # update used capacity
                    self.used_capacity = self.LOC[step+1]      
                
                charge_flow_rate = charge/(self.timestep*60)
                return(-charge_flow_rate) # return hydrogen absorbed [kg/s]
                
            else: # discharge H tank (this logic allows to back-calculate the LOC[0], it's useful for long term storage systems)
                
                if (self.used_capacity==self.max_capacity):  # the max_capacity has been reached, so LOC[step+1] can't become negative 
                       
                    discharge = min(-hyd,self.LOC[step]) # how much hydrogen can H tank supply?
                    self.LOC[step+1] = self.LOC[step]-discharge # discharge H tank
            
                else: # the max_capacity has not yet been reached, so LOC[step+1] may become negative and then the past LOC may be shifted upwards  
                                                      
                    discharge = min(-hyd,self.LOC[step]+self.max_capacity-self.used_capacity) # how much hydrogen can H tank supply?
                    self.LOC[step+1] = self.LOC[step]-discharge                                  # discharge H tank
                    if self.LOC[step+1] < 0:                                                  # if the level of charge has become negative
                        self.used_capacity += - self.LOC[step+1]                              # incrase the used capacity
                        self.LOC[:step+2] += - self.LOC[step+1]                                  # traslate the past LOC array
                
                discharge_flow_rate = discharge/(self.timestep*60)
                return(discharge_flow_rate) # return hydrogen supplied [kg/s]
            
        else:
            constant_demand = constant_demand*self.timestep*60       # Passing from kg/s to kg
            charge = hyd - constant_demand                           # how much hydrogen can H tank absorb?
            self.LOC[step+1] = self.LOC[step]+charge                 # charge H tank
            if step == (len(self.LOC)-2):                            # at the end of simulation. LOC array has self.simulation_hours + 1 values
                self.max_capacity   = max(self.LOC)+abs(min(self.LOC))  # [kg] max tank capacity
                self.shift          = abs(min(self.LOC))                # Hydrogen amount in storage at time 0
                self.LOC            = self.LOC + self.shift             # shifting the Level Of Charge curve to avoid negative minimum value (minimum is now at 0kg)
                                                                        # It is now possible to define how much H2 must be present in storage at the beginning of simulation. 
                self.tank_volume = round(self.max_capacity/self.density,2)   # [m^3] tank volume
            
            charge_flow_rate = charge/(self.timestep*60)  
            return(charge_flow_rate)
        
        
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