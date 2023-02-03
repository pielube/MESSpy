import numpy as np
import os
import sys 
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),os.path.pardir)))   # temorarily adding constants module path 
import constants as c

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
        
        self.pressure = parameters['pressure']          # H tank storage pressure
        self.LOC = np.zeros(simulation_hours+1)         # array H tank level of Charge 
        self.max_capacity = parameters['max capacity']  # H tank max capacity [kg]
        self.used_capacity = 0                          # H tank used capacity <= max_capacity [kg]      
        
    def use(self,h,hyd):
        """
        The H tank can supply or absorb hydrogen
     
        h: int hour to be simulated
        hyd: hydrogen requested (hyd<0) or provided (hyd>0) [kg]
      
        output : hydrogen supplied or absorbed that hour [kg]
        """
        
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