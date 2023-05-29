import numpy as np
import os
import sys 
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),os.path.pardir)))   # temorarily adding constants module path 
import constants as c

class O2_tank:    
    
    def __init__(self,parameters,simulation_hours):
        
        """
        Create a O2_tank object
    
        parameters : dictionary
            'max capacity': float [kg]
            'pressure': float [bar]
            'self discharge': ?
                      
        output : O2 tank object able to:
            supply or abrosrb oxygen .use(h,oxy)
            record the level of charge .LOC
            calculate its own volume (pressure) .volume(pressure)
        """
        self.pressure = parameters['pressure']          # H tank storage pressure
        self.LOC = np.zeros(simulation_hours+1)         # array H tank level of Charge 
        self.max_capacity = parameters['max capacity']  # H tank max capacity [kg]
        self.used_capacity = 0                          # H tank used capacity <= max_capacity [kg]      
        
    def use(self,h,oxy,constant_demand=False):
        """
        The O2 tank can supply or absorb oxygen
     
        h: int hour to be simulated
        oxy: oxygen requested (oxy<0) or provided (oxy>0) [kg]
      
        output : oxygen supplied or absorbed that hour [kg]
        """
        if self.max_capacity:
            if oxy >= 0:                                         # charge H tank
                
                charge = min(oxy,self.max_capacity-self.LOC[h])  # how much oxygen can H tank absorb?
                self.LOC[h+1] = self.LOC[h]+charge               # charge H tank
                
                if self.LOC[h+1] > self.used_capacity: # update used capacity
                    self.used_capacity = self.LOC[h+1]      
                
                return(-charge) # return oxygen absorbed
                
            else: # discharge H tank (this logic allows to back-calculate the LOC[0], it's useful for long term storage systems)
                
                if (self.used_capacity==self.max_capacity):  # the max_capacity has been reached, so LOC[h+1] can't become negative 
                       
                    discharge = min(-oxy,self.LOC[h]) # how much oxygen can H tank supply?
                    self.LOC[h+1] = self.LOC[h]-discharge # discharge H tank
            
                else: # the max_capacity has not yet been reached, so LOC[h+1] may become negative and then the past LOC may be shifted upwards  
                                                      
                    discharge = min(-oxy,self.LOC[h]+self.max_capacity-self.used_capacity) # how much oxygen can H tank supply?
                    self.LOC[h+1] = self.LOC[h]-discharge                                  # discharge H tank
                    if self.LOC[h+1] < 0:                                                  # if the level of charge has become negative
                        self.used_capacity += - self.LOC[h+1]                              # incrase the used capacity
                        self.LOC[:h+2] += - self.LOC[h+1]                                  # traslate the past LOC array
                        
                return(discharge) # return oxygen supplied
            
        else: 
            charge = oxy - constant_demand        # how much oxygen can H tank absorb?
            self.LOC[h+1] = self.LOC[h]+charge    # charge H tank
            if h == (len(self.LOC)-2):            # at the end of simulation. LOC array has self.simulation_hours + 1 values
                self.max_capacity   = max(self.LOC)+abs(min(self.LOC))  # [kg] max tank capacity
                self.shift          = abs(min(self.LOC))                # oxygen amount in storage at time 0
                self.LOC            = self.LOC + self.shift             # shifting the Level Of Charge curve to avoid negative minimum value (minimum is now at 0kg)
                                                                        # It is now possible to define how much H2 must be present in storage at the beginning of simulation. 
            return(charge)
    
    def sizing(self,htankmaxcapacity):
        """
        With this function the oxygen tank sizing is simplified as a direct consequence of hydrogne production and consumption
     
        htankmaxcapacity: hydrogen tank size defined in 'supply-led' mode
        
        output: oxygen tank size
        """
        constant = 9 # [-] units of oxygen produced per each unit of hydrogen
        self.max_capacity = htankmaxcapacity*constant
        
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