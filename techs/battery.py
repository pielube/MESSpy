import numpy as np
import os
import sys 
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),os.path.pardir)))   # temorarily adding constants module path 
import constants as c

class battery:    
    
    def __init__(self,parameters,simulation_hours):
        """
        Create a battery object
    
        parameters : dictionary
            'nominal capacity': float [kWh]
            'max E-rate': float charging and discharging rate related to maximum capacity [kW/kWh]
            'efficiency': float charge and discharge efficiency
            
            'ageing': bool true if aging has to be calculated
            'life cycles': int number of life cycles to reach the end of battery life
            'end life capacity': float maximum capacity left at end of life [%]
            
            'collective': int 0: no collective rules. 1: priority to csc and then charge or discharge the battery.
        
        output : battery object able to:
            supply or abrosrb electricity .use(h,e)
            record the level of charge .LOC
            take account of ageing .calculate_aging()   
        """
        
        self.nom_capacity = parameters['nominal capacity'] # battery early life max capacity [kWh]
        self.max_capacity = parameters['nominal capacity'] # battery max capacity [kWh]
        self.max_E_rate = parameters['max E-rate'] # battery max_E_rate [kW/kWh]
        self.eta = parameters['efficiency'] # float charge and discharge efficiency
        
        self.ageing = parameters['ageing'] # bool true if aging has to be calculated
        self.LC = parameters['life cycles'] # int number of life cycles to reach the end of battery life
        self.deg = self.max_capacity*(100-parameters['end life capacity'])/100 # float capacity that is going to be degradated in life cycles
        self.ageing_day = 7 # How often ageing has to bee calculated? [days]
        self.completed_cycles = 0 # float initialise completed_cycles, this parameter is usefull to calculate replacements
        self.replacements = [] # list initialise: h at which replecaments occur
        self.ageing_history = [[0],[self.max_capacity]] # list initialise ageing history. [completed_cycles],[max_capacity]
    
        self.LOC = np.zeros(simulation_hours+1) # array battery level of Charge 
        self.used_capacity = 0 # battery used capacity <= max_capacity [kWh]
      
        self.collective = parameters['collective'] # int 0: no collective rules. 1: priority to csc and then charge or discharge the battery.

    def use(self,h,e):
        """
        The battery can supply or absorb electricity
     
        h: int hour to be simulated
        e: electricity requested (e<0) or provided (e>0) [kWh]
      
        output : electricity supplied or absorbed that hour [kWh]
        """
        
        if self.ageing and (h/24%self.ageing_day == 0) and h!=0: # if aeging == True and it's time to calculate it
            self.calculate_ageing(h)
            
        if e >= 0: # charge battery
        
            if e > self.max_capacity*self.max_E_rate:
                e = self.max_capacity*self.max_E_rate
                              
            if e*self.eta < (self.max_capacity-self.LOC[h]): # if battery can't be full charged
                self.LOC[h+1] = self.LOC[h]+e*self.eta # charge 
            
            else: # if batter can be full charged
                self.LOC[h+1] = self.max_capacity
                e = (self.max_capacity-self.LOC[h]) / self.eta
            
            if self.LOC[h+1] > self.used_capacity: # update used capacity
                self.used_capacity = self.LOC[h+1] 
                   
            return(-e) # return electricity absorbed
        
            
        else: # discharge battery (this logic allows to back-calculate the LOC[0], it's useful for long term storage systems)
            
            e = e/self.eta # how much energy is really required
            if(self.used_capacity==self.nom_capacity):  # the nom_capacity has been reached, so LOC[h+1] can't become negative 
                   
                discharge = min(-e,self.LOC[h],self.max_capacity*self.max_E_rate) # how much electricity can battery supply?
                
                self.LOC[h+1] = self.LOC[h]-discharge # discharge battery
                
            else: # the max_capacity has not yet been reached, so LOC[h+1] may become negative and then the past LOC may be translated   
                                                  
                discharge = min(-e,self.LOC[h]+self.max_capacity-self.used_capacity,self.max_capacity*self.max_E_rate) # how much electricity can battery supply?
                self.LOC[h+1] = self.LOC[h]-discharge # discharge battery
                if self.LOC[h+1] < 0: # if the level of charge has become negative
                    self.used_capacity += - self.LOC[h+1] # incrase the used capacity
                    self.LOC[:h+2] += - self.LOC[h+1]  # traslate the past LOC array
            
            return(discharge*self.eta) # return electricity supplied
        
    def calculate_ageing(self,h):      
        
        # degradation (equivalent number of cycles, life cycles, end of life capacity)
        cycles = self.rainflow(h) # number of equivalent cycles completed
        self.max_capacity += - (cycles / self.LC) * self.deg      
        self.completed_cycles += cycles
        self.ageing_history[0].append(self.completed_cycles)
        self.ageing_history[1].append(self.max_capacity)

        
        if self.completed_cycles > self.LC: # replacement
            self.replacements.append(h)
            self.completed_cycles = 0
            self.max_capacity = self.nom_capacity
        
    def rainflow(self,h):
        
        #https://ieeexplore.ieee.org/document/7741532
        
        LOC = self.LOC[h-self.ageing_day*24:h] # part of the LOC whose contribution to aging is to be calculated
        
        # elimination of the plains 
        new=[LOC[0]] # initialise new LOC withouth plains
        for i in range(1,len(LOC)):
            if LOC[i]!=LOC[i-1]:
                new.append(LOC[i])
        LOC=new # new LOC withouth plains
        
        # elimination of what is not a pick or a valley
        def PoV(a,b,c):
            r=0
            if b>=a and b>=c: #peak
                r=1
            if b<=a and b<=c: #valley
                r=1
            return(r)     
        
        new=[LOC[0]] # initialise new LOC without pick or valley         
        for i in range(1,len(LOC)-1):
            r=PoV(LOC[i-1],LOC[i],LOC[i+1])
            if r==1:
                new.append(LOC[i])
        new.append(LOC[len(LOC)-1])
        LOC=new # new LOC without pick or valley
            
        # find half and full cicles and calculate their depth
        hc=[] #depth of half cycles
        fc=[] #depth of full cycles
        stop=0
        while(len(LOC)>3):
            if stop==1:
                break
                  
            for i in range(len(LOC)-2):    
                Rx=abs(LOC[i]-LOC[i+1])
                Ry=abs(LOC[i+1]-LOC[i+2])
                if Rx<=Ry: #half cycle finded
                    if i==0:
                        hc.append(Rx)
                        LOC.pop(i)
                        break
                if Rx<=Ry: #full cycle finded
                    fc.append(Rx)
                    LOC.pop(i)
                    LOC.pop(i)
                    break
                # if only half cicles remain and len(LOC) is still >3 we need break
                if i==len(LOC)-3:
                    stop=1
                    break
                
        # final cicles control
        for i in range(len(LOC)-1):
            hc.append(abs(LOC[i]-LOC[i+1]))
            
        # calculate the equivalent number of cycles
        n_cycles = 0        
        for c in hc:
            # depth of the cycle / depth of a complete cycle / 2 because it's an half cycle
            n_cycles += 0.5 * c / self.max_capacity            
        for c in fc:
            # depth of the cycle / depth of a complete cycle
            n_cycles += c / self.max_capacity 
        
        return(n_cycles)
    
    def tech_cost(self,tech_cost):
        """
        Parameters
        ----------
        tech_cost : dict
            'cost per unit': float [€/kWh]
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
        
        size = self.nom_capacity # kWH
        
        if tech_cost['cost per unit'] == 'default price correlation':
            C0 = 800 # €/kWh
            scale_factor = 0.8 # 0:1
            C = size * C0 **  scale_factor
        else:
            C = size * tech_cost['cost per unit']

        tech_cost['total cost'] = tech_cost.pop('cost per unit')
        tech_cost['total cost'] = C
        tech_cost['OeM'] = tech_cost['OeM'] *C /100 # €


        self.cost = tech_cost


###########################################################################################################################################################

if __name__ == "__main__":
    
    """
    Functional test
    """
    
    
    
    
    
    
    
    
    
    
    
    
    