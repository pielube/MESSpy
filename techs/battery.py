import numpy as np
import os
import sys 
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),os.path.pardir)))   # temorarily adding constants module path 
from core import constants as c

class battery:    
    
    def __init__(self,parameters):
        """
        Create a battery object
    
        parameters : dictionary
            'nominal capacity': float [kWh]
            'max charging power': maximum input power [kW]
            'max discharging power': maximum output power [kW]
            'charging efficiency': float
            'discharging efficiency': float
        
            'ageing': bool true if aging has to be calculated
            'life cycles': int number of life cycles to reach the end of battery life
            'end life capacity': float maximum capacity left at end of life [%]
            
            'collective': int 0: no collective rules. 1: priority to csc and then charge or discharge the battery.
        
        output : battery object able to:
            supply or abrosrb electricity .use(h,e)
            record the level of charge .LOC
            take account of ageing .calculate_aging()   
        """
        
        self.cost = False # will be updated with tec_cost()

        self.nom_capacity = parameters['nominal capacity']*c.kWh2kJ # battery early life max capacity [kJ]
        self.max_capacity = parameters['nominal capacity']*c.kWh2kJ # battery max capacity [kJ]

        self.etaC = parameters['charging efficiency'] # float 
        self.etaD = parameters['discharging efficiency'] # float  
        
        self.MpowerC = parameters['max charging power'] # float [kW]
        self.MpowerD = parameters['max discharging power'] # float [kW]
        
        self.ageing = parameters['ageing'] # bool true if aging has to be calculated
        self.LC = parameters['life cycles'] # int number of life cycles to reach the end of battery life
        self.deg = self.max_capacity*(100-parameters['end life capacity'])/100 # float capacity that is going to be degradated in life cycles
        self.ageing_day = 7 # How often ageing has to bee calculated? [days]
        self.completed_cycles = 0 # float initialise completed_cycles, this parameter is usefull to calculate replacements
        self.replacements = [] # list initialise: h at which replecaments occur
        self.ageing_history = [[0],[self.max_capacity]] # list initialise ageing history. [completed_cycles],[max_capacity]
    
        self.LOC = np.zeros(c.timestep_number+1) # array battery level of Charge 
        self.used_capacity = 0 # battery used capacity <= max_capacity [kWh]
      
        self.collective = parameters['collective'] # int 0: no collective rules. 1: priority to csc and then charge or discharge the battery.

    def use(self,step,p):
        """
        The battery can supply or absorb electricity
     
        step: int step to be simulated
        timestep: int selected time resolution for the simulation [min]
        p: power requested (p<0) or provided (p>0) [kW]
      
        output : electricity supplied or absorbed that step [kW]
        """

        if self.ageing and (step*c.timestep/60/24%self.ageing_day == 0) and step!=0: # if aeging == True and it's time to calculate it
            self.calculate_ageing(step)
            
        if p >= 0: # charge battery
        
            if p > self.MpowerC:
                p = self.MpowerC
                              
            if p*self.etaC*c.P2E < (self.max_capacity-self.LOC[step]): # if battery can't be full charged [kWh]
                self.LOC[step+1] = self.LOC[step]+p*self.etaC*c.P2E # charge [kWh]
            
            else: # if batter can be full charged
                self.LOC[step+1] = self.max_capacity
                p = ((self.max_capacity-self.LOC[step]) / self.etaC)/c.P2E
            
            if self.LOC[step+1] > self.used_capacity: # update used capacity
                self.used_capacity = self.LOC[step+1] 
                   
            return(-p) # return electricity absorbed
        
            
        else: # discharge battery (this logic allows to back-calculate the LOC[0], it's useful for long term storage systems)
            
            p = p/self.etaD # how much power is really required
            if(self.used_capacity==self.nom_capacity):  # the nom_capacity has been reached, so LOC[step+1] can't become negative 
                   
                discharge = min(-p,self.LOC[step]/c.P2E,self.max_capacity*self.MpowerD) # how much power can battery supply? [kW]         
                self.LOC[step+1] = self.LOC[step]-discharge*c.P2E # discharge battery
                
            else: # the max_capacity has not yet been reached, so LOC[step+1] may become negative and then the past LOC may be translated   
                                                  
                discharge = min(-p,(self.LOC[step]+self.max_capacity-self.used_capacity)/c.P2E,self.max_capacity*self.MpowerD) # how much power can battery supply?
                self.LOC[step+1] = self.LOC[step]-discharge*c.P2E # discharge battery
                if self.LOC[step+1] < 0: # if the level of charge has become negative
                    self.used_capacity += - self.LOC[step+1] # incrase the used capacity
                    self.LOC[:step+2] += - self.LOC[step+1]  # traslate the past LOC array
            
            return(discharge*self.etaD) # return electricity supplied
        
    def calculate_ageing(self,step):      
        
        # degradation (equivalent number of cycles, life cycles, end of life capacity)
        cycles = self.rainflow(step,c.timestep) # number of equivalent cycles completed
        self.max_capacity += - (cycles / self.LC) * self.deg      
        self.completed_cycles += cycles
        self.ageing_history[0].append(self.completed_cycles)
        self.ageing_history[1].append(self.max_capacity)

        
        if self.completed_cycles > self.LC: # replacement
            self.replacements.append(step)
            self.completed_cycles = 0
            self.max_capacity = self.nom_capacity
        
    def rainflow(self,step,timestep):
        
        #https://ieeexplore.ieee.org/document/7741532
        
        LOC = self.LOC[int(step-self.ageing_day*24*60/timestep):step] # part of the LOC whose contribution to aging is to be calculated
        
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
        
        size = self.nom_capacity/c.kWh2kJ # kWh
        
        if tech_cost['cost per unit'] == 'default price correlation':
            C0 = 800 # €/kWh
            scale_factor = 0.8 # 0:1
            C = C0 * size **  scale_factor
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
    
    
    
    
    
    
    
    
    
    
    
    
    