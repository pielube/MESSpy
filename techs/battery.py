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
            'charging efficiency': float efficiency of charging process
            'discharging efficiency': float efficiency of discharging process
            'depth of discharge': minimum SOC float
            'self discharge rate': hourly SOC loss for self-discharge float                                                        
        
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
        
        self.DoD = parameters['depth of discharge']  # float 
        
        self.self_discharge = parameters['self discharge rate'] / 100 # float [%]                                                                                        
        self.ageing = parameters['ageing'] # bool true if aging has to be calculated
        self.LC = parameters['life cycles'] # int number of life cycles to reach the end of battery life
        self.EOL = parameters['end life capacity'] # end of life capacity                                                                 
        self.deg = self.max_capacity*(100-self.EOL)/100 # float capacity that is going to be degradated in life cycles
        self.ageing_day = 7 # How often ageing has to bee calculated? [days]
        self.completed_cycles = 0 # float initialise completed_cycles, this parameter is usefull to calculate replacements
        self.replacements = [] # list initialise: h at which replecaments occur

        self.LOC = np.zeros(c.timestep_number+1) # array battery level of Charge 
        self.used_capacity = 0 # battery used capacity <= max_capacity [kWh]
      
        self.collective = parameters['collective'] # int 0: no collective rules. 1: priority to csc and then charge or discharge the battery.
       
        self.T = 52.5 # °C Operative temperature
        
        self.SOH = np.zeros(int(c.timestep_number/(self.ageing_day*60*24/c.timestep))+1) # % state of health
        self.SOH[0] = 1
        self.SOH_cal = np.zeros(len(self.SOH))
        self.SOH_cal[0] = 1
        self.SOH_cyc = np.zeros(len(self.SOH))
        self.SOH_cyc[0] = 1
        self.D_cal = np.zeros(len(self.SOH))
        self.D_cyc = np.zeros(len(self.SOH))
        self.D_tot = np.zeros(len(self.SOH))
        
        self.ageing_history = [[0],[self.max_capacity],[self.SOH],[self.SOH_cal],[self.SOH_cyc]] # list initialise ageing history. [completed_cycles],[max_capacity]            
        
    def use(self,step,p):
        """
        The battery can supply or absorb electricity
     
        step: int step to be simulated
        timestep: int selected time resolution for the simulation [min]
        p: power requested (p<0) or provided (p>0) [kW]
      
        output : electricity supplied or absorbed that step [kW]
        """
        #Apply self_discharge 

        self.LOC[step] = self.LOC[step]*(1-self.self_discharge)
        if self.ageing and (step*c.timestep/60/24%self.ageing_day == 0) and step!=0: # if aeging == True and it's time to calculate it
            self.calculate_ageing(step)
            
        if p >= 0: # charge battery
        
            if p > self.MpowerC:
                p = self.MpowerC
                              
            charge = p*self.etaC*c.P2E
            if charge < (self.max_capacity-self.LOC[step]): # if battery can't be full charged [kWh]
                self.LOC[step+1] = self.LOC[step]+charge # charge [kWh]
            
            else: # if batter can be full charged
                self.LOC[step+1] = self.max_capacity
                p = ((self.max_capacity-self.LOC[step]) / self.etaC)/c.P2E
            
            if self.LOC[step+1] > self.used_capacity: # update used capacity
                self.used_capacity = self.LOC[step+1] 
                   
            return(-p) # return electricity absorbed
        
            
        else: # discharge battery (this logic allows to back-calculate the LOC[0], it's useful for long term storage systems)
            
            p = p/self.etaD # how much power is really required
            
            min_LOC = self.max_capacity*self.DoD
            
            if self.LOC[step] + (self.nom_capacity-self.used_capacity) > min_LOC:
                
                if(self.used_capacity==self.nom_capacity):  # the nom_capacity has been reached, so LOC[step+1] can't become negative 
                       
                    discharge = min(-p,(self.LOC[step]-min_LOC)/c.P2E,self.max_capacity*self.MpowerD) # how much power can battery supply? [kW]         
                    self.LOC[step+1] = self.LOC[step]-discharge*c.P2E # discharge battery
                    
                else: # the max_capacity has not yet been reached, so LOC[step+1] may become negative and then the past LOC may be translated   
                                                      
                    discharge = min(-p,(self.LOC[step]-min_LOC+self.max_capacity-self.used_capacity)/c.P2E,self.max_capacity*self.MpowerD) # how much power can battery supply
                    self.LOC[step+1] = self.LOC[step]-discharge*c.P2E # discharge battery
                    if self.LOC[step+1] < min_LOC: # if the level of charge has become negative
                        self.used_capacity += (min_LOC - self.LOC[step+1]) # incrase the used capacity
                        self.LOC[:step+2] += (min_LOC - self.LOC[step+1])  # traslate the past LOC array
              
                return(discharge*self.etaD) # return electricity supplied
            
            else:
                #Battery is below minimum SOC, can't be discharged further
                self.LOC[step+1] = self.LOC[step]
                return 0                                            
        
    def calculate_ageing(self,step):      
        
        # degradation (equivalent number of cycles, life cycles, end of life capacity)
        cycles = self.rainflow(step,c.timestep) # number of equivalent cycles completed
                                                                  
        self.completed_cycles += cycles
        
        SOH_step = int(step/(self.ageing_day*60*24/c.timestep))
        
        k_ageing = 0.003273 + 0.00155*((self.T-52.5)/7.5) + 0.0000640*((self.T-52.5)/7.5)**2 + 0.00146*((self.LOC[step]/self.max_capacity-0.6)/0.2) + 0.000163*((self.LOC[step]/self.max_capacity-0.6)/0.2)**2 + 0.000687*((self.T-52.5)/7.5)*((self.LOC[step]/self.max_capacity-0.6)/0.2)
        alpha = 0.000323*c.NEPERO**(3586.3/(self.T+273.15))
        
        # delta_t = (self.ageing_day*60*24/c.timestep)
        
        D_cal = k_ageing*((2-self.SOH[SOH_step-1])**(-alpha))

        D_cyc = cycles/self.LC
        
        D_tot = D_cal + D_cyc
        
        self.SOH[SOH_step] = self.SOH[SOH_step-1] - D_tot*(1-self.EOL/100)
        
        self.SOH_cal[SOH_step] = self.SOH_cal[SOH_step-1] - D_cal*(1-self.EOL/100)
        
        self.SOH_cyc[SOH_step] = self.SOH_cyc[SOH_step-1] - D_cyc*(1-self.EOL/100)
        
        self.max_capacity = self.nom_capacity * (self.SOH[SOH_step])
        self.ageing_history[0].append(self.completed_cycles)
        self.ageing_history[1].append(self.max_capacity)
        self.ageing_history[2] = self.SOH
        self.ageing_history[3] = self.SOH_cal
        self.ageing_history[4] = self.SOH_cyc
        
        if self.SOH[SOH_step] <= self.EOL/100:
            self.replacements.append(step)
            self.completed_cycles = 0
            self.max_capacity = self.nom_capacity
            self.SOH[SOH_step] = 1                                  
        
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
    
    
    
    
    
    
    
    
    
    
    
    
    