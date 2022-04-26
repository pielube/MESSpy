import numpy as np

class battery:    
    
    def __init__(self,parameters,simulation_hours):
        """
        Create a battery object
    
        parameters : dictionary
            'max capacity': float [kWh]
            'E-rate': float charging and discharging rate related to maximum capacity [kW/kWh]
            'ageing': bool true if aging has to be calculated
            'collective': int 0: no collective rules. 1: charging collective rules. 2: charging and disharging collective rules.
            
        output : battery object able to:
            supply or abrosrb electricity .use(h,e)
            record the level of charge .LOC
            take account of ageing .calculate_aging()   
        """
                
        self.ageing = parameters['ageing'] # bool true if aging has to be calculated
        self.collective = parameters['collective'] # int 0: no collective rules. 1: charging collective rules. 2: charging and disharging collective rules.
        self.max_capacity = parameters['max capacity'] # battery max capacity [kWh]
        self.E_rate= parameters['E-rate'] # battery E-rate [kW/kWh]
        
        self.LOC = np.zeros(simulation_hours+1) # array battery level of Charge 
        self.used_capacity = 0 # battery used capacity <= max_capacity [kWh]
        
        self.time_used = 0
        
    def use(self,h,e):
        """
        The battery can supply or absorb electricity
     
        h: int hour to be simulated
        e: electricity requested (e<0) or provided (e>0) [kWh]
      
        output : electricity supplied or absorbed that hour [kWh]
        """
        
        self.time_used += 1
        
        if e >= 0: # charge battery
            
            charge = min(e,self.max_capacity-self.LOC[h],self.max_capacity*self.E_rate) # how much electricity can battery absorb?
            self.LOC[h+1] = self.LOC[h]+charge # charge battery
            
            if self.LOC[h+1] > self.used_capacity: # update used capacity
                self.used_capacity = self.LOC[h+1] 
                   
            return(-charge) # return electricity absorbed
        
            
        else: # discharge battery (this logic allows to back-calculate the LOC[0], it's useful for long term storage systems)
            
            if(self.used_capacity==self.max_capacity):  # the max_capacity has been reached, so LOC[h+1] can't become negative 
                   
                discharge = min(-e,self.LOC[h],self.max_capacity*self.E_rate) # how much electricity can battery supply?
                self.LOC[h+1] = self.LOC[h]-discharge # discharge battery
                
            else: # the max_capacity has not yet been reached, so LOC[h+1] may become negative and then the past LOC may be translated   
                                                  
                discharge = min(-e,self.LOC[h]+self.max_capacity-self.used_capacity,self.max_capacity*self.E_rate) # how much electricity can battery supply?
                self.LOC[h+1] = self.LOC[h]-discharge # discharge battery
                if self.LOC[h+1] < 0: # if the level of charge has become negative
                    self.used_capacity += - self.LOC[h+1] # incrase the used capacity
                    self.LOC[:h+2] += - self.LOC[h+1]  # traslate the past LOC array
            
            return(discharge) # return electricity supplied
    
        
        if self.ageing and h == 100: ### ogni quanto calcolo l'aging e aggiorno capacity???
            pass
            battery.calculate_aging()
        
    def calculate_ageing(self):
        pass
    
        # degradation 
        # corrosion
    


###########################################################################################################################################################

if __name__ == "__main__":
    
    """
    Functional test
    """
    
    inp_test = {'max capacity': 10,
                'E-rate': 0.5,
                'ageing': False,
                'collective': 0
                }
    
    b1 = battery(inp_test,simulation_hours=8760)   
    
    flow = [+5,-2,-15,+3,-1,-6,+20,0]
    #flow = [-2,-10, +2, +10, +1, -4, -20]
    for h in range(len(flow)):
        b1.use(h,flow[h])        
    
    print(b1.LOC[:15])
    
    
    
    
    
    
    
    
    
    
    
    