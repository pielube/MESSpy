import numpy as np

class battery:    
    
    def __init__(self,parameters,simulation_hours):
        """
        Create a battery object
    
        parameters : dictionary
            'max capacity': float [kWh]
            'ageing': bool true if aging has to be calculated
            'self discharge': ?
                      
        output : battery object able to:
            supply or abrosrb electricity .use(h,e)
            record the state of charge .SOC
            take account of ageing .calculate_aging()   
        """
                
        self.SOC = np.zeros(simulation_hours+1) # array battery State of Charge 
        self.ageing = parameters['ageing'] # bool true if aging has to be calculated
        self.max_capacity = parameters['max capacity'] # battery max capacity [kWh]
        self.used_capacity = 0 # battery used capacity <= max_capacity [kWh]
               
    def use(self,h,e):
        """
        The battery can supply or absorb electricity
     
        h: int hour to be simulated
        e: electricity requested (e<0) or provided (e>0) [kWh]
      
        output : electricity supplied or absorbed that hour [kWh]
        """
        
        if e >= 0: # charge battery
            
            charge = min(e,self.max_capacity-self.SOC[h]) # how much electricity can battery absorb?
            self.SOC[h+1] = self.SOC[h]+charge # charge battery
            
            if self.SOC[h+1] > self.used_capacity: # update used capacity
                self.used_capacity = self.SOC[h+1] 
                
            return(-charge) # return electricity absorbed
            
        else: # discharge battery (this logic allows to back-calculate the SOC[0], it's useful for long term storage systems)
            
            if(self.used_capacity==self.max_capacity):  # the max_capacity has been reached, so SOC[h+1] can't become negative 
                   
                discharge = min(-e,self.SOC[h]) # how much electricity can battery supply?
                self.SOC[h+1] = self.SOC[h]-discharge # discharge battery
                
            else: # the max_capacity has not yet been reached, so SOC[h+1] may become negative and then the past SOC may be translated   
                                                  
                discharge = min(-e,self.SOC[h]+self.max_capacity-self.used_capacity) # how much electricity can battery supply?
                self.SOC[h+1] = self.SOC[h]-discharge # discharge battery
                if self.SOC[h+1] < 0: # if the state of charge has become negative
                    self.used_capacity += - self.SOC[h+1] # incrase the used capacity
                    self.SOC[:h+2] += - self.SOC[h+1]  # traslate the past SOC array

            return(discharge) # return electricity supplied
        
        
# =============================================================================
#         if self.ageing and h == 100: ### ogni quanto calcolo l'aging e aggiorno capacity???
#             battery.calculate_aging()
#         
#     def calculate_ageing(self):
#         ### based on SOC self.max_capacity must decrease 
# =============================================================================
    


###########################################################################################################################################################

if __name__ == "__main__":
    
    """
    Functional test
    """
    
    inp_test = {'max capacity': 10,
                'ageing': False
                }
    
    b1 = battery(inp_test,simulation_hours=8760)   
    
    flow = [+5,-2,-15,+3,-1,-6,+20,0]
    #flow = [-2,-10, +2, +10, +1, -4, -20]
    for h in range(len(flow)):
        b1.use(h,flow[h])        
    
    
    
    
    
    
    
    
    
    
    
    