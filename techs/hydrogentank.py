import numpy as np
from CoolProp.CoolProp import PropsSI

class H_tank:    
    
    def __init__(self,parameters,simulation_hours):
        
        """
        Create an H tank object
    
        parameters : dictionary
            'max capacity': float [kg]
            'pressure': float [bar]
            'self dischare': ?
                      
        output : H tank object able to:
            supply or abrosrb hydrogen .use(h,hyd)
            record the state of charge .SOC
            calculate its own volume (pressure) .volume(pressure)
        """
        
        self.pressure = parameters['pressure'] # H tank storage pressure
        self.SOC = np.zeros(simulation_hours+1) # array H tank State of Charge 
        self.max_capacity = parameters['max capacity'] # H tank max capacity [kg]
        self.used_capacity = 0 # H tank used capacity <= max_capacity [kg]      
        
    def use(self,h,hyd):
        """
        The H tank can supply or absorb hydrogen
     
        h: int hour to be simulated
        hyd: hydrogen requested (hyd<0) or provided (hyd>0) [kg]
      
        output : hydrogen supplied or absorbed that hour [kg]
        """
        
        if hyd >= 0: # charge H tank
            
            charge = min(hyd,self.max_capacity-self.SOC[h]) # how much hydrogen can H tank absorb?
            self.SOC[h+1] = self.SOC[h]+charge # charge H tank
            
            if self.SOC[h+1] > self.used_capacity: # update used capacity
                self.used_capacity = self.SOC[h+1]      
            
            return(-charge) # return hydrogen absorbed
    
            
        else: # discharge H tank (this logic allows to back-calculate the SOC[0], it's useful for long term storage systems)
            
            if (self.used_capacity==self.max_capacity):  # the max_capacity has been reached, so SOC[h+1] can't become negative 
                   
                discharge = min(-hyd,self.SOC[h]) # how much hydrogen can H tank supply?
                self.SOC[h+1] = self.SOC[h]-discharge # discharge H tank
        
        
            else: # the max_capacity has not yet been reached, so SOC[h+1] may become negative and then the past SOC may be translated   
                                                  
                discharge = min(-hyd,self.SOC[h]+self.max_capacity-self.used_capacity) # how much hydrogen can H tank supply?
                self.SOC[h+1] = self.SOC[h]-discharge # discharge H tank
                if self.SOC[h+1] < 0: # if the state of charge has become negative
                    self.used_capacity += - self.SOC[h+1] # incrase the used capacity
                    self.SOC[:h+2] += - self.SOC[h+1]  # traslate the past SOC array
                    
            return(discharge) # return hydrogen supplied
        
    def SOC_volume(self):
        """
        [kg] --> [m3]
        function of H tank pressure
      
        output : SOC [m3]
        """
        
        SOC = self.SOC / PropsSI('D', 'T', 298.15, 'P', self.pressure*1e5, 'H2')  # H tank SOC [m3]    
        
        return(SOC)