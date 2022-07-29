import numpy as np

class intertialtank:    
    
    def __init__(self,parameters,simulation_hours):
        """
        Create a battery object
    
        parameters : dictionary
            
        
        output : thermal_tank object able to:
            
        """
        
        self.volume = parameters['volume'] # m3
        self.mass = self.volume * 1000 # kg
        self.cp = 4200 # J/kgK
        
        self.t_min_heat = parameters['t min heat'] # C°
        self.t_max_heat = parameters['t max heat'] # C°
        self.t_min_cool = parameters['t min cool'] # C°
        self.t_max_cool = parameters['t max cool'] # C°
        
        self.t = 25 # C°
        self.mode = "heat" # can be changed to "cool"
        
       
    def use(self,e):
        """
        The tank can supply or absorb thermal energy
     
        h: int hour to be simulated
        e: energy requested (e<0) or provided (e>0) [kWh]
      
        output : thermal energy supplied or absorbed that hour [kWh]
        
        """


        ## dispersione termica?!        
        ## t min_max ???
        
        # if e < 0 decrease tank temperature (hp switched on cooling or heat required)
        # if e > 0 incrase tank temperature (hp switched on heating or cool required)
        t = self.t + (e*3600*1000) / (self.mass*self.cp) # J / kg
        
        if self.mode == "heat":            
                
            if t < self.t_min_heat: # tank riched the minimum temperature (irradiation system temperature), so heat pump must switch on
                e = self.mass*self.cp*(self.t_min_heat-self.t) # heat that can be provided by the tank
                self.t = self.t_min_heat # minimum temperature riched
                
            if t > self.t_max_heat: # tank reached the maximum temperature (heat pump output temperature), so heat pump must switch off
                e = self.mass*self.cp*(self.t_max_heat-self.t)
                self.t = self.t_max_heat
                ### heat pump regulation?
                
            else:
                self.t = t
                
        
        if self.mode == "cool":
            
            if t < self.t_min_cool: # tank riched the minimum temperature (heat pumpt output temperature), so heat pump must switch off
                e = self.mass*self.cp*(self.t_min_cool-self.t)
                self.t = self.t_min_cool
                ### heat pump regulation?
                
            if t > self.t_max_cool: # tank reached the maximum temperature (irradiation system temperature), so heat pump must switch on
                e = self.mass*self.cp*(self.t_max_cool-self.t)
                self.t = self.t_max_cool
            
            else:
                self.t = t
            
        return(-e) 
    
    
            
            
            