import numpy as np
import pickle
import os
from location import location

class REC:
    
    def __init__(self,structure,simulation_years):
        """
        Create a Renewable Energy Comunity object composed of several locations (producers, consumers, prosumers)
    
        structure : dictionary (all the inputs are optional)
            'location_1_name': inputs required to create a location object (see Location.py)
            'location_2_name': 
                ...
            'location_n_name':
                        
        output : REC object able to:
            simulate the energy flows of each present locations .REC_simulation
            record REC energy balances .energy_balance (electricity, heat, cool, gas and hydrogen) 
        """
        
        self.simulation_hours = int(simulation_years*8760) # hourly timestep      
        
        self.energy_balance = {'electricity': {}, 'heat': {}, 'cool': {}, 'hydrogen': {}, 'gas': {}} # initialise energy balances dictionaries of each energy carrier
       
        self.locations = {} # initialise REC locations dictionary
        
        ### create location objects and add them to the REC locations dictionary
        for location_name in structure: # location_name are the keys of 'structure' dictionary and will be used as keys of REC 'locations' dictionary too
            self.locations[location_name] = location(structure[location_name],self.simulation_hours) # create location object and add it to REC 'locations' dictionary                
            
                  
    def REC_energy_simulation(self):
        """
        Simulate the REC every hour
        
        output :
            updating location energy balances
            updating REC energy balances
        """
        
        ### initialise REC electricity balances
        self.energy_balance['electricity']['from grid'] = np.zeros(self.simulation_hours) # array of collective self consumed energy
        self.energy_balance['electricity']['into grid'] = np.zeros(self.simulation_hours) # array of collective self consumed energy
        self.energy_balance['electricity']['collective self consumption'] = np.zeros(self.simulation_hours) # array of collective self consumed energy
        
        for h in range(self.simulation_hours): # h: hour to simulate from 0 to simulation_hours 
            for location_name in self.locations: # each locations 
                self.locations[location_name].loc_energy_simulation(h) # simulate a single location updating its energy balances
                
            ### solve electricity grid 
                if self.locations[location_name].energy_balance['electricity']['grid'][h] < 0:
                    self.energy_balance['electricity']['into grid'] += self.locations[location_name].energy_balance['electricity']['grid'][h]
                else:                                                     
                    self.energy_balance['electricity']['from grid'] += self.locations[location_name].energy_balance['electricity']['grid'][h]
                
            self.energy_balance['electricity']['collective self consumption'][h] = min(-self.energy_balance['electricity']['into grid'][h],self.energy_balance['electricity']['from grid'][h]) # calculate REC collective self consumption how regulation establishes      
    

    def save(self,simulation_name):
        """
        Save REC and each location energy balances
        
        simulationa_name : str 
        
        output: 
            balances/simulation_name.pkl
            soc/simulation_name.pkl
        
        """
        balances = {}
        SOC = {}
        balances['REC'] = self.energy_balance
        
        for location_name in self.locations:
            balances[location_name] = self.locations[location_name].energy_balance
            
            SOC[location_name] = {}
            
            tech_name = 'battery'
            if tech_name in self.locations[location_name].technologies:
                SOC[location_name][tech_name] = self.locations[location_name].technologies[tech_name].SOC
        
            tech_name = 'H tank'
            if tech_name in self.locations[location_name].technologies:
                SOC[location_name][tech_name] = self.locations[location_name].technologies[tech_name].SOC_volume()
        
        directory = './results'
        if not os.path.exists(directory):
            os.makedirs(directory)
        
        with open('results/balances_'+simulation_name+".pkl", 'wb') as f:
            pickle.dump(balances, f) 
            
        with open('results/SOC_'+simulation_name+".pkl", 'wb') as f:
            pickle.dump(SOC, f) 
            
            
    def reset(self):
        pass
        # azzera gli SOC, i bilanci non dovrebbe esserci bisogno..? 
        
    
        
   