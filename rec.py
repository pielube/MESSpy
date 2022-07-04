import numpy as np
import pickle
import os
from location import location

class REC:
    
    def __init__(self,structure,general,path):
        """
        Create a Renewable Energy Comunity object composed of several locations (producers, consumers, prosumers)
    
        structure : dictionary (all the inputs are optional)
            'location_1_name': inputs required to create a location object (see Location.py)
            'location_2_name': 
                ...
            'location_n_name':
                
        general : dictionary
            'simulation years': number of years to be simulated
            'latitude': float
            'longitude': float
            'reference year': int year [2005 - 2015] used for output data and to get data from PVGIS if TMY = False
            'TMY': bool true if data of TMY is to be used              
                        
        output : REC object able to:
            simulate the energy flows of each present locations .REC_simulation
            record REC energy balances .energy_balance (electricity, heat, cool, gas and hydrogen) 
        """

        self.locations = {} # initialise REC locations dictionary
        self.energy_balance = {'electricity': {}, 'heat': {}, 'cool': {}, 'hydrogen': {}, 'gas': {}} # initialise energy balances dictionaries of each energy carrier
        
        self.simulation_hours = int(general['simulation years']*8760) # hourly timestep  
        
        ### create location objects and add them to the REC locations dictionary
        for location_name in structure: # location_name are the keys of 'structure' dictionary and will be used as keys of REC 'locations' dictionary too
            self.locations[location_name] = location(structure[location_name],general,location_name,path) # create location object and add it to REC 'locations' dictionary                
            
                  
    def REC_energy_simulation(self):
        """
        Simulate the REC every hour
        
        output :
            updating location energy balances
            updating REC energy balances
        """
        
        ### initialise REC electricity balances
        self.energy_balance['electricity']['from grid'] = np.zeros(self.simulation_hours) # array of electricity withdrawn from the grid from the whole rec
        self.energy_balance['electricity']['into grid'] = np.zeros(self.simulation_hours) # array of electricity withdrawn from the grid
        self.energy_balance['electricity']['collective self consumption'] = np.zeros(self.simulation_hours) # array of collective self consumed electricity from the whole rec
        
        
        for h in range(self.simulation_hours): # h: hour to simulate from 0 to simulation_hours 
            for location_name in self.locations: # each locations 
                self.locations[location_name].loc_energy_simulation(h) # simulate a single location updating its energy balances
                
            ### solve electricity grid 
                if self.locations[location_name].energy_balance['electricity']['grid'][h] < 0:
                    self.energy_balance['electricity']['into grid'][h] += self.locations[location_name].energy_balance['electricity']['grid'][h] # electricity fed into the grid from the whole rec at hour h
                else:                                                     
                    self.energy_balance['electricity']['from grid'][h] += self.locations[location_name].energy_balance['electricity']['grid'][h] # electricity withdrawn from the grid the whole rec at hour h
                
                
            ### calculate collective self consumption and who contributed to it
            self.energy_balance['electricity']['collective self consumption'][h] = min(-self.energy_balance['electricity']['into grid'][h],self.energy_balance['electricity']['from grid'][h]) # calculate REC collective self consumption how regulation establishes      
            
            if self.energy_balance['electricity']['collective self consumption'][h] > 0:
                for location_name in self.locations:
                    if self.locations[location_name].energy_balance['electricity']['grid'][h] < 0: # contribution as producer
                        self.locations[location_name].energy_balance['electricity']['collective self consumption'][h] = - self.energy_balance['electricity']['collective self consumption'][h] * self.locations[location_name].energy_balance['electricity']['grid'][h] / self.energy_balance['electricity']['into grid'][h]
                    else: # contribution as consumer
                        self.locations[location_name].energy_balance['electricity']['collective self consumption'][h] = self.energy_balance['electricity']['collective self consumption'][h] * self.locations[location_name].energy_balance['electricity']['grid'][h] / self.energy_balance['electricity']['from grid'][h]


            ### solve smart batteries
            for location_name in self.locations:
                
                # battery.collective = 1: 
                # REC tels to location how mutch electricity can be absorbed or supplied by battery every hour, without decreasing the collective-self-consumption
                
                if 'battery' in self.locations[location_name].technologies and self.locations[location_name].technologies['battery'].collective == 1:
                    
                    # how much energy can be absorbed or supplied by the batteries cause it's not usefull for collective-self-consumption
                    E = - self.locations[location_name].energy_balance['electricity']['grid'][h] + self.locations[location_name].energy_balance['electricity']['collective self consumption'][h]
                      
                    self.locations[location_name].energy_balance['electricity']['battery'][h] = self.locations[location_name].technologies['battery'].use(h,E) # electricity absorbed(-) by battery
                    self.locations[location_name].energy_balance['electricity']['grid'][h] += - self.locations[location_name].energy_balance['electricity']['battery'][h] # update grid balance (locatiom)
                  
                    if self.locations[location_name].energy_balance['electricity']['battery'][h] < 0:
                        self.energy_balance['electricity']['into grid'][h] += - self.locations[location_name].energy_balance['electricity']['battery'][h] # update grid balance (rec)
                    else:
                        self.energy_balance['electricity']['from grid'][h] += - self.locations[location_name].energy_balance['electricity']['battery'][h] # update grid balance (rec)


    def save(self,simulation_name):
        """
        Save REC and each location energy balances
        
        simulationa_name : str 
        
        output: 
            balances/simulation_name.pkl
            LOC/simulation_name.pkl
        """
        
        balances = {}
        LOC = {}
        ageing = {}
        balances['REC'] = self.energy_balance
        
        for location_name in self.locations:
            balances[location_name] = self.locations[location_name].energy_balance
            
            LOC[location_name] = {}
            ageing[location_name] = {}
            
            tech_name = 'battery'
            if tech_name in self.locations[location_name].technologies:
                LOC[location_name][tech_name] = self.locations[location_name].technologies[tech_name].LOC
                if self.locations[location_name].technologies[tech_name].ageing:
                    ageing[location_name][tech_name] = [self.locations[location_name].technologies[tech_name].replacements,self.locations[location_name].technologies[tech_name].ageing_history]
                
            tech_name = 'H tank'
            if tech_name in self.locations[location_name].technologies:
                LOC[location_name][tech_name] = self.locations[location_name].technologies[tech_name].LOC
        
        directory = './results'
        if not os.path.exists(directory):
            os.makedirs(directory)
        
        with open('results/balances_'+simulation_name+".pkl", 'wb') as f:
            pickle.dump(balances, f) 
            
        with open('results/LOC_'+simulation_name+".pkl", 'wb') as f:
            pickle.dump(LOC, f) 
            
        with open('results/ageing_'+simulation_name+".pkl", 'wb') as f:
            pickle.dump(ageing, f) 
            
            
            
    def reset(self):
        """
        Use this function before a new simulation if you don't want to recreate the rec object
        
        output: initialise LOC
        """
        
        to_reset = ['battery','H tank'] # technologies for which the LOC must be reset
        for location_name in self.locations: # each location
            for tech_name in to_reset:
                if tech_name in self.locations[location_name].technologies: 
                    self.locations[location_name].technologies[tech_name].LOC = np.zeros(self.simulation_hours+1) # array level of Charge 
                    self.locations[location_name].technologies[tech_name].used_capacity = 0 # used capacity <= max_capacity   
    
        
   