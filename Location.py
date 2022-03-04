import numpy as np
import pandas as pd
from techs import PV,battery,H_tank,fuel_cell,electrolyzer
# from battery import battery
# from hydrogentank import H_tank
# from fuelcell import fuel_cell
# from electrolyzer import electrolyzer

class location:
    
    def __init__(self,system,simulation_hours):
        """
        Create a location object (producer, consumer or prosumer) 
    
        system: dictionary (all the inputs are optional)
            'demand'': dictionary
                'electricity': str 'file_name.csv' hourly time series of electricity demand 8760 values [kWh]
                'hydrogen': str 'file_name.csv' hourly time series of hydrogen demand 8760 values [kg/hr]
                'gas': str 'file_name.csv' hourly time series of gas demand 8760 values [kg/hr]
                'heat': str 'file_name.csv' hourly time series of heating demand 8760 values [kWh]
                'cool': str 'file_name.csv' hourly time series of cooling demand 8760 values [kWh]
            'PV': dictionary parameters needed to create PV object (see PV.py)
            'battery': dictionary parameters needed to create battery object (see Battery.py)
            'hydrogen demand': str 'file_name.csv' hourly time series of hydrogen demand [kg]
            'electrolyzer': dictionary parameters needed to create electrolyzer object (see electrolyzer.py)
            'H tank': dictionary parameters needed to create H_tank object (see H_tank.py)
            'fuel cell': dictionary parameters needed to create fuel_cell object (see fuel_cell.py)
            
        output : location object able to:
            simulate the energy flows of present technologies .loc_simulation
            record energy balances .energy_balance (electricity, heat, gas and hydrogen)
        """
        
        self.technologies = {} # initialise technologies dictionary
        self.energy_balance = {'electricity': {}, 'heat': {}, 'cool': {}, 'hydrogen': {}, 'gas': {}} # initialise energy balances dictionaries
        
        # create the objects of present technologies and add them to the technologies dictionary
        # initialise energy balance and add them to the energy_balance dictionary
        
        for carrier in system['demand']:            
            self.energy_balance[carrier]['demand'] = np.tile(pd.read_csv('Loads/'+system['demand'][carrier])['0'].to_numpy(),int(simulation_hours/8760)) # hourly energy carrier needed for the entire simulation
            self.energy_balance[carrier]['in'] = np.zeros(simulation_hours) # array energy carrier bought from the grid
  
        if 'PV' in system:
            self.technologies['PV'] = PV(system['PV'],simulation_hours) # PV object created and add to 'technologies' dictionary
            self.energy_balance['electricity']['production'] = np.zeros(simulation_hours) # array of electricity production
            self.energy_balance['electricity']['out'] = np.zeros(simulation_hours) # array of electricity feed in to the grid
            
        if 'battery' in system:
            self.technologies['battery'] = battery(system['battery'],simulation_hours) # battery object created and to 'technologies' dictionary
            
        if 'electrolyzer' in system:
            self.technologies['electrolyzer'] = electrolyzer(system['electrolyzer'],simulation_hours) # electrolyzer object created and to 'technologies' dictionary
            self.energy_balance['hydrogen']['production'] = np.zeros(simulation_hours) # array of hydrogen production
            self.energy_balance['hydrogen']['out'] = np.zeros(simulation_hours) # array of hydrogen feed in to the grid
            
        if 'fuel cell' in system:
            self.technologies['fuel cell'] = fuel_cell(system['fuel cell'],simulation_hours) # Fuel cell object created and to 'technologies' dictionary
            
        if 'H tank' in system:
            self.technologies['H tank'] = H_tank(system['H tank'],simulation_hours) # H tank object created and to 'technologies' dictionary
        
        ### determine the location type: it could be usefull to decide a simulation order: work in progress...
        if 'demand' in self.energy_balance['electricity']:
            self.location_type = 'consumer' # it consumes energy
            if 'production' in self.energy_balance['electricity']:            
                self.location_type = 'prosumer' # it produces and consumes energy
        else:  
            self.location_type = 'producer' # it produces energy
  
            
    def loc_energy_simulation(self,h):
        """
        Simulate the location
        
        input : int hour to simulate
    
        output : updating of location energy balances
        """
        
        EB = {'electricity': 0, 'heat': 0, 'hydrogen': 0, 'gas': 0} # initialise Energy Balances     
        
        for carrier in EB: # for each energy carrier
            if 'demand' in self.energy_balance[carrier]:                
                EB[carrier] += -self.energy_balance[carrier]['demand'][h] # energy balance update: - energy demand
             
        if 'PV' in self.technologies: 
            self.energy_balance['electricity']['production'][h] += self.technologies['PV'].use(h) # electricity produced from PV
            EB['electricity'] += self.energy_balance['electricity']['production'][h] # electricity balance update: + electricity produced from PV
        
        if 'battery' in self.technologies:
            EB['electricity'] +=  self.technologies['battery'].use(h,EB['electricity']) # electricity balance update: +- electricity absorbed or supplied by battery
            
        if 'electrolyzer' in self.technologies:
            if EB['electricity'] > 0:
                if 'H tank' in self.technologies: # if hydrogen is stored in a tank
                    storable_hydrogen = self.technologies['H tank'].max_capacity-self.technologies['H tank'].SOC[h-1] # the tank can't be full
                else: # if hydrogen is sold to the grid 
                    storable_hydrogen = 99999999999 # there are no limits f.i an hydrogen producer
                hyd,ele = self.technologies['electrolyzer'].use(h,EB['electricity'],storable_hydrogen)
                self.energy_balance['hydrogen']['production'][h] += hyd # array of hydrogen produced
                EB['hydrogen'] += hyd
                EB['electricity'] += ele
                
        if 'fuel cell' in self.technologies:
            if EB['electricity'] < 0:                
                hyd,ele = self.technologies['fuel cell'].use(h,EB['electricity'],self.technologies['H tank'].SOC[h-1])
                EB['hydrogen'] += hyd
                EB['electricity'] += ele
                
        if 'H tank' in self.technologies:
            EB['hydrogen'] += self.technologies['H tank'].use(h,EB['hydrogen'])
        
        for carrier in EB:
            if 'out' in self.energy_balance[carrier] and EB[carrier] > 0:
                self.energy_balance[carrier]['out'][h] += EB[carrier]
                # and from_grid = 0
                
            if 'in' in self.energy_balance[carrier] and EB[carrier] < 0:
                self.energy_balance[carrier]['in'][h] += -EB[carrier]
                # and to grid = 0
            
            
            
                
        
        

