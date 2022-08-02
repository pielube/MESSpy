import numpy as np
import pandas as pd
from techs import heatpump,boiler_el,boiler_ng,PV,battery,H_tank,fuel_cell,electrolyzer


class location:
    
    def __init__(self,system,general,location_name,path):
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
            
        general: dictionary 
            see rec.py
            
        output : location object able to:
            simulate the energy flows of present technologies .loc_simulation
            record energy balances .energy_balance (electricity, heat, gas and hydrogen)
        """
        
        self.name = location_name
        self.technologies = {} # initialise technologies dictionary
        self.energy_balance = {'electricity': {}, 'heat': {}, 'hydrogen': {}, 'gas': {}} # initialise energy balances dictionaries
       
        self.simulation_hours = int(general['simulation years']*8760) # hourly timestep     
        
        # create the objects of present technologies and add them to the technologies dictionary
        # initialise energy balance and add them to the energy_balance dictionary
        
        for carrier in self.energy_balance:
            if system['grid'][carrier]:
                self.energy_balance[carrier]['grid'] = np.zeros(self.simulation_hours) # array energy carrier bought from the grid (-) or feed into the grid (+)

            if carrier in system['demand']:            
                self.energy_balance[carrier]['demand'] = - np.tile(pd.read_csv(path+'/loads/'+system['demand'][carrier])['0'].to_numpy(),int(self.simulation_hours/8760)) # hourly energy carrier needed for the entire simulation

        if 'heatpump' in system:
            self.technologies['heatpump'] = heatpump(system['heatpump']) # heatpump object created and add to 'technologies' dictionary
            self.energy_balance['electricity']['heatpump'] = np.zeros(self.simulation_hours) # array heatpump electricity balance
            self.energy_balance['heat']['heatpump'] = np.zeros(self.simulation_hours) # array heatpump heat balance
            self.energy_balance['heat']['inertial tank'] = np.zeros(self.simulation_hours) # array inertial tank heat balance
            
        if 'boiler_el' in system:
            self.technologies['boiler_el'] = boiler_el(system['boiler_el']) # boiler_el object created and add to 'technologies' dictionary
            self.energy_balance['electricity']['boiler_el'] = np.zeros(self.simulation_hours) # array boiler_el electricity balance
            self.energy_balance['heat']['boiler_el'] = np.zeros(self.simulation_hours) # array boiler_el heat balance
            
        if 'boiler_ng' in system:
            self.technologies['boiler_ng'] = boiler_ng(system['boiler_ng']) # boiler_ng object created and add to 'technologies' dictionary
            self.energy_balance['gas']['boiler_ng'] = np.zeros(self.simulation_hours) # array boiler_ng gas balance
            self.energy_balance['heat']['boiler_ng'] = np.zeros(self.simulation_hours) # array boiler_ng heat balance 

        if 'PV' in system:
            self.technologies['PV'] = PV(system['PV'],general,self.simulation_hours,self.name,path) # PV object created and add to 'technologies' dictionary
            self.energy_balance['electricity']['PV'] = np.zeros(self.simulation_hours) # array PV electricity balance
            
        if 'battery' in system:
            self.technologies['battery'] = battery(system['battery'],self.simulation_hours) # battery object created and to 'technologies' dictionary
            self.energy_balance['electricity']['battery'] = np.zeros(self.simulation_hours) # array battery electricity balance
                           
        if 'electrolyzer' in system:
            self.technologies['electrolyzer'] = electrolyzer(system['electrolyzer'],self.simulation_hours) # electrolyzer object created and to 'technologies' dictionary
            self.energy_balance['electricity']['electrolyzer'] = np.zeros(self.simulation_hours) # array electrolyzer electricity balance
            self.energy_balance['hydrogen']['electrolyzer'] = np.zeros(self.simulation_hours) # array electrolyzer hydrogen balance
            
        if 'fuel cell' in system:
            self.technologies['fuel cell'] = fuel_cell(system['fuel cell'],self.simulation_hours) # Fuel cell object created and to 'technologies' dictionary
            self.energy_balance['electricity']['fuel cell'] = np.zeros(self.simulation_hours) # array fuel cell electricity balance
            self.energy_balance['hydrogen']['fuel cell'] = np.zeros(self.simulation_hours) # array fuel cell hydrogen balance
            
        if 'H tank' in system:
            self.technologies['H tank'] = H_tank(system['H tank'],self.simulation_hours) # H tank object created and to 'technologies' dictionary
            self.energy_balance['hydrogen']['H tank'] = np.zeros(self.simulation_hours) # array H tank hydrogen balance

        self.energy_balance['electricity']['collective self consumption'] = np.zeros(self.simulation_hours) # array contribution to collective-self-consumption as producer (-) or as consumer (+)

    def loc_energy_simulation(self,h,weather):
        """
        Simulate the location
        
        input : int hour to simulate
    
        output : updating of location energy balances
        """
        
        EB = {'electricity': 0, 'heat': 0, 'hydrogen': 0, 'gas': 0} # initialise Energy Balances     
        
        for carrier in EB: # for each energy carrier
            if 'demand' in self.energy_balance[carrier]:                
                EB[carrier] += self.energy_balance[carrier]['demand'][h] # energy balance update: energy demand(-) # n.b cooling demand (+)
                        
        if 'heatpump' in self.technologies:         
            self.energy_balance['electricity']['heatpump'][h], self.energy_balance['heat']['inertialtank'][h], self.energy_balance['heat']['heatpump'][h] = self.technologies['heatpump'].use(weather['temp_air'][h],EB['heat']) 
            EB['electricity'] += self.energy_balance['electricity']['heatpump'][h] # electricity absorbed by heatpump
            EB['heat'] += self.energy_balance['heat']['inertialtank'][h] # heat or cool supplied by inertialtank
            EB['heat'] += self.energy_balance['heat']['heatpump'][h]  # heat or cool supplied by heatpump
            
        if 'boiler_el' in self.technologies: 
            self.energy_balance['electricity']['boiler_el'][h], self.energy_balance['heat']['boiler_el'][h] = self.technologies['boiler_el'].use(EB['heat'],1) # el consumed and heat produced from boiler_el
            EB['electricity'] += self.energy_balance['electricity']['boiler_el'][h] # elecricity balance update: - electricity consumed by boiler_el
            EB['heat'] += self.energy_balance['heat']['boiler_el'][h] # heat balance update: + heat produced by boiler_el

        if 'boiler_ng' in self.technologies: 
            self.energy_balance['gas']['boiler_ng'][h], self.energy_balance['heat']['boiler_ng'][h] = self.technologies['boiler_ng'].use(EB['heat'],1) # ng consumed and heat produced from boiler_ng
            EB['gas'] += self.energy_balance['electricity']['boiler_ng'][h] # gas balance update: - gas consumed by boiler_ng
            EB['heat'] += self.energy_balance['heat']['boiler_ng'][h] # heat balance update: + heat produced by boiler_ng

        if 'PV' in self.technologies: 
            self.energy_balance['electricity']['PV'][h] = self.technologies['PV'].use(h) # electricity produced from PV
            EB['electricity'] += self.energy_balance['electricity']['PV'][h] # elecricity balance update: + electricity produced from PV
            
        if 'heatpump' in self.technologies:
            self.energy_balance['electricity']['heatpump'][h], self.energy_balance['heat']['inertialtank'][h], self.energy_balance['heat']['heatpump'][h] = self.technologies['heatpump'].use(weather['temp_air'][h],EB['heat']) 
            EB['electricity'] += self.energy_balance['electricity']['heatpump'][h] # electricity absorbed by heatpump

        if 'battery' in self.technologies:
            if self.technologies['battery'].collective == 0: 
                self.energy_balance['electricity']['battery'][h] = self.technologies['battery'].use(h,EB['electricity']) # electricity absorbed(-) or supplied(+) by battery
                EB['electricity'] += self.energy_balance['electricity']['battery'][h]  # electricity balance update: +- electricity absorbed or supplied by battery
                
        if 'electrolyzer' in self.technologies:
            if EB['electricity'] > 0:
                if 'H tank' in self.technologies: # if hydrogen is stored in a tank
                    storable_hydrogen = self.technologies['H tank'].max_capacity-self.technologies['H tank'].LOC[h] # the tank can't be full
                else: # if hydrogen is sold to the grid 
                    storable_hydrogen = 99999999999 # there are no limits, f.i an hydrogen producer
                    
                self.energy_balance['hydrogen']['electrolyzer'][h], self.energy_balance['electricity']['electrolyzer'][h] = self.technologies['electrolyzer'].use(h,EB['electricity'],storable_hydrogen) # hydrogen supplied by electrolyzer(+) and electricity absorbed(-) 
                EB['hydrogen'] += self.energy_balance['hydrogen']['electrolyzer'][h]
                EB['electricity'] += self.energy_balance['electricity']['electrolyzer'][h]
                
        if 'fuel cell' in self.technologies:
            if EB['electricity'] < 0:      
                if 'H tank' in self.technologies: # if hydrogen is stored in a tank
                    available_hyd = self.technologies['H tank'].LOC[h] + self.technologies['H tank'].max_capacity - self.technologies['H tank'].used_capacity
                else: # if hydrogen is purchased
                    available_hyd = 99999999999 # there are no limits
                self.energy_balance['hydrogen']['fuel cell'][h], self.energy_balance['electricity']['fuel cell'][h] =self.technologies['fuel cell'].use(h,EB['electricity'],available_hyd) # hydrogen absorbeed by fuel cell(-) and electricity supplied(+) 
                EB['hydrogen'] += self.energy_balance['hydrogen']['fuel cell'][h]
                EB['electricity'] += self.energy_balance['electricity']['fuel cell'][h]
                
        if 'H tank' in self.technologies:
            self.energy_balance['hydrogen']['H tank'][h] = self.technologies['H tank'].use(h,EB['hydrogen'])
            EB['hydrogen'] += self.energy_balance['hydrogen']['H tank'][h]
        
        for carrier in EB:    
            if 'grid' in self.energy_balance[carrier]:
                self.energy_balance[carrier]['grid'][h] = - EB[carrier] # energy from grid(+) or into grid(-)

            
            
            
                
        
        

