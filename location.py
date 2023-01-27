import numpy as np
import pandas as pd
from techs import heatpump,boiler_el,boiler_ng,boiler_h2,PV,wind,battery,H_tank,fuel_cell,electrolyzer,inverter,chp_gt


class location:
    
    def __init__(self,system,general,location_name,path,check,rec_name):
        """
        Create a location object (producer, consumer or prosumer) 
    
        system: dictionary (all the inputs are optional)
            'demand'': dictionary
                'electricity':   str 'file_name.csv' hourly time series of electricity demand 8760 values [kWh]
                'heating water': str 'file_name.csv' hourly time series of heating and dhw demand 8760 values [kWh]
                'cooling water': str 'file_name.csv' hourly time series of cooling demand 8760 values [kWh]
                'hydrogen':      str 'file_name.csv' hourly time series of hydrogen demand 8760 values [kg/h]
                'gas':           str 'file_name.csv' hourly time series of gas demand 8760 values [kWh]
                'process steam': str 'file_name.csv' hourly time series of process steam demand 8760 values [kg/h]
            'PV':           dictionary parameters needed to create PV object (see PV.py)
            'inverter':     dictionary parameters needed to create inverter object (see inverter.py)
            'wind':         dictionary parameters needed to create wind object (see wind.py)
            'battery':      dictionary parameters needed to create battery object (see Battery.py)
            'electrolyzer': dictionary parameters needed to create electrolyzer object (see electrolyzer.py)
            'H tank':       dictionary parameters needed to create H_tank object (see H_tank.py)
            'heatpump':     dictionary parameters needed to create heat pump object (see heatpump.py)
            'boiler_ng':    dictionary parameters needed to create fuel cell object (see boiler.py)
            'boiler_el':    dictionary parameters needed to create fuel cell object (see boiler.py)
            'boiler_h2':    dictionary parameters needed to create fuel cell object (see boiler.py)
            'chp_gt':       dicitonary parameters needed to create a chp object based on gas turbine technoology (see chp_gt.py)
            
            
        general: dictionary 
            see rec.py
            
        output : location object able to:
            simulate the energy flows of present technologies .loc_simulation
            record energy balances .energy_balance (electricity, heating water, cooling water, gas, hydrogen)
        """

        self.system = system
        self.name = location_name
        self.technologies = {} # initialise technologies dictionary
        self.energy_balance = {'electricity': {}, 'heating water': {}, 'cooling water': {}, 'hydrogen': {}, 'gas': {}, 'process steam': {}} # initialise energy balances dictionaries

        self.simulation_hours = int(general['simulation years']*8760) # hourly timestep     
        
        # create the objects of present technologies and add them to the technologies dictionary
        # initialise energy balance and add them to the energy_balance dictionary
        
        for carrier in self.energy_balance:
            if carrier in system['grid'] and system['grid'][carrier]:
                self.energy_balance[carrier]['grid'] = np.zeros(self.simulation_hours) # array energy carrier bought from the grid (-) or feed into the grid (+)

            if carrier in system['demand']:
                if carrier == 'hydrogen' or carrier == 'process steam':
                    self.energy_balance[carrier]['demand'] = - np.tile(pd.read_csv(path+'/loads/'+system['demand'][carrier])['kg'].to_numpy(),int(self.simulation_hours/8760))   # hourly energy carrier needed for the entire simulation
                # if carrier == 'process steam':
                    # self.energy_balance[carrier]['demand'] = - np.tile(pd.read_csv(path+'/loads/'+system['demand'][carrier])['kg'].to_numpy(),int(self.simulation_hours/8760))   # hourly energy carrier needed for the entire simulation
                else:                                                                                                                                                                                                          
                    self.energy_balance[carrier]['demand'] = - np.tile(pd.read_csv(path+'/loads/'+system['demand'][carrier])['kWh'].to_numpy(),int(self.simulation_hours/8760))  # hourly energy carrier needed for the entire simulation

        if 'chp_gt' in system:
            self.technologies['chp_gt'] = chp_gt(system['chp_gt'],self.simulation_hours) # chp_gt object created and added to 'technologies' dictionary
            self.energy_balance['process steam']['chp_gt'] = np.zeros(self.simulation_hours) # array chp_gt process steam balance 
            self.energy_balance['electricity']['chp_gt'] = np.zeros(self.simulation_hours) # array chp_gt electricity balance
            self.energy_balance['hydrogen']['chp_gt'] = np.zeros(self.simulation_hours) # array chp_gt process hydrogen balance 
        
        if 'heatpump' in system:
            self.technologies['heatpump'] = heatpump(system['heatpump'],self.simulation_hours) # heatpump object created and add to 'technologies' dictionary
            self.energy_balance['electricity']['heatpump'] = np.zeros(self.simulation_hours) # array heatpump electricity balance
            self.energy_balance['heating water']['heatpump'] = np.zeros(self.simulation_hours) # array heatpump heat balance
            self.energy_balance['heating water']['inertial TES'] = np.zeros(self.simulation_hours) # array inertial tank heat balance
                
        if 'boiler_el' in self.system:
            self.technologies['boiler_el'] = boiler_el(self.system['boiler_el']) # boiler_el object created and add to 'technologies' dictionary
            self.energy_balance['electricity']['boiler_el'] = np.zeros(self.simulation_hours) # array boiler_el electricity balance
            self.energy_balance['heating water']['boiler_el'] = np.zeros(self.simulation_hours) # array boiler_el heat balance
               
        if 'boiler_ng' in self.system:
            self.technologies['boiler_ng'] = boiler_ng(self.system['boiler_ng'])  # boiler_ng object created and add to 'technologies' dictionary
            self.energy_balance['gas']['boiler_ng'] = np.zeros(self.simulation_hours) # array boiler_ng gas balance
            self.energy_balance['heating water']['boiler_ng'] = np.zeros(self.simulation_hours) # array boiler_ng heat balance 

        if 'PV' in self.system:
            self.technologies['PV'] = PV(self.system['PV'],general,self.simulation_hours,self.name,path,check,rec_name) # PV object created and add to 'technologies' dictionary
            self.energy_balance['electricity']['PV'] = np.zeros(self.simulation_hours) # array PV electricity balance
           
        if 'inverter' in self.system:
            self.technologies['inverter'] = inverter(self.system['inverter'],self.simulation_hours) # inverter object created and add to 'technologies' dictionary
            self.energy_balance['electricity']['inverter'] = np.zeros(self.simulation_hours) # array inverter electricity balance
            
        if 'wind' in self.system:
            self.technologies['wind'] = wind(self.system['wind']) # wind object created and add to 'technologies' dictionary
            self.energy_balance['electricity']['wind'] = np.zeros(self.simulation_hours) # array wind electricity balance 
           
        if 'battery' in self.system:
            self.technologies['battery'] = battery(self.system['battery'],self.simulation_hours) # battery object created and to 'technologies' dictionary
            self.energy_balance['electricity']['battery'] = np.zeros(self.simulation_hours) # array battery electricity balance
                           
        if 'electrolyzer' in self.system:
            self.technologies['electrolyzer'] = electrolyzer(self.system['electrolyzer'],self.simulation_hours) # electrolyzer object created and to 'technologies' dictionary
            self.energy_balance['electricity']['electrolyzer'] = np.zeros(self.simulation_hours) # array electrolyzer electricity balance
            self.energy_balance['hydrogen']['electrolyzer'] = np.zeros(self.simulation_hours) # array electrolyzer hydrogen balance
            
        if 'fuel cell' in self.system:
            self.technologies['fuel cell'] = fuel_cell(self.system['fuel cell'],self.simulation_hours) # Fuel cell object created and to 'technologies' dictionary
            self.energy_balance['electricity']['fuel cell'] = np.zeros(self.simulation_hours)     # array fuel cell electricity balance
            self.energy_balance['hydrogen']['fuel cell'] = np.zeros(self.simulation_hours)        # array fuel cell hydrogen balance
            self.energy_balance['heating water']['fuel cell']=np.zeros(self.simulation_hours)              #array fuel cell heat balance used
            
        if 'H tank' in self.system:
            self.technologies['H tank'] = H_tank(self.system['H tank'],self.simulation_hours) # H tank object created and to 'technologies' dictionary
            self.energy_balance['hydrogen']['H tank'] = np.zeros(self.simulation_hours)  # array H tank hydrogen balance

        if 'boiler_h2' in self.system:
            self.technologies['boiler_h2'] = boiler_h2(self.system['boiler_h2'])                 # boiler_h2 object created and added to 'technologies' dictionary
            self.energy_balance['hydrogen']['boiler_h2'] = np.zeros(self.simulation_hours)  # array boiler_h2 gas balance
            self.energy_balance['heating water']['boiler_h2'] = np.zeros(self.simulation_hours)      # array boiler_h2 heat balance 

        self.energy_balance['electricity']['collective self consumption'] = np.zeros(self.simulation_hours) # array contribution to collective-self-consumption as producer (-) or as consumer (+)
        self.energy_balance['heating water']['collective self consumption'] = np.zeros(self.simulation_hours) # array contribution to collective-self-consumption as producer (-) or as consumer (+)---heat----mio!!!
        self.energy_balance['process steam']['collective self consumption'] = np.zeros(self.simulation_hours) # array contribution to collective-self-consumption as producer (-) or as consumer (+)---heat----mio!!!

    def loc_energy_simulation(self,h,weather):
        """
        Simulate the location
        
        input : int hour to simulate
    
        output : updating of location energy balances
        """
        
        eb = {'electricity': 0, 'heating water': 0, 'cooling water': 0, 'hydrogen': 0, 'gas': 0, 'process steam': 0} # initialise Energy Balances     
        
        
        ### production from renewables (pv, wind, inverter)
        if 'PV' in self.technologies: 
            self.energy_balance['electricity']['PV'][h] = self.technologies['PV'].use(h) # electricity produced from PV
            eb['electricity'] += self.energy_balance['electricity']['PV'][h] # elecricity balance update: + electricity produced from PV

        if 'wind' in self.technologies: 
            self.energy_balance['electricity']['wind'][h] = self.technologies['wind'].use(h,weather['wind_speed'][h]) # electricity produced from wind
            eb['electricity'] += self.energy_balance['electricity']['wind'][h] # elecricity balance update: + electricity produced from wind
    
        if 'inverter' in self.technologies:
            self.energy_balance['electricity']['inverter'][h] = self.technologies['inverter'].use(h,eb['electricity']) # electricity lost in conversion by the inverter
            eb['electricity'] += self.energy_balance['electricity']['inverter'][h] # electricity balance update: - electricity lost in conversion by the invertert
     
        ### demand
        for carrier in eb: # for each energy carrier
            if 'demand' in self.energy_balance[carrier]:                
                eb[carrier] += self.energy_balance[carrier]['demand'][h] # energy balance update: energy demand(-) # n.b cooling demand (+)              
       
        ### production from cogeneration system (chp_gt)
        if 'chp_gt' in self.technologies:
            if 'H tank' in self.technologies:   # if hydrogen is stored in a tank
                available_hyd = self.technologies['H tank'].LOC[h] + self.technologies['H tank'].max_capacity - self.technologies['H tank'].used_capacity                                                                  
            else:                               # if hydrogen is sold to the grid 
                available_hyd = 99999999999 # there are no limits, f.i an hydrogen producer

            use = self.technologies['chp_gt'].use(h,weather['temp_air'][h],eb['process steam'])     # saving fuel cell working parameters for the current timeframe
            self.energy_balance['process steam']['chp_gt'][h] = use[0]   # produced steam (+)
            self.energy_balance['electricity']['chp_gt'][h] =   use[1]   # produced electricity (+)
            self.energy_balance['hydrogen']['chp_gt'][h] =      use[2]   # hydrogen required by chp system to run (-)  

            eb['hydrogen'] += self.energy_balance['hydrogen']['chp_gt'][h]            
            eb['process steam'] += self.energy_balance['process steam']['chp_gt'][h]    
            eb['electricity'] += self.energy_balance['electricity']['chp_gt'][h]

        ### other techologies (boiler_el, boiler_ng, heatpump, battery, electrolyzer, fuel cell, boiler_H2, H tank)
                    
        if 'boiler_el' in self.technologies: 
            self.energy_balance['electricity']['boiler_el'][h], self.energy_balance['heating water']['boiler_el'][h] = self.technologies['boiler_el'].use(eb['heating water'],1) # el consumed and heat produced from boiler_el
            eb['electricity'] += self.energy_balance['electricity']['boiler_el'][h] # elecricity balance update: - electricity consumed by boiler_el
            eb['heating water'] += self.energy_balance['heating water']['boiler_el'][h] # heat balance update: + heat produced by boiler_el

        if 'boiler_ng' in self.technologies: 
            self.energy_balance['gas']['boiler_ng'][h], self.energy_balance['heating water']['boiler_ng'][h] = self.technologies['boiler_ng'].use(eb['heating water'],1) # ng consumed and heat produced from boiler_ng
            eb['gas'] += self.energy_balance['gas']['boiler_ng'][h] # gas balance update: - gas consumed by boiler_ng
            eb['heating water'] += self.energy_balance['heating water']['boiler_ng'][h] # heat balance update: + heat produced by boiler_ng
    
        if 'heatpump' in self.technologies:     
            self.energy_balance['electricity']['heatpump'][h], self.energy_balance['heating water']['heatpump'][h], self.energy_balance['heating water']['inertial TES'][h] = self.technologies['heatpump'].use(weather['temp_air'][h],eb['heating water'],eb['electricity'],h) 
         
            eb['electricity'] += self.energy_balance['electricity']['heatpump'][h] # electricity absorbed by heatpump
            self.energy_balance['electricity']['demand'][h] += self.energy_balance['electricity']['heatpump'][h] # add heatpump demand to 'electricity demand'
            eb['heating water'] += self.energy_balance['heating water']['inertial TES'][h] + self.energy_balance['electricity']['heatpump'][h] # heat or cool supplied by HP or inertial TES
    
        if 'battery' in self.technologies:
            if self.technologies['battery'].collective == 0: 
                self.energy_balance['electricity']['battery'][h] = self.technologies['battery'].use(h,eb['electricity']) # electricity absorbed(-) or supplied(+) by battery
                eb['electricity'] += self.energy_balance['electricity']['battery'][h]  # electricity balance update: +- electricity absorbed or supplied by battery
       
        if 'electrolyzer' in self.technologies:
            if eb['electricity'] > 0:
                if 'H tank' in self.technologies:   # if hydrogen is stored in a tank
                    storable_hydrogen = self.technologies['H tank'].max_capacity-self.technologies['H tank'].LOC[h] # the tank can't be full
                    if storable_hydrogen>self.technologies['H tank'].max_capacity*0.00001:
                        self.energy_balance['hydrogen']['electrolyzer'][h], self.energy_balance['electricity']['electrolyzer'][h] = self.technologies['electrolyzer'].use(h,eb['electricity'],storable_hydrogen) # hydrogen supplied by electrolyzer(+) and electricity absorbed(-) 
                        eb['hydrogen']    += self.energy_balance['hydrogen']['electrolyzer'][h]
                        eb['electricity'] += self.energy_balance['electricity']['electrolyzer'][h]                                                                  
                else:                               # if hydrogen is sold to the grid 
                    storable_hydrogen = 99999999999 # there are no limits, f.i an hydrogen producer
                    
                    self.energy_balance['hydrogen']['electrolyzer'][h],self.energy_balance['electricity']['electrolyzer'][h] = self.technologies['electrolyzer'].use(h,eb['electricity'],storable_hydrogen)      # [:2] # hydrogen supplied by electrolyzer(+) # electricity absorbed by the electorlyzer(-) 
                  
                    eb['hydrogen']    += self.energy_balance['hydrogen']['electrolyzer'][h]
                    eb['electricity'] += self.energy_balance['electricity']['electrolyzer'][h]
        
                
        if 'fuel cell' in self.technologies:
            if eb['electricity'] < 0:      
                if 'H tank' in self.technologies: # if hydrogen is stored in a tank
                    available_hyd = self.technologies['H tank'].LOC[h] + self.technologies['H tank'].max_capacity - self.technologies['H tank'].used_capacity
                else:                             # if hydrogen is purchased
                    available_hyd = 99999999999   # there are no limits

                use = self.technologies['fuel cell'].use(h,eb['electricity'],available_hyd)     # saving fuel cell working parameters for the current timeframe
                self.energy_balance['hydrogen']['fuel cell'][h] =    use[0] # hydrogen absorbed by fuel cell(-)
                self.energy_balance['electricity']['fuel cell'][h] = use[1] # electricity supplied(+) 

                if use[2] < -eb['heating water']: #all of the heat producted by FC is used      
                    self.energy_balance['heating water']['fuel cell'][h]=use[2] # heat loss of fuel cell
                else:
                    self.energy_balance['heating water']['fuel cell'][h]=-eb['heating water'] # heat loss of fuel cell- demand
                    # self.energy_balance['heating water']['fuel cell']['unused'][h]=self.technologies['fuel cell'].use(h,eb['electricity'],available_hyd)[2]+eb['heating water']
                if 'grid' in self.energy_balance['heating water']:
                    if use[2] > -eb['heating water']:      
                        self.energy_balance['heating water']['grid'][h]=use[2]+eb['heating water']
                    else:
                        self.energy_balance['heating water']['grid'][h]=0
          
                eb['hydrogen'] += self.energy_balance['hydrogen']['fuel cell'][h]
                eb['electricity'] += self.energy_balance['electricity']['fuel cell'][h]
                eb['heating water'] += self.energy_balance['heating water']['fuel cell'][h] 
                
        if 'boiler_h2' in self.technologies: 
            if 'H tank' in self.technologies: # if hydrogen is stored in a tank
                available_hyd = self.technologies['H tank'].LOC[h] + self.technologies['H tank'].max_capacity - self.technologies['H tank'].used_capacity-self.energy_balance['hydrogen']['fuel cell'][h]
            else: # if hydrogen is purchased
                available_hyd = 99999999999 # there are no limits
            self.energy_balance['hydrogen']['boiler_h2'][h], self.energy_balance['heating water']['boiler_h2'][h] = self.technologies['boiler_h2'].use(eb['heating water'],available_hyd,1)[1:3] #h2 consumed from boiler_h2 and heat produced by boiler_h2
            eb['hydrogen'] += self.energy_balance['hydrogen']['boiler_h2'][h] # hydrogen balance update: - hydrogen consumed by boiler_h2
            eb['heating water'] += self.energy_balance['heating water']['boiler_h2'][h] # heat balance update: + heat produced by boiler_h2
                
        if 'H tank' in self.technologies:
            self.energy_balance['hydrogen']['H tank'][h] = self.technologies['H tank'].use(h,eb['hydrogen'])
            eb['hydrogen'] += self.energy_balance['hydrogen']['H tank'][h]
        
        for carrier in eb:    
            if 'grid' in self.energy_balance[carrier]:
                if carrier=='heating water':
                    continue
                else:
                    self.energy_balance[carrier]['grid'][h] = - eb[carrier] # energy from grid(+) or into grid(-) 
                                                                  
                
        
        

