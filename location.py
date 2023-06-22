import numpy as np
import pandas as pd
from techs import (
                    heatpump,
                    boiler_el,
                    boiler_ng,
                    boiler_h2,
                    PV,
                    wind,
                    battery,
                    H_tank,
                    HPH_tank,
                    O2_tank,
                    fuel_cell,
                    electrolyzer,
                    inverter,
                    chp_gt,
                    Chp,
                    Absorber,
                    hydrogen_compressor,
                    Compressor
                  )

class location:
    
    def __init__(self,system,general,location_name,path,check,file_structure,file_general):
        """
        Create a location object (producer, consumer or prosumer) 
    
        system: dictionary (all the inputs are optional)
            'demand'': dictionary
                'electricity':              str 'file_name.csv' hourly time series of electricity demand 8760 values [kWh]
                'heating water':            str 'file_name.csv' hourly time series of heating and dhw demand 8760 values [kWh]
                'cooling water':            str 'file_name.csv' hourly time series of cooling demand 8760 values [kWh]
                'hydrogen':                 str 'file_name.csv' hourly time series of hydrogen demand 8760 values [kg/h]
                'HP hydrogen':              str 'file_name.csv' hourly time series of High-Pressure hydrogen demand 8760 values [kg/h]
                'gas':                      str 'file_name.csv' hourly time series of gas demand 8760 values [kWh]
                'process steam':            str 'file_name.csv' hourly time series of process steam demand 8760 values [kg/h]
                'process heat':             str 'file_name.csv' hourly time series of process heat demand 8760 values [kWh]
                'process hot water':        str 'file_name.csv' hourly time series of process hot water demand 8760 values [kWh]
                'process cold water':       str 'file_name.csv' hourly time series of process cold water demand (absorber, 7-12 °C) 8760 values [kWh]
                'process chilled water':    str 'file_name.csv' hourly time series of process chilled water demand (absorber, 1-5 °C) 8760 values [kWh]
            'PV':                       dictionary parameters needed to create PV object (see PV.py)
            'inverter':                 dictionary parameters needed to create inverter object (see inverter.py)
            'wind':                     dictionary parameters needed to create wind object (see wind.py)
            'battery':                  dictionary parameters needed to create battery object (see Battery.py)
            'electrolyzer':             dictionary parameters needed to create electrolyzer object (see electrolyzer.py)
            'H tank':                   dictionary parameters needed to create H_tank object (see H_tank.py)
            'O2 tank':                  dictionary parameters needed to create O2_tank object (see O2_tank.py)
            'HPH tank':                 dictionary parameters needed to create High Pressure H_tank object (see H_tank.py)
            'heatpump':                 dictionary parameters needed to create heat pump object (see heatpump.py)
            'boiler_ng':                dictionary parameters needed to create fuel cell object (see boiler.py)
            'boiler_el':                dictionary parameters needed to create fuel cell object (see boiler.py)
            'boiler_h2':                dictionary parameters needed to create fuel cell object (see boiler.py)
            'chp_gt':                   dicitonary parameters needed to create a chp object based on gas turbine technoology (see chp_gt.py)
            'hydrogen_compressor':      dicitonary parameters needed to create a mhhc object (see mhhc compressor.py)
            'mechanical compressor':    dicitonary parameters needed to create a mechanical object (see compressor.py)
        general: dictionary 
            see rec.py
            
        output : location object able to:
            simulate the energy flows of present technologies .loc_simulation
            record energy balances .energy_balance (electricity, heating water, cooling water, gas, hydrogen)
        """
        
        self.system = dict(sorted(system.items(), key=lambda item: item[1]['priority'])) # ordered by priority
        self.name = location_name
        self.technologies = {}                                  # initialise technologies dictionary
        self.energy_balance = {'electricity'            : {},   # initialise energy balances dictionaries
                               'heating water'          : {}, 
                               'cooling water'          : {}, 
                               'hydrogen'               : {},
                               'LP hydrogen'            : {},
                               'HP hydrogen'            : {},
                               'water'                  : {},
                               'oxygen'                 : {},
                               'gas'                    : {}, 
                               'process steam'          : {},
                               'process heat'           : {},
                               'process hot water'      : {},
                               'process cold water'     : {},
                               'process chilled water'  : {}} 

        self.simulation_hours = int(general['simulation years']*8760) # hourly timestep     
        
        # create the objects of present technologies and add them to the technologies dictionary
        # initialise energy balance and add them to the energy_balance dictionary
        
        for carrier in self.energy_balance:
            
            if f"{carrier} grid" in self.system:
                self.energy_balance[carrier]['grid'] = np.zeros(self.simulation_hours) # array energy carrier bought from the grid (-) or feed into the grid (+)

            if f"{carrier} demand" in self.system:
                if carrier in ['hydrogen','HP hydrogen','process steam']:
                    self.energy_balance[carrier]['demand'] = - np.tile(pd.read_csv(path+'/loads/'+system[f"{carrier} demand"]['serie'])['kg'].to_numpy(),int(self.simulation_hours/8760))   # hourly energy carrier needed for the entire simulation
                    if carrier == 'hydrogen' or carrier == 'HP hydrogen':
                        self.hydrogen_demand = carrier  # demand can be defined as 'hydrogen demand' or 'HP hydrogen demand' depending on the required delivery pressure
                else:                                                                                                                                                                                                          
                    self.energy_balance[carrier]['demand'] = - np.tile(pd.read_csv(path+'/loads/'+system[f"{carrier} demand"]['serie'])['kWh'].to_numpy(),int(self.simulation_hours/8760))  # hourly energy carrier needed for the entire simulation

        if 'chp_gt' in self.system:
            self.technologies['chp_gt'] = chp_gt(system['chp_gt'],self.simulation_hours) # chp_gt object created and added to 'technologies' dictionary
            self.energy_balance['process steam']['chp_gt'] = np.zeros(self.simulation_hours) # array chp_gt process steam balance 
            self.energy_balance['electricity']['chp_gt'] = np.zeros(self.simulation_hours) # array chp_gt electricity balance
            self.energy_balance['hydrogen']['chp_gt'] = np.zeros(self.simulation_hours) # array chp_gt process hydrogen balance 
        
        if 'chp' in self.system:
            self.technologies['chp'] = Chp(system['chp'],self.simulation_hours) # chp object created and added to 'technologies' dictionary
            self.energy_balance[self.technologies['chp'].th_out]['chp'] = np.zeros(self.simulation_hours) # array chp thermal output balance (process steam/hot water)
            self.energy_balance['electricity']['chp'] = np.zeros(self.simulation_hours) # array chp electricity balance
            self.energy_balance[self.technologies['chp'].fuel]['chp'] = np.zeros(self.simulation_hours) # array chp fuel consumption balance
            self.energy_balance['process heat']['chp'] = np.zeros(self.simulation_hours) # array chp process heat balance
            self.energy_balance['process hot water']['chp'] = np.zeros(self.simulation_hours) # array chp process hot water balance
            self.energy_balance['process cold water']['chp'] = np.zeros(self.simulation_hours) # array chp process cold water balance
       
        if 'absorber' in self.system:
            self.technologies['absorber'] = Absorber(system['absorber'],self.simulation_hours) # absorber object created and added to 'technologies' dictionary
            self.energy_balance['process heat']['absorber'] = np.zeros(self.simulation_hours) # array absorber process steam balance 
            self.energy_balance['process hot water']['absorber'] = np.zeros(self.simulation_hours) # array absorber process steam balance 
            self.energy_balance['process cold water']['absorber'] = np.zeros(self.simulation_hours) # array absorber process steam balance 
            
        if 'heatpump' in self.system:
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
            self.technologies['PV'] = PV(self.system['PV'],general,self.simulation_hours,self.name,path,check,file_structure,file_general) # PV object created and add to 'technologies' dictionary
            self.energy_balance['electricity']['PV'] = np.zeros(self.simulation_hours) # array PV electricity balance
           
        if 'inverter' in self.system:
            self.technologies['inverter'] = inverter(self.system['inverter'],self.simulation_hours) # inverter object created and add to 'technologies' dictionary
            self.energy_balance['electricity']['inverter'] = np.zeros(self.simulation_hours)        # array inverter electricity balance
            
        if 'wind' in self.system:
            self.technologies['wind'] = wind(self.system['wind'],self.simulation_hours,path)    # wind object created and add to 'technologies' dictionary
            self.energy_balance['electricity']['wind'] = np.zeros(self.simulation_hours)        # array wind electricity balance 
           
        if 'battery' in self.system:
            self.technologies['battery'] = battery(self.system['battery'],self.simulation_hours)    # battery object created and to 'technologies' dictionary
            self.energy_balance['electricity']['battery'] = np.zeros(self.simulation_hours)         # array battery electricity balance
                           
        if 'electrolyzer' in self.system:
            self.technologies['electrolyzer'] = electrolyzer(self.system['electrolyzer'],self.simulation_hours) # electrolyzer object created and to 'technologies' dictionary
            self.energy_balance['electricity']['electrolyzer']              = np.zeros(self.simulation_hours) # array electrolyzer electricity balance
            self.energy_balance['oxygen']['electrolyzer']                   = np.zeros(self.simulation_hours) # array electrolyzer oxygen balance
            self.energy_balance['water']['electrolyzer']                    = np.zeros(self.simulation_hours) # array electrolyzer water balance
            self.energy_balance['hydrogen']['electrolyzer']                 = np.zeros(self.simulation_hours) # array electrolyzer hydrogen balance
            if self.technologies['electrolyzer'].strategy == "full-time" and not self.system["electricity grid"]["draw"]:
                raise ValueError(f"Full-time electrolyzers operation considered without electricity grid connection in {self.name} location.\n\
                Options to fix the problem: \n\
                    (a) - Insert electricity grid withdrawal in studycase.json\n\
                    (b) - Change electrolyzers strategy in studycase.json")
    
        if 'fuel cell' in self.system:
            self.technologies['fuel cell'] = fuel_cell(self.system['fuel cell'],self.simulation_hours) # Fuel cell object created and to 'technologies' dictionary
            self.energy_balance['electricity']['fuel cell']     = np.zeros(self.simulation_hours)     # array fuel cell electricity balance
            self.energy_balance['hydrogen']['fuel cell']        = np.zeros(self.simulation_hours)     # array fuel cell hydrogen balance
            self.energy_balance['heating water']['fuel cell']   = np.zeros(self.simulation_hours)     # array fuel cell heat balance used
            
        if 'hydrogen compressor' in self.system:
            self.technologies['hydrogen compressor'] = hydrogen_compressor(self.system['hydrogen compressor'],self.simulation_hours) # MHHC compressor object created and to 'technologies' dictionary
            self.energy_balance['hydrogen']['hydrogen compressor']  = np.zeros(self.simulation_hours)     # array hydrogen compressor hydrogen compressed
            self.energy_balance['gas']['hydrogen compressor']       = np.zeros(self.simulation_hours)     # array hydrogen compressor heating water balanced used
        
        if 'H tank' and 'HPH tank' in self.system: 
            self.technologies['H tank'] = H_tank(self.system['H tank'],self.simulation_hours) # H tank object created and to 'technologies' dictionary
            self.energy_balance['hydrogen']['H tank'] = np.zeros(self.simulation_hours)         # array H tank hydrogen balance
            self.technologies['HPH tank'] = HPH_tank(self.system['HPH tank'],self.simulation_hours)   # HPH tank object created and to 'technologies' dictionary
            self.energy_balance['HP hydrogen']['HPH tank'] = np.zeros(self.simulation_hours)        # array HPH tank hydrogen balance
            
            self.tank_stream = {'H tank':'hydrogen',        # dictionary assigning different hydrogen streams to different storage technologies - necessary for loc_energy_simulation
                                'HPH tank':'HP hydrogen'}
        
        if 'O2 tank' in self.system: 
            self.technologies['O2 tank'] = O2_tank(self.system['O2 tank'],self.simulation_hours) # LPH tank object created and to 'technologies' dictionary
            self.energy_balance['oxygen']['O2 tank'] = np.zeros(self.simulation_hours)         # array LPH tank hydrogen balance
        
            
        if 'mechanical compressor' in self.system:
            if "electrolyzer" not in self.system:
                raise ValueError(f"Electorlyzer not present in the {self.name} location. The model as it is considers hydrogen compression\n\
                                    only when hydrogen is produced in situ. Compressor is directly connected either to \n\
                                    electrolyzer or a buffer tank, check priorities in studycase.json.\n\
                                    Options to fix the problem: \n\
                    (a) - Insert electrolyzer technology in studycase.json\n")
            if 'O2 tank' not in self.system:
                maxflowrate = self.technologies['electrolyzer'].maxh2prod_stack
            else:
                maxflowrate =   self.technologies['electrolyzer'].maxh2prod_stack +\
                                self.technologies['electrolyzer'].maxh2prod_stack*self.technologies['electrolyzer'].oxy  # maximum hydrogen mass flow rate + maximum oxygen mass flow rate
            self.technologies['mechanical compressor'] = Compressor(self.system['mechanical compressor'],self.simulation_hours,maxflowrate = maxflowrate) # H tank object created and to 'technologies' dictionary
            self.energy_balance['electricity']['mechanical compressor']    = np.zeros(self.simulation_hours) # array H tank hydrogen balance
            self.energy_balance['hydrogen']['mechanical compressor']       = np.zeros(self.simulation_hours) # array of hydrogen flow entering the mechanical compressor from LPH tank
            self.energy_balance['HP hydrogen']['mechanical compressor']    = np.zeros(self.simulation_hours) # array of compressed hydrogen flow sent toward HPH tank
            self.energy_balance['cooling water']['mechanical compressor']  = np.zeros(self.simulation_hours) # array of water flow to be fed to the refrigeration system 
 
        if 'H tank' in self.system and not 'HPH tank' in self.system:
            if 'hydrogen demand' in self.system or 'HP hydrogen demand' in self.system:
                if self.system[self.hydrogen_demand+' demand']['strategy'] == 'supply-led' and self.system['H tank']['max capacity'] != False:
                    raise ValueError(f"Adjust {self.name} location system in studycase.json. When the system is operated in supply-led mode, H tank size cannot be defined in advance.\n\
            Options to fix the problem: \n\
            (a) - Insert false for 'max capacity' among H tank parameters in studycase.json\n\
            (b) - Switch to 'demand-led' in 'hydrogen-demand'('strategy')\
            ")
            self.technologies['H tank'] = H_tank(self.system['H tank'],self.simulation_hours)   # H tank object created and to 'technologies' dictionary
            self.energy_balance['hydrogen']['H tank'] = np.zeros(self.simulation_hours)         # array H tank hydrogen balance
            
            self.tank_stream = {'H tank':'hydrogen'}     # dictionary assigning hydrogen stream to H tank storage technologies - necessary for loc_energy_simulation
        
        if 'boiler_h2' in self.system:
            self.technologies['boiler_h2'] = boiler_h2(self.system['boiler_h2'])                    # boiler_h2 object created and added to 'technologies' dictionary
            self.energy_balance['hydrogen']['boiler_h2'] = np.zeros(self.simulation_hours)          # array boiler_h2 gas balance
            self.energy_balance['heating water']['boiler_h2'] = np.zeros(self.simulation_hours)     # array boiler_h2 heat balance 

        self.energy_balance['electricity']['collective self consumption']   = np.zeros(self.simulation_hours) # array contribution to collective-self-consumption as producer (-) or as consumer (+)
        self.energy_balance['heating water']['collective self consumption'] = np.zeros(self.simulation_hours) # array contribution to collective-self-consumption as producer (-) or as consumer (+)---heat----mio!!!
        self.energy_balance['process steam']['collective self consumption'] = np.zeros(self.simulation_hours) # array contribution to collective-self-consumption as producer (-) or as consumer (+)---heat----mio!!!
             
    def loc_energy_simulation(self,h,weather):
        """
        Simulate the location
        
        input : int hour to simulate
    
        output : updating of location energy balances
        """
        
        eb = {}
        for carrier in self.energy_balance:
            eb[carrier] = 0 # initialise hourly Energy Balances     
            
        for tech_name in self.system: # (which is ordered py priority)
        
            if tech_name == 'PV': 
                self.energy_balance['electricity']['PV'][h] = self.technologies['PV'].use(h) # electricity produced from PV
                eb['electricity'] += self.energy_balance['electricity']['PV'][h] # elecricity balance update: + electricity produced from PV
    
            if tech_name == 'wind': 
                self.energy_balance['electricity']['wind'][h] = self.technologies['wind'].use(h,weather['wind_speed'][h]) # electricity produced from wind
                eb['electricity'] += self.energy_balance['electricity']['wind'][h] # elecricity balance update: + electricity produced from wind
                        
            if tech_name == 'boiler_el': 
                self.energy_balance['electricity']['boiler_el'][h], self.energy_balance['heating water']['boiler_el'][h] = self.technologies['boiler_el'].use(eb['heating water'],1) # el consumed and heat produced from boiler_el
                eb['electricity'] += self.energy_balance['electricity']['boiler_el'][h] # elecricity balance update: - electricity consumed by boiler_el
                eb['heating water'] += self.energy_balance['heating water']['boiler_el'][h] # heat balance update: + heat produced by boiler_el
    
            if tech_name == 'boiler_ng': 
                self.energy_balance['gas']['boiler_ng'][h], self.energy_balance['heating water']['boiler_ng'][h] = self.technologies['boiler_ng'].use(eb['heating water'],1) # ng consumed and heat produced from boiler_ng
                eb['gas'] += self.energy_balance['gas']['boiler_ng'][h] # gas balance update: - gas consumed by boiler_ng
                eb['heating water'] += self.energy_balance['heating water']['boiler_ng'][h] # heat balance update: + heat produced by boiler_ng
        
            if tech_name == 'heatpump':     
                self.energy_balance['electricity']['heatpump'][h], self.energy_balance['heating water']['heatpump'][h], self.energy_balance['heating water']['inertial TES'][h] = self.technologies['heatpump'].use(weather['temp_air'][h],eb['heating water'],eb['electricity'],h) 
             
                eb['electricity'] += self.energy_balance['electricity']['heatpump'][h] # electricity absorbed by heatpump
                self.energy_balance['electricity']['demand'][h] += self.energy_balance['electricity']['heatpump'][h] # add heatpump demand to 'electricity demand'
                eb['heating water'] += self.energy_balance['heating water']['inertial TES'][h] + self.energy_balance['electricity']['heatpump'][h] # heat or cool supplied by HP or inertial TES
        
            if tech_name == 'battery':
                if self.technologies['battery'].collective == 0: 
                    self.energy_balance['electricity']['battery'][h] = self.technologies['battery'].use(h,eb['electricity']) # electricity absorbed(-) or supplied(+) by battery
                    eb['electricity'] += self.energy_balance['electricity']['battery'][h]  # electricity balance update: +- electricity absorbed or supplied by battery
           
            if tech_name == 'chp_gt':
                if "hydrogen grid" in self.system and self.system["hydrogen grid"]["draw"]: # hydrogen can be withdranw from an hydrogen grid
                    available_hyd = 9999999999999999999 
                elif 'H tank' in self.system:   # only hydrogen inside H tank can be used
                    available_hyd = self.technologies['H tank'].LOC[h] + self.technologies['H tank'].max_capacity - self.technologies['H tank'].used_capacity                                                                  
                else:
                    available_hyd = max(0,eb['hydrogen']) # hydrogen is produced in the same timestep by an electrolyzer with a higher priority than chp
                if available_hyd > 0:
                    use = self.technologies['chp_gt'].use(h,weather['temp_air'][h],eb['process steam'],available_hyd)     # saving chp_gt working parameters for the current timeframe
                    self.energy_balance['process steam']['chp_gt'][h] = use[0]   # produced steam (+)
                    self.energy_balance['electricity']['chp_gt'][h] =   use[1]   # produced electricity (+)
                    self.energy_balance['hydrogen']['chp_gt'][h] =      use[2]   # hydrogen required by chp system to run (-)  
    
                eb['hydrogen'] += self.energy_balance['hydrogen']['chp_gt'][h]            
                eb['process steam'] += self.energy_balance['process steam']['chp_gt'][h]    
                eb['electricity'] += self.energy_balance['electricity']['chp_gt'][h] 
           
            if tech_name == 'chp':
                strategy    = self.technologies['chp'].strategy     # thermal-load follow or electric-load follow 
                coproduct   = self.technologies['chp'].coproduct    # process co-product depending on the  approache chosen above 
                if self.system['chp']['Fuel'] == 'hydrogen':
                    if "hydrogen grid" in self.system and self.system["hydrogen grid"]["draw"]: # hydrogen can be withdranw from an hydrogen grid
                        available_hyd = 9999999999999999999 
                    elif 'H tank' in self.system:   # only hydrogen inside H tank can be used
                        available_hyd = self.technologies['H tank'].LOC[h] + self.technologies['H tank'].max_capacity - self.technologies['H tank'].used_capacity                                                                  
                    else:
                        available_hyd = max(0,eb['hydrogen']) # hydrogen is produced in the same timestep by an electrolyzer with a higher priority than chp
                    use = self.technologies['chp'].use(h,weather['temp_air'][h],eb[strategy], eb[coproduct], available_hyd)     # saving chp working parameters for the current timeframe
                else:
                    use = self.technologies['chp'].use(h,weather['temp_air'][h],eb[strategy], eb[coproduct])                    # saving chp working parameters for the current timeframe
                
                self.energy_balance[self.technologies['chp'].th_out]['chp'][h]  = use[0]   # produced thermal output (+) (steam/hot water)
                self.energy_balance['electricity']['chp'][h]                    = use[1]   # produced electricity (+)
                self.energy_balance[self.technologies['chp'].fuel]['chp'][h]    = use[2]   # fuel required by chp system to run (-)  
                self.energy_balance['process heat']['chp'][h]                   = use[3]   # process heat produced by chp system (+)  
                self.energy_balance['process hot water']['chp'][h]              = use[4]   # process heat produced by chp system (+)  

                eb[self.technologies['chp'].fuel] += self.energy_balance[self.technologies['chp'].fuel]['chp'][h]            
                eb[self.technologies['chp'].th_out] += self.energy_balance[self.technologies['chp'].th_out]['chp'][h]    
                eb['electricity'] += self.energy_balance['electricity']['chp'][h] 
                eb['process heat'] += self.energy_balance['process heat']['chp'][h] 
                
# =============================================================================
#                 # if q_th != q_th_chp :     # WIP:  logica in cui si tiene conto di avere priorità per la domanda termica dell'assorbitore               \
#                                             #       si interroga l'assorbitore con la domanda frigorifera dell'utenza, tramite funzione inversa          \
#                                             #       si risale al calore che sarebbe necessario avere in input all'assorbitore (q_in =q_frigo/COP)        \
#                                             #       e si somma alla domanda termica con cui si va ad interrogare il chp in questo if tech_name == 'chp'  \
#                                             #       A questo punto si verifica quanto calore è in grado di generare il CHP e si gestiscono i flussi      \     
#                                             #       qui con comandi specifici o con una funzione chp.trigeneration() che suggeriva Marco. 
# =============================================================================

            if tech_name == 'absorber':  
                # self.energy_balance['process cold water']['absorber'][h] = self.technologies['absorber'].use(h, eb[self.technologies['chp'].th_out])  # cold energy produced via the absorption cycle (+) - here accounting for\
                self.energy_balance['process cold water']['absorber'][h] = self.technologies['absorber'].use(h, eb['process heat'])  # cold energy produced via the absorption cycle (+)
            
            if tech_name == 'electrolyzer':
                
                if self.technologies['electrolyzer'].strategy == 'hydrogen-first': # electrolyzer activated when renewable energy is available
                    if eb['electricity'] > 0: # electrolyzer activated only when renewable energy is available
                        if "hydrogen grid" in self.system and self.system["hydrogen grid"]["feed"]: # hydrogen can be fed into an hydrogen grid
                            producible_hyd = 9999999999999999999 
                        # elif 'H tank' in self.system and self.system['hydrogen demand']['strategy'] != 'supply-led':   # hydrogen can only be stored into an H tank 
                        elif 'H tank' in self.system and 'HPH tank' not in self.system and self.system[self.hydrogen_demand+' demand']['strategy'] != 'supply-led':   # hydrogen can only be stored into an H tank 
                            producible_hyd  = self.technologies['H tank'].max_capacity-self.technologies['H tank'].LOC[h] # the tank can't be full
                            if producible_hyd < self.technologies['H tank'].max_capacity*0.00001: # to avoid unnecessary iteration
                                producible_hyd = 0
                        elif 'H tank' in self.system and 'HPH tank' not in self.system and self.system[self.hydrogen_demand+' demand']['strategy'] == 'supply-led':   # hydrogen can only be stored into an H tank 
                            producible_hyd = 9999999999999999999   # electrolyzer can produce continuously as the storage capacity is infinite. Tank is dimensioned at the end of simulation
                        elif 'H tank' in self.system and 'HPH tank' in self.system:
                            # storable_hydrogen_lp    = self.technologies['LPH tank'].max_capacity-self.technologies['LPH tank'].LOC[h] # the tank can't be full
                            producible_hyd   = self.technologies['H tank'].max_capacity-self.technologies['H tank'].LOC[h] # the tank can't be full                        
                        else:
                            producible_hyd = max(0,-eb['hydrogen']) # hydrogen is consumed by a technology which have higher priority than tank
                        if producible_hyd > 0:
                            self.energy_balance['hydrogen']['electrolyzer'][h],   \
                            self.energy_balance['electricity']['electrolyzer'][h],\
                            self.energy_balance['oxygen']['electrolyzer'][h],     \
                            self.energy_balance['water']['electrolyzer'][h]        = self.technologies['electrolyzer'].use(h,eb['electricity'],producible_hyd)      # [:2] # hydrogen supplied by electrolyzer(+) # electricity absorbed by the electorlyzer(-) 
                            
                            eb['hydrogen']      += self.energy_balance['hydrogen']['electrolyzer'][h]
                            eb['electricity']   += self.energy_balance['electricity']['electrolyzer'][h]
                            eb['oxygen']        += self.energy_balance['oxygen']['electrolyzer'][h]
                            eb['water']         += self.energy_balance['water']['electrolyzer'][h]
                
                elif self.technologies['electrolyzer'].strategy == 'full-time': # electrolyzer working continuously at each time step of the simulation
                    if "electricity grid" in self.system and self.system["electricity grid"]["draw"]:  # to assure full-time operation the system must be connected to the grid
                        # if self.system[self.hydrogen_demand+' demand']['strategy'] == 'supply-led':    
                        producible_hyd = 9999999999999999999
                        self.energy_balance['hydrogen']['electrolyzer'][h],     \
                        self.energy_balance['electricity']['electrolyzer'][h],  \
                        self.energy_balance['oxygen']['electrolyzer'][h],       \
                        self.energy_balance['water']['electrolyzer'][h]         = self.technologies['electrolyzer'].use(h,eb['electricity'],producible_hyd)      # [:2] # hydrogen supplied by electrolyzer(+) # electricity absorbed by the electorlyzer(-) 
                    
                        eb['hydrogen']      += self.energy_balance['hydrogen']['electrolyzer'][h]
                        eb['electricity']   += self.energy_balance['electricity']['electrolyzer'][h]
                        eb['oxygen']        += self.energy_balance['oxygen']['electrolyzer'][h]
                        eb['water']         += self.energy_balance['water']['electrolyzer'][h]
                        
                        # if 'H tank' in self.system: 
                        #     self.energy_balance['hydrogen']['H tank'][h] = self.technologies['H tank'].use(h,eb['hydrogen'])
                        #     eb['hydrogen'] += self.energy_balance['hydrogen']['H tank'][h]
                            
                    # elif self.system[self.hydrogen_demand+' demand']['strategy'] == 'demand-led':     # only availabke with a minimum constant demand
                        
                        
                    #     if 'H tank' in self.system and 'HPH tank' not in self.system:
                    #           producible_hyd  = self.technologies['H tank'].max_capacity-self.technologies['H tank'].LOC[h] + abs(eb['hydrogen']) # the tank can't be full
                    #           if producible_hyd < self.technologies['H tank'].max_capacity*0.00001: # to avoid unnecessary iteration
                    #               producible_hyd = 0
                    #           if producible_hyd > 0:
                    #               self.energy_balance['hydrogen']['electrolyzer'][h],   \
                    #               self.energy_balance['electricity']['electrolyzer'][h],\
                    #               self.energy_balance['oxygen']['electrolyzer'][h],     \
                    #               self.energy_balance['water']['electrolyzer'][h]        = self.technologies['electrolyzer'].useh2(h,producible_hyd)      # [:2] # hydrogen supplied by electrolyzer(+) # electricity absorbed by the electorlyzer(-) 
                               
                    #               eb['hydrogen']      += self.energy_balance['hydrogen']['electrolyzer'][h]
                    #               eb['electricity']   += self.energy_balance['electricity']['electrolyzer'][h]
                    #               eb['oxygen']        += self.energy_balance['oxygen']['electrolyzer'][h]
                    #               eb['water']         += self.energy_balance['water']['electrolyzer'][h]
                            
                             
   
                if h == (self.simulation_hours - 1) and ('hydrogen demand' in self.system or 'HP hydrogen demand' in self.system):
                    if self.system[self.hydrogen_demand+' demand']['strategy'] == 'supply-led':  # activates only at the final step of simulation
                        self.constant_flow = sum(self.energy_balance['hydrogen']['electrolyzer'])/self.simulation_hours # [kg/h] constant hydrogen output based on the total production
                    
            if tech_name == 'hydrogen compressor':   #!!! WIP to be modified by Andrea
                if self.energy_balance['hydrogen']['electrolyzer'][h] > 0:
                    if "hydrogen grid" in self.system and self.system["hydrogen grid"]["feed"]: # hydrogen can be fed into an hydrogen grid
                        storable_hyd = 9999999999999999999 
                    elif 'H tank' in self.system:   # hydrogen can only be stored into an H tank 
                        storable_hydrogen = self.technologies['H tank'].max_capacity-self.technologies['H tank'].LOC[h] # the tank can't be full
                    if storable_hydrogen>self.technologies['H tank'].max_capacity*0.00001:
                        self.energy_balance['hydrogen']['hydrogen compressor'][h], self.energy_balance['gas']['hydrogen compressor'][h] = self.technologies['hydrogen compressor'].use(h,self.energy_balance['hydrogen']['electrolyzer'][h],storable_hydrogen) # hydrogen compressed by the compressor (+) and heat requested to make it work expressed as heating water need (-) 
                        eb['gas'] += self.energy_balance['gas']['hydrogen compressor'][h]
                        # eb['hydrogen']=...self.energy_balance['hydrogen']['hydrogen compressor'][h]?? come ne tengo conto di quanto comprimo? in linea teorica ne dovrei sempre comprimere esattamente quanto me ne entra perchè il controllo sullo sotrable hydrogen lho gia fatto nell'elettrolizzatore'
            
            if tech_name == 'mechanical compressor':   
                if 'HPH tank' not in self.system:
                    if 'O2 tank' not in self.system:
                        if "electricity grid" in self.system and self.system["electricity grid"]["draw"]:
                            self.energy_balance['hydrogen']['mechanical compressor'][h], \
                            self.energy_balance['electricity']['mechanical compressor'][h] = self.technologies['mechanical compressor'].use(h,massflowrate= self.energy_balance['hydrogen']['electrolyzer'][h])[:2] # hydrogen compressed by the compressor (+) and electricity consumption (-) 
                            
                            eb['electricity']   += self.energy_balance['electricity']['mechanical compressor'][h]
                    
                        elif "electricity grid" not in self.system or self.system["electricity grid"]["draw"] == False:   # if the system is configurated as fully off-grid, relying only on RES production
                            if self.energy_balance['hydrogen']['electrolyzer'][h] > 0 :  # if hydrogen has been produced by the electrolyzer and electricity is available in the system
                                a = self.technologies['mechanical compressor'].use(h,massflowrate= self.energy_balance['hydrogen']['electrolyzer'][h])[1] # [kWh] compressor energy consumption for a certain h2 mass flow rate
                                if abs(a) <= eb['electricity']:     # there is enough renewable electricity to power the compressor 
                                    self.energy_balance['hydrogen']['mechanical compressor'][h],    \
                                    self.energy_balance['electricity']['mechanical compressor'][h], \
                                    self.energy_balance['cooling water']['mechanical compressor'][h]= self.technologies['mechanical compressor'].use(h,massflowrate= self.energy_balance['hydrogen']['electrolyzer'][h]) # hydrogen compressed by the compressor (+) and electricity consumption (-) 
                                    
                                    eb['electricity']   += self.energy_balance['electricity']['mechanical compressor'][h]
                                    
                                elif abs(a) > eb['electricity']:    # if available electricity in the system is not enough to power the compression system - enter the loop to reallocate the energy among the components
                                    a1  = 1     # % of available electricity fed to the electrolyzer
                                    a11 = 0     # % of available electricity fed to the compressor
                                    en  = eb['electricity'] + abs(self.energy_balance['electricity']['electrolyzer'][h]) # [kWh] electric energy available at time h before entering the electorlyzer
                                    el  = self.energy_balance['electricity']['electrolyzer'][h]
                                    hy  = self.energy_balance['hydrogen']['electrolyzer'][h]
                                    ox  = self.energy_balance['oxygen']['electrolyzer'][h]
                                    wa  = self.energy_balance['water']['electrolyzer'][h]
                                    
                                    # Iteration parameters
                                    i   = 0             # initializing iteration count
                                    maxiter = 10000     # max number of iterations allowed
                                    abs_err = 0.00001   # absolute error allowed
                                    
                                    while a1 >= 0:       # while loop necessary to iterate in the redistribution of renewable electricity to satisfy both electrolyzer and compressor demand
                                        hydrogen_ele,  \
                                        electricity_ele = self.technologies['electrolyzer'].use(h,a1*en,producible_hyd)[:2]  # [kg] of produced H2 and [kWh] of consumed electricity for the given energy input  
                                        a = -self.technologies['mechanical compressor'].use(h,massflowrate= hydrogen_ele)[1] # [kWh] compressor energy consumption for a certain h2 mass flow rate
                                        b1 = a/en
                                        a11 = 1-b1
                                        i += 1      # updating iteration count
                                    
                                        if abs(a1-a11) < abs_err or i > maxiter:    # strict tolerance for convergence 
                                            break
                                        else: 
                                            a1=a11    
                                    
                                    # Electorlyzer balances update and overwriting
                                    self.energy_balance['hydrogen']['electrolyzer'][h],   \
                                    self.energy_balance['electricity']['electrolyzer'][h],\
                                    self.energy_balance['oxygen']['electrolyzer'][h],     \
                                    self.energy_balance['water']['electrolyzer'][h]        = self.technologies['electrolyzer'].use(h,a1*en,producible_hyd)      # [:2] # hydrogen supplied by electrolyzer(+) # electricity absorbed by the electorlyzer(-) 
                                    
                                    eb['hydrogen']      += self.energy_balance['hydrogen']['electrolyzer'][h]    - hy
                                    eb['electricity']   += self.energy_balance['electricity']['electrolyzer'][h] - el
                                    eb['oxygen']        += self.energy_balance['oxygen']['electrolyzer'][h]      - ox
                                    eb['water']         += self.energy_balance['water']['electrolyzer'][h]       + wa
    
                                    # Compressor balances update and overwriting
                                    self.energy_balance['hydrogen']['mechanical compressor'][h],    \
                                    self.energy_balance['electricity']['mechanical compressor'][h], \
                                    self.energy_balance['cooling water']['mechanical compressor'][h]   = self.technologies['mechanical compressor'].use(h,massflowrate= self.energy_balance['hydrogen']['electrolyzer'][h]) # hydrogen compressed by the compressor (+) and electricity consumption (-) 
    
                                    eb['electricity']   += self.energy_balance['electricity']['mechanical compressor'][h]
    
                                
                                else:  # if no hydrogen has been produced at time h
                                    self.energy_balance['hydrogen']['mechanical compressor'][h]     = 0
                                    self.energy_balance['electricity']['mechanical compressor'][h]  = 0
                                    
                    elif 'O2 tank' in self.system:   # simplified approach for oxygen compression. To be updated
                        massflow_tot = (self.energy_balance['hydrogen']['electrolyzer'][h])+(self.energy_balance['oxygen']['electrolyzer'][h])
                        if "electricity grid" in self.system and self.system["electricity grid"]["draw"]:
                            self.energy_balance['hydrogen']['mechanical compressor'][h], \
                            self.energy_balance['electricity']['mechanical compressor'][h] = self.technologies['mechanical compressor'].use(h,massflowrate = massflow_tot )[:2] # hydrogen compressed by the compressor (+) and electricity consumption (-) 
                            
                            eb['electricity']   += self.energy_balance['electricity']['mechanical compressor'][h]
                    
                        elif "electricity grid" not in self.system or self.system["electricity grid"]["draw"] == False:   # if the system is configurated as fully off-grid, relying only on RES production
                            if self.energy_balance['hydrogen']['electrolyzer'][h] > 0 :  # if hydrogen has been produced by the electrolyzer and electricity is available in the system
                                a = self.technologies['mechanical compressor'].use(h,massflowrate = massflow_tot)[1] # [kWh] compressor energy consumption for a certain h2 mass flow rate
                                if abs(a) <= eb['electricity']:     # there is enough renewable electricity to power the compressor 
                                    self.energy_balance['hydrogen']['mechanical compressor'][h],    \
                                    self.energy_balance['electricity']['mechanical compressor'][h], \
                                    self.energy_balance['cooling water']['mechanical compressor'][h]= self.technologies['mechanical compressor'].use(h,massflowrate= massflow_tot) # hydrogen compressed by the compressor (+) and electricity consumption (-) 
                                    
                                    eb['electricity']   += self.energy_balance['electricity']['mechanical compressor'][h]
                                    
                                elif abs(a) > eb['electricity']:    # if available electricity in the system is not enough to power the compression system - enter the loop to reallocate the energy among the components
                                    a1  = 1     # % of available electricity fed to the electrolyzer
                                    a11 = 0     # % of available electricity fed to the compressor
                                    en  = eb['electricity'] + abs(self.energy_balance['electricity']['electrolyzer'][h]) # [kWh] electric energy available at time h before entering the electorlyzer
                                    el  = self.energy_balance['electricity']['electrolyzer'][h]
                                    hy  = self.energy_balance['hydrogen']['electrolyzer'][h]
                                    ox  = self.energy_balance['oxygen']['electrolyzer'][h]
                                    wa  = self.energy_balance['water']['electrolyzer'][h]
                                    
                                    # Iteration parameters
                                    i   = 0             # initializing iteration count
                                    maxiter = 10000     # max number of iterations allowed
                                    abs_err = 0.00001   # absolute error allowed
                                    
                                    while a1 >= 0:       # while loop necessary to iterate in the redistribution of renewable electricity to satisfy both electrolyzer and compressor demand
                                        hydrogen_ele,  \
                                        electricity_ele = self.technologies['electrolyzer'].use(h,a1*en,producible_hyd)[:2]  # [kg] of produced H2 and [kWh] of consumed electricity for the given energy input  
                                        a = -self.technologies['mechanical compressor'].use(h,massflowrate= hydrogen_ele + hydrogen_ele*7.93)[1] # [kWh] compressor energy consumption for a certain h2 mass flow rate
                                        b1 = a/en
                                        a11 = 1-b1
                                        i += 1      # updating iteration count
                                    
                                        if abs(a1-a11) < abs_err or i > maxiter:    # strict tolerance for convergence 
                                            break
                                        else: 
                                            a1=a11    
                                    
                                    # Electorlyzer balances update and overwriting
                                    self.energy_balance['hydrogen']['electrolyzer'][h],   \
                                    self.energy_balance['electricity']['electrolyzer'][h],\
                                    self.energy_balance['oxygen']['electrolyzer'][h],     \
                                    self.energy_balance['water']['electrolyzer'][h]        = self.technologies['electrolyzer'].use(h,a1*en,producible_hyd)      # [:2] # hydrogen supplied by electrolyzer(+) # electricity absorbed by the electorlyzer(-) 
                                    
                                    eb['hydrogen']      += self.energy_balance['hydrogen']['electrolyzer'][h]    - hy
                                    eb['electricity']   += self.energy_balance['electricity']['electrolyzer'][h] - el
                                    eb['oxygen']        += self.energy_balance['oxygen']['electrolyzer'][h]      - ox
                                    eb['water']         += self.energy_balance['water']['electrolyzer'][h]       + wa
    
                                    # Compressor balances update and overwriting
                                    self.energy_balance['hydrogen']['mechanical compressor'][h],    \
                                    self.energy_balance['electricity']['mechanical compressor'][h], \
                                    self.energy_balance['cooling water']['mechanical compressor'][h]   = self.technologies['mechanical compressor'].use(h,massflowrate= self.energy_balance['hydrogen']['electrolyzer'][h] + self.energy_balance['oxygen']['electrolyzer'][h]) # hydrogen compressed by the compressor (+) and electricity consumption (-) 
    
                                    eb['electricity']   += self.energy_balance['electricity']['mechanical compressor'][h]
    
                                
                                else:  # if no hydrogen has been produced at time h
                                    self.energy_balance['hydrogen']['mechanical compressor'][h]     = 0
                                    self.energy_balance['electricity']['mechanical compressor'][h]  = 0
                        
                             
                
                if 'H tank' in self.system and 'HPH tank' in self.system:
                    # self.energy_balance['hydrogen']['H tank'][h] = self.technologies['H tank'].use(h,eb['hydrogen'])
                    # eb['hydrogen'] += self.energy_balance['hydrogen']['H tank'][h]
                    available_hyd_lp = self.technologies['H tank'].LOC[h] + self.technologies['H tank'].max_capacity - self.technologies['H tank'].used_capacity
                    storable_hydrogen_hp = self.technologies['HPH tank'].max_capacity-self.technologies['HPH tank'].LOC[h]
                    
                    if self.technologies['HPH tank'].LOC[h] == self.technologies['HPH tank'].max_capacity:  # if High-Pressure-Tank is full 
                        self.energy_balance['HP hydrogen']['mechanical compressor'][h]      = 0     
                        self.energy_balance['electricity']['mechanical compressor'][h]      = 0    
                        self.energy_balance['cooling water']['mechanical compressor'][h]    = 0
                        self.energy_balance['hydrogen']['mechanical compressor'][h]         = 0
                        
                        eb['HP hydrogen']    += 0
                        eb['electricity']    += 0   # compressor not working
                        eb['hydrogen']       += 0   # compressor not working
                    else:  # if there is enough room available in the High Pressure Tank, the compressor is activated
                        self.energy_balance['HP hydrogen']['mechanical compressor'][h],     \
                        self.energy_balance['electricity']['mechanical compressor'][h],     \
                        self.energy_balance['cooling water']['mechanical compressor'][h]    = self.technologies['mechanical compressor'].use(h, available_hyd_lp=available_hyd_lp ,storable_hydrogen_hp=storable_hydrogen_hp) # hydrogen supplied by H tank (+) and electricity absorbed(-) 
                        self.energy_balance['hydrogen']['mechanical compressor'][h] = - self.energy_balance['HP hydrogen']['mechanical compressor'][h]
                        
                        eb['HP hydrogen'] += self.energy_balance['HP hydrogen']['mechanical compressor'][h]
                        eb['hydrogen']    += self.energy_balance['hydrogen']['mechanical compressor'][h]
                        eb['electricity'] += self.energy_balance['electricity']['mechanical compressor'][h]

            if tech_name == 'fuel cell':
                if eb['electricity'] < 0: #? this condition must be solved if you want to produce electricity to be fed into the gird
                    if "hydrogen grid" in self.system and self.system["hydrogen grid"]["draw"]: # hydrogen can be withdranw from an hydrogen grid
                        available_hyd = 9999999999999999999 
                    elif 'H tank' in self.system:   # only hydrogen inside H tank can be used
                        available_hyd = self.technologies['H tank'].LOC[h] + self.technologies['H tank'].max_capacity - self.technologies['H tank'].used_capacity                                                                  
                    else:
                        available_hyd = max(0,eb['hydrogen']) # hydrogen is produced by an electrolyzer with a higher priority than fc
                    if available_hyd > 0:
                        use = self.technologies['fuel cell'].use(h,eb['electricity'],available_hyd)     # saving fuel cell working parameters for the current timeframe
                        self.energy_balance['hydrogen']['fuel cell'][h] =    use[0] # hydrogen absorbed by fuel cell(-)
                        self.energy_balance['electricity']['fuel cell'][h] = use[1] # electricity supplied(+) 
        
                        if use[2] < -eb['heating water']: #all of the heat producted by FC is used      
                            self.energy_balance['heating water']['fuel cell'][h]=use[2] # heat loss of fuel cell
                        else:
                            self.energy_balance['heating water']['fuel cell'][h]=-eb['heating water'] # heat loss of fuel cell- demand

                        eb['hydrogen'] += self.energy_balance['hydrogen']['fuel cell'][h]
                        eb['electricity'] += self.energy_balance['electricity']['fuel cell'][h]
                        eb['heating water'] += self.energy_balance['heating water']['fuel cell'][h] 
                    
            if tech_name == 'boiler_h2':                                                                                                                                   
                if "hydrogen grid" in self.system and self.system["hydrogen grid"]["draw"]: # hydrogen can be withdranw from an hydrogen grid
                    available_hyd = 9999999999999999999 
                elif 'H tank' in self.system:   # only hydrogen inside H tank can be used
                    available_hyd = self.technologies['H tank'].LOC[h] + self.technologies['H tank'].max_capacity - self.technologies['H tank'].used_capacity + eb['hydrogen']
                                         
                self.energy_balance['hydrogen']['boiler_h2'][h], self.energy_balance['gas']['boiler_h2'][h] = self.technologies['boiler_h2'].use(eb['gas'],available_hyd,1)[1:3] #h2 consumed from boiler_h2 and heat produced by boiler_h2
                eb['hydrogen'] += self.energy_balance['hydrogen']['boiler_h2'][h] # hydrogen balance update: - hydrogen consumed by boiler_h2
                eb['gas'] += self.energy_balance['gas']['boiler_h2'][h] # heat balance update: + heat produced by boiler_h2
           
            #!!!ANDREA HA MESSO QUESTO COME PROCESSO NEL LOCATION, QUELLO SOTTO A COSA é DOVUTO   WIP to be modified by Andrea                                                                                                                                                          
            if tech_name == 'boiler_h2': 
                if eb['electricity'] < 0: #? this condition must be solved if you want to produce electricity to fed into the gird
                    if "hydrogen grid" in self.system and self.system["hydrogen grid"]["draw"]: # hydrogen can be withdranw from an hydrogen grid
                        available_hyd = 9999999999999999999 
                    elif 'H tank' in self.system:   # only hydrogen inside H tank can be used
                        available_hyd = self.technologies['H tank'].LOC[h] + self.technologies['H tank'].max_capacity - self.technologies['H tank'].used_capacity                                                                  
                    else:
                        available_hyd = max(0,eb['hydrogen']) # hydrogen is produced by an electrolyzer which have higher priority than boiler_h2
                    if available_hyd > 0:
                        self.energy_balance['hydrogen']['boiler_h2'][h], self.energy_balance['heating water']['boiler_h2'][h] = self.technologies['boiler_h2'].use(eb['heating water'],available_hyd,1)[1:3] #h2 consumed from boiler_h2 and heat produced by boiler_h2
                        eb['hydrogen'] += self.energy_balance['hydrogen']['boiler_h2'][h] # hydrogen balance update: - hydrogen consumed by boiler_h2
                        eb['heating water'] += self.energy_balance['heating water']['boiler_h2'][h] # heat balance update: + heat produced by boiler_h2
                    
            # if tech_name in ['H tank','HPH tank']:     #versione buona 
            #     if self.system['hydrogen demand']['strategy'] != 'supply-led':
            #         self.energy_balance[self.tank_stream[tech_name]][tech_name][h] = self.technologies[tech_name].use(h,eb[self.tank_stream[tech_name]])
            #         eb[self.tank_stream[tech_name]] += self.energy_balance[self.tank_stream[tech_name]][tech_name][h]
            #     elif self.system['hydrogen demand']['strategy'] == 'supply-led' and h == (self.simulation_hours - 1):
            #         prod = self.energy_balance['hydrogen']['electrolyzer']
            #         for h in range(self.simulation_hours):
            #             self.energy_balance[self.tank_stream[tech_name]][tech_name][h] = self.technologies[tech_name].use(h,prod[h],constant_demand=self.constant_flow )                        
            #     else:
            #         pass
        
            if tech_name == 'H tank':
                if 'HPH tank' not in self.system and ('hydrogen demand' in self.system or 'HP hydrogen demand' in self.system):
                    if self.system[self.hydrogen_demand+' demand']['strategy'] == 'demand-led':
                        self.energy_balance['hydrogen']['H tank'][h] = self.technologies['H tank'].use(h,eb['hydrogen'])
                        eb['hydrogen'] += self.energy_balance['hydrogen']['H tank'][h]
                    elif self.system[self.hydrogen_demand+' demand']['strategy'] == 'supply-led' and h == (self.simulation_hours - 1):
                        prod = self.energy_balance['hydrogen']['electrolyzer']
                        for h in range(self.simulation_hours):
                            self.energy_balance['hydrogen']['H tank'][h] = self.technologies['H tank'].use(h,prod[h],constant_demand=self.constant_flow )                        
                    # else:
                        # pass
                else:
                    self.energy_balance['hydrogen']['H tank'][h] = self.technologies['H tank'].use(h,eb['hydrogen'])
                    eb['hydrogen'] += self.energy_balance['hydrogen']['H tank'][h]
            
            if tech_name == 'HPH tank':
                self.energy_balance['HP hydrogen']['HPH tank'][h] = self.technologies['HPH tank'].use(h,eb['HP hydrogen'])
                eb['HP hydrogen'] += self.energy_balance['HP hydrogen']['HPH tank'][h]
                
            if tech_name == 'O2 tank':
                if 'oxygen demand' in self.system and self.system['oxygen demand']['strategy'] != 'supply-led':
                    self.energy_balance['oxygen']['O2 tank'][h] = self.technologies['O2 tank'].use(h,eb['oxygen'])
                    eb['oxygen'] += self.energy_balance['oxygen']['O2 tank'][h]
                elif self.system['hydrogen demand']['strategy'] == 'supply-led' and h == (self.simulation_hours - 1):
                    self.technologies['O2 tank'].sizing(self.technologies['H tank'].max_capacity)
                else:
                    pass
            
            if tech_name == 'inverter':
                self.energy_balance['electricity']['inverter'][h] = self.technologies['inverter'].use(h,eb['electricity']) # electricity lost in conversion by the inverter
                eb['electricity'] += self.energy_balance['electricity']['inverter'][h] # electricity balance update: - electricity lost in conversion by the invertert

            ### demand and grid   
            for carrier in eb: # for each energy carrier
                if tech_name == f"{carrier} demand":                
                    eb[carrier] += self.energy_balance[carrier]['demand'][h]    # energy balance update: energy demand(-)  
                if tech_name == f"{carrier} grid":
                    if eb[carrier] > 0 and self.system[f"{carrier} grid"]['feed'] or eb[carrier] < 0 and self.system[f"{carrier} grid"]['draw']:
                        self.energy_balance[carrier]['grid'][h] = - eb[carrier] # energy from grid(+) or into grid(-) 
                        eb[carrier] += self.energy_balance[carrier]['grid'][h]  # elecricity balance update                                                                                  
            
        
        

