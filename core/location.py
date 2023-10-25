import numpy as np
import pandas as pd
from techs import (heatpump, boiler_el, boiler_ng, boiler_h2, PV, wind, battery, H_tank, HPH_tank, O2_tank, fuel_cell, electrolyzer, inverter, chp_gt, Chp, Absorber, hydrogen_compressor, Compressor)
from core import constants as c

class location:
    
    def __init__(self,system,location_name,path,check,file_structure,file_general):
        """
        Create a location object (producer, consumer or prosumer) 
    
        system: dictionary (all the inputs are optional)
            'demand'': dictionary
                'electricity':              str 'file_name.csv' time series of electricity demand  [kW]
                'heating water':            str 'file_name.csv' time series of heating and dhw demand  [kW]
                'cooling water':            str 'file_name.csv' time series of cooling demand  [kW]
                'hydrogen':                 str 'file_name.csv' time series of hydrogen demand  [kg/s]
                'HP hydrogen':              str 'file_name.csv' time series of High-Pressure hydrogen demand  [kg/s]
                'gas':                      str 'file_name.csv' time series of gas demand  [kW]
                'process steam':            str 'file_name.csv' time series of process steam demand  [kg/s]
                'process heat':             str 'file_name.csv' time series of process heat demand  [kW]
                'process hot water':        str 'file_name.csv' time series of process hot water demand  [kW]
                'process cold water':       str 'file_name.csv' time series of process cold water demand (absorber, 7-12 °C)  [kW]
                'process chilled water':    str 'file_name.csv' time series of process chilled water demand (absorber, 1-5 °C)  [kW]
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
            
        output : location object able to:
            simulate the energy flows of present technologies .loc_simulation
            record power balances (electricity, heating water, cooling water, gas, hydrogen)
        """
        
        self.system = dict(sorted(system.items(), key=lambda item: item[1]['priority'])) # ordered by priority
        self.name = location_name
        self.technologies = {}                                 # initialise technologies dictionary
        self.power_balance = {'electricity'            : {},   # initialise power balances dictionaries
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

        # create the objects of present technologies and add them to the technologies dictionary
        # initialise power balance and add them to the power_balance dictionary
        
        for carrier in self.power_balance: 
            
            if f"{carrier} grid" in self.system: 
                self.power_balance[carrier]['grid'] = np.zeros(c.timestep_number) # array energy carrier bought from the grid (-) or feed into the grid (+)

            if f"{carrier} demand" in self.system: # read and check demand series
                
                if carrier == 'hydrogen' or carrier == 'HP hydrogen':
                    self.hydrogen_demand = carrier  # demand can be defined as 'hydrogen demand' or 'HP hydrogen demand' depending on the required delivery pressure
                    if self.system[carrier+' demand']['strategy'] == 'supply-led' and self.system[carrier+' demand']['serie'] != False:
                        raise ValueError(f"Warning in {self.name} location: supply-led strategy is not consistent with providing a demand series.\n\
                        Options to fix the problem: \n\
                            (a) - Insert 'false' at {carrier} demand 'serie' in studycase.json\n\
                            (b) - Change 'strategy' to 'demand-led' in studycase.json")
                            
                if carrier in ['hydrogen','HP hydrogen','process steam'] and self.system[carrier+' demand']['serie'] != False:
                    self.power_balance[carrier]['demand'] = - pd.read_csv(path+'/loads/'+system[f"{carrier} demand"]['serie'])['kg/s'].to_numpy() 
                    
                elif carrier in ['hydrogen','HP hydrogen','process steam'] and not self.system[carrier+' demand']['serie']:
                    if self.system[carrier+' demand']['strategy'] == 'supply-led' and self.system[carrier+' demand']['serie'] == False:
                        self.power_balance[carrier]['demand'] =  np.zeros(c.timestep_number)
                        
                else:     
                     self.power_balance[carrier]['demand'] = - pd.read_csv(path+'/loads/'+system[f"{carrier} demand"]['serie'])['kW'].to_numpy() 
                     
                ### check demand serie length
                if len(self.power_balance[carrier]['demand']) == c.timestep_number: 
                    pass 
                
                elif len(self.power_balance[carrier]['demand'])<c.timestep_number:
                                    
                    if c.timestep_number%len(self.power_balance[carrier]['demand']) == 0:
                        self.power_balance[carrier]['demand'] = np.tile(self.power_balance[carrier]['demand'],int(c.timestep_number/len(self.power_balance[carrier]['demand']))  ) 
                else:
                    raise ValueError(f"Warning! Checks the length of the {carrier} demand series you input for {self.name}. In general.json it checks the timestep and the simulation length")

        if 'chp_gt' in self.system:
            self.technologies['chp_gt'] = chp_gt(system['chp_gt'],c.timestep_number) # chp_gt object created and added to 'technologies' dictionary
            self.power_balance['process steam']['chp_gt'] = np.zeros(c.timestep_number) # array chp_gt process steam balance 
            self.power_balance['electricity']['chp_gt'] = np.zeros(c.timestep_number) # array chp_gt electricity balance
            self.power_balance['hydrogen']['chp_gt'] = np.zeros(c.timestep_number) # array chp_gt process hydrogen balance 
        
        if 'chp' in self.system:
            self.technologies['chp'] = Chp(system['chp'],c.timestep_number) # chp object created and added to 'technologies' dictionary
            self.power_balance[self.technologies['chp'].th_out]['chp'] = np.zeros(c.timestep_number) # array chp thermal output balance (process steam/hot water)
            self.power_balance['electricity']['chp'] = np.zeros(c.timestep_number) # array chp electricity balance
            self.power_balance[self.technologies['chp'].fuel]['chp'] = np.zeros(c.timestep_number) # array chp fuel consumption balance
            self.power_balance['process heat']['chp'] = np.zeros(c.timestep_number) # array chp process heat balance
            self.power_balance['process hot water']['chp'] = np.zeros(c.timestep_number) # array chp process hot water balance
            self.power_balance['process cold water']['chp'] = np.zeros(c.timestep_number) # array chp process cold water balance
       
        if 'absorber' in self.system:
            self.technologies['absorber'] = Absorber(system['absorber'],c.timestep_number) # absorber object created and added to 'technologies' dictionary
            self.power_balance['process heat']['absorber'] = np.zeros(c.timestep_number) # array absorber process steam balance 
            self.power_balance['process hot water']['absorber'] = np.zeros(c.timestep_number) # array absorber process steam balance 
            self.power_balance['process cold water']['absorber'] = np.zeros(c.timestep_number) # array absorber process steam balance 
            
        if 'heatpump' in self.system:
            self.technologies['heatpump'] = heatpump(system['heatpump'],c.timestep_number) # heatpump object created and add to 'technologies' dictionary
            self.power_balance['electricity']['heatpump'] = np.zeros(c.timestep_number) # array heatpump electricity balance
            self.power_balance['heating water']['heatpump'] = np.zeros(c.timestep_number) # array heatpump heat balance
            self.power_balance['heating water']['inertial TES'] = np.zeros(c.timestep_number) # array inertial tank heat balance
                
        if 'boiler_el' in self.system:
            self.technologies['boiler_el'] = boiler_el(self.system['boiler_el']) # boiler_el object created and add to 'technologies' dictionary
            self.power_balance['electricity']['boiler_el'] = np.zeros(c.timestep_number) # array boiler_el electricity balance
            self.power_balance['heating water']['boiler_el'] = np.zeros(c.timestep_number) # array boiler_el heat balance
               
        if 'boiler_ng' in self.system:
            self.technologies['boiler_ng'] = boiler_ng(self.system['boiler_ng'])  # boiler_ng object created and add to 'technologies' dictionary
            self.power_balance['gas']['boiler_ng'] = np.zeros(c.timestep_number) # array boiler_ng gas balance
            self.power_balance['heating water']['boiler_ng'] = np.zeros(c.timestep_number) # array boiler_ng heat balance 

        if 'PV' in self.system:
            self.technologies['PV'] = PV(self.system['PV'],self.name,path,check,file_structure,file_general) # PV object created and add to 'technologies' dictionary
            self.power_balance['electricity']['PV'] = np.zeros(c.timestep_number) # array PV electricity balance
           
        if 'inverter' in self.system:
            self.technologies['inverter'] = inverter(self.system['inverter'],c.timestep_number) # inverter object created and add to 'technologies' dictionary
            self.power_balance['electricity']['inverter'] = np.zeros(c.timestep_number)        # array inverter electricity balance
            
        if 'wind' in self.system:
            self.technologies['wind'] = wind(self.system['wind'],c.timestep_number,c.timestep,path)    # wind object created and add to 'technologies' dictionary
            self.power_balance['electricity']['wind'] = np.zeros(c.timestep_number)        # array wind electricity balance 
           
        if 'battery' in self.system:
            self.technologies['battery'] = battery(self.system['battery'])    # battery object created and to 'technologies' dictionary
            self.power_balance['electricity']['battery'] = np.zeros(c.timestep_number)         # array battery electricity balance
                           
        if 'electrolyzer' in self.system:
            self.technologies['electrolyzer'] = electrolyzer(self.system['electrolyzer'],c.timestep_number) # electrolyzer object created and to 'technologies' dictionary
            self.power_balance['electricity']['electrolyzer']              = np.zeros(c.timestep_number) # array electrolyzer electricity balance
            self.power_balance['oxygen']['electrolyzer']                   = np.zeros(c.timestep_number) # array electrolyzer oxygen balance
            self.power_balance['water']['electrolyzer']                    = np.zeros(c.timestep_number) # array electrolyzer water balance
            self.power_balance['hydrogen']['electrolyzer']                 = np.zeros(c.timestep_number) # array electrolyzer hydrogen balance
            if 'hydrogen demand' not in self.system and 'hydrogen grid' not in self.system and self.system['electrolyzer']['only_renewables'] == False :
                raise ValueError(f"Electrolyzers operation considered only for renewable energy long term storage in {self.name} location. It can be powered only by renewables\n\
                                 Options to fix the problem: \n\
                                     (a) - Change electrolyzer 'only_renewables' parameter to 'true' in studycase.json\\")
            if self.technologies['electrolyzer'].strategy == "full-time" and (not self.system["electricity grid"]["draw"] or self.system['electrolyzer']['only_renewables']):
                raise ValueError(f"Full-time electrolyzers operation considered without electricity grid connection in {self.name} location.\n\
                Options to fix the problem: \n\
                    (a) - Insert electricity grid withdrawal in studycase.json\n\
                    (b) - Change electrolyzers strategy in studycase.json\n\
                    (c) - Check electrolyzers for 'only_renewables' value in studycase.json to be 'false'\\ ")
            if self.technologies['electrolyzer'].strategy == "full-time" and self.system[self.hydrogen_demand+' demand']['strategy'] == 'supply-led':  
                if 'mechanical compressor' in self.system and 'tank' in self.system: # when electrolyzers working continuously in supply-mode there is no need for storage
                    raise ValueError(f"Full-time electrolyzers operation considered in supply-led mode in {self.name} location. Compression and storage must not be considered\n\
                    Options to fix the problem: \n\
                        (a) - Remove 'mechanical compressor' and 'H tank' from studycase.json\\")
            if self.technologies['electrolyzer'].strategy == "full-time" and self.system[self.hydrogen_demand+' demand']['strategy'] == 'demand-led':
                raise ValueError(f"Warning in {self.name} location: full-time electrolyzers operation not consistent in demand-led mode.\n\
                Feasible combinations: \n\
                    (a) - hydrogen demand:'supply-led' & electrolyzer strategy:'full-time' \n\
                    (b) - hydrogen demand:'demand-led' & electrolyzer strategy:'hydrogen-first' ")
            if self.system['electrolyzer']['only_renewables'] == False and self.system["electricity grid"]["draw"] == False:
                raise ValueError(f"If 'only_renewables' strategy for electrolyzers operation is set 'false', grid connection in {self.name} location must be allowed.\n\
                Options to fix the problem: \n\
                    (a) - Insert electricity grid withdrawal in studycase.json\n\
                    (b) - Change electrolyzers 'only_renewables' strategy in studycase.json\n\\ ")
                
        if 'fuel cell' in self.system:
            self.technologies['fuel cell'] = fuel_cell(self.system['fuel cell'],c.timestep_number,c.timestep) # Fuel cell object created and to 'technologies' dictionary
            self.power_balance['electricity']['fuel cell']     = np.zeros(c.timestep_number)     # array fuel cell electricity balance
            self.power_balance['hydrogen']['fuel cell']        = np.zeros(c.timestep_number)     # array fuel cell hydrogen balance
            self.power_balance['heating water']['fuel cell']   = np.zeros(c.timestep_number)     # array fuel cell heat balance used
            
        if 'hydrogen compressor' in self.system:
            self.technologies['hydrogen compressor'] = hydrogen_compressor(self.system['hydrogen compressor'],c.timestep_number) # MHHC compressor object created and to 'technologies' dictionary
            self.power_balance['hydrogen']['hydrogen compressor']  = np.zeros(c.timestep_number)     # array hydrogen compressor hydrogen compressed
            self.power_balance['gas']['hydrogen compressor']       = np.zeros(c.timestep_number)     # array hydrogen compressor heating water balanced used
        
        if 'H tank' and 'HPH tank' in self.system: 
            self.technologies['H tank'] = H_tank(self.system['H tank'],c.timestep_number) # H tank object created and to 'technologies' dictionary
            self.power_balance['hydrogen']['H tank'] = np.zeros(c.timestep_number)         # array H tank hydrogen balance
            self.technologies['HPH tank'] = HPH_tank(self.system['HPH tank'],c.timestep_number)   # HPH tank object created and to 'technologies' dictionary
            self.power_balance['HP hydrogen']['HPH tank'] = np.zeros(c.timestep_number)        # array HPH tank hydrogen balance
            
            self.tank_stream = {'H tank':'hydrogen',        # dictionary assigning different hydrogen streams to different storage technologies - necessary for loc_energy_simulation
                                'HPH tank':'HP hydrogen'}
        
        if 'O2 tank' in self.system: 
            self.technologies['O2 tank'] = O2_tank(self.system['O2 tank'],c.timestep_number) # LPH tank object created and to 'technologies' dictionary
            self.power_balance['oxygen']['O2 tank'] = np.zeros(c.timestep_number)         # array LPH tank hydrogen balance
        
            
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
            self.technologies['mechanical compressor'] = Compressor(self.system['mechanical compressor'],c.timestep_number,maxflowrate = maxflowrate) # H tank object created and to 'technologies' dictionary
            self.power_balance['electricity']['mechanical compressor']    = np.zeros(c.timestep_number) # array H tank hydrogen balance
            self.power_balance['hydrogen']['mechanical compressor']       = np.zeros(c.timestep_number) # array of hydrogen flow entering the mechanical compressor from LPH tank
            self.power_balance['HP hydrogen']['mechanical compressor']    = np.zeros(c.timestep_number) # array of compressed hydrogen flow sent toward HPH tank
            self.power_balance['cooling water']['mechanical compressor']  = np.zeros(c.timestep_number) # array of water flow to be fed to the refrigeration system 
 
        if 'H tank' in self.system and not 'HPH tank' in self.system:
            if 'hydrogen demand' in self.system or 'HP hydrogen demand' in self.system:
                if self.system[self.hydrogen_demand+' demand']['strategy'] == 'supply-led' and self.system['H tank']['max capacity'] != False:
                    raise ValueError(f"Adjust {self.name} location system in studycase.json. When the system is operated in supply-led mode, H tank size cannot be defined in advance.\n\
            Options to fix the problem: \n\
            (a) - Insert false for 'max capacity' among H tank parameters in studycase.json\n\
            (b) - Switch to 'demand-led' in 'hydrogen-demand'('strategy')\
            ")
            self.technologies['H tank'] = H_tank(self.system['H tank'],c.timestep_number)   # H tank object created and to 'technologies' dictionary
            self.power_balance['hydrogen']['H tank'] = np.zeros(c.timestep_number)         # array H tank hydrogen balance
            
            self.tank_stream = {'H tank':'hydrogen'}     # dictionary assigning hydrogen stream to H tank storage technologies - necessary for loc_energy_simulation
        
        if 'boiler_h2' in self.system:
            self.technologies['boiler_h2'] = boiler_h2(self.system['boiler_h2'])                    # boiler_h2 object created and added to 'technologies' dictionary
            self.power_balance['hydrogen']['boiler_h2'] = np.zeros(c.timestep_number)          # array boiler_h2 gas balance
            self.power_balance['heating water']['boiler_h2'] = np.zeros(c.timestep_number)     # array boiler_h2 heat balance 

        self.power_balance['electricity']['collective self consumption']   = np.zeros(c.timestep_number) # array contribution to collective-self-consumption as producer (-) or as consumer (+)
        #self.power_balance['heating water']['collective self consumption'] = np.zeros(c.timestep_number) # array contribution to collective-self-consumption as producer (-) or as consumer (+)---heat----mio!!!
        #self.power_balance['process steam']['collective self consumption'] = np.zeros(c.timestep_number) # array contribution to collective-self-consumption as producer (-) or as consumer (+)---heat----mio!!!
             
    def loc_power_simulation(self,step,weather):
        """
        Simulate the location
        
        input : int step to simulate
    
        output : updating of location power balances
        """
        
        pb = {} # power balance [kW] [kg/s]
        
        for carrier in self.power_balance:
           pb[carrier] = 0 # initialise power balances 
            
        for tech_name in self.system: # (which is ordered py priority)
        
            if tech_name == 'PV': 
                self.power_balance['electricity']['PV'][step] = self.technologies['PV'].use(step) # electricity produced from PV
                pb['electricity'] += self.power_balance['electricity']['PV'][step] # elecricity balance update: + electricity produced from PV
    
            if tech_name == 'wind': 
                self.power_balance['electricity']['wind'][step] = self.technologies['wind'].use(step,weather['wind_speed'][step]) # electricity produced from wind
                pb['electricity'] += self.power_balance['electricity']['wind'][step] # elecricity balance update: + electricity produced from wind
                        
            if tech_name == 'boiler_el': 
                self.power_balance['electricity']['boiler_el'][step], self.power_balance['heating water']['boiler_el'][step] = self.technologies['boiler_el'].use(pb['heating water']) # el consumed and heat produced from boiler_el
                pb['electricity'] += self.power_balance['electricity']['boiler_el'][step] # elecricity balance update: - electricity consumed by boiler_el
                pb['heating water'] += self.power_balance['heating water']['boiler_el'][step] # heat balance update: + heat produced by boiler_el
    
            if tech_name == 'boiler_ng': 
                self.power_balance['gas']['boiler_ng'][step], self.power_balance['heating water']['boiler_ng'][step] = self.technologies['boiler_ng'].use(pb['heating water']) # ng consumed and heat produced from boiler_ng
                pb['gas'] += self.power_balance['gas']['boiler_ng'][step] # gas balance update: - gas consumed by boiler_ng
                pb['heating water'] += self.power_balance['heating water']['boiler_ng'][step] # heat balance update: + heat produced by boiler_ng
        
            if tech_name == 'heatpump':     
                self.power_balance['electricity']['heatpump'][step], self.power_balance['heating water']['heatpump'][step], self.power_balance['heating water']['inertial TES'][step] = self.technologies['heatpump'].use(weather['temp_air'][step],pb['heating water'],pb['electricity'],step,self.P2E) 
             
                pb['electricity'] += self.power_balance['electricity']['heatpump'][step] # electricity absorbed by heatpump
                self.power_balance['electricity']['demand'][step] += self.power_balance['electricity']['heatpump'][step] # add heatpump demand to 'electricity demand'
                pb['heating water'] += self.power_balance['heating water']['inertial TES'][step] + self.power_balance['electricity']['heatpump'][step] # heat or cool supplied by HP or inertial TES
        
            if tech_name == 'battery':
                if self.technologies['battery'].collective == 0: 
                    self.power_balance['electricity']['battery'][step] = self.technologies['battery'].use(step,pb['electricity']) # electricity absorbed(-) or supplied(+) by battery
                    pb['electricity'] += self.power_balance['electricity']['battery'][step]  # electricity balance update: +- electricity absorbed or supplied by battery
           
            if tech_name == 'chp_gt':
                if "hydrogen grid" in self.system and self.system["hydrogen grid"]["draw"]: # hydrogen can be withdranw from an hydrogen grid
                    available_hyd = 9999999999999999999 
                elif 'H tank' in self.system:   # only hydrogen inside H tank can be used
                    available_hyd = self.technologies['H tank'].LOC[step] + self.technologies['H tank'].max_capacity - self.technologies['H tank'].used_capacity                                                                  
                else:
                    available_hyd = max(0,pb['hydrogen']) # hydrogen is produced in the same timestep by an electrolyzer with a higher priority than chp
                if available_hyd > 0:
                    use = self.technologies['chp_gt'].use(step,weather['temp_air'][step],pb['process steam'],available_hyd)     # saving chp_gt working parameters for the current timeframe
                    self.power_balance['process steam']['chp_gt'][step] = use[0]   # produced steam (+)
                    self.power_balance['electricity']['chp_gt'][step] =   use[1]   # produced electricity (+)
                    self.power_balance['hydrogen']['chp_gt'][step] =      use[2]   # hydrogen required by chp system to run (-)  
    
                pb['hydrogen'] += self.power_balance['hydrogen']['chp_gt'][step]            
                pb['process steam'] += self.power_balance['process steam']['chp_gt'][step]    
                pb['electricity'] += self.power_balance['electricity']['chp_gt'][step] 
           
            if tech_name == 'chp':
                strategy    = self.technologies['chp'].strategy     # thermal-load follow or electric-load follow 
                coproduct   = self.technologies['chp'].coproduct    # process co-product depending on the  approache chosen above 
                if self.system['chp']['Fuel'] == 'hydrogen':
                    if "hydrogen grid" in self.system and self.system["hydrogen grid"]["draw"]: # hydrogen can be withdranw from an hydrogen grid
                        available_hyd = 9999999999999999999 
                    elif 'H tank' in self.system:   # only hydrogen inside H tank can be used
                        available_hyd = self.technologies['H tank'].LOC[step] + self.technologies['H tank'].max_capacity - self.technologies['H tank'].used_capacity                                                                  
                    else:
                        available_hyd = max(0,pb['hydrogen']) # hydrogen is produced in the same timestep by an electrolyzer with a higher priority than chp
                    use = self.technologies['chp'].use(step,weather['temp_air'][step],pb[strategy],pb[coproduct], available_hyd)     # saving chp working parameters for the current timeframe
                else:
                    use = self.technologies['chp'].use(step,weather['temp_air'][step],pb[strategy],pb[coproduct])                    # saving chp working parameters for the current timeframe
                
                self.power_balance[self.technologies['chp'].th_out]['chp'][step]  = use[0]   # produced thermal output (+) (steam/hot water)
                self.power_balance['electricity']['chp'][step]                    = use[1]   # produced electricity (+)
                self.power_balance[self.technologies['chp'].fuel]['chp'][step]    = use[2]   # fuel required by chp system to run (-)  
                self.power_balance['process heat']['chp'][step]                   = use[3]   # process heat produced by chp system (+)  
                self.power_balance['process hot water']['chp'][step]              = use[4]   # process heat produced by chp system (+)  

                pb[self.technologies['chp'].fuel] += self.power_balance[self.technologies['chp'].fuel]['chp'][step]            
                pb[self.technologies['chp'].th_out] += self.power_balance[self.technologies['chp'].th_out]['chp'][step]    
                pb['electricity'] += self.power_balance['electricity']['chp'][step] 
                pb['process heat'] += self.power_balance['process heat']['chp'][step]             

            if tech_name == 'absorber':  
                self.power_balance['process cold water']['absorber'][step],self.power_balance['process heat']['absorber'][step] = self.technologies['absorber'].use(step,pb['process heat'])  # cold energy produced via the absorption cycle (+)
                pb['process heat'] += self.power_balance['process heat']['absorber'][step]             
                pb['process cold water'] += self.power_balance['process cold water']['absorber'][step]
                
            if tech_name == 'electrolyzer':
                
                if self.technologies['electrolyzer'].strategy == 'hydrogen-first' and self.technologies['electrolyzer'].only_renewables == True: # electrolyzer activated when renewable energy is available
                    if pb['electricity'] > 0: # electrolyzer activated only when renewable energy is available
                        if "hydrogen grid" in self.system and self.system["hydrogen grid"]["feed"] and 'H tank' not in self.system: # hydrogen can be fed into an hydrogen grid
                            producible_hyd = 9999999999999999999
                        elif 'hydrogen demand' not in self.system: # hydrogen-energy-storage configuration, only renewable energy is stored in the form of hydrogen to be converted back into electricity via fuel cell 
                            producible_hyd  = self.technologies['H tank'].max_capacity-self.technologies['H tank'].LOC[step] # the tank can't be full
                            if producible_hyd < self.technologies['H tank'].max_capacity*0.00001: # to avoid unnecessary iteration
                                producible_hyd = 0
                        elif 'H tank' in self.system and 'HPH tank' not in self.system and self.system[self.hydrogen_demand+' demand']['strategy'] == 'demand-led':   # hydrogen can only be stored into an H tank 
                            producible_hyd  = self.technologies['H tank'].max_capacity-self.technologies['H tank'].LOC[step] + (-pb['hydrogen']) # the tank can't be full
                            if producible_hyd < self.technologies['H tank'].max_capacity*0.00001:                           # to avoid unnecessary iteration
                                producible_hyd = 0
                        elif 'H tank' in self.system and 'HPH tank' not in self.system and self.system[self.hydrogen_demand+' demand']['strategy'] == 'supply-led':   # hydrogen can only be stored into an H tank 
                            producible_hyd = 9999999999999999999   # electrolyzer can produce continuously as the storage capacity is infinite. Tank is dimensioned at the end of simulation
                        elif 'H tank' in self.system and 'HPH tank' in self.system:
                            # storable_hydrogen_lp    = self.technologies['LPH tank'].max_capacity-self.technologies['LPH tank'].LOC[step] # the tank can't be full
                            producible_hyd   = self.technologies['H tank'].max_capacity-self.technologies['H tank'].LOC[step] # the tank can't be full                        
                        else:
                            producible_hyd = max(0,-pb['hydrogen']) # hydrogen is consumed by a technology with a higher priority than tank
                        if producible_hyd > 0:
                            self.power_balance['hydrogen']['electrolyzer'][step],   \
                            self.power_balance['electricity']['electrolyzer'][step],\
                            self.power_balance['oxygen']['electrolyzer'][step],     \
                            self.power_balance['water']['electrolyzer'][step]        = self.technologies['electrolyzer'].use(step,c.timestep,storable_hydrogen=producible_hyd,e=pb['electricity'])      # [:2] # hydrogen supplied by electrolyzer(+) # electricity absorbed by the electorlyzer(-) 
                            
                            pb['hydrogen']      += self.power_balance['hydrogen']['electrolyzer'][step]
                            pb['electricity']   += self.power_balance['electricity']['electrolyzer'][step]
                            pb['oxygen']        += self.power_balance['oxygen']['electrolyzer'][step]
                            pb['water']         += self.power_balance['water']['electrolyzer'][step]
                            
                elif self.technologies['electrolyzer'].strategy == 'hydrogen-first' and self.technologies['electrolyzer'].only_renewables == False: # electrolyzer working both with energy from renewables and from grid, but giving precedence to electricity from renewables
                    
                    if 'hydrogen grid' in self.system and self.system["hydrogen grid"]["feed"] and 'H tank' not in self.system: # hydrogen can be fed into an hydrogen grid
                        producible_hyd = 9999999999999999999
                    elif 'H tank' in self.system and 'HPH tank' not in self.system and self.system[self.hydrogen_demand+' demand']['strategy'] == 'demand-led':   # hydrogen can only be stored into an H tank 
                        producible_hyd  = self.technologies['H tank'].max_capacity-self.technologies['H tank'].LOC[step] + (-pb['hydrogen']) # the tank can't be full
                        if producible_hyd < self.technologies['H tank'].max_capacity*0.00001: # to avoid unnecessary iteration
                            producible_hyd = 0
                    elif 'H tank' in self.system and 'HPH tank' not in self.system and self.system[self.hydrogen_demand+' demand']['strategy'] == 'supply-led':   # hydrogen can only be stored into an H tank 
                        producible_hyd = 9999999999999999999   # electrolyzer can produce continuously as the storage capacity is infinite. Tank is dimensioned at the end of simulation
                    elif 'H tank' in self.system and 'HPH tank' in self.system:
                        producible_hyd   = self.technologies['H tank'].max_capacity-self.technologies['H tank'].LOC[step] # the tank can't be full                        
                    else:
                        producible_hyd = max(0,-pb['hydrogen']) # hydrogen is consumed by a technology which have higher priority than tank
                    if producible_hyd > 0:
                        self.power_balance['hydrogen']['electrolyzer'][step],   \
                        self.power_balance['electricity']['electrolyzer'][step],\
                        self.power_balance['oxygen']['electrolyzer'][step],     \
                        self.power_balance['water']['electrolyzer'][step]        = self.technologies['electrolyzer'].use(step,storable_hydrogen=producible_hyd,e=pb['electricity'])      # [:2] # hydrogen supplied by electrolyzer(+) # electricity absorbed by the electorlyzer(-) 
                            
                    available_hyd = self.technologies['H tank'].LOC[step] + self.technologies['H tank'].max_capacity - self.technologies['H tank'].used_capacity
                    if available_hyd + self.power_balance['hydrogen']['electrolyzer'][step] < -pb['hydrogen']:     # hydrogen produced from electrolyzer with only renewables and the tank are not sufficient to cover the hydrogen demand in this hour --> need for grid interaction
                        hyd_from_ele = (-pb['hydrogen']) - available_hyd     # The electrolyzer must generate only the quantity of hydrogen the tank can't cover (thus using even electricity from grid) 
                        self.power_balance['hydrogen']['electrolyzer'][step],   \
                        self.power_balance['electricity']['electrolyzer'][step],\
                        self.power_balance['oxygen']['electrolyzer'][step],     \
                        self.power_balance['water']['electrolyzer'][step]        = self.technologies['electrolyzer'].use(step,hydrog=hyd_from_ele)      # [:2] # hydrogen supplied by electrolyzer(+) # electricity absorbed by the electorlyzer(-) 
                        
                        pb['hydrogen']      += self.power_balance['hydrogen']['electrolyzer'][step]
                        pb['electricity']   += self.power_balance['electricity']['electrolyzer'][step]
                        pb['oxygen']        += self.power_balance['oxygen']['electrolyzer'][step]
                        pb['water']         += self.power_balance['water']['electrolyzer'][step]

                    else:
                        pb['hydrogen']      += self.power_balance['hydrogen']['electrolyzer'][step]
                        pb['electricity']   += self.power_balance['electricity']['electrolyzer'][step]
                        pb['oxygen']        += self.power_balance['oxygen']['electrolyzer'][step]
                        pb['water']         += self.power_balance['water']['electrolyzer'][step]
                
                elif self.technologies['electrolyzer'].strategy == 'full-time': # electrolyzer working continuously at each time step of the simulation
                    if "electricity grid" in self.system and self.system["electricity grid"]["draw"]:  # to assure full-time operation the system must be connected to the grid
                        # if self.system[self.hydrogen_demand+' demand']['strategy'] == 'supply-led':    
                        producible_hyd = 9999999999999999999
                        self.power_balance['hydrogen']['electrolyzer'][step],     \
                        self.power_balance['electricity']['electrolyzer'][step],  \
                        self.power_balance['oxygen']['electrolyzer'][step],       \
                        self.power_balance['water']['electrolyzer'][step]         = self.technologies['electrolyzer'].use(step,storable_hydrogen=producible_hyd,e=pb['electricity'])      # [:2] # hydrogen supplied by electrolyzer(+) # electricity absorbed by the electorlyzer(-) 
                    
                        pb['hydrogen']      += self.power_balance['hydrogen']['electrolyzer'][step]
                        pb['electricity']   += self.power_balance['electricity']['electrolyzer'][step]
                        pb['oxygen']        += self.power_balance['oxygen']['electrolyzer'][step]
                        pb['water']         += self.power_balance['water']['electrolyzer'][step]
                        
   
                if step == (c.timestep_number - 1) and ('hydrogen demand' in self.system or 'HP hydrogen demand' in self.system):
                    if self.system[self.hydrogen_demand+' demand']['strategy'] == 'supply-led':  # activates only at the final step of simulation
                        self.constant_flow = sum(self.power_balance['hydrogen']['electrolyzer'])/c.timestep_number # [kg/s] constant hydrogen output based on the total production
                    
            if tech_name == 'hydrogen compressor':   #!!! WIP to be modified by Andrea
                if self.power_balance['hydrogen']['electrolyzer'][step] > 0:
                    if "hydrogen grid" in self.system and self.system["hydrogen grid"]["feed"]: # hydrogen can be fed into an hydrogen grid
                        storable_hyd = 9999999999999999999 
                    elif 'H tank' in self.system:   # hydrogen can only be stored into an H tank 
                        storable_hydrogen = self.technologies['H tank'].max_capacity-self.technologies['H tank'].LOC[step] # the tank can't be full
                    if storable_hydrogen>self.technologies['H tank'].max_capacity*0.00001:
                        self.power_balance['hydrogen']['hydrogen compressor'][step], self.power_balance['gas']['hydrogen compressor'][step] = self.technologies['hydrogen compressor'].use(step,c.timestep,self.power_balance['hydrogen']['electrolyzer'][step],storable_hydrogen) # hydrogen compressed by the compressor (+) and heat requested to make it work expressed as heating water need (-) 
                        pb['gas'] += self.power_balance['gas']['hydrogen compressor'][step]
                        #pb['hydrogen']=...self.power_balance['hydrogen']['hydrogen compressor'][step]?? come ne tengo conto di quanto comprimo? in linea teorica ne dovrei sempre comprimere esattamente quanto me ne entra perchè il controllo sullo sotrable hydrogen lho gia fatto nell'elettrolizzatore'
            
            if tech_name == 'mechanical compressor':   
                if 'HPH tank' not in self.system:
                    if 'O2 tank' not in self.system:
                        if "electricity grid" in self.system and self.system["electricity grid"]["draw"]:
                            if 'hydrogen demand' in self.system:
                                massflow = max(0,self.power_balance['hydrogen']['electrolyzer'][step] + self.power_balance['hydrogen']['demand'][step])
                            else:
                                massflow = max(0,self.power_balance['hydrogen']['electrolyzer'][step])
                            self.power_balance['hydrogen']['mechanical compressor'][step], \
                            self.power_balance['electricity']['mechanical compressor'][step] = self.technologies['mechanical compressor'].use(step,massflowrate= massflow)[:2] # hydrogen compressed by the compressor (+) and electricity consumption (-) 
                            
                            pb['electricity']   += self.power_balance['electricity']['mechanical compressor'][step]
                    
                        elif "electricity grid" not in self.system or self.system["electricity grid"]["draw"] == False:   # if the system is configurated as fully off-grid, relying only on RES production
                            if self.power_balance['hydrogen']['electrolyzer'][step] > 0 :  # if hydrogen has been produced by the electrolyzer and electricity is available in the system
                                if 'hydrogen demand' in self.system:
                                    demand = self.power_balance['hydrogen']['demand'][step] # hydrogen demand at timestep h
                                    massflow = max(0,self.power_balance['hydrogen']['electrolyzer'][step] + demand )  # hydrogen mass flow rate to be compressed and stored in hydrogen tank
                                else:           # no hydrogen demand
                                    demand = 0  # hydrogen demand at timestep h
                                    massflow = self.power_balance['hydrogen']['electrolyzer'][step]  # in case no hydrogen demand is present all produced hydrogen is compressed and flows through the hydrogne tank
                                a = self.technologies['mechanical compressor'].use(step,massflowrate= massflow)[1] # [kW] compressor energy consumption for a certain h2 mass flow rate
                                if abs(a) <=pb['electricity']:     # there is enough renewable electricity to power the compressor 
                                    self.power_balance['hydrogen']['mechanical compressor'][step],    \
                                    self.power_balance['electricity']['mechanical compressor'][step], \
                                    self.power_balance['cooling water']['mechanical compressor'][step]= self.technologies['mechanical compressor'].use(step,massflowrate= massflow) # hydrogen compressed by the compressor (+) and electricity consumption (-) 
                                    
                                    pb['electricity']   += self.power_balance['electricity']['mechanical compressor'][step]
                                    
                                elif abs(a) >pb['electricity']:    # if available electricity in the system is not enough to power the compression system - enter the loop to reallocate the energy among the components
                                    a1  = 1     # % of available electricity fed to the electrolyzer
                                    a11 = 0     # % of available electricity fed to the compressor
                                    en  =pb['electricity'] + abs(self.power_balance['electricity']['electrolyzer'][step]) # [kW] electric energy available at time h before entering the electorlyzer
                                    el  = self.power_balance['electricity']['electrolyzer'][step]
                                    hy  = self.power_balance['hydrogen']['electrolyzer'][step]
                                    ox  = self.power_balance['oxygen']['electrolyzer'][step]
                                    wa  = self.power_balance['water']['electrolyzer'][step]
                                    
                                    # Iteration parameters
                                    i   = 0             # initializing iteration count
                                    maxiter = 10000     # max number of iterations allowed
                                    abs_err = 0.00001   # absolute error allowed
                                    
                                    while a1 >= 0:       # while loop necessary to iterate in the redistribution of renewable electricity to satisfy both electrolyzer and compressor demand
                                        hydrogen_ele,  \
                                        electricity_ele = self.technologies['electrolyzer'].use(step,storable_hydrogen=producible_hyd,e=a1*en)[:2]  # [kg] of produced H2 and [kW] of consumed electricity for the given energy input  
                                        massflow = max(0, hydrogen_ele + demand)
                                        a = -self.technologies['mechanical compressor'].use(step,massflowrate= massflow)[1] # [kW] compressor energy consumption for a certain h2 mass flow rate
                                        b1 = a/en
                                        a11 = 1-b1
                                        i += 1      # updating iteration count
                                    
                                        if abs(a1-a11) < abs_err or i > maxiter:    # strict tolerance for convergence 
                                            break
                                        else: 
                                            a1=a11  
                                            
                                    # Electorlyzer balances update and overwriting
                                    self.power_balance['hydrogen']['electrolyzer'][step],   \
                                    self.power_balance['electricity']['electrolyzer'][step],\
                                    self.power_balance['oxygen']['electrolyzer'][step],     \
                                    self.power_balance['water']['electrolyzer'][step]        = self.technologies['electrolyzer'].use(step,storable_hydrogen=producible_hyd,e=a1*en)      # [:2] # hydrogen supplied by electrolyzer(+) # electricity absorbed by the electorlyzer(-) 
                                    
                                    pb['hydrogen']      += self.power_balance['hydrogen']['electrolyzer'][step]    - hy
                                    pb['electricity']   += self.power_balance['electricity']['electrolyzer'][step] - el
                                    pb['oxygen']        += self.power_balance['oxygen']['electrolyzer'][step]      - ox
                                    pb['water']         += self.power_balance['water']['electrolyzer'][step]       + wa
    
                                    # Compressor balances update and overwriting
                                    self.power_balance['hydrogen']['mechanical compressor'][step],    \
                                    self.power_balance['electricity']['mechanical compressor'][step], \
                                    self.power_balance['cooling water']['mechanical compressor'][step]   = self.technologies['mechanical compressor'].use(step,massflowrate= massflow) # hydrogen compressed by the compressor (+) and electricity consumption (-) 
    
                                    pb['electricity']   += self.power_balance['electricity']['mechanical compressor'][step]
    
                                else:  # if no hydrogen has been produced at time h
                                    self.power_balance['hydrogen']['mechanical compressor'][step]     = 0
                                    self.power_balance['electricity']['mechanical compressor'][step]  = 0
                                    
                    elif 'O2 tank' in self.system:   # simplified approach for oxygen compression. To be updated
                        massflow_tot = (self.power_balance['hydrogen']['electrolyzer'][step])+(self.power_balance['oxygen']['electrolyzer'][step])
                        if "electricity grid" in self.system and self.system["electricity grid"]["draw"]:
                            self.power_balance['hydrogen']['mechanical compressor'][step], \
                            self.power_balance['electricity']['mechanical compressor'][step] = self.technologies['mechanical compressor'].use(step,massflowrate = massflow_tot )[:2] # hydrogen compressed by the compressor (+) and electricity consumption (-) 
                            
                            pb['electricity']   += self.power_balance['electricity']['mechanical compressor'][step]
                    
                        elif "electricity grid" not in self.system or self.system["electricity grid"]["draw"] == False:   # if the system is configurated as fully off-grid, relying only on RES production
                            if self.power_balance['hydrogen']['electrolyzer'][step] > 0 :  # if hydrogen has been produced by the electrolyzer and electricity is available in the system
                                a = self.technologies['mechanical compressor'].use(step,massflowrate = massflow_tot)[1] # [kW] compressor energy consumption for a certain h2 mass flow rate
                                if abs(a) <=pb['electricity']:     # there is enough renewable electricity to power the compressor 
                                    self.power_balance['hydrogen']['mechanical compressor'][step],    \
                                    self.power_balance['electricity']['mechanical compressor'][step], \
                                    self.power_balance['cooling water']['mechanical compressor'][step]= self.technologies['mechanical compressor'].use(step,massflowrate= massflow_tot) # hydrogen compressed by the compressor (+) and electricity consumption (-) 
                                    
                                    pb['electricity']   += self.power_balance['electricity']['mechanical compressor'][step]
                                    
                                elif abs(a) >pb['electricity']:    # if available electricity in the system is not enough to power the compression system - enter the loop to reallocate the energy among the components
                                    a1  = 1     # % of available electricity fed to the electrolyzer
                                    a11 = 0     # % of available electricity fed to the compressor
                                    en  =pb['electricity'] + abs(self.power_balance['electricity']['electrolyzer'][step]) # [kW] electric energy available at time h before entering the electorlyzer
                                    el  = self.power_balance['electricity']['electrolyzer'][step]
                                    hy  = self.power_balance['hydrogen']['electrolyzer'][step]
                                    ox  = self.power_balance['oxygen']['electrolyzer'][step]
                                    wa  = self.power_balance['water']['electrolyzer'][step]
                                    
                                    # Iteration parameters
                                    i   = 0             # initializing iteration count
                                    maxiter = 10000     # max number of iterations allowed
                                    abs_err = 0.00001   # absolute error allowed
                                    
                                    while a1 >= 0:       # while loop necessary to iterate in the redistribution of renewable electricity to satisfy both electrolyzer and compressor demand
                                        hydrogen_ele,  \
                                        electricity_ele = self.technologies['electrolyzer'].use(step,a1*en,producible_hyd)[:2]  # [kg] of produced H2 and [kW] of consumed electricity for the given energy input  
                                        a = -self.technologies['mechanical compressor'].use(step,massflowrate= hydrogen_ele + hydrogen_ele*7.93)[1] # [kW] compressor energy consumption for a certain h2 mass flow rate
                                        b1 = a/en
                                        a11 = 1-b1
                                        i += 1      # updating iteration count
                                    
                                        if abs(a1-a11) < abs_err or i > maxiter:    # strict tolerance for convergence 
                                            break
                                        else: 
                                            a1=a11    
                                    
                                    # Electorlyzer balances update and overwriting
                                    self.power_balance['hydrogen']['electrolyzer'][step],   \
                                    self.power_balance['electricity']['electrolyzer'][step],\
                                    self.power_balance['oxygen']['electrolyzer'][step],     \
                                    self.power_balance['water']['electrolyzer'][step]        = self.technologies['electrolyzer'].use(step,a1*en,producible_hyd)      # [:2] # hydrogen supplied by electrolyzer(+) # electricity absorbed by the electorlyzer(-) 
                                    
                                    pb['hydrogen']      += self.power_balance['hydrogen']['electrolyzer'][step]    - hy
                                    pb['electricity']   += self.power_balance['electricity']['electrolyzer'][step] - el
                                    pb['oxygen']        += self.power_balance['oxygen']['electrolyzer'][step]      - ox
                                    pb['water']         += self.power_balance['water']['electrolyzer'][step]       + wa
    
                                    # Compressor balances update and overwriting
                                    self.power_balance['hydrogen']['mechanical compressor'][step],    \
                                    self.power_balance['electricity']['mechanical compressor'][step], \
                                    self.power_balance['cooling water']['mechanical compressor'][step]   = self.technologies['mechanical compressor'].use(step,massflowrate= self.power_balance['hydrogen']['electrolyzer'][step] + self.power_balance['oxygen']['electrolyzer'][step]) # hydrogen compressed by the compressor (+) and electricity consumption (-) 
    
                                    pb['electricity']   += self.power_balance['electricity']['mechanical compressor'][step]
    
                                
                                else:  # if no hydrogen has been produced at time h
                                    self.power_balance['hydrogen']['mechanical compressor'][step]     = 0
                                    self.power_balance['electricity']['mechanical compressor'][step]  = 0
                        
                             
                
                if 'H tank' in self.system and 'HPH tank' in self.system:
                    # self.power_balance['hydrogen']['H tank'][step] = self.technologies['H tank'].use(h,pb['hydrogen'])
                    #pb['hydrogen'] += self.power_balance['hydrogen']['H tank'][step]
                    available_hyd_lp = self.technologies['H tank'].LOC[step] + self.technologies['H tank'].max_capacity - self.technologies['H tank'].used_capacity
                    storable_hydrogen_hp = self.technologies['HPH tank'].max_capacity-self.technologies['HPH tank'].LOC[step]
                    
                    if self.technologies['HPH tank'].LOC[step] == self.technologies['HPH tank'].max_capacity:  # if High-Pressure-Tank is full 
                        self.power_balance['HP hydrogen']['mechanical compressor'][step]      = 0     
                        self.power_balance['electricity']['mechanical compressor'][step]      = 0    
                        self.power_balance['cooling water']['mechanical compressor'][step]    = 0
                        self.power_balance['hydrogen']['mechanical compressor'][step]         = 0
                        
                        pb['HP hydrogen']    += 0
                        pb['electricity']    += 0   # compressor not working
                        pb['hydrogen']       += 0   # compressor not working
                    else:  # if there is enough room available in the High Pressure Tank, the compressor is activated
                        self.power_balance['HP hydrogen']['mechanical compressor'][step],     \
                        self.power_balance['electricity']['mechanical compressor'][step],     \
                        self.power_balance['cooling water']['mechanical compressor'][step]    = self.technologies['mechanical compressor'].use(step, available_hyd_lp=available_hyd_lp ,storable_hydrogen_hp=storable_hydrogen_hp) # hydrogen supplied by H tank (+) and electricity absorbed(-) 
                        self.power_balance['hydrogen']['mechanical compressor'][step] = - self.power_balance['HP hydrogen']['mechanical compressor'][step]
                        
                        pb['HP hydrogen'] += self.power_balance['HP hydrogen']['mechanical compressor'][step]
                        pb['hydrogen']    += self.power_balance['hydrogen']['mechanical compressor'][step]
                        pb['electricity'] += self.power_balance['electricity']['mechanical compressor'][step]

            if tech_name == 'fuel cell':
                if pb['electricity'] < 0: #? this condition must be solved if you want to produce electricity to be fed into the gird
                    if "hydrogen grid" in self.system and self.system["hydrogen grid"]["draw"]: # hydrogen can be withdranw from an hydrogen grid
                        available_hyd = 9999999999999999999 
                    elif 'H tank' in self.system:   # only hydrogen inside H tank can be used
                        available_hyd = self.technologies['H tank'].LOC[step] + self.technologies['H tank'].max_capacity - self.technologies['H tank'].used_capacity                                                                  
                    else:
                        available_hyd = max(0,pb['hydrogen']) # hydrogen is produced by an electrolyzer with a higher priority than fc
                    if available_hyd > 0:
                        use = self.technologies['fuel cell'].use(step,self.P2E,pb['electricity'],available_hyd)     # saving fuel cell working parameters for the current timeframe
                        self.power_balance['hydrogen']['fuel cell'][step] =    use[0] # hydrogen absorbed by fuel cell(-)
                        self.power_balance['electricity']['fuel cell'][step] = use[1] # electricity supplied(+) 
        
                        if use[2] < -pb['heating water']: #all of the heat producted by FC is used      
                            self.power_balance['heating water']['fuel cell'][step]=use[2] # heat loss of fuel cell
                        else:
                            self.power_balance['heating water']['fuel cell'][step]=-pb['heating water'] # heat loss of fuel cell- demand

                        pb['hydrogen'] += self.power_balance['hydrogen']['fuel cell'][step]
                        pb['electricity'] += self.power_balance['electricity']['fuel cell'][step]
                        pb['heating water'] += self.power_balance['heating water']['fuel cell'][step] 
                    
            if tech_name == 'boiler_h2':                                                                                                                                   
                if "hydrogen grid" in self.system and self.system["hydrogen grid"]["draw"]: # hydrogen can be withdranw from an hydrogen grid
                    available_hyd = 9999999999999999999 
                elif 'H tank' in self.system:   # only hydrogen inside H tank can be used
                    available_hyd = self.technologies['H tank'].LOC[step] + self.technologies['H tank'].max_capacity - self.technologies['H tank'].used_capacity +pb['hydrogen']
                                         
                self.power_balance['hydrogen']['boiler_h2'][step], self.power_balance['gas']['boiler_h2'][step] = self.technologies['boiler_h2'].use(pb['gas'],available_hyd,1)[1:3] #h2 consumed from boiler_h2 and heat produced by boiler_h2
                pb['hydrogen'] += self.power_balance['hydrogen']['boiler_h2'][step] # hydrogen balance update: - hydrogen consumed by boiler_h2
                pb['gas'] += self.power_balance['gas']['boiler_h2'][step] # heat balance update: + heat produced by boiler_h2
           
            #!!!ANDREA HA MESSO QUESTO COME PROCESSO NEL LOCATION, QUELLO SOTTO A COSA é DOVUTO   WIP to be modified by Andrea                                                                                                                                                          
            if tech_name == 'boiler_h2': 
                if pb['electricity'] < 0: #? this condition must be solved if you want to produce electricity to fed into the gird
                    if "hydrogen grid" in self.system and self.system["hydrogen grid"]["draw"]: # hydrogen can be withdranw from an hydrogen grid
                        available_hyd = 9999999999999999999 
                    elif 'H tank' in self.system:   # only hydrogen inside H tank can be used
                        available_hyd = self.technologies['H tank'].LOC[step] + self.technologies['H tank'].max_capacity - self.technologies['H tank'].used_capacity                                                                  
                    else:
                        available_hyd = max(0,pb['hydrogen']) # hydrogen is produced by an electrolyzer which have higher priority than boiler_h2
                    if available_hyd > 0:
                        self.power_balance['hydrogen']['boiler_h2'][step], self.power_balance['heating water']['boiler_h2'][step] = self.technologies['boiler_h2'].use(pb['heating water'],available_hyd,1)[1:3] #h2 consumed from boiler_h2 and heat produced by boiler_h2
                        pb['hydrogen'] += self.power_balance['hydrogen']['boiler_h2'][step] # hydrogen balance update: - hydrogen consumed by boiler_h2
                        pb['heating water'] += self.power_balance['heating water']['boiler_h2'][step] # heat balance update: + heat produced by boiler_h2
                    
            # if tech_name in ['H tank','HPH tank']:     #versione buona 
            #     if self.system['hydrogen demand']['strategy'] != 'supply-led':
            #         self.power_balance[self.tank_stream[tech_name]][tech_name][step] = self.technologies[tech_name].use(h,pb[self.tank_stream[tech_name]])
            #        pb[self.tank_stream[tech_name]] += self.power_balance[self.tank_stream[tech_name]][tech_name][step]
            #     elif self.system['hydrogen demand']['strategy'] == 'supply-led' and h == (c.timestep_number - 1):
            #         prod = self.power_balance['hydrogen']['electrolyzer']
            #         for h in range(c.timestep_number):
            #             self.power_balance[self.tank_stream[tech_name]][tech_name][step] = self.technologies[tech_name].use(h,prod[step],constant_demand=self.constant_flow )                        
            #     else:
            #         pass
        
            if tech_name == 'H tank':
                if 'HPH tank' not in self.system and ('hydrogen demand' in self.system or 'HP hydrogen demand' in self.system):
                    if self.system[self.hydrogen_demand+' demand']['strategy'] == 'demand-led':
                        self.power_balance['hydrogen']['H tank'][step] = self.technologies['H tank'].use(step,pb['hydrogen'])
                        pb['hydrogen'] += self.power_balance['hydrogen']['H tank'][step]
                    elif self.system[self.hydrogen_demand+' demand']['strategy'] == 'supply-led' and step == (c.timestep_number - 1):
                        prod = self.power_balance['hydrogen']['electrolyzer']
                        for step in range(c.timestep_number):
                            self.power_balance['hydrogen']['H tank'][step] = self.technologies['H tank'].use(step,prod[step],constant_demand=self.constant_flow )                        
                    # else:
                        # pass
                else:
                    self.power_balance['hydrogen']['H tank'][step] = self.technologies['H tank'].use(step,pb['hydrogen'])
                    pb['hydrogen'] += self.power_balance['hydrogen']['H tank'][step]
            
            if tech_name == 'HPH tank':
                self.power_balance['HP hydrogen']['HPH tank'][step] = self.technologies['HPH tank'].use(step,pb['HP hydrogen'])
                pb['HP hydrogen'] += self.power_balance['HP hydrogen']['HPH tank'][step]
                
            if tech_name == 'O2 tank':
                if 'oxygen demand' in self.system and self.system['oxygen demand']['strategy'] != 'supply-led':
                    self.power_balance['oxygen']['O2 tank'][step] = self.technologies['O2 tank'].use(step,pb['oxygen'])
                    pb['oxygen'] += self.power_balance['oxygen']['O2 tank'][step]
                elif self.system['hydrogen demand']['strategy'] == 'supply-led' and step == (c.timestep_number - 1):
                    self.technologies['O2 tank'].sizing(self.technologies['H tank'].max_capacity)
                else:
                    pass
            
            if tech_name == 'inverter':
                self.power_balance['electricity']['inverter'][step] = self.technologies['inverter'].use(step,pb['electricity']) # electricity lost in conversion by the inverter
                pb['electricity'] += self.power_balance['electricity']['inverter'][step] # electricity balance update: - electricity lost in conversion by the invertert

            ### demand and grid   
            for carrier in pb: # for each energy carrier
                if tech_name == f"{carrier} demand":                
                    pb[carrier] += self.power_balance[carrier]['demand'][step]    # power balance update: energy demand(-)  
                    
                if tech_name == f"{carrier} grid":
                    if pb[carrier] > 0 and self.system[f"{carrier} grid"]['feed'] or pb[carrier] < 0 and self.system[f"{carrier} grid"]['draw']:
                        self.power_balance[carrier]['grid'][step] = -pb[carrier] # energy from grid(+) or into grid(-) 
                        pb[carrier] += self.power_balance[carrier]['grid'][step]  # electricity balance update      

#%%            
        ### Global check on power balances at the end of every timestep
        for carrier in pb:
            if carrier == 'heating water':
                continue
            tol = 0.0001  # [-] tolerance on error
            if pb[carrier] != 0:
                maxvalues = []
                for arr in self.power_balance[carrier]:
                    maxvalues.append(max(self.power_balance[carrier][arr]))
                m = max(maxvalues)
                if abs(pb[carrier]) > abs(m*tol):
                    if pb[carrier] >0:  sign = 'positive'
                    else:               sign = 'negative'
                    raise ValueError(f'Warning: {carrier} balance at the end of timestep {step} shows {sign} value of {round(pb[carrier],2)} \n\
                    It means there is an overproduction not fed to grid or demand is not satisfied.\n\
                    Options to fix the problem: \n\
                        (a) - Include {carrier} grid[\'draw\']: true if negative or {carrier} grid[\'feed\']: true if positive in studycase.json \n\
                        (b) - Vary components size or demand series in studycase.json')
