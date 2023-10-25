import numpy as np
import pickle
import csv
import os
import pandas as pd
import pvlib #https://github.com/pvlib
from core import location
from core import constants as c

class REC:
    
    def __init__(self,structure,general,file_structure,file_general,path):
        """
        Create a Renewable Energy Comunity object composed of several locations (producers, consumers, prosumers)
    
        structure : dictionary (all the inputs are optional)
            'location_1_name': inputs required to create a location object (see Location.py)
            'location_2_name': 
                ...
            'location_n_name':
                
        general : dictionary
            'simulation years': number of years to be simulated
            'timestep': time step in minutes
            'latitude': float
            'longitude': float    
            'UTC time zone': int 0,1,2 [UTC] es. Italy is in UTC+1 time zone EUROPEAN DATABASE
            'DST': boolean, Daily saving time (fusorario)
            'weather': if "TMY" weather database based on typical meteorological year is used
                if "filename.csv" a different database can be used (upload it in input/weather)
                in this case 'latitude' and 'longitude' are ignored
                        
        output : REC object able to:
            simulate the power flows of each present locations .REC_simulation
            record REC power balances (electricity, heat, gas and hydrogen) 
        
        If general.json is the same of the previous simulation, neither the meteorological data nor the PV has to be updated, 
        otherwise they are downloaded from PVgis considering the typical meteorological year.
        
        """

        ##############################################################################################
        ### Check if new input files have to been downloaded from PV gis 
        check = True # Used to check if TMY have to been downloaded from PVgis or the old one can be used
        check_pv = True # Used to check if PV_production series have to been downloaded from PVgis or the old one can be used
        directory = './previous_simulation'
        if not os.path.exists(directory): os.makedirs(directory)

        if os.path.exists(f"previous_simulation/{file_general}.pkl"):
            with open(f"previous_simulation/{file_general}.pkl", 'rb') as f: ps_general = pickle.load(f) # previous simulation general
            par_to_check = ['latitude','longitude','UTC time zone','DST']
            for par in par_to_check:
                if ps_general[par] != general[par]:
                    check = False          
        else:
            check = False 
        if os.path.exists(f"previous_simulation/{file_general}_{file_structure}.pkl"):
            with open(f"previous_simulation/{file_general}_{file_structure}.pkl", 'rb') as f: ps_general = pickle.load(f) # previous simulation general
            par_to_check = ['latitude','longitude','UTC time zone','DST']
            for par in par_to_check:
                if ps_general[par] != general[par]:
                    check_pv = False          
        else:
            check_pv = False
        self.weather = self.weather_generation(general,path,check,file_general) # check if metereological data have to been downloaded from PVgis or has already been done in a previous simulation
        if check == False:
            with open(f"previous_simulation/{file_general}.pkl", 'wb') as f: pickle.dump(general, f)
        if check_pv == False:
            with open(f"previous_simulation/{file_general}_{file_structure}.pkl", 'wb') as f: pickle.dump(general, f)
        ##############################################################################################

        ### Add global variables to constants.py as c, to make them known to the other modules.
        c.timestep = general['timestep'] # timestep [min]
        c.timestep_number = int( general['simulation years']* 365*24*60 / c.timestep ) # number of timestep [#]

        c.simulation_years = general['simulation years']
        c.P2E = c.timestep*60 # conversion factor from kW to kJ or from kg/s to kg
        c.latitude = general['latitude']
        c.longitude = general["longitude"]
        c.UTC = general["UTC time zone"] # int 0,1,2 [UTC] es. Italy is in UTC+1 time zone EUROPEAN DATABASE
        c.DST = general["DST"] # boolean, Daily saving time (fusorario)

        self.locations = {} # initialise REC locations dictionary
        self.power_balance = {'electricity': {}, 'heating water': {}, 'cooling water': {}, 'hydrogen': {}, 'gas': {}, 'process steam': {}} # initialise power balances dictionaries
     
        ### create location objects and add them to the REC locations dictionary
        for location_name in structure: # location_name are the keys of 'structure' dictionary and will be used as keys of REC 'locations' dictionary too
            self.locations[location_name] = location.location(structure[location_name],location_name,path,check_pv,file_structure,file_general) # create location object and add it to REC 'locations' dictionary                
                     

    def REC_power_simulation(self):
        """
        Simulate the REC every hour
        
        output :
            updating location power balances
            updating REC power balances
        """
        
        ### initialise REC electricity balances
        self.power_balance['electricity']['from grid'] = np.zeros(c.timestep_number) # array of electricity withdrawn from the grid from the whole rec
        self.power_balance['electricity']['into grid'] = np.zeros(c.timestep_number) # array of electricity withdrawn from the grid
        self.power_balance['electricity']['collective self consumption'] = np.zeros(c.timestep_number) # array of collective self consumed electricity from the whole rec
        self.count = []
        
        ### simulation core
        for step in range(c.timestep_number): # step to simulate
            for location_name in self.locations: # each locations 
                self.locations[location_name].loc_power_simulation(step,self.weather) # simulate a single location updating its power balances
                
            ### solve electricity grid 
                if 'grid' in self.locations[location_name].power_balance['electricity']:
                    if self.locations[location_name].power_balance['electricity']['grid'][step] < 0:
                        self.power_balance['electricity']['into grid'][step] += self.locations[location_name].power_balance['electricity']['grid'][step] # electricity fed into the grid from the whole rec at step step
                    else:                                                     
                        self.power_balance['electricity']['from grid'][step] += self.locations[location_name].power_balance['electricity']['grid'][step] # electricity withdrawn from the grid the whole rec at step step
                
               
            ###################################################################################################################################
            ### solve smart heatpumps (REC_surplus == True) (only available with timestep == 60)
            HPs_available = []
            TESs_temperature = []            
            if - self.power_balance['electricity']['into grid'][step] - self.power_balance['electricity']['from grid'][step] > 0: # if there is surplus
            
                # find the HPs available
                for location_name in self.locations:
                    if 'heatpump' in self.locations[location_name].technologies:
                            if self.locations[location_name].technologies['heatpump'].REC_surplus:
                                if c.timestep != 60:
                                    raise ValueError("Warning! Heatpumps with strategy REC_suruplus == True only work with timestep == 60 ")
                                if self.locations[location_name].technologies['heatpump'].mode == 1:
                                    if self.locations[location_name].technologies['heatpump'].satisfaction_story[step] in [0,1,2,3]:   
                                        HPs_available.append(location_name)
                                        TESs_temperature.append(self.locations[location_name].technologies['heatpump'].i_TES_t)
                            
                # order locations according to iTES temperature 
                HPs_available = [HPs_available for _,HPs_available in sorted(zip(TESs_temperature,HPs_available))]
             
            while - self.power_balance['electricity']['into grid'][step] - self.power_balance['electricity']['from grid'][step] > 0: # while there is surplus
                surplus = - self.power_balance['electricity']['into grid'][step] - self.power_balance['electricity']['from grid'][step]
                if HPs_available == []:
                        break
                location_name = HPs_available[0]
                HPs_available = HPs_available[1:]
 
                # clean balance
                self.locations[location_name].power_balance['electricity']['demand'][step] += - self.locations[location_name].power_balance['electricity']['heatpump'][step]
                self.locations[location_name].power_balance['electricity']['grid'][step] += self.locations[location_name].power_balance['electricity']['heatpump'][step]
                self.power_balance['electricity']['from grid'][step] += self.locations[location_name].power_balance['electricity']['heatpump'][step] # electricity fed into the grid from the whole rec at hour h

                # resimulate 
                self.locations[location_name].technologies['heatpump'].i_TES_t = self.locations[location_name].technologies['heatpump'].i_TES_story[step]
                self.locations[location_name].power_balance['electricity']['heatpump'][step], self.locations[location_name].power_balance['heat']['heatpump'][step], self.locations[location_name].power_balance['heat']['inertial TES'][step] = self.locations[location_name].technologies['heatpump'].use(self.weather['temp_air'][step],self.locations[location_name].power_balance['heat']['demand'][step],surplus,step) 
  
                # updata balance
                self.locations[location_name].power_balance['electricity']['demand'][step] += self.locations[location_name].power_balance['electricity']['heatpump'][step]
                self.locations[location_name].power_balance['electricity']['grid'][step] += - self.locations[location_name].power_balance['electricity']['heatpump'][step]
                self.power_balance['electricity']['from grid'][step] += - self.locations[location_name].power_balance['electricity']['heatpump'][step] # electricity fed into the grid from the whole rec at hour h
             
            ###################################################################################################################################
            ### calculate collective self consumption and who contributed to it
            self.power_balance['electricity']['collective self consumption'][step] = min(-self.power_balance['electricity']['into grid'][step],self.power_balance['electricity']['from grid'][step]) # calculate REC collective self consumption how regulation establishes      
            
            if self.power_balance['electricity']['collective self consumption'][step] > 0:
                for location_name in self.locations:
                    if self.locations[location_name].power_balance['electricity']['grid'][step] < 0: # contribution as producer
                        self.locations[location_name].power_balance['electricity']['collective self consumption'][step] = - self.power_balance['electricity']['collective self consumption'][step] * self.locations[location_name].power_balance['electricity']['grid'][step] / self.power_balance['electricity']['into grid'][step]
                    else: # contribution as consumer
                        self.locations[location_name].power_balance['electricity']['collective self consumption'][step] = self.power_balance['electricity']['collective self consumption'][step] * self.locations[location_name].power_balance['electricity']['grid'][step] / self.power_balance['electricity']['from grid'][step]

            ###################################################################################################################################
            ### solve smart batteries (only available with timestep == 60)
            for location_name in self.locations:
                
                # battery.collective = 1: 
                # REC tels to location how mutch electricity can be absorbed or supplied by battery every hour, without decreasing the collective-self-consumption
                
                if 'battery' in self.locations[location_name].technologies and self.locations[location_name].technologies['battery'].collective == 1:
                    
                    if c.timestep != 60:
                        raise ValueError("Warning! Batteries with strategy collective == 1 only work with timestep == 60 ")

                    # how much energy can be absorbed or supplied by the batteries cause it's not usefull for collective-self-consumption
                    E = - self.locations[location_name].power_balance['electricity']['grid'][step] + self.locations[location_name].power_balance['electricity']['collective self consumption'][step]
                      
                    self.locations[location_name].power_balance['electricity']['battery'][step] = self.locations[location_name].technologies['battery'].use(step,E) # electricity absorbed(-) by battery
                    self.locations[location_name].power_balance['electricity']['grid'][step] += - self.locations[location_name].power_balance['electricity']['battery'][step] # update grid balance (locatiom)
                  
                    if self.locations[location_name].power_balance['electricity']['battery'][step] < 0:
                        self.power_balance['electricity']['into grid'][step] += - self.locations[location_name].power_balance['electricity']['battery'][step] # update grid balance (rec)
                    else:
                        self.power_balance['electricity']['from grid'][step] += - self.locations[location_name].power_balance['electricity']['battery'][step] # update grid balance (rec)
                

    def save(self,simulation_name,f,sep=';',dec=','):
        """
        Save REC and each location power balances
        
        simulationa_name : str 
        f: 'csv' or 'pkl'
        using sep and dec you can choose the separator and decima of the .csv format
        
        output: 
            balances/simulation_name.pkl
            LOC/simulation_name.pkl
        """
        
        balances = {}
        LOC = {}
        ageing = {}
        electrolyzer = {}                 
        balances['REC'] = self.power_balance
        parameters = {}        
        tech_cost = {}
        
        for location_name in self.locations:
            balances[location_name] = self.locations[location_name].power_balance
            # parameters[location_name] = self.locations[location_name].tech_param                                                                                
            
            LOC[location_name] = {}
            ageing[location_name] = {}
            electrolyzer[location_name] = {}    
            parameters[location_name] = {}                            
            
            tech_name = 'battery'
            if tech_name in self.locations[location_name].technologies:
                LOC[location_name][tech_name] = self.locations[location_name].technologies[tech_name].LOC
                if self.locations[location_name].technologies[tech_name].ageing:
                    ageing[location_name][tech_name] = [self.locations[location_name].technologies[tech_name].replacements,self.locations[location_name].technologies[tech_name].ageing_history]
                
            tech_name = 'H tank'
            if tech_name in self.locations[location_name].technologies:
                LOC[location_name][tech_name] = self.locations[location_name].technologies[tech_name].LOC
                
            tech_name = 'heatpump'
            if tech_name in self.locations[location_name].technologies:
                LOC[location_name]['inertial TES'] = self.locations[location_name].technologies[tech_name].i_TES_story
        
            tech_name = 'electrolyzer'
            if tech_name in self.locations[location_name].technologies:
                parameters[location_name][tech_name] = {}      
                if self.locations[location_name].technologies['electrolyzer'].model == 'PEM General':
                    parameters[location_name][tech_name]['efficiency'] = self.locations[location_name].technologies[tech_name].EFF
                    parameters[location_name][tech_name]['hourly capacity factor'] =   ((-balances[location_name]['electricity'][tech_name])/             \
                                                                                (self.locations[location_name].technologies[tech_name].MaxPowerStack))*100
                    parameters[location_name][tech_name]['capacity factor'] =   ((-balances[location_name]['electricity'][tech_name].sum())/             \
                                                                                (self.locations[location_name].technologies[tech_name].MaxPowerStack*self.timestep_number))*100
                    
            tech_name = 'fuel cell'
            if tech_name in self.locations[location_name].technologies:
                if self.locations[location_name].technologies['fuel cell'].model == 'PEM General':
                    parameters[location_name][tech_name] = {}      
                    parameters[location_name][tech_name]['cell voltage'] = self.locations[location_name].technologies[tech_name].VOLT
                    parameters[location_name][tech_name]['current density'] = self.locations[location_name].technologies[tech_name].CURR_DENS
                if self.locations[location_name].technologies['fuel cell'].model == 'SOFC':
                    parameters[location_name][tech_name] = {}
                    parameters[location_name][tech_name]['efficiency'] = self.locations[location_name].technologies[tech_name].EFF
                    parameters[location_name][tech_name]['efficiency last module'] = self.locations[location_name].technologies[tech_name].EFF_last_module
                    parameters[location_name][tech_name]['number modules used'] = self.locations[location_name].technologies[tech_name].n_modules_used
            
            tech_name = 'hydrogen compressor'
            if tech_name in self.locations[location_name].technologies:
                parameters[location_name][tech_name] = {}      
                parameters[location_name][tech_name]['compressor number used'] = self.locations[location_name].technologies[tech_name].n_compressors_used                                                                                           
            
            tech_cost[location_name] = {}
            for tech_name in self.locations[location_name].technologies:
                if self.locations[location_name].technologies[tech_name].cost and (not hasattr(self.locations[location_name].technologies[tech_name], 'property') or self.locations[location_name].technologies[tech_name].property) :
                    tech_cost[location_name][tech_name] = self.locations[location_name].technologies[tech_name].cost
        
        
        if f == 'pkl':
            directory = './results'
            if not os.path.exists(directory): os.makedirs(directory)
            directory = './results/pkl'
            if not os.path.exists(directory): os.makedirs(directory)
            with open('results/pkl/balances_'+simulation_name+".pkl", 'wb') as f: pickle.dump(balances, f)
            with open('results/pkl/tech_params_'+simulation_name+".pkl", 'wb') as f: pickle.dump(parameters, f)
            with open('results/pkl/LOC_'+simulation_name+".pkl", 'wb') as f: pickle.dump(LOC, f)             
            with open('results/pkl/ageing_'+simulation_name+".pkl", 'wb') as f: pickle.dump(ageing, f)   
            with open('results/pkl/tech_cost_'+simulation_name+".pkl", 'wb') as f: pickle.dump(tech_cost, f)   
            
        if f == 'csv':
            directory = './results'
            if not os.path.exists(directory): os.makedirs(directory)
            directory = './results/csv'
            if not os.path.exists(directory): os.makedirs(directory)
            balances_csv = {}
            
            for loc_name in balances:
                for carrier in balances[loc_name]:
                    for tech in balances[loc_name][carrier]:
                        key = f"{loc_name} - {carrier} - {tech}"
                        balances_csv[key] = balances[loc_name][carrier][tech]
            
            df = pd.DataFrame(balances_csv)
            df = df.round(4)
            df.to_csv('results/csv/balances_'+simulation_name+'.csv',index=False,sep=sep,decimal=dec)
            
        
    def weather_generation(self,general,path,check,file_general):
        """
        
        If the meteorological data have not already been downloaded and saved in a previous simulation, 
        then they are downloaded from PVgis considering the typical meteorological year.

        Parameters
        ----------
        general : see REC __init___

        Returns
        -------
        previous_simulation/files.csv

        """                        

        if check and os.path.exists(f"{path}/weather/TMY_{file_general}.csv"): # if the prevoius weather series can be used
            weather = pd.read_csv(f"{path}/weather/TMY_{file_general}.csv")
        
        else: # if new weather data must be downoladed from PV gis
            print('Downolading typical metereological year data from PVGIS for '+file_general)   
                            
            latitude = general['latitude']
            longitude = general['longitude']

            weather = pvlib.iotools.get_pvgis_tmy(latitude, longitude, map_variables=True)[0]
            
            # time zone correction
            if general['UTC time zone'] > 0:
                refindex = weather.index
                refindex = refindex.shift(general['UTC time zone']*60,'T')
                weather.index = refindex
                
                we2 = pd.DataFrame(data=weather[-general['UTC time zone']:], index=None, columns=weather.columns)
                weather = weather[:-general['UTC time zone']]
                
                reindex = weather.index[:general['UTC time zone']]
                reindex = reindex.shift(-general['UTC time zone']*60,'T')
                we2.index = reindex   
                
                weather = pd.concat([we2,weather])
                weather['Local time']=weather.index
                weather.set_index('Local time',inplace=True)

            # Daily saving time (DST) correction 
            # Is CEST (Central European Summertime) observed? if yes it means that State is applying DST
            # DST lasts between last sunday of march at 00:00:00+UTC+1 and last sunday of october at 00:00:00+UTC+2
            # For example in Italy DST in 2022 starts in March 27th at 02:00:00 and finishes in October 30th at 03:00:00
            if general['DST']==True:

                zzz_in = weather[weather.index.month==3]
                zzz_in = zzz_in[zzz_in.index.weekday==6]
                zzz_in = zzz_in[zzz_in.index.hour==1+general['UTC time zone']]
                zzz_in = pd.Series(zzz_in.index).unique()[-1]
              
                zzz_end = weather[weather.index.month==10]
                zzz_end = zzz_end[zzz_end.index.weekday==6]
                zzz_end = zzz_end[zzz_end.index.hour==1+general['UTC time zone']]
                zzz_end = pd.Series(zzz_end.index).unique()[-1]
                
                weather.loc[zzz_in:zzz_end] = weather.loc[zzz_in:zzz_end].shift(60,'T')
                weather = weather.interpolate(method='linear')

                weather['Local time - DST'] = weather.index
                weather.set_index('Local time - DST',inplace=True)  
                                             
            weather.to_csv(f"{path}/weather/TMY_{file_general}.csv") 
            
        weather = pd.DataFrame(np.repeat(weather.values, 4, axis=0), columns=weather.columns)

        return(weather)
   
    def tech_cost(self,tech_cost):
        for location_name in self.locations:
            for tech_name in self.locations[location_name].technologies:
                self.locations[location_name].technologies[tech_name].tech_cost(tech_cost[tech_name])
                