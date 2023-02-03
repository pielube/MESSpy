import numpy as np
import pickle
import random
import os
import pandas as pd
import pvlib #https://github.com/pvlib
from location import location

class REC:
    
    def __init__(self,structure,general,path,rec_name):
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
            'time zone': int 0,1,2 [UTC] es. Italy is in UTC+1 time zone EUROPEAN DATABASE
            'weather': if "TMY" weather database based on typical meteorological year is used
                if "filename.csv" a different database can be used (upload it in input/weather)
                in this case 'latitude' and 'longitude' are ignored
                        
        output : REC object able to:
            simulate the energy flows of each present locations .REC_simulation
            record REC energy balances .energy_balance (electricity, heat, gas and hydrogen) 
        """
        
        
        """
        
        If general.json is the same of the previous simulation, neither the meteorological data nor the PV has to be updated, 
        otherwise they are downloaded from PVgis considering the typical meteorological year.
        
        """
        
        check = True # True if no general parameters are changed from the old simulation
        
        directory = './previous_simulation'
        if not os.path.exists(directory):
            os.makedirs(directory)
           
        if os.path.exists('previous_simulation/general_'+rec_name+'.pkl'):
            with open('previous_simulation/general_'+rec_name+'.pkl', 'rb') as f:
                ps_general = pickle.load(f) # previous simulation general
            par_to_check = ['latitude','longitude','UTC time zone','DST']
            for par in par_to_check:
                if ps_general[par] != general[par]:
                    check = False  
        else:
            check = False


        self.weather = self.weather_generation(general,path,check) # check if metereological data have to been downloaded from PVgis or has already been done in a previous simulation

        self.locations = {} # initialise REC locations dictionary
        self.energy_balance = {'electricity': {}, 'heating water': {}, 'cooling water': {}, 'hydrogen': {}, 'gas': {}, 'process steam': {}} # initialise energy balances dictionaries
                                                                                                        
        
        self.simulation_hours = int(general['simulation years']*8760) # hourly timestep  
        
        ### create location objects and add them to the REC locations dictionary
        for location_name in structure: # location_name are the keys of 'structure' dictionary and will be used as keys of REC 'locations' dictionary too
            self.locations[location_name] = location(structure[location_name],general,location_name,path,check,rec_name) # create location object and add it to REC 'locations' dictionary                
                     
            
        if check == False:
            with open('previous_simulation/general_'+rec_name+'.pkl', 'wb') as f:
                pickle.dump(general, f)
            
            
            
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
        self.count = []
        
        ### simulation core
        for h in range(self.simulation_hours): # h: hour to simulate from 0 to simulation_hours 
            for location_name in self.locations: # each locations 
                self.locations[location_name].loc_energy_simulation(h,self.weather) # simulate a single location updating its energy balances
                
            ### solve electricity grid 
                if 'grid' in self.locations[location_name].energy_balance['electricity']:
                    if self.locations[location_name].energy_balance['electricity']['grid'][h] < 0:
                        self.energy_balance['electricity']['into grid'][h] += self.locations[location_name].energy_balance['electricity']['grid'][h] # electricity fed into the grid from the whole rec at hour h
                    else:                                                     
                        self.energy_balance['electricity']['from grid'][h] += self.locations[location_name].energy_balance['electricity']['grid'][h] # electricity withdrawn from the grid the whole rec at hour h
                
                
            ### solve smart heatpumps (REC_surplus == True)
            HPs_available = []
            TESs_temperature = []            
            if - self.energy_balance['electricity']['into grid'][h] - self.energy_balance['electricity']['from grid'][h] > 0: # if there is surplus
            
                # find the HPs available
                for location_name in self.locations:
                    if 'heatpump' in self.locations[location_name].technologies:
                            if self.locations[location_name].technologies['heatpump'].REC_surplus:
                                if self.locations[location_name].technologies['heatpump'].mode == 1:
                                    if self.locations[location_name].technologies['heatpump'].satisfaction_story[h] in [0,1,2,3]:   
                                        HPs_available.append(location_name)
                                        TESs_temperature.append(self.locations[location_name].technologies['heatpump'].i_TES_t)
                            
                # order locations according to iTES temperature 
                HPs_available = [HPs_available for _,HPs_available in sorted(zip(TESs_temperature,HPs_available))]
             
            while - self.energy_balance['electricity']['into grid'][h] - self.energy_balance['electricity']['from grid'][h] > 0: # while there is surplus
                surplus = - self.energy_balance['electricity']['into grid'][h] - self.energy_balance['electricity']['from grid'][h]
                if HPs_available == []:
                        break
                location_name = HPs_available[0]
                HPs_available = HPs_available[1:]
 
                # clean balance
                self.locations[location_name].energy_balance['electricity']['demand'][h] += - self.locations[location_name].energy_balance['electricity']['heatpump'][h]
                self.locations[location_name].energy_balance['electricity']['grid'][h] += self.locations[location_name].energy_balance['electricity']['heatpump'][h]
                self.energy_balance['electricity']['from grid'][h] += self.locations[location_name].energy_balance['electricity']['heatpump'][h] # electricity fed into the grid from the whole rec at hour h

                # resimulate 
                self.locations[location_name].technologies['heatpump'].i_TES_t = self.locations[location_name].technologies['heatpump'].i_TES_story[h]
                self.locations[location_name].energy_balance['electricity']['heatpump'][h], self.locations[location_name].energy_balance['heat']['heatpump'][h], self.locations[location_name].energy_balance['heat']['inertial TES'][h] = self.locations[location_name].technologies['heatpump'].use(self.weather['temp_air'][h],self.locations[location_name].energy_balance['heat']['demand'][h],surplus,h) 
  
                # updata balance
                self.locations[location_name].energy_balance['electricity']['demand'][h] += self.locations[location_name].energy_balance['electricity']['heatpump'][h]
                self.locations[location_name].energy_balance['electricity']['grid'][h] += - self.locations[location_name].energy_balance['electricity']['heatpump'][h]
                self.energy_balance['electricity']['from grid'][h] += - self.locations[location_name].energy_balance['electricity']['heatpump'][h] # electricity fed into the grid from the whole rec at hour h
             
                
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
        electrolyzer = {}                 
        balances['REC'] = self.energy_balance
        parameters = {}        
        tech_cost = {}
        
        for location_name in self.locations:
            balances[location_name] = self.locations[location_name].energy_balance
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
                    
            tech_name = 'fuel cell'
            if tech_name in self.locations[location_name].technologies:
                if self.locations[location_name].technologies['fuel cell'].model == 'PEM General':
                    parameters[location_name][tech_name] = {}      
                    parameters[location_name][tech_name]['cell voltage'] = self.locations[location_name].technologies[tech_name].VOLT
                    parameters[location_name][tech_name]['current density'] = self.locations[location_name].technologies[tech_name].CURR_DENS
            
            tech_cost[location_name] = {}
            for tech_name in self.locations[location_name].technologies:
                tech_cost[location_name][tech_name] = self.locations[location_name].technologies[tech_name].cost
        
        directory = './results'
        if not os.path.exists(directory): os.makedirs(directory)
        
        with open('results/balances_'+simulation_name+".pkl", 'wb') as f: pickle.dump(balances, f)      
        with open('results/tech_params_'+simulation_name+".pkl", 'wb') as f: pickle.dump(parameters, f)
        with open('results/LOC_'+simulation_name+".pkl", 'wb') as f: pickle.dump(LOC, f)             
        with open('results/ageing_'+simulation_name+".pkl", 'wb') as f: pickle.dump(ageing, f)   
        with open('results/tech_cost_'+simulation_name+".pkl", 'wb') as f: pickle.dump(tech_cost, f)   
        
       
    def weather_generation(self,general,path,check):
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
                                
        if check and os.path.exists(path+'/weather/weather_TMY.csv'): # if the prevoius weather series can be used
            weather = pd.read_csv(path+'/weather/weather_TMY.csv')
        
        else: # if new weather data must be downoladed from PV gis
            print('Downolading typical metereological year data from PVGIS')   
                            
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
            weather.to_csv(path+'/weather/weather_TMY.csv')   

        return(weather)
   
    def tech_cost(self,tech_cost):
        for location_name in self.locations:
            for tech_name in self.locations[location_name].technologies:
                self.locations[location_name].technologies[tech_name].tech_cost(tech_cost[tech_name])
                