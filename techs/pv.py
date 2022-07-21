import pvlib #https://github.com/pvlib
import pandas as pd
import numpy as np
import os
import pickle

class PV:    
    
    def __init__(self,parameters,general,simulation_hours,location_name,path):
        """
        Create a PV object based on PV production taken from PVGIS data 
    
        parameters : dictionary
            'peakP': float peak DC power [kWp] 
            'losses': float losses in cables, power inverters, dirt (sometimes snow), over the years loss of power [%]
            'tilt':  float surface tilt [deg]
            'azimuth': float azimuth angle 0 = south, 180 = north [deg]  
            'serie': if "TMY" production serie based on typical meteorological year is used
                if "filename.csv" a different serie can be used (upload it in input/production)
                "filename.csv" must be the hourly time series of PV production 8760 values [Wh]
                in this case 'peakP', 'azimuth' and 'tilt' are ignored
            
        general: dictironary
            see rec.py

        output : PV object able to:
            produce electricity .use(h)
        """
        
        if parameters['serie'] == "TMY":
            ### If PV serie have already been downloaded and saved as file.csv, this file is used
            ### Otherwise new serie is downloaded from PVgis (type meteorological year)
            
            check = True # True if no PV parameters are changed from the old simulation
            
            directory = './previous_simulation'
            if not os.path.exists(directory):
                os.makedirs(directory)
           
            if os.path.exists('previous_simulation/parameters_'+location_name+'.pkl'):
                with open('previous_simulation/parameters_'+location_name+'.pkl', 'rb') as f:
                    ps_parameters = pickle.load(f) # previous simulation location parameters
                par_to_check = ['tilt','azimuth','losses']
                for par in par_to_check:
                    if ps_parameters[par] != parameters[par]:
                        check = False
            else:
                check = False
          
            if os.path.exists('previous_simulation/general.pkl'):
                with open('previous_simulation/general.pkl', 'rb') as f:
                    ps_general = pickle.load(f) # previous simulation general
                par_to_check = ['latitude','longitude']
                for par in par_to_check:
                    if ps_general[par] != general[par]:
                        check = False  
            else:
                check = False
                                    
            name_serie = location_name + '_PV_TMY.csv'
            if check and os.path.exists(path+'/production/'+name_serie): # if the prevoius pv serie can be used
                pv = pd.read_csv(path+'/production/'+name_serie)['P'].to_numpy()
            
            else: # if a new pv serie must be downoladed from PV gis
                print('downolading a new PV serie from PVgis for '+location_name)   
                
                # save new parameters in previous_simulation
                with open('previous_simulation/parameters_'+location_name+'.pkl', 'wb') as f:
                    pickle.dump(parameters, f) 
                    
                with open('previous_simulation/general.pkl', 'wb') as f:
                    pickle.dump(general, f)               
                    
                latitude = general['latitude']
                longitude = general['longitude']
    
                losses = parameters['losses']
                tilt = parameters['tilt']
                azimuth = parameters['azimuth']
                
                weather = pvlib.iotools.get_pvgis_tmy(latitude, longitude, map_variables=True)[0]
                refindex = weather.index
                refindex = refindex.shift(10,'T')
    
                # Actual production calculation (extract all available data points)
                res = pvlib.iotools.get_pvgis_hourly(latitude,longitude,surface_tilt=tilt,surface_azimuth=azimuth,pvcalculation=True,peakpower=1,loss=losses)
                
                # Index to select TMY relevant data points
                pv = res[0]['P'] 
                pv = pv[refindex]
            
                series_frame = pd.DataFrame(pv)
                series_frame.to_csv(path+'/production/'+name_serie)
            
            peakP = parameters['peakP']
            pv = pv * peakP/1000        
            # electricity produced every hour in the reference_year [kWh]
            self.production = np.tile(pv,int(simulation_hours/8760))
            # electricity produced every hour for the entire simulation [kWh]
            
        else:
            # read production serie if "TMY" not have to be used
            pv = pd.read_csv(path+'/production/'+parameters['serie'])['P'].to_numpy()
            pv = pv * (1-parameters['losses'])          # add losses
            pv = pv/1000                             # Wh -> kWh
            self.production = np.tile(pv,int(simulation_hours/8760))
            
            
    def use(self,h):
        """
        Produce electricity
        
        h : int hour to be simulated
    
        output : float electricity produced that hour [kWh]    
        """
        
        return(self.production[h])

