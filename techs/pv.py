import pvlib #https://github.com/pvlib
import pandas as pd
import numpy as np
import os
import pickle
import os
import sys 
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),os.path.pardir)))   # temorarily adding constants module path 
import constants as c

class PV:    
    
    def __init__(self,parameters,general,simulation_hours,location_name,path,check,rec_name):
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
            
            # check = True # True if no PV parameters are changed from the old simulation
            
            directory = './previous_simulation'
            if not os.path.exists(directory):
                os.makedirs(directory)
           
            if os.path.exists('previous_simulation/parameters_'+rec_name+location_name+'.pkl'):
                with open('previous_simulation/parameters_'+rec_name+location_name+'.pkl', 'rb') as f:
                    ps_parameters = pickle.load(f) # previous simulation location parameters
                par_to_check = ['tilt','azimuth','losses']
                for par in par_to_check:
                    if ps_parameters[par] != parameters[par]:
                        check = False
            else:
                check = False
          
            # if os.path.exists('previous_simulation/general.pkl'):
            #     with open('previous_simulation/general.pkl', 'rb') as f:
            #         ps_general = pickle.load(f) # previous simulation general
            #     par_to_check = ['latitude','longitude','UTC time zone','DST']
            #     for par in par_to_check:
            #         if ps_general[par] != general[par]:
            #             check = False  
            # else:
            #     check = False
                                    
            name_serie = rec_name + location_name + '_PV_TMY.csv'
            if check and os.path.exists(path+'/production/'+name_serie): # if the prevoius pv serie can be used
                pv = pd.read_csv(path+'/production/'+name_serie)['P'].to_numpy()
            
            else: # if a new pv serie must be downoladed from PV gis
                print('Downolading a new PV serie from PVgis for '+rec_name+location_name)   
                    
                latitude = general['latitude']
                longitude = general['longitude']
    
                losses = parameters['losses']
                tilt = parameters['tilt']
                azimuth = parameters['azimuth']
                
                weather = pvlib.iotools.get_pvgis_tmy(latitude, longitude, map_variables=True)[0]
                
                # Actual production calculation (extract all available data points)
                res = pvlib.iotools.get_pvgis_hourly(latitude,longitude,surface_tilt=tilt,surface_azimuth=azimuth,pvcalculation=True,peakpower=1,loss=losses,optimalangles=False)
                
                
                # Index to select TMY relevant data points
                pv = res[0]['P']
                refindex = weather.index
                shift_minutes = int(str(pv.index[0])[14:16])
                refindex = refindex.shift(shift_minutes,'T')
                pv = pv[refindex]
                pv = pd.DataFrame(pv)
                
                # time zone correction
                if general['UTC time zone'] > 0:
                    pv_index = pv.index
                    pv_index = pv_index.shift(general['UTC time zone']*60,'T')
                    pv.index = pv_index
                    
                    pv2 = pd.DataFrame(data=pv[-general['UTC time zone']:], index=None, columns=pv.columns)
                    pv = pv[:-general['UTC time zone']]
                    
                    reindex = pv.index[:general['UTC time zone']]
                    reindex = reindex.shift(-general['UTC time zone']*60,'T')
                    pv2.index = reindex  
                    pv = pd.concat([pv2,pv])
                    
                    pv['Local time']=pv.index
                    pv.set_index('Local time',inplace=True)
                    
                # Daily saving time (DST) correction
                # Is CEST (Central European Summertime) observed? if yes it means that State is applying DST
                # DST lasts between last sunday of march at 00:00:00+UTC+1 and last sunday of october at 00:00:00+UTC+2
                # For example in Italy DST in 2022 starts in March 27th at 02:00:00 and finishes in October 30th at 03:00:00
                if general['DST']==True:
                
                    zzz_in=pv[pv.index.month==3]
                    zzz_in=zzz_in[zzz_in.index.weekday==6]
                    zzz_in=zzz_in[zzz_in.index.hour==1+general['UTC time zone']]
                    zzz_in = pd.Series(zzz_in.index).unique()[-1]
                  
                    zzz_end=pv[pv.index.month==10]
                    zzz_end=zzz_end[zzz_end.index.weekday==6]
                    zzz_end=zzz_end[zzz_end.index.hour==1+general['UTC time zone']]
                    zzz_end = pd.Series(zzz_end.index).unique()[-1]
                    
                    pv.loc[zzz_in:zzz_end] = pv.loc[zzz_in:zzz_end].shift(60,'T')
                    pv=pv.interpolate(method='linear')
                
                    pv['Local time - DST']=pv.index
                    pv.set_index('Local time - DST',inplace=True)
                
                # save series .csv
                pv.to_csv(path+'/production/'+name_serie)
            
                # save new parameters in previous_simulation
                with open('previous_simulation/parameters_'+rec_name+location_name+'.pkl', 'wb') as f:
                    pickle.dump(parameters, f) 
                    
                # with open('previous_simulation/general.pkl', 'wb') as f:
                #     pickle.dump(general, f)             
            
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

