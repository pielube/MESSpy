import pvlib #https://github.com/pvlib
import pandas as pd
import numpy as np
import os
import pickle

class PV:    
    
    def __init__(self,parameters,general,simulation_hours,location_name):
        """
        Create a PV object based on PV production taken from PVGIS data 
    
        parameters : dictionary
            'peakP': float peak DC power [kWp] 
            'losses': float losses in cables, power inverters, dirt (sometimes snow), over the years loss of power [%]
            'tilt':  float surface tilt [deg]
            'azimuth': float azimuth angle 0 = south, 180 = north [deg]  
            
        general: dictironary
            see rec.py

        output : PV object able to:
            produce electricity .use(h)
        """
        
        
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
      
        if os.path.exists('previous_simulation/general_'+location_name+'.pkl'):
            with open('previous_simulation/general_'+location_name+'.pkl', 'rb') as f:
                ps_general = pickle.load(f) # previous simulation general
            par_to_check = ['latitude','longitude','reference year','TMY']
            for par in par_to_check:
                if ps_general[par] != general[par]:
                    check = False  
        else:
            check = False
                                
        name_serie = location_name + 'PV.csv'
        if check and os.path.exists('previous_simulation/'+name_serie): # if the prevoius pv serie can be used
            pv = pd.read_csv('previous_simulation/'+name_serie)['P'].to_numpy()
        
        else: # if a new pv serie must be downoladed from PV gis
            print('downolading a new PV serie from PVgis for '+location_name)   
            
            # save new parameters in previous_simulation
            with open('previous_simulation/parameters_'+location_name+'.pkl', 'wb') as f:
                pickle.dump(parameters, f) 
                
            with open('previous_simulation/general_'+location_name+'.pkl', 'wb') as f:
                pickle.dump(general, f)               
                
            latitude = general['latitude']
            longitude = general['longitude']
            tmybool = general['TMY']
            year = general['reference year']
            
            losses = parameters['losses']
            tilt = parameters['tilt']
            azimuth = parameters['azimuth']
            
            index60min = pd.date_range(start=str(year)+'-01-01 00:00:00',end=str(year)+'-12-31 23:00:00',freq='60T')
            
            if tmybool:
                weather = pvlib.iotools.get_pvgis_tmy(latitude, longitude, map_variables=True)[0]
                refindex = weather.index
                refindex = refindex.shift(10,'T')
            else:
                refindex = pd.date_range(start=str(year)+'-01-01 00:00:00',end=str(year)+'-12-31 23:00:00',freq='60T',tz='utc')
                refindex = refindex.shift(10,'T')
            
            # Actual production calculation (extract all available data points)
            res = pvlib.iotools.get_pvgis_hourly(latitude,longitude,surface_tilt=tilt,surface_azimuth=azimuth,pvcalculation=True,peakpower=1,loss=losses)
            
            # Index to select TMY relevant data points
            pv = res[0]['P'] 
            pv = pv[refindex]
            pv.index = index60min
    
            series_frame = pd.DataFrame(pv)
            series_frame.to_csv('previous_simulation/'+name_serie)
        
        peakP = parameters['peakP']
        pv = pv * peakP/1000        
        # electricity produced every hour in the reference_year [kWh]
        self.production = np.tile(pv,int(simulation_hours/8760))
        # electricity produced every hour for the entire simulation [kWh]
        
    def use(self,h):
        """
        Produce electricity
        
        h : int hour to be simulated
    
        output : float electricity produced that hour [kWh]    
        """
        
        return(self.production[h])




###########################################################################################################################################################

if __name__ == "__main__":
    
    """
    Functional test
    """
    
    inp_test = {'latitude': 50.6,
                'longitude': 5.6,
                'Ppeak': 5,
                'losses': 0.1,
                'tilt': 30,
                'azimuth': 0,
                'reference year':2016,
                'TMY': False}
    
    pv1 = PV(inp_test,2,1)
    
    for h in range(8760):
        print(pv1.use(h))
    print(sum(pv1.production))
    