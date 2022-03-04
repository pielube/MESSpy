import pvlib #https://github.com/pvlib
import pandas as pd
import numpy as np

class PV:    
    
    def __init__(self,parameters,simulation_hours):
        """
        Create a PV object based on PV production taken from PVGIS data 
    
        parameters : dictionary
            'latitude': float
            'longitude': float
            'peakP': float peak DC power [kWp]
            'reference year': int year [2005 - 2015] used for output data and to get data from PVGIS if TMY = False 
            'losses': float losses in cables, power inverters, dirt (sometimes snow), over the years loss of power [%]
            'tilt':  float surface tilt [deg]
            'azimuth': float azimuth angle 0 = south, 180 = north [deg]
            'TMY': bool true if data of TMY is to be used    

        output : PV object able to:
            produce electricity .use(h)
        """
                
        latitude = parameters['latitude']
        longitude = parameters['longitude']
        peakP = parameters['peakP']
        year = parameters['reference year']
        losses = parameters['losses']
        tilt = parameters['tilt']
        azimuth = parameters['azimuth']
        tmybool = parameters['TMY']
        
        index60min = pd.date_range(start=str(year)+'-01-01 00:00:00',end=str(year)+'-12-31 23:00:00',freq='60T')
        
        if tmybool:
            weather = pvlib.iotools.get_pvgis_tmy(latitude, longitude, map_variables=True)[0]
            refindex = weather.index
            refindex = refindex.shift(10,'T')
        else:
            refindex = pd.date_range(start=str(year)+'-01-01 00:00:00',end=str(year)+'-12-31 23:00:00',freq='60T',tz='utc')
            refindex = refindex.shift(10,'T')
        
        # Actual production calculation (extract all available data points)
        res = pvlib.iotools.get_pvgis_hourly(latitude,longitude,surface_tilt=tilt,surface_azimuth=azimuth,pvcalculation=True,peakpower=peakP/1000,loss=losses)
        
        # Index to select TMY relevant data points
        pv = res[0]['P'] 
        pv = pv[refindex]
        pv.index = index60min
        
        # electricity produced every hour in the reference_year [kWh]
        self.production = np.tile(pv,int(simulation_hours/8760))
        # electricity produced every hour for the entire simulation [kWh]
        
        self.energy_balance = {'electricity': {'out': self.production}}
        
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
    