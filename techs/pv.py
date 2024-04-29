import pvlib #https://github.com/pvlib
import pandas as pd
import numpy as np
import os
import pickle
import sys 
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),os.path.pardir)))   # temorarily adding constants module path 
from core import constants as c

class PV:    
    
    def __init__(self,parameters,location_name,path,check,file_structure,file_general):
        """
        Create a PV object based on PV production taken from PVGIS data 
    
        parameters : dictionary
            'peakP': float peak DC power [kWp] 
            'losses': float losses in cables, power inverters, dirt (sometimes snow), over the years loss of power [%]
            'tilt':  float surface tilt [deg]
            'azimuth': float azimuth angle 0 = south, 180 = north [deg]  
            'serie': if "TMY" production serie based on typical meteorological year is used
                if INT [2005-2016] a serie of the specific year is used
                if "filename.csv" a different serie can be used (upload it in input/production)
                "filename.csv" must be the hourly time series of PV production 8760 values [Wh]
                in this case 'peakP', 'azimuth' and 'tilt' are ignored
            'Max field width': optional: float >0
            'Max field length': optional: float >0
            
        general: dictironary
            see rec.py

        output : PV object able to:
            produce electricity .use(h)
        """
        
        self.cost = False # will be updated with tec_cost()
        if 'owned' in parameters:
            self.property = parameters['owned']         # bool value to take into account if the plant is owned or only electricity purchase is considered. 
                                                        # Main impact on economic assessment and key parameters

        if parameters['serie'] == "TMY" or type(parameters['serie']) == int:
            ### If PV serie have already been downloaded and saved as file.csv, this file is used
            ### Otherwise new serie is downloaded from PVgis (type meteorological year)
            
            # check = True # True if no PV parameters are changed from the old simulation
            
            directory = './previous_simulation'
            if not os.path.exists(directory): os.makedirs(directory)
           
            if os.path.exists(f"previous_simulation/{file_structure}_{location_name}.pkl"):
                with open(f"previous_simulation/{file_structure}_{location_name}.pkl", 'rb') as f: ps_parameters = pickle.load(f) # previous simulation location parameters
                par_to_check = ['tilt','azimuth','losses','serie']
                for par in par_to_check:
                    if ps_parameters[par] != parameters[par]:
                        check = False
                        
            else:
                check = False
                                    
            name_serie = f"PV_{parameters['serie']}_{location_name}_{file_general}_{file_structure}.csv"
            if check and os.path.exists(path+'/production/'+name_serie): # if the prevoius pv serie can be used
                pv = pd.read_csv(path+'/production/'+name_serie)['P'].to_numpy()
            
            else: # if a new pv serie must be downoladed from PV gis
                print(f"Downolading a new PV serie from PVgis for {location_name}_{file_general}_{file_structure}") 
                    
                losses = parameters['losses']
                tilt = parameters['tilt']
                azimuth = parameters['azimuth']
                
                
                if parameters['serie'] == 'TMY':
                    weather = pvlib.iotools.get_pvgis_tmy(c.latitude, c.longitude, map_variables=True)[0]
                    # Actual production calculation (extract all available data points)
                    res = pvlib.iotools.get_pvgis_hourly(c.latitude,c.longitude,surface_tilt=tilt,surface_azimuth=azimuth,pvcalculation=True,peakpower=1,loss=losses,optimalangles=False)
                    # Index to select TMY relevant data points
                    pv = res[0]['P']
                    refindex = weather.index
                    shift_minutes = int(str(pv.index[0])[14:16])
                    refindex = refindex.shift(shift_minutes,'min')
                    pv = pv[refindex]
                    
                else: # INT
                    year = parameters['serie']
                    res = pvlib.iotools.get_pvgis_hourly(c.latitude,c.longitude,start=year,end=year,surface_tilt=tilt,surface_azimuth=azimuth,pvcalculation=True,peakpower=1,loss=losses,optimalangles=False)
                    pv = res[0]['P']
                
                pv = pd.DataFrame(pv)
                
                # time zone correction
                if c.UTC > 0:
                    pv_index = pv.index
                    pv_index = pv_index.shift(c.UTC*60,'min')
                    pv.index = pv_index
                    
                    pv2 = pd.DataFrame(data=pv[-c.UTC:], index=None, columns=pv.columns)
                    pv = pv[:-c.UTC]
                    
                    reindex = pv.index[:c.UTC]
                    reindex = reindex.shift(-c.UTC*60,'min')
                    pv2.index = reindex  
                    pv = pd.concat([pv2,pv])
                    
                    pv['Local time']=pv.index
                    pv.set_index('Local time',inplace=True)
                    
                # Daily saving time (DST) correction
                # Is CEST (Central European Summertime) observed? if yes it means that State is applying DST
                # DST lasts between last sunday of march at 00:00:00+UTC+1 and last sunday of october at 00:00:00+UTC+2
                # For example in Italy DST in 2022 starts in March 27th at 02:00:00 and finishes in October 30th at 03:00:00
                if c.DST==True:
                
                    zzz_in=pv[pv.index.month==3]
                    zzz_in=zzz_in[zzz_in.index.weekday==6]
                    zzz_in=zzz_in[zzz_in.index.hour==1+c.UTC]
                    zzz_in = pd.Series(zzz_in.index).unique()[-1]
                  
                    zzz_end=pv[pv.index.month==10]
                    zzz_end=zzz_end[zzz_end.index.weekday==6]
                    zzz_end=zzz_end[zzz_end.index.hour==1+c.UTC]
                    zzz_end = pd.Series(zzz_end.index).unique()[-1]
                    
                    pv.loc[zzz_in:zzz_end] = pv.loc[zzz_in:zzz_end].shift(60,'min')
                    pv=pv.interpolate(method='linear')
                
                    pv['Local time - DST']=pv.index
                    pv.set_index('Local time - DST',inplace=True)
                
                # save series .csv
                pv.to_csv(path+'/production/'+name_serie)
            
                # save new parameters in previous_simulation
                with open(f"previous_simulation/{file_structure}_{location_name}.pkl", 'wb') as f:
                    pickle.dump(parameters, f)            
            
            self.peakP = parameters['peakP']
            pv = pv * self.peakP/1000
            
            # electricity produced every hour in the reference_year [kWh] [kW]
            self.production = np.tile(pv,int(c.timestep_number*c.timestep/60/8760)) # from 1 year to simlation length years
        
            # from hourly to timestep
            if c.timestep < 60:
                self.production =  np.repeat(self.production, 60/c.timestep) # [kW] creating a production series alligned with selected timestep 
        else:
            # read a specific production serie expressed as kW/kWpeak
            self.peakP = parameters['peakP']
            pv = pd.read_csv(path+'/production/'+parameters['serie'])['P'].to_numpy()
            pv = pv * (1-parameters['losses']/100)      # add losses if to be added
            pv = pv*self.peakP                          # kWh
            self.production = np.tile(pv,c.timestep_number*c.timestep/60/8760)
            if len(self.production != c.timestep_number):
                raise ValueError(f"Warning! Checks the length and timestep of the PV production you input for {location_name}.")


    def use(self,step):
        """
        Produce electricity
        
        step : int step to be simulated
    
        output : float electricity produced that hour [kWh]    
        """
        
        return(self.production[step])

    def tech_cost(self,tech_cost):
        """
        Parameters
        ----------
        tech_cost : dict
            'cost per unit': float [€/kW]
            'OeM': float, percentage on initial investment [%]
            'refud': dict
                'rate': float, percentage of initial investment which will be rimbursed [%]
                'years': int, years for reimbursment
            'replacement': dict
                'rate': float, replacement cost as a percentage of the initial investment [%]
                'years': int, after how many years it will be replaced

        Returns
        -------
        self.cost: dict
            'total cost': float [€]
            'OeM': float, percentage on initial investment [%]
            'refud': dict
                'rate': float, percentage of initial investment which will be rimbursed [%]
                'years': int, years for reimbursment
            'replacement': dict
                'rate': float, replacement cost as a percentage of the initial investment [%]
                'years': int, after how many years it will be replaced
        """
        
        tech_cost = {key: value for key, value in tech_cost.items()}
        
        size = self.peakP
        
        if tech_cost['cost per unit'] == 'default price correlation':
            C0 = 1200 # €/kW
            scale_factor = 0.9 # 0:1
            C = C0 * size ** scale_factor
        else:
            C = size * tech_cost['cost per unit']

        tech_cost['total cost'] = tech_cost.pop('cost per unit')
        tech_cost['total cost'] = C
        tech_cost['OeM'] = tech_cost['OeM'] *C /100 # €

        self.cost = tech_cost    
        
        
        
###########################################################################################################################################################

if __name__ == "__main__":
    
    """
    Functional test
    """
    
    
    
    
    