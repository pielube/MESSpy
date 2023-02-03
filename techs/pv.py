import pvlib #https://github.com/pvlib
import pandas as pd
import numpy as np
import os
import pickle
import os
import sys 
import math as m
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
            'Max field width': optional: float >0
            'Max field length': optional: float >0
            
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
            
            self.peakP = parameters['peakP']
            pv = pv * self.peakP/1000
            
            # electricity produced every hour in the reference_year [kWh]
            self.production = np.tile(pv,int(simulation_hours/8760))
            # electricity produced every hour for the entire simulation [kWh]
            
            
        else:
            # read production serie if "TMY" not have to be used
            pv = pd.read_csv(path+'/production/'+parameters['serie'])['P'].to_numpy()
            pv = pv * (1-parameters['losses'])          # add losses
            pv = pv/1000                             # Wh -> kWh
            self.production = np.tile(pv,int(simulation_hours/8760))
        
        if 'Max field width' and 'Max field length' in parameters:   
            # Covered Surface calculation
            # Considering modules having size 1m x 1,5m and power 250W it means that, for each kWp, four modules are needed, thus covering 6 m2 of surface (4m width and 1,5m high)
            # Knowing the maximum field width by input, the number of panel rows can be found as follows
            module_width = 1 #m
            module_height = 1.5 #m
            field_ideal_width = self.peakP*4*module_width #m                                                  # each kWp is 4m width
            max_field_width = parameters ['Max field width'] #m
            self.n_rows = field_ideal_width/max_field_width
            latitude = general['latitude'] #°
            sun_angle_winter_solstice_rad = (90-(latitude+23.5))*m.pi/180 #rad                           # https://it.science19.com/how-to-calculate-winter-solstice-sun-angle-1948
            tilt_rad = (parameters['tilt']*m.pi)/180 #rad
            rows_distance =  module_height*(m.sin(tilt_rad)/m.tan(sun_angle_winter_solstice_rad))  #m    # d = [sen(β)]/[tg(δ)]·L  https://www.ediltecnico.it/85611/come-calcolare-la-distanza-minima-di-installazione-di-file-di-pannelli-fotovoltaici/
            row_length = module_height*m.cos(tilt_rad)+rows_distance #m
            
            if self.n_rows >= 1:
                self.surfac_cov = max_field_width*row_length*(int(self.n_rows)-1)+max_field_width*module_height*m.cos(tilt_rad)+(self.n_rows-int(self.n_rows))*max_field_width*row_length #m2
                if self.n_rows-int(self.n_rows) != 0:
                    field_length = row_length*int(self.n_rows)+module_height*m.cos(tilt_rad) #m
                else:
                    field_length = row_length*(int(self.n_rows)-1)+module_height*m.cos(tilt_rad) #m
                self.surface_cov_rectangular = field_length*max_field_width #m2
                
            if 0 < self.n_rows < 1:
                self.surface = self.n_rows*max_field_width*module_height*m.cos(tilt_rad)
                field_length = module_height*m.cos(tilt_rad)
                self.surface_cov_rectangular = self.surface
                
            max_field_length = parameters ['Max field length'] #m
            if field_length > max_field_length:
                print('Warning!! The surface covered by the '+rec_name+location_name+' PV field is higher than the maximum available one')
            
        
    def use(self,h):
        """
        Produce electricity
        
        h : int hour to be simulated
    
        output : float electricity produced that hour [kWh]    
        """
        
        return(self.production[h])

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
            C = size * C0 **  scale_factor
        else:
            C = size * tech_cost['cost per unit']

        tech_cost['total cost'] = tech_cost.pop('cost per unit')
        tech_cost['total cost'] = C
        tech_cost['OeM'] = tech_cost['OeM'] *C /100 # €

        self.cost = tech_cost    