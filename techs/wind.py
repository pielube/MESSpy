import numpy as np
import pandas as pd
import math
import warnings
import os
import sys 
import pvlib
import pickle    
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),os.path.pardir)))   # temorarily adding constants module path 
from core import constants as c
import matplotlib.pyplot as plt

class wind:    
    
    def __init__(self, parameters, location_name, path, check, file_structure, file_general):
        """
        Create a wind object based on the specified model
    
        parameters : dictionary
            'model': str type of model to be used for wind
                    - power curve -> based on https://doi.org/10.1016/j.est.2021.103893.
                                    If a series is given, [in kW/kWpeak], wind production data can be retrieved from https://www.renewables.ninja/.
                                    If a series is given, only one parameter needs to be defined: 'Npower'.                                                                 
                    - betz -> simple model based on Betz theory
                    - detailed -> more detailed model based on
                    Saint-Drenan, Yves-Marie, et al. 
                    "A parametric model for wind turbine power curves incorporating environmental conditions." 
                    Renewable Energy 157 (2020): 754-768.
                              
            'area': float swept area [m2] e.g. 39.6 m^2 (Aircon 10/10 kW)
            'efficiency': float total efficiency = Betz*efficiency [-] default: 0.45 (ca. 0.593*0.76, i.e. Betz*efficiency)
            'Npower': float rated power [kW] # NOTE: useless for 'betz' and 'detailed' methods
            'WSrated': float rated wind speed [m/s] e.g. 11.0 m/s (Aircon 10/10 kW)
            'WScutin': float cut in wind speed [m/s] e.g.  2.5 m/s (Aircon 10/10 kW)
            'WScutoff': float cut off wind speed [m/s] e.g. 32.0 m/s (Aircon 10/10 kW)
            
            'omega_min': [rpm] OPTIONAL default = from eq.
            'omega_max': [rpm] OPTIONAL default = from eq.
            'beta': [°] e.g. 0°
            'cp_max': [-] OPTIONAL default = 0.44; values from 0.4 to 0.5
            'idx': [-] e.g. 5; values from 0 to 5
            
            'z_i': [m] wind turbine height
            'z_hub': [m] hub height, e.g. 30 m (Aircon 10/10 kW)
            'alpha': [-] Hellman or shear coefficient, values from 0 to 0.4, default 0.4
            'Vu': [°/m] Veer coefficient, values from 0 to 0.75 °/m
            'Nbands': [-] Number of horizontal bands
            
            'serie': if "TMY" production series based on typical meteorological year, or a specific year, or a custom CSV file
            'ageing': bool, if wind turbine performance degradation must be considered
            'degradation_factor': float, [%] of performance loss every year
        output : wind object able to:
            produce electricity .use(step)
        """
        self.parameters         = parameters                                                                                          
        self.model              = self.parameters['model']      # [-]  wind technology model                                                     
        
        if self.model not in ['betz','detailed','power curve']:
            raise ValueError("Warning: selected model for wind techonology is not among the available options.\n\
            Option to fix the problem: \n\
                (a) - In studycase.py file choose one among 'betz','detailed' and 'power curve' models.")

        self.cost               = False                         # will be updated with tech_cost()
        self.property           = self.parameters['owned']      # bool value to take into account if the plant is owned or only electricity purchase is considered. It only impacts impact on economic assessment and key parameters
        self.vw_ci              = self.parameters['WScutin']    # Cut-in wind speed (m/s)
        self.vw_r               = self.parameters['WSrated']    # Rated wind speed (m/s)                                                        
        self.vw_co              = self.parameters['WScutoff']   # Cut-out wind speed (m/s)
        
        if self.model == 'power curve':            
            self.Npower             = self.parameters['Npower']     # [kW] wind technology nominal power
            self.z_i                = self.parameters['z_i']        # Height of the wind turbine (m)
            self.href               = 10                            # Reference height for wind speed data from PVGIS (m)
            self.alpha              = self.parameters.get('alpha', 0.4)    # Exponent law coefficient (default to 0.14)                        
            
        if self.model == 'betz':
            self.z_i                = self.parameters['z_i']    # Height of the wind turbine (m)
            self.href               = 10                        # Reference height for wind speed data from PVGIS (m)
            self.alpha              = self.parameters.get('alpha', 0.14)    # Exponent law coefficient (default to 0.14)        
            self.rho                     = c.AIRSDENSITY             # [kg/m3] air density assumed constant. Must be upgraded to be a function of external weather conditions. 
            self.area               = self.parameters['area']
            self.efficiency         = self.parameters['efficiency']
            # Nominal power calculation                
            self.Npower = 0.5*self.rho*self.area*self.vw_r**3*self.efficiency/1000
            print(f"The size of the wind turbine is: {round(self.Npower,2)} kW")
                    
        if self.model == 'detailed':   
            self.z_i                = self.parameters['z_i']    # Height of the wind turbine (m)
            self.alpha              = self.parameters.get('alpha', 0.14)    # Exponent law coefficient (default to 0.14) 
            self.rho                     = c.AIRSDENSITY             # [kg/m3] air density assumed constant. Must be upgraded to be a function of external weather conditions. 
            self.area               = self.parameters['area']
            self.Nbands             = self.parameters['Nbands']
            self.cp_max                  = self.parameters.get('cp_max', 0.44)        # Consider 0.44 if not specified
            self.beta               = self.parameters.get('beta', 0)        # Consider 0.44 if not specified
            self.idx                = self.parameters.get('idx', 5)        # Consider 5 if not specified
            self.z_hub              = self.parameters.get('z_hub', 30)        # Consider 30m if not specified
            self.Vu                 = self.parameters.get('Vu', 0.5)        # Consider 0.5°/m if not specified
            # Nominal power calculation
            nominal_power_coefficent = self.cpfunc(self.vw_r,self.area,self.beta,self.idx,self.cp_max) 
            self.Npower = 0.5*self.rho*self.area*self.vw_r**3*nominal_power_coefficent/1000
            print(f"The size of the wind turbine is: {round(self.Npower,2)} kW")
            
        self.timestep           = c.timestep            # [min] simulation timestep 
        self.timestep_number    = c.timestep_number     # [-] number of timesteps considered in the simulation
        self.simulation_years   = c.simulation_years    # [-] simulation years
        self.latitude           = c.latitude
        self.longitude          = c.longitude

        
        if self.parameters['serie'] == "TMY" or type(self.parameters['serie']) == int:
            ### If wind serie has already been downloaded and saved as file.csv, this file is used
            ### Otherwise new serie is downloaded from PVGIS
            
            # check = True # True if no wind parameters are changed from the old simulation
                
            # Directory for storing previous simulation data
            directory = './previous_simulation'
            if not os.path.exists(directory):
                os.makedirs(directory)
    
            # Checking if the previous simulation exists
            if os.path.exists(f"{directory}/wind_{file_structure}_{location_name}.pkl"):
                with open(f"{directory}/wind_{file_structure}_{location_name}.pkl", 'rb') as f:
                    ps_parameters = pickle.load(f)  # Load previous simulation parameters
                    par_to_check = ['model', 'WScutin', 'WSrated', 'WScutoff', 'z_i', 'alpha', 'serie', 'area', 'efficiency', 'Nbands', 'cp_max', 'beta', 'idx', 'z_hub', 'Vu']
                    for par in par_to_check:
                        if par in ps_parameters and par in self.parameters:  # Some parameters haven't to be defined
                            if ps_parameters[par] != self.parameters[par]:
                                check = False
    
            else:
                check = False

            name_serie = f"Wind_{self.parameters['serie']}_{location_name}_{file_general}_{file_structure}.csv"
            if check and os.path.exists(path + '/production/' + name_serie):  # If previous wind series can be used
                wind_data = pd.read_csv(path + '/production/' + name_serie)['P'].to_numpy()
                
                if self.model == 'power curve':
                    wind_data = wind_data * self.Npower
    
            else:  # Download new wind data from PVGIS
                print(f"Downloading new wind series from PVGIS for {location_name}_{file_general}_{file_structure}")
                # Retrieve wind speed data from PVGIS based on selected 'serie'
                if self.parameters['serie'] == "TMY":
                    weather = pvlib.iotools.get_pvgis_tmy(self.latitude, self.longitude, map_variables=True)[0]
                    res = pvlib.iotools.get_pvgis_hourly(self.latitude, self.longitude)
                    wind_speed_data = res[0]['wind_speed']  # Wind speed at 10m height
                    refindex = weather.index
                    shift_minutes = int(str(wind_speed_data.index[0])[14:16])
                    refindex = refindex.shift(shift_minutes,'min')
                    wind_speed_data = wind_speed_data[refindex]
           
                else:  # If specific year 
                    year = self.parameters['serie']
                    res = pvlib.iotools.get_pvgis_hourly(self.latitude, self.longitude, start=year, end=year)
                    wind_speed_data = res[0]['wind_speed']

                # Remove 29th of february if present
                wind_speed_data = wind_speed_data[~((wind_speed_data.index.month == 2) & (wind_speed_data.index.day == 29))]
                
                wind_speed_data = pd.DataFrame(wind_speed_data)

                # Time correction (UTC, DST)
                if c.UTC > 0:
                    wind_speed_data_index = wind_speed_data.index
                    wind_speed_data_index = wind_speed_data_index.shift(c.UTC * 60, 'min')
                    wind_speed_data.index = wind_speed_data_index
        
                    wind_speed_data2 = pd.DataFrame(data=wind_speed_data[-c.UTC:], index=None, columns=wind_speed_data.columns)
                    wind_speed_data = wind_speed_data[:-c.UTC]
                    
                    reindex = wind_speed_data.index[:c.UTC]
                    reindex = reindex.shift(-c.UTC*60,'min')
                    wind_speed_data2.index = reindex  
                    wind_speed_data = pd.concat([wind_speed_data2,wind_speed_data])
                    
                    wind_speed_data['Local time']=wind_speed_data.index
                    wind_speed_data.set_index('Local time',inplace=True)

                # Daily saving time (DST) correction
                # Is CEST (Central European Summertime) observed? if yes it means that State is applying DST
                # DST lasts between last sunday of march at 00:00:00+UTC+1 and last sunday of october at 00:00:00+UTC+2
                # For example in Italy DST in 2022 starts in March 27th at 02:00:00 and finishes in October 30th at 03:00:00
                if c.DST==True:
                
                    zzz_in=wind_speed_data[wind_speed_data.index.month==3]
                    zzz_in=zzz_in[zzz_in.index.weekday==6]
                    zzz_in=zzz_in[zzz_in.index.hour==1+c.UTC]
                    zzz_in = pd.Series(zzz_in.index).unique()[-1]
                  
                    zzz_end=wind_speed_data[wind_speed_data.index.month==10]
                    zzz_end=zzz_end[zzz_end.index.weekday==6]
                    zzz_end=zzz_end[zzz_end.index.hour==1+c.UTC]
                    zzz_end = pd.Series(zzz_end.index).unique()[-1]
                    
                    wind_speed_data.loc[zzz_in:zzz_end] = wind_speed_data.loc[zzz_in:zzz_end].shift(60,'min')
                    wind_speed_data=wind_speed_data.interpolate(method='linear')
                
                    wind_speed_data['Local time - DST']=wind_speed_data.index
                    wind_speed_data.set_index('Local time - DST',inplace=True)
   
    
                if self.model == 'power curve':  # https://doi.org/10.1016/j.est.2021.103893
                    # Step 2: Correct wind speed for the turbine's height using the power law
                    corrected_wind_speed = wind_speed_data * (self.z_i / self.href) ** self.alpha
                    # Step 3: Calculate the power output based on the corrected wind speed
                    hourly_power_output_1kW = []
                    for vw in corrected_wind_speed['wind_speed']:
                        if vw <= self.vw_ci or vw >= self.vw_co:
                            hourly_power_output_1kW.append(0)
                        elif self.vw_ci < vw <= self.vw_r:
                            # Power output between cut-in and rated wind speed (cubic interpolation)
                            power = 1 * ((vw ** 3 - self.vw_ci ** 3) / (self.vw_r ** 3 - self.vw_ci ** 3))
                            hourly_power_output_1kW.append(power)
                        elif self.vw_r < vw <= self.vw_co:
                            # Power output at rated wind speed
                            hourly_power_output_1kW.append(1)
                        else:
                            hourly_power_output_1kW.append(0)
        
                    wind_data = pd.DataFrame(hourly_power_output_1kW, columns=['P'])
                    wind_speed_index = corrected_wind_speed.index
                    wind_data.set_index(wind_speed_index, inplace=True)
                    # save series .csv Saving power per kW installed
                    wind_data.to_csv(path + '/production/' + name_serie)
                    # Save new parameters in previous_simulation
                    wind_data = np.array(wind_data['P'])
                    with open(f"{directory}/wind_{file_structure}_{location_name}.pkl", 'wb') as f:
                        pickle.dump(self.parameters, f)

                    wind_data = wind_data * self.Npower
                 
            
                if self.model == 'betz':
                    
                    # Step 2: Correct wind speed for the turbine's height using the power law
                    corrected_wind_speed = wind_speed_data * (self.z_i / self.href) ** self.alpha

                    # Calculate the power output based on the corrected wind speed
                    hourly_power_output = []
                    for vw in corrected_wind_speed['wind_speed']:
                        if vw < self.vw_ci:
                            vw = 0
                            hourly_power_output.append(0)
                        elif self.vw_ci < vw <= self.vw_r:
                            vw = vw
                            hourly_power_output.append(0.5*self.rho*self.area*vw**3*self.efficiency/1000) # [kW] power generation in the considered timestep
                        elif self.vw_r < vw < self.vw_co:
                            vw = self.vw_r
                            hourly_power_output.append(0.5*self.rho*self.area*vw**3*self.efficiency/1000) # [kW] power generation in the considered timestep
                        elif vw > self.vw_co:
                            vw = 0
                            hourly_power_output.append(0)
                     
                    wind_data = pd.DataFrame(hourly_power_output, columns=['P'])
                    wind_speed_index = corrected_wind_speed.index
                    wind_data.set_index(wind_speed_index, inplace=True)
                    # save series .csv In betz model it is saved the total power produced by the turbine (not the power per kW)
                    wind_data.to_csv(path + '/production/' + name_serie)
                    wind_data = np.array(wind_data['P'])
                    # Save new parameters in previous_simulation
                    with open(f"{directory}/wind_{file_structure}_{location_name}.pkl", 'wb') as f:
                        pickle.dump(self.parameters, f)


                if self.model == 'detailed':

                    # Step 2: Correct wind speed for the turbine's height using the power law
                    corrected_wind_speed = []
                    for vw in wind_speed_data['wind_speed']:
                        corrected_wind_speed.append(self.eqspeed(vw,self.z_i,self.z_hub,self.alpha,self.area,self.Vu,self.Nbands))                    
                    corrected_wind_speed = pd.DataFrame(corrected_wind_speed,columns=['wind_speed'])
                    wind_speed_index = wind_speed_data.index
                    corrected_wind_speed.set_index(wind_speed_index, inplace=True)

                    # Calculate the power output based on the corrected wind speed
                    hourly_power_output = []
                    for vw in corrected_wind_speed['wind_speed']:
                        if vw < self.vw_ci:
                            vw = 0                            
                            hourly_power_output.append(0)
                        elif self.vw_ci < vw <= self.vw_r:
                            vw = vw
                            powercoeff = self.cpfunc(vw,self.area,self.beta,self.idx,self.cp_max)
                            power_output = 0.5*self.rho*self.area*vw**3*powercoeff/1000
                            if power_output < 0:
                                power_output = 0
                            hourly_power_output.append(power_output) # [kW] power generation in the considered timestep
                        elif self.vw_r < vw < self.vw_co:
                            vw = self.vw_r
                            power_output = 0.5*self.rho*self.area*vw**3*powercoeff/1000
                            if power_output < 0:
                                power_output = 0
                            hourly_power_output.append(power_output) # [kW] power generation in the considered timestep
                        elif vw > self.vw_co:
                            vw = 0
                            hourly_power_output.append(0)

                    wind_data = pd.DataFrame(hourly_power_output, columns=['P'])
                    wind_speed_index = corrected_wind_speed.index
                    wind_data.set_index(wind_speed_index, inplace=True)
                    # save series .csv In betz model it is saved the total power produced by the turbine (not the power per kW)
                    wind_data.to_csv(path + '/production/' + name_serie)
                    wind_data = np.array(wind_data['P'])
                    # Save new parameters in previous_simulation
                    with open(f"{directory}/wind_{file_structure}_{location_name}.pkl", 'wb') as f:
                        pickle.dump(self.parameters, f)
             
        else:
            # read a specific production serie expressed as kW/kWpeak
            wind_data = pd.read_csv(path+'/production/'+self.parameters['serie'])['P'].to_numpy()
            wind_data = wind_data * self.Npower                          # kWh
            self.production = np.tile(wind_data,int(c.timestep_number*c.timestep/60/8760))
            if len(self.production) != c.timestep_number:
                raise ValueError(f"Warning! Checks the length and timestep of the wind production you input for {location_name}.")

        # Aging calculation
        if parameters['ageing'] == False:
            # electricity produced every hour in the reference_year [kWh] [kW]
            self.production = np.tile(wind_data,int(c.timestep_number*c.timestep/60/8760)) # from 1 year to simlation length years
        else:
            self.degradation = parameters['degradation factor']
            n_years = int(c.timestep_number*c.timestep/60/8760)
            annual_ts_number = 60*8760/c.timestep
            self.production = np.tile(wind_data,n_years)   # no degradation
            for i in range(1,n_years+1):
                self.production[int(i*annual_ts_number):int((i+1)*annual_ts_number)] *= ((1-self.degradation/100)**i)  # apply degradation
        # from hourly to timestep
        if c.timestep < 60 and (self.parameters['serie'] == "TMY" or type(self.parameters['serie']) == int):
            self.production =  np.repeat(self.production, 60/c.timestep) # [kW] creating a production series alligned with selected timestep 

        
    def use(self,step):
        """
        Produce electricity
        
        step : int step to be simulated
        windspeed: float wind speed [m/s]
    
        output : float electricity produced that timestep [kW]
    
        """

        power_output = self.production[step]               
        return(power_output)


    def cpfunc(self,ws,area,beta,idx,cp_max):
        
        c1 = [0.73, 0.5, 0.5176, 0.77, 0.5,  0.22]
        c2 = [151., 116.,116.,151.,116.,120.]
        c3 = [0.58,0.4,0.4,0.,0.,0.4]
        c4 = [0.0,0.0,0.0,0.0,0.4,0.]
        c5 = [0.002,0.,0.,0.,0.,0.]        
        x =  [2.14,0.0,0.0,0.0,0.0,0.0]
        c6 = [13.2,5.0,5.0,13.65,5.0,5.0]
        c7 = [18.4,21.,21.,18.4,21.,12.5]
        c8 = [0.0,0.0,0.006795,0.0,0.0,0.0]
        c9 = [-0.02,0.089,0.089,0.0,0.08,0.08]
        c10 = [0.03,0.035,0.035,0.0,0.035,0.035]
        
        diam = (4*area/math.pi)**(1/2)
        
        # cp_max and lambda_opt for beta = 0
        lambdas_0 = np.arange(0.1,12.,0.1)
        cps_0 = c1[idx]*((c2[idx]/(1.0/((1.0/lambdas_0)-c10[idx])))-c6[idx])*np.exp(-c7[idx]/(1.0/((1.0/lambdas_0)-c10[idx])))+c8[idx]*lambdas_0 # Eq. 2(I) + 2(II) (beta=0)
        cp_max_0 = np.max(cps_0)
        lambda_opt_0 = lambdas_0[np.argmax(cps_0)]
                        
        # if omega_min and omega_max are not given as inputs:     
        if 'omega_min' in self.parameters:
            omega_min = self.parameters['omega_min']
        else:
            omega_min = 1046.558*diam**(-1.0911)*2*math.pi/60
            # omega_min = 188.8*diam**(-0.7081)*2*math.pi/60 # alternative equation proposed by Niccolò Baldi  
        
        if 'omega_max' in self.parameters:
            omega_max = self.parameters['omega_max']
        else:
            omega_max = 705.406*diam**(-0.8349)*2*math.pi/60
            # omega_max = 793.7*diam**(-0.8504)*2*math.pi/60 # alternative equation proposed by Niccolò Baldi

        omega = np.minimum(omega_max, np.maximum(omega_min, lambda_opt_0/(diam/2)*ws)) # Eq. 5, all speeds in [rad/s]

        # cp and lambda for the actual beta
        lambdaparam = 0.
        cp = 0
        if not(ws == 0.):
            lambdaparam = omega*(diam/2)/ws # Eq. 3
            cp = c1[idx]*(c2[idx]/(1.0/((lambdaparam+c9[idx]*beta)**(-1)-c10[idx]*((beta**3)+1.0)**(-1)))-c3[idx]*beta-c4[idx]*(1.0/((lambdaparam+c9[idx]*beta)**(-1)-c10[idx]*((beta**3)+1.0)**(-1)))*beta-c5[idx]*beta**x[idx]-c6[idx])*np.exp(-c7[idx]/(1.0/((lambdaparam+c9[idx]*beta)**(-1)-c10[idx]*((beta**3)+1.0)**(-1))))+c8[idx]*lambdaparam # Eq. 2 (I) + (II)
        
        # scaling cps wrt cp_max
        cp = cp_max/cp_max_0*cp

        return cp
    
    def hbandareas(self,area,Nbands):
        
        radius = (area/math.pi)**(1/2)
        
        area_bands = np.zeros(Nbands)
        area_tmp = 0.
        
        for i in range(int(Nbands/2)):
            
            dist_chord_center = radius-i*(2*radius/ Nbands)

            if Nbands % 2 == 0:
                theta = 2*math.acos(dist_chord_center/radius)        
                area_bands[i] = theta*180/math.pi/360*math.pi*radius**2-(radius*dist_chord_center*(math.sin(theta/2)))-area_tmp  # formula per il calcolo dell'area di un segmento circolare a due basi  
                area_bands[Nbands-1-i] = area_bands[i]
                area_tmp += area_bands[i]
            else:
                if dist_chord_center>0:
                    theta = 2*math.acos(dist_chord_center/radius)        
                    area_bands[i] = theta*180/math.pi/360*math.pi*radius**2-(radius*dist_chord_center*(math.sin(theta/2)))-area_tmp # formula per il calcolo dell'area di un segmento circolare a due basi  
                    area_bands[Nbands-1-i] = area_bands[i]
                    area_tmp += area_bands[i]
                else:
                    area_bands[i] = area-2*(sum(area_bands[0:i-1]))

        return area_bands
    
    def windshear(self,ws,z_hub,z_i,alpha):
        
        ws_shear_i = ws*(z_i/z_hub)**alpha # Eq. 7 wind speed corrected considering vertical wind profile for each h band
        
        return ws_shear_i
    
    def windveer(self,z_hub,z_i,Vu):
        
        deltaphi_i = Vu*(z_i-z_hub)  # Eq. 8
        
        return deltaphi_i
    
    def eqspeed(self,ws,z_i,z_hub,alpha,area,Vu,Nbands):
        
        diam = (4*area/math.pi)**(1/2)
        h = np.zeros(Nbands)
        ws_eq = 0.
        
        for i in range(Nbands):
            h[i] = z_hub-diam/2+diam/(Nbands*2)+i*diam/Nbands # height of the barycenter of the various h bands at which compute U_i e DeltaPhi_i
            
        area_bands = self.hbandareas(area,Nbands)
        
        for i in range(Nbands):
            ws_shear_i = self.windshear(ws,z_hub,z_i,alpha)
            deltaphi_i = self.windveer(z_hub,z_i,Vu)
            ws_eq = ws_eq+(area_bands[i]/area*(ws_shear_i*math.cos(deltaphi_i*2*math.pi/360))**3) # Eq. 9 (senza radice cubica)
           
        ws_eq = ws_eq**(1./3.)
        
        return ws_eq

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

        size = self.Npower # [kW] wind techology rated power
         
        if tech_cost['cost per unit'] == 'default price correlation':
            C0 = 1270 # €/kW
            scale_factor = 0.8 # 0:1
            C = C0 * size **  scale_factor
        else:
            C = size * tech_cost['cost per unit']

        tech_cost['total cost'] = tech_cost.pop('cost per unit')
        tech_cost['total cost'] = C
        tech_cost['OeM'] = tech_cost['OeM'] *C /100 # €


        self.cost = tech_cost

#%%##########################################################################################

if __name__ == "__main__":
    
    """
    Functional test
    """

