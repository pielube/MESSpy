import numpy as np
import pandas as pd
import math
import warnings
import os
import sys 
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),os.path.pardir)))   # temorarily adding constants module path 
from core import constants as c
import matplotlib.pyplot as plt

class wind:    
    
    def __init__(self,parameters,path=False,timestep_number=False,timestep=False,simulation_years=False):
        """
        Create a wind object based on the specified model
    
        parameters : dictionary
            'model': str type of model to be used for wind
                    betz -> simple model based on Betz theory
                    detailed -> more detailed model based on
                    Saint-Drenan, Yves-Marie, et al. 
                    "A parametric model for wind turbine power curves incorporating environmental conditions." 
                    Renewable Energy 157 (2020): 754-768.
                    simple -> wind production data retrieved from https://www.renewables.ninja/.
                              When using this model, only one more parameter needs 
                              to be defined: 'Npower'.
                              
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
            
            'z_i': [m] wind turbine height, (?)
            'z_hub': [m] hub height, e.g. 30 m (Aircon 10/10 kW)
            'alpha': [-] Hellman or shear coefficient, values from 0 to 0.4
            'Vu': [°/m] Veer coefficient, values from 0 to 0.75 °/m
            'Nbands': [-] Number of horizontal bands
              
        timestep_number : int number of timesteps considered in the simulation

        general: dictionary
            see rec.py

        output : wind object able to:
            produce electricity .use(step)
        """
        self.parameters         = parameters
        self.Npower             = self.parameters['Npower']     # [kW] wind technology nominal power
        self.model              = self.parameters['model']      # [-]  wind technology model
        self.cost               = False                         # will be updated with tec_cost()
        self.property           = self.parameters['owned']      # bool value to take into account if the plant is owned or only electricity purchase is considered. It only impacts impact on economic assessment and key parameters
        if __name__ == "__main__":                          # if code is being executed from chp_gt.py script
            self.timestep           = timestep              # [min] simulation timestep if launched from main
            self.timestep_number    = timestep_number       # [-] number of timesteps considered in the simulation
            self.simulation_years   = 1                     # [-] time horizon for the considered data series - important when reading production data from a file. For the sake of simplicity, in functional test, only one year of power production data is considered. Regardless of the selected timestep
        else:
            self.timestep           = c.timestep            # [min] simulation timestep if script is launched from wind.py
            self.timestep_number    = c.timestep_number     # [-] number of timesteps considered in the simulation
            self.simulation_years   = c.simulation_years 
        
        if self.model not in ['betz','detailed','simple']:
            raise ValueError("Warning: selected model for wind techonology is not among the available options.\n\
            Option to fix the problem: \n\
                (a) - In studycase.py file choose one among 'betz','detailed' and 'simple' models.")

        if self.model == 'simple':
            if self.simulation_years != 1: 
                warning_message = f"Warning: In current case study {self.simulation_years} have been specified. \n\
                                    Model 'simple' is selected for wind power: power production data series is provided as input data \n\
                                    If analysis aim is to include wind power production and variation over more than one year, a multi-year data set must be provided. \n\
                                    Check consistency."
                warnings.warn(warning_message, UserWarning)
            self.series = self.parameters["series"]  # selected dataset for wind hourly production
                           
            'Data Extraction - Wind power generation'
            main = False                        # check
            if __name__ == "__main__":          # if code is being executed from wind.py file
                directory = r'../input_dev'     # check for input folder
                if os.path.exists(directory):   # if 'input_dev' exists
                    pass
                else:                           # if not check in 'input_test'
                    directory = r'../input_test'
                os.chdir(directory +'/production')          
            else: 
                os.chdir(f"{path}\\production") # if code is being executed from main
                main = True
            
            self.wind_prod  = (pd.read_csv(self.series,usecols=["kW"]).values).reshape(-1,)             # [-] power production. Importing wind production data series for the selected location. Expressed as ratio kWprod/1KWrated
            if self.timestep < 60 and len(self.wind_prod) != self.timestep_number:                      # check if selected timestep is smaller than hourly and if the provided production series (csv file) is already in the desired form (number of timestep)
                self.wind_prod = np.repeat(self.wind_prod, 60/self.timestep)                            # [-] creating a production series alligned with selected timestep 
            self.prod_1kw   = np.tile(self.wind_prod,int(self.timestep_number*self.timestep/60/8760))   # [-] creating the production series needed for the entire simulation in case 1 typical year for production is considered. From 1 year to the duration of the simulation
            self.wprod      = self.prod_1kw*self.Npower                                                 # [kW] wind power production data series
            
            if main ==  True:                 # if code is being executed from main, change directory back to main
                os.chdir(r'../..')

        
    def use(self,step,ws_input):
        """
        Produce electricity
        
        step : int step to be simulated
        windspeed: float wind speed [m/s]
    
        output : float electricity produced that timestep [kW]
    
        """
        if self.model == 'simple':
            power_output = self.wprod[step]
        
        else:
            ws_turbine=ws_input
            rho=1.225 # [kg/m3] air density assumed constant. Must be upgraded to be a function of external weather conditions. 
            
            if self.model == 'betz':
                            
                if ws_turbine<self.parameters['WScutin']:
                    ws_turbine = 0
                elif self.parameters['WSrated']<ws_turbine<self.parameters['WScutout']:
                    ws_turbine = self.parameters['WSrated']
                elif ws_turbine>self.parameters['WScutout']:
                   ws_turbine = 0
                
                power_output = 0.5*rho*self.parameters['area']*ws_turbine**3*self.parameters['efficiency']/1000 # [kW] power generation in the considered timestep
                
            elif self.model == 'detailed':
                
                if ws_turbine<self.parameters['WScutin']:
                    ws_turbine = 0.
                elif self.parameters['WSrated']<ws_turbine<self.parameters['WScutout']:
                    ws_turbine = self.parameters['WSrated']
                elif ws_turbine>self.parameters['WScutout']:
                   ws_turbine = 0.
                
                cp_max = 0.44
                
                if 'cp_max' in self.parameters:
                    cp_max = self.parameters['cp_max']
                    
                powercoeff = self.cpfunc(ws_turbine,self.parameters['area'],self.parameters['beta'],self.parameters['idx'],cp_max=cp_max)
                wseq = self.eqspeed(ws_turbine,self.parameters['z_i'],self.parameters['z_hub'],self.parameters['alpha'],self.parameters['area'],self.parameters['Vu'],self.parameters['Nbands'])
                
                power_output = 0.5*rho*self.parameters['area']*wseq**3*powercoeff/1000 # [kW] power generation in the considered timestep
                
        return(power_output)


    def cpfunc(self,ws,area,beta,idx,cp_max=0.44):
        
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
    
    inp_test = {  "model"   : "simple",
                  "series"  : "15windproduction.csv",
                  "Npower"  : 2000,
                  "owned"   : True}
    
    timestep        = 15                # [min] selected timestep for the simulation
    timestep_number = 60/timestep*8760  # [-] number of steps for the simulation                 
    
    wind = wind(inp_test,timestep_number=timestep_number,timestep=timestep)         # creating wind turbine/wind park object
    
    if inp_test['model'] == 'simple':
        fig, ax = plt.subplots(figsize=(10,5),dpi=600)
        ax.plot(wind.wprod,linewidth=1)
        ax.set_ylabel('Produced power [kW]')
        ax.set_xlabel('Time [year]')
        ax.set_title('Wind production over time')
        xticks = list(np.linspace(0, len(wind.wprod), 13).astype(int))
        xticklabels = ['            Jan','            Feb','             Mar','            Apr','            May','             Jun','           Jul','             Aug','           Sep','          Oct','          Nov','           Dec','']
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)
        ax.grid(alpha=0.3)
        
    
