import CoolProp
from CoolProp.CoolProp import PropsSI
from CoolProp.Plots import PropertyPlot
from scipy.interpolate import interp1d
import warnings
import math as m
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import os
import sys 
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),os.path.pardir)))   # temporarily adding constants module path 
from core import constants as c

class Compressor:
    
    def __init__(self,parameters,timestep_number,timestep=False,maxflowrate_ele=False):
        """
        Create a compressor object
    
        parameters : dictionary
            'compressor model'  :'simple_compressor',                       # specific consumption calculation based on interpolation on a given dataset 
                                 'normal_compressor',                       # single-stage compressor without refrigeration - detailed thermodynamic calculations
                                 'compressor_with_refrigeration',           # single-stage compressor with refrigeration - detailed thermodynamic calculations
                                 'multistage_compressor_with_refrigeration' # multi-stage compressor with refrigeration - detailed thermodynamic calculations
            'P_out'                 : max pressure [bar]          
            'P_in'                  : inlet pressure [bar]
            'T_in'                  : inlet Temperature [K]
            'Npower'                : power [kW]                    # Nominal power of compressor, to be provided if 'flow_rate' is not provided 
            'flow_rate'             : hydrogen flow rate [kg/s]     # Nominal mass flow rate of compressor, to be provided if 'Npower' is not provided 
            'pressure losses IC'    : pressure losses in heat exchangers [%]
            'T_IC'                  : Temperature of intercooler [K]
            'n_stages'              : number of compression stages
            'only_renewables'       : operational strategy. Working with only renewable energy or not, boolean value
        
        timestep_number : int number of timesteps considered in the simulation [-]
        timestep        : int number of minutes considered as temporal resolution [min]
        maxflowrate     : float maximum mass flow rate produced by the upstream electrolysis system [kg/s]
        
        """
        
        self.model              = parameters['compressor model']
        self.fluid              = parameters['fluid']                   # [-] type of selected fluid 
        self.only_renewables    = parameters['only_renewables']
        self.P_in               = parameters['P_in']
        self.P_out              = parameters['P_out']
        self.hyd                = np.zeros(timestep_number)             # [kg/s] hydrogen flow rate
        self.T_in               = parameters['T_in']                    # [°C] fluid inlet temperature
        self.conversion         = c.kWh2kJ                              # [kJ/kWh]
        self.eta_motor          = 0.95                                  # [-] assumed efficiency of electric motor driving the compressor https://transitionaccelerator.ca/wp-content/uploads/2023/04/TA-Technical-Brief-1.1_TEEA-Hydrogen-Compression_PUBLISHED.pdf
        if parameters.get('Npower') and parameters.get('Nflow_rate'):
            raise ValueError("Warning: redundancy issue. Both nominal power and nominal flow rate have been defined as input.\n\
            Option to fix the problem: \n\
                (a) - Choose only one as input parameter")
        if 'Npower' in parameters       :  self.Npower      = parameters['Npower']      # [kW] nominal power of compressor
        if 'Nflow_rate' in parameters   :  self.Nflowrate   = parameters['Nflow_rate']  # [kg/s]  nominal flow rate that can be processed by the compressor
        if maxflowrate_ele:             # if code is being executed from main maxflowrate_ele is a given parameter: maximum hydrogen producible by the electrolyzers. Compressor functioning only in conjuction with electorlyzers. 
            if parameters.get('Npower') and parameters.get('Npower') < maxflowrate_ele:
                warning_message = f"Warning: The specified flow rate is lower than maximum flow rate producible by electorlysers. \n\
                                    Compressor 'flow_rate' parameter should be increased to at least {maxflowrate_ele:.5f} kg/s in studycase file."
                warnings.warn(warning_message, UserWarning)
                
        if __name__ == "__main__":                                      # if code is being executed from compressor.py script
            self.timestep   = timestep                                  # [min]       simulation timestep if launched from current script
        else:
            self.timestep   = c.timestep                                # [min]       simulation timestep if launched from main
 
        if self.model == 'simple_compressor':
            '''
            Specific power consumption data
            Ref. FCH Study - https://hsweb.hs.uni-hamburg.de/projects/star-formation/hydrogen/P2H_Full_Study_FCHJU.pdf
            
            Power consumption is estimated from adiabatic compression model with 
            an efficiency of 50%. This efficiency considers the efficiency of 
            electrical power transformation and auxiliary systems such as the 
            cooling circuit.
            
            '''
            if self.fluid != 'Hydrogen':  
                raise ValueError(f"Warning: {self.fluid} has been selected as working fluid. 'simple_compressor' model is only valid for operation with hydrogen.\n\
                Option to fix the problem: \n\
                    (a) - Insert 'Hydrogen' as working fluid when defining compressor parameters")
            if parameters.get('Npower') and parameters.get('Nflow_rate'): # if both parameters exist and 
                raise ValueError("Warning: if 'simple_compressor' model has been selected neither the size nor the nominal flow rate need to be defined. Sizing results as a consequence of maximum hydrogen production capacity of electrolizers.\n\
                Option to fix the problem: \n\
                    (a) - Either remove 'Npower' and 'flow_rate' keys from parameters or assign false/0 value")
                   
            self.Nflowrate = maxflowrate_ele  # [kg/s] maximum flow rate that can be processed by the compressor - size is automatically defined. Must be defined to initialise compressor model, check already present in location.py.
            
            # Original set of data from which linear approximation is derived
            powerconsumption = {
                                'Pressure in [barg]'            : [0, 0, 15, 15, 30, 30, 60, 60],
                                'Pressure out [barg]'           : [200, 500, 200, 500, 200, 500, 200, 500],
                                'Pressure difference [bar]'     : [200, 500, 185, 485, 170, 470, 140, 440],
                                'Energy consumption [kWh/kg]'   : [5, 6.3, 2.4, 3.5, 1.7, 2.7, 1.1, 2]                
                                }
            # Slope and intercetp values calculated on 'power consumption' dataset
            coefficients    = {
                                'Pressure in [barg]': [0,15,30,60],
                                'Slope'             : [0.00433,0.00367,0.0033,0.003],
                                'Intercept'         : [4.1333,1.7216,1.1333,0.68]
                                }
            # creating a DataFrame to better handle the data
            df = pd.DataFrame(coefficients)
            
            slope       = df.loc[df['Pressure in [barg]'] == self.P_in, 'Slope'].iloc[0]
            intercept   = df.loc[df['Pressure in [barg]'] == self.P_in, 'Intercept'].iloc[0]
            
            self.en_cons     = (intercept+slope*(self.P_out-self.P_in))*self.conversion     # [kJ/kgH2] specific energy consumption of the simple compressor model
        
        else:  # if any other model is selected

            # initialisation of variables - thermodynamic points of compression process
            self.P_points   = [] # [bar] Pressure
            self.T_points   = [] # [K] Temperature
            self.h_points   = [] # [kJ/kg] Enthalpy
            self.s_points   = [] # [kJ/kgK] Entropy
            self.rho_points = [] # [kg/m^3] Density
                        
            # thermodynamic and chemical properties
            self.hydrogen_HHV   = c.HHVH2*1000                  # [kJ/kg]
            self.hydrogen_LHV   = c.LHVH2*1000                  # [kJ/kg]
            self.CP_air         = c.CP_AIR                      # [kJ/kgK]
            self.CV_air         = c.CV_AIR                      # [kJ/kgK]
            self.CP_H2          = c.CP_H2                       # [kJ/kg*K] at T=25 °C e P=1 atm
            self.CV_H2          = c.CV_H2                       # [kJ/kg*K] at T=25 °C e P=1 atm
            self.R_univ         = c.R_UNIVERSAL                 # [J/mol*K] Molar ideal gas constant
            self.M_H2           = c.H2MOLMASS*1000              # [kg/kmol]                   
            self.M_air          = c.AIRMOLMASS*1000             # [kg/kmol]
            self.R_H2           = self.R_univ/(self.M_H2/1000)  # [J/kgK]
            self.R_air          = self.R_univ/(self.M_air/1000) # [J/kgK]
            
            # calculation of polytropic exponents hydrogen and air 
            self.beta_max       = self.P_out/self.P_in                      # [-] compression ratio
            self.eta_pol        = 0.70                                      # [-] Ref: The TechnoEconomics of Hydrogen Compression 
            self.gamma          = self.CP_H2/self.CV_H2                     # [-]
            self.gamma_air      = self.CP_air/self.CV_air                   # [-]
            self.epsilon        = (self.gamma-1)/self.gamma                 # [-]
            self.epsilon_air    = (self.gamma_air-1)/self.gamma_air         # [-]
            self.exp_pol_H2     = (self.eta_pol/self.epsilon)/(self.eta_pol/self.epsilon-1)         # [-]
            self.exp_pol_air    = (self.eta_pol/self.epsilon_air)/(self.eta_pol/self.epsilon_air-1) # [-]
            self.omega          = (self.exp_pol_H2-1)/self.exp_pol_H2       # [-]
            self.omega_air      = (self.exp_pol_air-1)/self.exp_pol_air     # [-]
            
            ####################################################################################
            # Single-stage compressor operation without refrigeration
            if self.model == 'normal_compressor': 
                self.n_stages   = 1
                self.delta_P    = 0      # [-] pressure losses 
                
                self.P_target = np.zeros(self.n_stages)
                
                for i in range(self.n_stages*2):          # point 0 is the compressor input, point 1 is the compressor output
    
                    self.P_points.append(self.P_in)
                    self.T_points.append(self.T_in)
                    self.h_points.append(PropsSI('H', 'T', self.T_points[i], 'P', self.P_points[i]*100000, self.fluid)/1000)    # [kJ/kg]
                    self.s_points.append(PropsSI('S', 'T', self.T_points[i], 'P', self.P_points[i]*100000, self.fluid)/1000)    # [kJ/kgK]
                    self.rho_points.append(PropsSI('D', 'T', self.T_points[i], 'P', self.P_points[i]*100000, self.fluid))       # [kg/m^3]
                
                self.comp_power_list    = []
                self.comp_lav_spec      = []          
                self.comp_beta_targ     = self.beta_max**(1/self.n_stages)
                                               
                self.comp_power_list.append(0.)
                self.comp_lav_spec.append(0.)
                self.P_target = self.P_in*self.comp_beta_targ**self.n_stages
                P_in    = self.P_points[0]
                h_in    = self.h_points[0]
                s_in    = self.s_points[0]
                T_in    = self.T_points[0]
                rho_in  = self.rho_points[0]
                
                # compression from 0 to 1
                P_out       = self.comp_beta_targ*P_in
                s_out_iso   = s_in
                h_out_iso   = PropsSI('H', 'P', P_out*100000, 'S', s_out_iso*1000, self.fluid)/1000     # [kJ/kg]
                self.T_is   = PropsSI('T', 'P', P_out*100000, 'S', s_out_iso*1000, self.fluid)
                self.eta_is = (self.comp_beta_targ**self.epsilon - 1)/(self.comp_beta_targ**(self.epsilon/self.eta_pol)-1)
                h_out       = (h_out_iso - h_in) / self.eta_is + h_in # [kJ/kg]
                
                self.T_points[1]    = PropsSI('T', 'P', P_out*100000, 'H', h_out*1000, self.fluid)
                self.P_points[1]    = P_out
                self.h_points[1]    = h_out             # [kJ/kg]
                self.s_points[1]    = PropsSI('S', 'P', P_out*100000, 'H', h_out*1000, self.fluid)/1000
                self.rho_points[1]  = PropsSI('D', 'P', P_out*100000, 'H', h_out*1000, self.fluid)
                
                self.comp_lav_spec[0] = (h_out-h_in)    # [kJ/kg] specific work of compression given in and out conditions
                
                if self.Nflowrate:      # if nominal mass flow rate is defined as input - compressor nominal power is defined as a consequence
                    self.Npower         = (self.Nflowrate*self.comp_lav_spec[0])/self.eta_motor     # [kW] compressor nominal power - simple assumption
                    
                elif self.Npower:       # if nominal power is defined as input - compressor mass flow rate is defined as a consequence
                    self.Nflowrate      = self.Npower*self.eta_motor/self.comp_lav_spec[0]          # [kg/s] mass flow rate - simple assumption
                    
                else:                   # if neither the size nor the nominal flow rate have been defined to initialise the compressor object
                    self.Nflowrate      = maxflowrate_ele
                    self.Npower         = (self.Nflowrate*self.comp_lav_spec[0])/self.eta_motor     # [kW] compressor nominal power - simple assumption

            ################################################################################################
            # Single-stage compressor operation with refrigeration 
            if self.model == 'compressor_with_refrigeration': 
                self.n_stages   = 1                                 # [-] single stage
                self.T_IC       = parameters['T_IC']                # [K] fluid temperature after refrigeration
                self.delta_P    = parameters['pressure losses IC']  # [%] 0-1 pressure losses linked to refrigeration
                self.P_target   = np.zeros(self.n_stages)
                
                if self.beta_max > 4:
                    warning_message = f"Warning: Desired compression ratio {round(self.beta_max,2)} is considered too high for a single stage compression. \n\
                Consider adopting 'multistage_compressor_with_refrigeration' model instead ."
                    warnings.warn(warning_message, UserWarning)
                
                for i in range(self.n_stages*2+1):           # point 0 is the compressor input, point 1 is the compressor output
    
                    self.P_points.append(self.P_in)
                    self.T_points.append(self.T_in)
                    self.h_points.append(PropsSI('H', 'T', self.T_points[i], 'P', self.P_points[i]*100000, self.fluid)/1000)
                    self.s_points.append(PropsSI('S', 'T', self.T_points[i], 'P', self.P_points[i]*100000, self.fluid)/1000)
                    self.rho_points.append(PropsSI('D', 'T', self.T_points[i], 'P', self.P_points[i]*100000, self.fluid))
                
                self.comp_power_list    = []
                self.comp_lav_spec      = []          
                self.IC_power_list      = []
                self.comp_beta_targ     = self.beta_max**(1/self.n_stages)       
                                        
                self.comp_lav_spec.append(0.)
                self.comp_power_list.append(0.)
                self.IC_power_list.append(0.)
                self.P_target = self.P_in*self.comp_beta_targ**self.n_stages
                
                P_in    = self.P_points[0]
                h_in    = self.h_points[0]
                s_in    = self.s_points[0]
                T_in    = self.T_points[0]
                rho_in  = self.rho_points[0]
                
                # compression 0 - 1
                P_out_ic        = self.P_target
                P_out           = (P_out_ic)/(1-self.delta_P)
                self.beta_new   = P_out/P_in
                s_out_iso       = s_in
                h_out_iso       = PropsSI('H', 'P', P_out*100000, 'S', s_out_iso*1000, self.fluid)/1000
                self.T_is       = PropsSI('T', 'P', P_out*100000, 'S', s_out_iso*1000, self.fluid)
                self.eta_is     = (self.comp_beta_targ**self.epsilon - 1)/(self.comp_beta_targ**(self.epsilon/self.eta_pol)-1)
                h_out           = (h_out_iso-h_in)/self.eta_is+h_in
                
                self.T_points[1]    = PropsSI('T', 'P', P_out*100000, 'H', h_out*1000, self.fluid)
                self.P_points[1]    = P_out
                self.h_points[1]    = h_out
                self.s_points[1]    = PropsSI('S', 'P', P_out*100000, 'H', h_out*1000, self.fluid)/1000
                self.rho_points[1]  = PropsSI('D', 'P', P_out*100000, 'H', h_out*1000, self.fluid)
                
                # interrefrigeration 1 - 2
                self.P_points[2]    = P_out_ic
                self.T_points[2]    = self.T_IC
                self.h_points[2]    = PropsSI('H', 'P', P_out_ic*100000, 'T', self.T_IC, self.fluid)/1000
                self.s_points[2]    = PropsSI('S', 'P', P_out_ic*100000, 'T', self.T_IC, self.fluid)/1000
                self.rho_points[2]  = PropsSI('D', 'P', P_out_ic*100000, 'T', self.T_IC, self.fluid)
                
                h_out_ic = PropsSI('H', 'P', P_out_ic*100000, 'T', self.T_points[2], self.fluid)/1000  # [kJ/kg]
                
                self.comp_lav_spec[0]   = h_out-h_in      # [kJ/kg] specific work of compression given in and out conditions
                self.delta_H            = h_out-h_out_ic
                
                # Alternative calculation for compressor power (without IC at last compressor) Ref:The TechnoEconomics of Hydrogen Compression TECHNICAL BRIEF
                beta_stage  = self.beta_new              # [-] diaphragm compressor
                eta_is      = self.eta_is                # [-]
                T_out       = self.T_in*(1+((self.P_out/self.P_in)**(self.epsilon)-1)/eta_is)   # [K]
                P_avg       = (self.P_out+self.P_in)/2                                          # [bar]
                T_avg       = (T_out+self.T_in)/2                                               # [K] 
                Z           = PropsSI('Z', 'T', T_avg, 'P', P_avg*100000,self.fluid)            # [-] Compressibility Factor
                MM          = PropsSI('MOLAR_MASS',self.fluid)                                  # [kg/mol]
                a           = (1/self.epsilon)*(Z/eta_is)*self.T_in*self.R_univ*((self.P_out/self.P_in)**(self.epsilon)-1)
                    
                if self.Nflowrate:      # if flow rate is defined as input - compressor nominal power is defined as a consequence
                    q_m                 = self.Nflowrate/MM                 # [mol/s] molar flow rate
                    Power               = q_m*a/1000                        # [kW] https://transitionaccelerator.ca/wp-content/uploads/2023/04/TA-Technical-Brief-1.1_TEEA-Hydrogen-Compression_PUBLISHED.pdf
                    self.Npower         = m.ceil(Power/self.eta_motor)      # [kW] compressor nominal power 
                    self.IC_power_list  = self.Nflowrate*self.delta_H       # [kW] Total heat to be removed by the cooling system
                elif self.Npower:       # if nominal power is defined as input - compressor mass flow rate is defined as a consequence
                    Power               = self.Npower*self.eta_motor        # [kW]
                    q_m                 = Power/a*1000                      # [mol/s] molar flow rate
                    self.Nflowrate      = round(q_m*MM,6)                   # [kg/s]
                    self.IC_power_list  = self.Nflowrate*self.delta_H       # [kW] Total heat to be removed by the cooling system
                else:                   # if neither the size nor the nominal flow rate have been defined to initialise the compressor object
                    self.Nflowrate      = maxflowrate_ele
                    q_m                 = self.Nflowrate/MM                 # [mol/s] molar flow rate
                    Power               = q_m*a/1000                        # [kW] https://transitionaccelerator.ca/wp-content/uploads/2023/04/TA-Technical-Brief-1.1_TEEA-Hydrogen-Compression_PUBLISHED.pdf
                    self.Npower         = m.ceil(Power/self.eta_motor)      # [kW] compressor nominal power 
                    self.IC_power_list  = self.Nflowrate*self.delta_H       # [kW] Total heat to be removed by the cooling system
    
            ########################################################################################à
            # Multistage compressor operation with refrigeration
            if parameters['compressor model'] == 'multistage_compressor_with_refrigeration':
                
                self.T_IC           = parameters['T_IC']
                self.delta_P        = parameters['pressure losses IC']  
                self.n_stages       = parameters['n_stages']  
                self.T_is_points    = []
                self.P_target       = np.zeros(self.n_stages)
                self.eta_is         = np.zeros(self.n_stages)
                self.beta_new       = np.zeros(self.n_stages)
                self.delta_H        = np.zeros(self.n_stages)
                
                for i in range(self.n_stages*2 + 1):            # point 0 is the input to the compressor, point 1 is the outlet of the compressor
    
                    self.P_points.append(self.P_in)
                    self.T_points.append(self.T_in)
                    self.h_points.append(PropsSI('H', 'T', self.T_points[i], 'P', self.P_points[i]*100000, self.fluid)/1000)
                    self.s_points.append(PropsSI('S', 'T', self.T_points[i], 'P', self.P_points[i]*100000, self.fluid)/1000)
                    self.rho_points.append(PropsSI('D', 'T', self.T_points[i], 'P', self.P_points[i]*100000, self.fluid))
                
                self.comp_lav_spec      = []     
                self.comp_power_list    = []
                self.IC_power_list      = []
                self.comp_beta_targ     = np.power(self.beta_max, 1 / self.n_stages)
                
                # P_ref       = 2
                # T_ref       = 298.15
                # deltaT_app  = 5
                                            
                for i in range(self.n_stages):
                                               
                    self.comp_lav_spec.append(0.)
                    self.comp_power_list.append(0.)
                    self.IC_power_list.append(0.)
                    self.P_target[i]    = self.P_in*self.comp_beta_targ**(i+1)
                    self.delta_H        = np.zeros(self.n_stages)
                    P_in                = self.P_points[i*2]
                    h_in                = self.h_points[i*2]
                    s_in                = self.s_points[i*2]
                    T_in                = self.T_points[i*2]
                    rho_in              = self.rho_points[i*2]  
                    
                    # compression stage i
                    P_out_ic            = self.P_target[i]
                    P_out               = (P_out_ic)/(1-self.delta_P)
                    self.beta_new[i]    = P_out/P_in
                    s_out_iso           = s_in
                    h_out_iso           = PropsSI('H', 'P', P_out*100000, 'S', s_out_iso*1000, self.fluid)/1000
                    self.T_is_points.append(PropsSI('T', 'P', P_out*100000, 'S', s_out_iso*1000, self.fluid))
                    self.eta_is[i]      = (self.beta_new[i]**self.epsilon - 1)/(self.beta_new[i]**(self.epsilon/self.eta_pol)-1)
                    h_out               = (h_out_iso - h_in) / self.eta_is[i] + h_in
                    
                    self.T_points[i*2 + 1]      = PropsSI('T', 'P', P_out*100000, 'H', h_out*1000, self.fluid)
                    self.P_points[i*2 + 1]      = P_out
                    self.h_points[i*2 + 1]      = h_out
                    self.s_points[i*2 + 1]      = PropsSI('S', 'P', P_out*100000, 'H', h_out*1000, self.fluid)/1000
                    self.rho_points[i*2 + 1]    = PropsSI('D', 'P', P_out*100000, 'H', h_out*1000, self.fluid)
                    
                    # interrefrigeration stage i
                    self.P_points[i*2 + 2]      = P_out_ic
                    self.T_points[i*2 + 2]      = self.T_IC
                    self.h_points[i*2 + 2]      = PropsSI('H', 'P', P_out_ic*100000, 'T', self.T_IC, self.fluid)/1000
                    self.s_points[i*2 + 2]      = PropsSI('S', 'P', P_out_ic*100000, 'T', self.T_IC, self.fluid)/1000
                    self.rho_points[i*2 + 2]    = PropsSI('D', 'P', P_out_ic*100000, 'T', self.T_IC, self.fluid)
                    
                    h_out_ic = PropsSI('H', 'P', P_out_ic*100000, 'T', self.T_IC, self.fluid)/1000
                    self.delta_H[i]             = h_out - h_out_ic
                    self.comp_lav_spec[i]       = h_out - h_in
                    self.delta_H                = sum(self.delta_H)
                
                # Alternative calculation for compressor power (without IC at last compressor) Ref:The TechnoEconomics of Hydrogen Compression TECHNICAL BRIEF
                beta_stage  = self.beta_new[0]              # [-] diaphragm compressor
                eta_is      = self.eta_is[0]                # [-]
                N           = round(m.log(self.P_out/self.P_in)/(m.log(beta_stage)))            # [-] Number of compressor stages
                T_out       = self.T_in*(1+((self.P_out/self.P_in)**(self.epsilon/N)-1)/eta_is) # [K]
                P_avg       = (self.P_out+self.P_in)/2                                          # [bar]
                T_avg       = (T_out+self.T_in)/2                                               # [K] 
                Z           = PropsSI('Z', 'T', T_avg, 'P', P_avg*100000,self.fluid)            # [-] Compressibility Factor
                MM          = PropsSI('MOLAR_MASS',self.fluid)                                  # [kg/mol]
                a           = N*(1/self.epsilon)*(Z/eta_is)*self.T_in*self.R_univ*((self.P_out/self.P_in)**(self.epsilon/N)-1)
                    
                if self.Nflowrate:      # if flow rate is defined as input - compressor nominal power is defined as a consequence
                    q_m                 = self.Nflowrate/MM             # [mol/s] molar flow rate
                    Power               = q_m*a/1000                    # [kW] https://transitionaccelerator.ca/wp-content/uploads/2023/04/TA-Technical-Brief-1.1_TEEA-Hydrogen-Compression_PUBLISHED.pdf
                    self.Npower         = m.ceil(Power/self.eta_motor)  # [kW] compressor nominal power 
                    self.IC_power_list  = self.Nflowrate*self.delta_H   # [kW] Total heat to be removed by the cooling system
                elif self.Npower:       # if nominal power is defined as input - compressor mass flow rate is defined as a consequence
                    Power               = self.Npower*self.eta_motor    # [kW]
                    q_m                 = Power/a*1000                  # [mol/s] molar flow rate
                    self.Nflowrate      = round(q_m*MM,6)               # [kg/s]
                    self.IC_power_list  = self.Nflowrate*self.delta_H   # [kW] Total heat to be removed by the cooling system
                else:                   # if neither the size nor the nominal flow rate have been defined to initialise the compressor object
                    self.Nflowrate      = maxflowrate_ele
                    q_m                 = self.Nflowrate/MM             # [mol/s] molar flow rate
                    Power               = q_m*a/1000                    # [kW] https://transitionaccelerator.ca/wp-content/uploads/2023/04/TA-Technical-Brief-1.1_TEEA-Hydrogen-Compression_PUBLISHED.pdf
                    self.Npower         = m.ceil(Power/self.eta_motor)  # [kW] compressor nominal power 
                    self.IC_power_list  = self.Nflowrate*self.delta_H   # [kW] Total heat to be removed by the cooling system
    ##########################################################################################################################################################
                    
                # # Hydrogen Refueling Station (HRS) application - isoentalpic transformation, control from T to dispenser to check if temperature limits are met
                # Pout_dispenser = 350  # [bar]
                # self.P_points.append(Pout_dispenser) #Pout_dispenser bar, serbatoio trenio
                # self.h_points.append(h_out_ic)
                # self.T_points.append(PropsSI('T', 'P', Pout_dispenser*100000, 'H', h_out_ic*1000, self.fluid))
                # self.s_points.append(PropsSI('S', 'P', Pout_dispenser*100000, 'H', h_out_ic*1000, self.fluid)/1000)
                # self.rho_points.append(PropsSI('D', 'P', Pout_dispenser*100000, 'H', h_out_ic*1000, self.fluid))
                
    ################################################################################################################################################################
                    
    def thermodynamics_points(self):
            
        """
        Representation of thermodynamic points
        
        """
        sim_step= 1
        if self.model == 'normal compressor':
            print('Simple adiabatic transformation')
        else:
            points = self.n_stages*2+1

            #Parameters:
            # fluid_name (string or AbstractState) – The name of the fluid to be plotted or a state instance
            # graph_type (string) – The graph type to be plotted, like “PH” or “TS”
            # axis (matplotlib.pyplot.gca(), Optional) – The current axis system to be plotted to. Default: create a new axis system
            # fig (matplotlib.pyplot.figure(), Optional) – The current figure to be plotted to. Default: create a new figure
            # unit_system (string, ['EUR','KSI','SI']) – Select the units used for the plotting. ‘EUR’ is bar, kJ, C; ‘KSI’ is kPa, kJ, K; ‘SI’ is Pa, J, K
            # tp_limits (string, ['NONE','DEF','ACHP','ORC']) – Select the limits in T and p.
            
            comp                = Compressor(inp_test, sim_step)
            p_values_round      = [m.ceil(val) for val in comp.P_points]           # rounding pressure values 
            unique_vals,indices = np.unique(p_values_round,return_index=True)
            max_y               = 300 # [°C] y-axis limit

            plt.rcParams['figure.dpi'] = 1000
            plot = PropertyPlot('Hydrogen', 'Ts', unit_system='EUR')
            plot.set_axis_limits([0, 50, -250, max_y])
            plot.calc_isolines(CoolProp.iQ, num=11)
            plot.calc_isolines(CoolProp.iP, iso_range=unique_vals, num=len(unique_vals), rounding=True)
            plot.draw()
            for i,p_val in enumerate(unique_vals):
                plt.text(comp.s_points[indices[i]], comp.T_points[indices[i]]-268.15, f'{p_val:.0f} bar', rotation=60, ha='right', va='center', color='dimgrey', fontsize=8, fontstyle='italic')

            # Functioning points for compressor operations
            for i in range(points-1):
                if i+1 < len(comp.s_points) and i+1 < len(comp.T_points):
                    x1, y1 = comp.s_points[i], comp.T_points[i] - 273.15
                    x2, y2 = comp.s_points[i+1], comp.T_points[i+1] - 273.15
                # x1, y1 = comp.s_points[i], comp.T_points[i] - 273.15
                # x2, y2 = comp.s_points[i+1], comp.T_points[i+1] - 273.15
                plt.plot(x1, y1, marker="o", color="tab:red",markeredgecolor="black")
                plt.annotate(i + 1, (x1 + 0.5, y1 + 10))
                # arrow parameters
                dx = x2 - x1
                dy = y2 - y1
                # arrow_length = np.sqrt(dx**2 + dy**2)
                plt.arrow(x1, y1, dx, dy, head_width=1, head_length=1, fc='tab:blue', ec='tab:blue')

            plt.plot(comp.s_points[points - 1], comp.T_points[points-1] - 273.15, marker="o", color="tab:red",markeredgecolor="black")
            plt.annotate(points, (comp.s_points[points-1] + 0.5, comp.T_points[points-1] - 273.15 + 10))
            plt.grid(alpha=0.2,zorder=-1)
            plt.title('Ts Diagram - multistage hydrogen compression')
            plt.show()
                                   
    def parametric_stage(self):
                
            P_max = np.linspace(70,875,24)
            n_stages = np.linspace(1,4,4)
            P_in = 35
            beta = P_max/P_in
            stadi = {'monostadio': [],
                     '2 stadi': [],
                     '3 stadi': [],
                     '4 stadi': []
                     }
            
            W_c = {'monostadio': [],
                   '2 stadi': [],
                   '3 stadi': [],
                   '4 stadi': []
                   }
    
            delta_h = {'monostadio': [],
                   '2 stadi': [],
                   '3 stadi': [],
                   '4 stadi': []
                   }
            
            LHV = {'monostadio': [],
                   '2 stadi': [],
                   '3 stadi': [],
                   '4 stadi': []
                   }
            
            inp_mono = {'P_max' : P_max[0],
                        'P_in': P_in,
                        'compressor model': 'normal compressor',
                        'flow_rate': 0.05,
                        'fluid': 'Hydrogen',
                        'T_in': 298.15                        
                        }
            
           
            inp_multi = {'P_max' : P_max[0],
                        'P_in': P_in,
                        'compressor model': 'multistage compressor with refrigeration',
                        'flow_rate': 0.05,
                        'fluid': 'Hydrogen',
                        'T_in': 298.15,
                        'perdita di pressione IC': 0.02,
                        'T_IC': 308.15,
                        'n_stages': int(n_stages[1])
                        }
            
            for k in sorted(stadi):
                for n in range(len(P_max)):
                    if k == 'monostadio':
                        inp_mono['P_max'] = P_max[n]
                        a = Compressor(inp_mono, 1)
                        stadi[k].append(a)
                        W_c[k].append(a.power)
                                           
                    else:
                        inp_multi['P_max'] = P_max[n]
                        c = Compressor(inp_multi, 1)
                        stadi[k].append(c)
                        W_c[k].append(c.power)
                        
                inp_multi['n_stages'] += 1
                
            for k in sorted(delta_h):
                for n in range(len(P_max)):
                    delta_h[k].append(W_c[k][n]/0.05)
                    LHV[k].append(delta_h[k][n]/self.hydrogen_LHV*100)
    
            plt.figure(dpi=1000)    
            color = ['orange', 'blue', 'pink', 'turquoise']
            color_2 = ['red', 'blue', 'green', 'turquoise']
    
            a = 0
                                     
            for k in sorted(stadi):
                plt.scatter(beta, W_c[k], s=20, color = color[a] ,edgecolors='k',zorder=3, label = str(k))
                # plt.plot(beta, W_c[k], ls = (0,(3,1,1,1)), color = color_2[a], label = str(k))
                a +=1  
    
            plt.legend()        
            plt.xlabel('$\\beta$ [-]')
            plt.ylabel('W$_{c}$ [kW]')
            plt.title('W$_{c}$ vs $\\beta$')
            plt.grid()
            plt.show()
                    
            b = 0
            plt.figure(dpi=1000)    
    
            for k in sorted(stadi):
                plt.scatter(beta, delta_h[k], s=20, color = color[b] ,edgecolors='k',zorder=3, label = str(k))
                # plt.plot(beta, delta_h[k], ls = (0,(3,1,1,1)), color = color_2[b], label = str(k))
                b +=1  
                
            plt.legend()        
            plt.xlabel('$\\beta$ [-]')
            plt.ylabel('L$_{spec}$ [kJ/kg]')
            plt.title('L$_{spec}$ vs $\\beta$')
            plt.grid()
            plt.show()
            
            c = 0
            plt.figure(dpi=1000)   
            for k in sorted(stadi):
                plt.scatter(beta, LHV[k], s=20, color = color[c],edgecolors='k',zorder=3, label = str(k))
                # plt.plot(beta, LHV[k], ls = (0,(3,1,1,1)), color = color_2[c], label = str(k))
                c +=1  
                
            plt.legend()        
            plt.xlabel('$\\beta$ [-]')
            plt.ylabel('%LHV [kJ/kg]')
            plt.title('%LHV vs $\\beta$')
            plt.grid()
    
            plt.show()
                
    def fluid_vs_air(self, fluid = 'Hydrogen'):
        ################### normal compressor #########################
        flow_rate_H2 = self.maxflowrate
        flow_rate = flow_rate_H2*self.M_air/self.M_H2
        
        beta = np.linspace(2,25,24)
        P_0 = 1.01325
        T_0 = 298.15
        h_0_air = 0
        s_0_air = 0
        P_in = 35 #bar
        T_in = 298.15 #K
        
        prop_air = {'P': [],
                    'T': [],
                    'h': [],
                    's': [],
                    'rho': [],
                    'eta_is': [],
                    'Delta H': [],
                    'Npower': np.zeros(len(beta))}
        
        fluid = 'Hydrogen'
        
        prop_fluid = {'P': [],
                    'T': [],
                    'h': [],
                    's': [],
                    'rho': [],
                    'eta_is': [],
                    'Delta H': [],
                    'Npower': np.zeros(len(beta))}
        
        s_in_air = self.CP_air*m.log(T_in/T_0)+self.R_air*m.log(P_in/P_0)
        s_in_fluid = PropsSI('S', 'T', T_in, 'P', P_in*100000, fluid)/1000
        h_in_air = self.CP_air*(T_in-T_0)
        h_in_fluid = PropsSI('H', 'T', T_in, 'P', P_in*100000, fluid)/1000
        
        eta_pol = 0.7 
        
        eta_is_air = np.zeros(len(beta))
        eta_is_fluid = np.zeros(len(beta))
        
        self.power_fluid_corr = np.zeros(len(beta))
        self.fract_power = np.zeros(len(beta))
        
        for i in range(len(beta)):
            prop_air['P'].append(P_in*beta[i])
            prop_air['eta_is'].append((beta[i]**self.epsilon_air - 1)/(beta[i]**(self.epsilon_air/eta_pol)-1))
            T_is_a = T_in*beta[i]**self.epsilon_air
            h_is_a = self.CP_air*(T_is_a-T_0)
            prop_air['h'].append((h_is_a - h_in_air) / prop_air['eta_is'][i] + h_in_air)
            prop_air['T'].append(T_in+(prop_air['h'][i]-h_0_air)/self.CP_air)
            prop_air['s'].append(self.CP_air*m.log(prop_air['T'][i]/T_0)+self.R_air*m.log(prop_air['P'][i]/P_0))
            prop_air['Delta H'].append(prop_air['h'][i] - h_in_air)
            prop_air['Npower'][i]=((prop_air['h'][i] - h_in_air)*flow_rate)
            
            
            prop_fluid['P'].append(P_in*beta[i])
            prop_fluid['eta_is'].append((beta[i]**self.epsilon - 1)/(beta[i]**(self.epsilon/eta_pol)-1))
            h_is_f = PropsSI('H', 'P', prop_fluid['P'][i]*100000, 'S', s_in_fluid*1000, fluid)/1000
            prop_fluid['h'].append((h_is_f - h_in_fluid) / prop_fluid['eta_is'][i] + h_in_fluid)
            prop_fluid['Delta H'].append(prop_fluid['h'][i] - h_in_fluid)
            prop_fluid['Power'][i] = ((prop_fluid['h'][i] - h_in_fluid)*flow_rate)
            
            self.fract_power[i] =  prop_fluid['Power'][i]/ prop_air['Power'][i]
            
            self.power_fluid_corr[i] = flow_rate_H2*(prop_fluid['h'][i] - h_in_fluid)
        
################################################################################################################################################

        # air flow rate and equivalent compression ratio
        
        # here a hydrogen compressor with defined beta and flow rate is considered, 
        # using air allows to use the equivalent air flow rate and the correct beta in the compressor
       
        prop_air_2 = {'P': [],
                      'T': [],
                      'h': [],
                      's': [],
                      'rho': [],
                      'eta_is': [],
                      'Delta H': [],
                      'Power': np.zeros(len(beta))}

        beta_air = np.linspace(2,25,24)    
        cost = np.zeros(len(beta_air))

        beta_fluid = np.zeros(len(beta_air))
        self.delta_beta = np.zeros(len(beta_air))
        beta_corr = np.zeros(len(beta_air))
        n_stages = np.zeros(len(beta_air))
        
        costi_corr  = {'Z_1': np.zeros(len(beta_air)),
                       'Z_2': np.zeros(len(beta_air))}

        costi       = {'Z_1': np.zeros(len(beta_air)),
                       'Z_2': np.zeros(len(beta_air))}    

        for n in range(len(beta_air)):
            cost[n] = ((beta_air[n]**(self.epsilon_air)-1)*self.CP_air*T_in/0.85)
            beta_fluid[n] = (1+0.85/(self.CP_H2*T_in)*(cost[n]))**(1/self.epsilon)
            self.delta_beta[n] = (beta_air[n])/(beta_fluid[n])
            beta_corr[n] = self.delta_beta[n]*beta_air[n]
            n_stages[n] = m.log(beta_air[n])/m.log(beta_fluid[n])
        
        
        # Thermodynamic points considering air with beta correction
        for i in range(len(beta)):
            prop_air_2['P'].append(P_in*beta_corr[i])
            prop_air_2['eta_is'].append((beta_corr[i]**self.epsilon_air - 1)/(beta_corr[i]**(self.epsilon_air/eta_pol)-1))
            T_is_a = T_in*beta_corr[i]**self.epsilon_air
            h_is_a = self.CP_air*(T_is_a-T_0)
            prop_air_2['h'].append((h_is_a - h_in_air) / prop_air_2['eta_is'][i] + h_in_air)
            prop_air_2['T'].append(T_in+(prop_air_2['h'][i]-h_0_air)/self.CP_air)
            prop_air_2['s'].append(self.CP_air*m.log(prop_air_2['T'][i]/T_0)+self.R_air*m.log(prop_air_2['P'][i]/P_0))
            prop_air_2['Delta H'].append(prop_air_2['h'][i] - h_in_air)
            prop_air_2['Power'][i]=((prop_air_2['h'][i] - h_in_air)*flow_rate)
            

        for n in range(len(beta_corr)):
            # m air con beta corr #
            costi_corr['Z_1'][n] = (39.5*flow_rate*beta_corr[n])/(0.9-0.85)*(m.log(beta_corr[n]))
            costi_corr['Z_2'][n] = 44.71*flow_rate/(0.95-0.85)*beta_corr[n]*m.log(beta_corr[n])
            # m air con beta air (2,3,4....) #
            costi['Z_1'][n] = (39.5*flow_rate*beta_air[n])/(0.9-0.85)*(m.log(beta_air[n]))
            costi['Z_2'][n] = 44.71*flow_rate/(0.95-0.85)*beta_air[n]*m.log(beta_air[n])
        
###########################################################################################################################################################
        
        # costs expressed as a function of mass flow rate and integration of compression ratio (beta)
        
        costi_air = {'Z_1': np.zeros(len(beta)), 'Z_2': np.zeros(len(beta)), 'Z_3': np.zeros(len(beta)), 'Z_4': np.zeros(len(beta)),
                      'Z_5': np.zeros(len(beta)),'Z_6': np.zeros(len(beta)), 'Z_7': np.zeros(len(beta)), 'Z_8': np.zeros(len(beta)),
                      'Z_9': np.zeros(len(beta))
                      }
        
        costi_fluid = {'Z_1': np.zeros(len(beta)),'Z_2': np.zeros(len(beta)),'Z_3': np.zeros(len(beta)), 'Z_4': np.zeros(len(beta)),
                      'Z_5': np.zeros(len(beta)), 'Z_6': np.zeros(len(beta)),'Z_7': np.zeros(len(beta)), 'Z_8': np.zeros(len(beta)),
                      'Z_9': np.zeros(len(beta))
                      }
        
        #correzione della portata#
        costi_fluid_corr = {'Z_1': np.zeros(len(beta)), 'Z_2': np.zeros(len(beta)), 'Z_3': np.zeros(len(beta)), 'Z_4': np.zeros(len(beta)),
                            'Z_5': np.zeros(len(beta)), 'Z_6': np.zeros(len(beta)), 'Z_7': np.zeros(len(beta)), 'Z_8': np.zeros(len(beta)),
                            'Z_9': np.zeros(len(beta))
                            }
        
        a = 0.04147
        b = 454.8
        c = 181000
        
        for i in range(len(beta)):
            
            costi_air['Z_1'][i] = 91562*(prop_air['Power'][i]/455)**0.67             # Ref: "Development of Cost Correlations for the Economic Assessment of Power Plant Equipment"
            costi_air['Z_2'][i] = 9642.2*prop_air['Power'][i]**0.46                  # //
            costi_air['Z_3'][i] = 10167.5*prop_air['Power'][i]**0.46                 # //
            costi_air['Z_4'][i] = m.log10(prop_air['Power'][i])  + a*prop_air['Power'][i]**2 + b+prop_air['Power'][i] + c # //
            costi_air['Z_5'][i] = 25421.72*prop_air['Power'][i]**0.61                #AN ECONOMIC ANALYSIS OF THREE HYDROGEN LIQUEFACTION SYSTEMS#
            costi_air['Z_6'][i] = 40038*prop_air['Power'][i]**0.6038                 #Techno-economic modelling and analysis of hydrogen fuelling stations#
            costi_air['Z_7'][i] = 3083.3*prop_air['Power'][i]**0.8355*2              #The TechnoEconomics of Hydrogen Compression TECHNICAL BRIEF, pipeline compressor
            costi_air['Z_8'][i] = 63684.6*prop_air['Power'][i]**0.4603*1.3           #The TechnoEconomics of Hydrogen Compression TECHNICAL BRIEF, HRS 350 bar
            costi_air['Z_9'][i] = 62209.9*prop_air['Power'][i]**0.6038*1.3           #The TechnoEconomics of Hydrogen Compression TECHNICAL BRIEF, HRS 700 bar
            
            costi_fluid['Z_1'][i] = 91562*(prop_fluid['Power'][i]/455)**0.67      
            costi_fluid['Z_2'][i] = 9642.2*prop_fluid['Power'][i]**0.46                
            costi_fluid['Z_3'][i] = 10167.5*prop_fluid['Power'][i]**0.46  
            costi_fluid['Z_4'][i] = m.log10(prop_fluid['Power'][i])  + a*prop_fluid['Power'][i]**2 + b+prop_fluid['Power'][i] + c
            costi_fluid['Z_5'][i] = 25421.72*prop_fluid['Power'][i]**0.61               
            costi_fluid['Z_6'][i] = 40038*prop_fluid['Power'][i]**0.6038                 
            costi_fluid['Z_7'][i] = 3083.3*prop_fluid['Power'][i]**0.8355*2            
            costi_fluid['Z_8'][i] = 63684.6*prop_fluid['Power'][i]**0.4603*1.3           
            costi_fluid['Z_9'][i] = 62209.9*prop_fluid['Power'][i]**0.6038*1.3   
            
            costi_fluid_corr['Z_1'][i] = 91562*(self.power_fluid_corr[i]/455)**0.67   
            costi_fluid_corr['Z_2'][i] = 9642.2*self.power_fluid_corr[i]**0.46                     
            costi_fluid_corr['Z_3'][i] = 10167.5*self.power_fluid_corr[i]**0.46  
            costi_fluid_corr['Z_4'][i] = m.log10(self.power_fluid_corr[i])  + a*self.power_fluid_corr[i]**2 + b+self.power_fluid_corr[i] + c
            costi_fluid_corr['Z_5'][i] = 25421.72*self.power_fluid_corr[i]**0.61               
            costi_fluid_corr['Z_6'][i] = 40038*self.power_fluid_corr[i]**0.6038                 
            costi_fluid_corr['Z_7'][i] = 3083.3*self.power_fluid_corr[i]**0.8355*2            
            costi_fluid_corr['Z_8'][i] = 63684.6*self.power_fluid_corr[i]**0.4603*1.3           
            costi_fluid_corr['Z_9'][i] = 62209.9*self.power_fluid_corr[i]**0.6038*1.3   
        
        plt.figure(dpi=1000)    
        color = ['orange', 'lightgreen', 'blue','brown', 'red', 
                  'pink', 'salmon', 'indigo', 'turquoise', 'gold', 'darkkhaki'
                  ]
        label = ['91.562(kW/455)$^{0,67}$, Sadeghi et al.',
                  '9.642,2(kW)$^{0,46}$, Parikhani et al.',
                  '10.167,5(kW)$^{0,46}$, Rezayanet al.',
                  'log(kW) + a(kW)$^{2}$ + b(kW) + c, Shamoushaki et al.',
                  '25.422(kW)$^{0,61}$, Syed et al.', 
                  '40.038(kW)$^{0,61}$, Blazquez-Diaz ',
                  '6.186,6(kW)$^{0,8355}$, HDSAM',
                  '82.790(kW)$^{0.4603}$, HDSAM',
                  '80.873(kW)$^{0,6038}$, HDSAM', 
                  'Saghafifar - Gadalla',
                  'Roosen et al'
                  ]
        a = 0
        b = 0
        c = 0
##############################################################################################################################################################
        # Costs - Using the same beta#
        
        for n in sorted(costi_air):
            if a <= 2:
                plt.plot(beta_air, costi_air[n], marker = 'o', color = color[a], zorder=3 , label = label[a])
                a +=1 
            elif a>=4 and a < 9:
                plt.plot(beta_air, costi_air[n], marker = '^', color = color[a], zorder=3, label = label[a])
                a +=1 
            else:
                a +=1
        
        plt.plot(beta_air, costi['Z_1'], marker = 's', color = color[9], zorder=3, label = label[9])
        plt.plot(beta_air, costi['Z_2'], marker = 's', color = color[10], zorder=3, label = label[10])
        plt.legend(fontsize=8, bbox_to_anchor=(1.04, 1), borderaxespad=0)
        plt.xlabel('$\\beta$ [-]')
        plt.ylabel('Costi [€]')
        plt.title('$\\beta$ vs Costi, m'+'\u0307'+'$_{air}$= '+str(round(flow_rate,4)))
        plt.grid()
        plt.show()
        
        plt.figure(dpi=1000)    
        for n in sorted(costi_fluid):
            if b <= 2:
                plt.plot(beta_air, costi_fluid[n], marker = 'o', color = color[b], zorder=3, label = label[b])
                b +=1 
            elif b>=4:
                plt.plot(beta_air, costi_fluid[n], marker = '^', color = color[b], zorder=3, label = label[b])
                b +=1
            else:
                b +=1
                
        plt.plot(beta_air, costi_corr['Z_1']*self.M_air/self.M_H2, marker = 's', color = color[9], zorder=3, label = label[9])
        plt.plot(beta_air, costi_corr['Z_2']*self.M_air/self.M_H2, marker = 's', color = color[10], zorder=3, label = label[10])
        plt.legend(fontsize=8, bbox_to_anchor=(1.04, 1), borderaxespad=0)
        plt.xlabel('$\\beta$ [-]')
        plt.ylabel('Costi [€]')
        plt.title('$\\beta$ vs Costi, m'+'\u0307'+'$_{H2}$= '+str(flow_rate))
        plt.grid()
        plt.show()
        
        plt.figure(dpi=1000)    
        for n in sorted(costi_air):
            if c <= 2:
                plt.plot(beta_air, costi_fluid_corr[n], marker = 'o', color = color[c], zorder=3, label = label[c])
                c +=1 
            elif c>=4 and c < 9:
                plt.plot(beta_air, costi_fluid_corr[n], marker = '^', color = color[c], zorder=3, label = label[c])
                c +=1 
            else:
                c +=1
        plt.plot(beta_air, costi_corr['Z_1'], marker = 's', color = color[9], zorder=3, label = label[9])
        plt.plot(beta_air, costi_corr['Z_2'], marker = 's', color = color[10], zorder=3, label = label[10])
        plt.legend(fontsize=8, bbox_to_anchor=(1.04, 1), borderaxespad=0)
        plt.xlabel('$\\beta$ [-]')
        plt.ylabel('Costi [€]')
        plt.title('$\\beta$ vs Costi, m'+'\u0307'+'$_{H2}$= '+str(round(flow_rate_H2,4)))
        plt.grid()
        plt.show()
                    
############################################################################################################################################################################################################    


    def use(self,step,available_hyd_lp=False,storable_hydrogen_hp=False,massflowrate=False):
        """
        Compressor object absorbs electricity and works on fluid
        
        storable_hydrogen_hp    : float storable hydrogen H tank max_capacity - SOC[h-1] [kg] or maximum absorbable production if there is no tank
        available_hyd_lp        : float available hydrogen H tank SOC[h-1] [kg]
        massflowrate            : float hydrogen flow rate to be processed in the timestep [kg/s]
        step                    : int step to be simulated [-]

        output : 
        float hydrogen compressed in the timestep [kg/s]    
        float power absorbed that hour [kW]
        float cooling needs [kW]
        """
        if self.model == 'simple_compressor':
            
            self.hyd[step]  = massflowrate                  # [kg/s] hydrogen mass flowrate in the considered timestep
            p_absorbed      = self.en_cons*massflowrate     # [kW] power consumption inimestep
            t_absorbed      = 0                             # [kW] cooling not considered in this model.  Power consumption is estimated from adiabatic compression model with 
                                                            #                                             an efficiency of 50%. This efficiency considers the efficiency of 
                                                            #                                             electrical power transformation and auxiliary systems such as the 
                                                            #                                             cooling circuit.
            
            return(self.hyd[step],-p_absorbed,t_absorbed)
        
        if self.model == "normal_compressor":
            
            self.hyd[step]  = self.flow_rate                # [kg/s] hydrogen mass flowrate in the considered timestep
            p_absorbed      = self.Npower                   # [kW] energia del compressore in un'ora di utilizzo
            t_absorbed      = 0                             # [kW] cooling not considered in this model
            
            if self.hyd[step] > available_hyd_lp:           # if the is not enough space in the H tank to store the hydrogen (H tank is nearly full)
                self.hyd[step]  = 0                         # partial load operations not allowed
                p_absorbed      = 0
                t_absorbed      = 0
            
            return(self.hyd[step],-p_absorbed,-t_absorbed)
                
        elif self.model == 'multistage_compressor_with_refrigeration' or self.model == 'compressor_with_refrigeration': 
            
            if massflowrate == False:                   # compressor working in between a buffer and an HPH storage. Only full load operations allowed. 
                self.hyd[step]  = self.Nflowrate        # [kg/s] hydrogen mass flow rate
                t_absorbed      = self.IC_power_list    # [kW] of required cooling
                p_absorbed      = self.Npower           # [kW] absorbed power 
                    
                if self.hyd[step] > available_hyd_lp:   # if the is not enough space in the H tank to store the hydrogen (H tank is nearly full)
                    self.hyd[step]  = 0
                    p_absorbed      = 0
                    t_absorbed      = 0
                    
                elif self.hyd[step] > storable_hydrogen_hp:     # if mass flow rate is higher than the storable amount of hydrogen
                    self.hyd[step]  = 0
                    p_absorbed      = 0
                    t_absorbed      = 0
                            
                return(self.hyd[step],-p_absorbed,-t_absorbed)
            
            else:   # working with a single storage pressure level. Streamlining of partial load functioning - Can be upgraded
                self.hyd[step]  = massflowrate
                p_absorbed      = (massflowrate*sum(self.comp_lav_spec))/self.eta_motor
                t_absorbed      = massflowrate*self.delta_H
                
                return(self.hyd[step],-p_absorbed,-t_absorbed)
                
    @property
    def comp_power(self):

        return sum(self.comp_power_list)
    
    def comp_lav(self):
        
        return sum(self.comp_lav_spec)
    
    def thermal_power(self):

        return sum(self.IC_power_list)
    
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
            'refund': dict
                'rate': float, percentage of initial investment which will be rimbursed [%]
                'years': int, years for reimbursment
            'replacement': dict
                'rate': float, replacement cost as a percentage of the initial investment [%]
                'years': int, after how many years it will be replaced
        """
        compressor_costs = {                                                                # Ref: https://hsweb.hs.uni-hamburg.de/projects/star-formation/hydrogen/P2H_Full_Study_FCHJU.pdf
                            'Gas Flow [kg/h]'               : [1,10,100,1000,10000],        # dataset to be updated 
                            'Compressed Tanks [€/(kg/h)]'   : [63460,29007,13259,6016,2770]
                            }
        
        costdf = pd.DataFrame(compressor_costs)
        
        x = costdf['Gas Flow [kg/h]']
        y = costdf['Compressed Tanks [€/(kg/h)]']
        
        f = interp1d(x, y,fill_value= 'extrapolate')
        
        hourly_massflow = self.Nflowrate*3600       # [kg/h] maxflowrate defined in kg/s --> *3600 to obtain kg/h
        specific_cost = f(hourly_massflow)          # [€/(kg/h)] 
        
        tech_cost = {key: value for key, value in tech_cost.items()}

        if self.model == 'simple_compressor':
            exchange_rate   = 0.91                                           # [2015USD/2023€]  exchange rate between USD and €
            correlation1    = (51901*(hourly_massflow**0.65))*exchange_rate # Ref: https://www.sciencedirect.com/science/article/pii/S0360319919330022?via%3Dihub
            correlation     = correlation1
        elif self.model != 'simple_compressor':
            size = self.Npower
            IF              = 1.3                                            # [-] Installation Factor
            correlation2    = 63684.6*(self.Npower**0.4603)*IF               # Ref: https://transitionaccelerator.ca/wp-content/uploads/2021/10/TA-Briefs-1.2-The-Techno-Economics-of-Hydrogen-Compression-FINALPDF.pdf
            correlation     = correlation2
        
        if tech_cost['cost per unit'] == 'default price correlation':
            C = correlation                         # [€] 
        else:
            C = size * tech_cost['cost per unit']   # [€]

        tech_cost['total cost'] = tech_cost.pop('cost per unit')
        tech_cost['total cost'] = C
        tech_cost['OeM'] = tech_cost['OeM'] *C /100 # €

        self.cost = tech_cost

###########################################################################################################################################################################

if __name__ == "__main__":
    
    from matplotlib.patches import Patch
    
    'Simple Compressor test'

    inp_test = {'P_out'             : 525,
                'P_in'              : 30,
                'compressor model'  : 'simple_compressor',
                'Nflow_rate'        : 0.03,
                'fluid'             : 'Hydrogen',
                'T_in'              : 343.15,
                'pressure losses IC': 0.00,
                'T_IC'              : 308.15,
                'n_stages'          : 3,
                'only_renewables'   : False}

    sim_steps   = 50      # [-] number of steps to be considered for the simulation - usually a time horizon of 1 year minimum is considered
    timestep    = 60      # [min] selected timestep for the simulation
        
    P_in    = np.array([0,15,30,60])
    P_out   = np.arange(100,701,50)
    en_cons = {}   # [kWh/kg]
    
    for i in P_in:
        inp_test['P_in'] = i
        pout = []
        en_cons[i] = pout
        for k in P_out:
            inp_test['P_out'] = k
            comp  = Compressor(inp_test,sim_steps,timestep=timestep)  # creating compressor object
            pout.append(round(comp.en_cons/3600,2))
    
    fig = plt.figure(dpi=600)
    for key, values in en_cons.items():
        plt.plot(P_out-key, values, label=key)

    plt.xlabel('P$_{diff}$ (P$_{out}$-P$_{in}$) [bar]')
    plt.ylabel('Energy consumption [kWh/kg$_{H2}$]')
    plt.legend(title='P$_{inlet}$ [barg]')
    plt.grid(alpha=0.3)
    plt.title('Simple compressor model - Specific consumption')
    plt.show()

    'Multistage Compressor With refrigeration test'
    
    inp_test['compressor model'] = 'multistage_compressor_with_refrigeration'
    P_in    = np.array([1,15,30,60])
    power_cons  = {} # # [kWh/kg] compressor power consumption
    
    for i in P_in:
        inp_test['P_in'] = i
        pout = []
        power_cons[i] = pout
        for k in P_out:
            inp_test['P_out'] = k
            comp  = Compressor(inp_test,sim_steps,timestep=timestep)  # creating compressor object
            power = sum(comp.comp_lav_spec)/comp.eta_motor/3600
            pout.append(round(power,2))
    
    fig = plt.figure(dpi=600)
    for key, values in power_cons.items():
        plt.plot(P_out-key, values, label=key)

    plt.xlabel('P$_{diff}$ (P$_{out}$-P$_{in}$) [bar]')
    plt.ylabel('Energy consumption [kWh/kg$_{H2}$]')
    plt.legend(title='P$_{inlet}$ [bar]')
    plt.grid(alpha=0.3)
    plt.title('Multistage compressor model - Specific consumption')
    plt.show()
    
    'Power consumption - FLow rate'
    
    flow        = np.linspace(0,inp_test['Nflow_rate'],sim_steps)    # [kg/s] hydrogen mass flow rate 
    power_cons  = [] # [kW] compressor power consumption
    cooling     = [] # [kW] thermal power need
    
    for i,val in enumerate(flow):
        a,b = comp.use(i,massflowrate=val)[1:]
        power_cons.append(-a)
        cooling.append(-b)
     
    fig, ax = plt.subplots(dpi=600)
    ax.plot(flow*3600,power_cons,color='tab:green',zorder=3)
    ax.set_xlabel('Hydrogen mass flow rate [kg/h]')
    ax.set_ylabel('Power consumption [kW]')
    ax.grid(alpha=0.3, zorder=-1)
    ax.set_title(f"Consumption vs flow rate (Nominal Power = {comp.Npower} kW)")
    
    'Specific work for different fluids'
    
    fluid = ['Hydrogen','Helium','Methane']    
    P_out   = np.arange(2,801,10)
    power_cons  = {}  # [kWh/kg] compressor power consumption
    inp_test['P_in'] = 1
    colors = ['#1f77b4', '#7fcdbb',  '#fc8d59']
    
    for name in fluid:
        inp_test['fluid'] = name
        pout = []
        power_cons[name] = pout
        for k in P_out:
            inp_test['P_out'] = k
            comp  = Compressor(inp_test,sim_steps,timestep=timestep)  # creating compressor object
            power = sum(comp.comp_lav_spec)/comp.eta_motor/3600
            pout.append(round(power,2))
    
    fig = plt.figure(dpi=600)
    i = 0
    for key, values in power_cons.items():
        plt.plot(P_out, values, label=key, color=colors[i])
        i += 1

    plt.xlabel('Final pressure [bar]')
    plt.ylabel('Energy consumption [kWh/kg$_{H2}$]')
    plt.legend(title='Fluid')
    plt.grid(alpha=0.3)
    plt.title('Multistage compressor model - Specific consumption')
    plt.show()
    
    'Energy consumption vs HHV - Hydrogen'
    
    inp_test['fluid'] = 'Hydrogen'
    pout    = []            # [kJ/kg] power cons
    hhv     = c.HHVH2*1000  # [kJ/kg] heating value
    for k in P_out:
        inp_test['P_out'] = k
        comp  = Compressor(inp_test,sim_steps,timestep=timestep)    # creating compressor object
        power = round(sum(comp.comp_lav_spec)/comp.eta_motor,2)     # [kJ/kg] 
        pout.append(round(power/hhv*100,2))
    
    fig = plt.figure(dpi=600)
    plt.plot(P_out, pout, color=colors[0])
    plt.xlabel('Final pressure [bar]')
    plt.ylabel('Compression energy as HHV % [-]')
    plt.grid(alpha=0.3)
    plt.title('Energy costs vs HHV - Hydrogen')
    plt.show()
    
    
    'Test - Thermodynamic transformation'
    
    inp_test = {'P_out'             : 525,
                'P_in'              : 30,
                'compressor model'  : 'multistage_compressor_with_refrigeration',
                'Nflow_rate'        : 0.03,
                'fluid'             : 'Hydrogen',
                'T_in'              : 343.15,
                'pressure losses IC': 0.00,
                'T_IC'              : 308.15,
                'n_stages'          : 3,
                'only_renewables'   : False}

    comp  = Compressor(inp_test,sim_steps,timestep=timestep)    # creating compressor object
    comp.thermodynamics_points()

  
    


