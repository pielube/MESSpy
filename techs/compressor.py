# import CoolProp
from CoolProp.CoolProp import PropsSI
from CoolProp.Plots import PropertyPlot
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from scipy.interpolate import interp1d
import math as m
import pandas as pd
import numpy
import matplotlib.pyplot as plt
import os
import sys 
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),os.path.pardir)))   # temorarily adding constants module path 
import constants as c

class Compressor:
    
    def __init__(self, parameters, simulation_hours, maxflowrate=False):      
        
        """
        Create a compressor object
    
        parameters : dictionary
            'compressor model' : 'normal compression', 'compression with refrigeration', 'multistade compression with refrigeration'
            'P_max': max pressure [bar]          
            'P_in': inlet pressure [bar]
            'T_in' : inlet Temperature
            'P_start' : optional, needed for initial mass calculation
            'Power' : optional, power [kW]
            'flow_rate' : hydrogen's flow rate [kg/h]
            'pressure losses' : pressure losses in heat exchangers [%]
            'n_stages' : number of compression stages
            'T_IC': Temperature of intercooler [K]
            
        """
        
        self.hyd = numpy.zeros(simulation_hours)
        self.timestep = 1
        if 'Power' in parameters:  self.Npower = parameters['Power']  # [kW] nominal power of compressor
        if maxflowrate: self.maxflowrate = maxflowrate            # [kg/h] maximum flow rate that can be processed by the compressor - size is automatically defined
        
        if parameters['compressor model'] == 'simple compressor':
            self.model = parameters['compressor model']

            '''
            POWER CONSUMPTION DATA  
            Ref. FCH Study - https://hsweb.hs.uni-hamburg.de/projects/star-formation/hydrogen/P2H_Full_Study_FCHJU.pdf
            
            Power consumption is estimated from adiabatic compression model with 
            an efficiency of 50%. This efficiency considers the efficiency of 
            electrical power transformation and auxiliary systems such as the 
            cooling circuit.
            
            '''
            # Original set of data from which linear approximation is derived
            powerconsumption = {
                                'Pressure in [barg]': [0, 0, 15, 15, 30, 30, 60, 60],
                                'Pressure out [barg]': [200, 500, 200, 500, 200, 500, 200, 500],
                                'Pressure difference [bar]': [200, 500, 185, 485, 170, 470, 140, 440],
                                'Energy consumption [kWh/kg]': [5, 6.3, 2.4, 3.5, 1.7, 2.7, 1.1, 2]                
                                }
            # Slope and intercetp values calculated on 'poerconsumption' dataset
            coefficients    = {
                                'Pressure in [barg]': [0,15,30,60],
                                'Slope': [0.00433,0.00367,0.0033,0.003],
                                'Intercept': [4.13,1.72,1.13,0.68]
                                }
            # creating a DataFrame to better handle the data
            df = pd.DataFrame(coefficients)
            
            self.P_in   = parameters['P_in']    # [bar] compressor inlet pressure
            self.P_out  = parameters['P_out']   # [bar] compressor outlet pressure
            slope       = df.loc[df['Pressure in [barg]'] == self.P_in, 'Slope'].iloc[0]
            intercept   = df.loc[df['Pressure in [barg]'] == self.P_in, 'Intercept'].iloc[0]
            
            self.en_cons     = intercept + slope*(self.P_out-self.P_in) # [kWh/kgH2] specific energy consumption of the simple compressor model
        
        else:  # if any other model il selected
            # Parameters #
            self.P_max  = parameters['P_out']
            self.P_in   = parameters['P_in']
            self.T_in   = parameters['T_in']
            self.fluid  = parameters['fluid']
            ####################################################
            
            # initialisation of variables
            self.P_points   = []
            self.T_points   = []
            self.h_points   = []
            self.s_points   = []
            self.rho_points = []
            
            ####################################################
            
            # thermodynamic and chemical properties
            self.hydrogen_HHV = c.HHVH2*1000    # [kJ/kg]
            self.hydrogen_LHV = c.LHVH2*1000    # [kJ/kg]
            self.CP_air = c.CP_AIR              # [kJ/kgK]
            self.CV_air = c.CV_AIR              # [kJ/kgK]
            self.CP_H2  = c.CP_H2               # [kJ/kg*K] at T=25 °C e P=1 atm
            self.CV_H2  = c.CV_H2               # [kJ/kg*K] at T=25 °C e P=1 atm
            self.R_univ = c.R_UNIVERSAL         # [J/kmol*K]
            self.M_H2   = c.H2MOLMASS*1000      # [kg/kmol]
            self.M_air  = c.AIRMOLMASS          # [kg/kmol]
            self.R_H2 = self.R_univ/self.M_H2   # [J/kgK]
            self.R_air = self.R_univ/self.M_air # [J/kgK]
            
            #####################################################
            
            # calculation of polytropic H2 and air exponents
            self.beta_max = self.P_max/self.P_in
            self.eta_pol = 0.70                     # Ref: The TechnoEconomics of Hydrogen Compression 
            self.gamma = self.CP_H2/self.CV_H2
            self.gamma_air = self.CP_air/self.CV_air
            self.epsilon = (self.gamma-1)/self.gamma
            self.epsilon_air = (self.gamma_air-1)/self.gamma_air
            self.exp_pol_H2 = (self.eta_pol/self.epsilon)/(self.eta_pol/self.epsilon-1)
            self.exp_pol_air = (self.eta_pol/self.epsilon_air)/(self.eta_pol/self.epsilon_air-1)
            self.omega = (self.exp_pol_H2-1)/self.exp_pol_H2
            self.omega_air = (self.exp_pol_air-1)/self.exp_pol_air
            
            ####################################################################################
            
            # Single-stage compressor operation without refrigeration
                
            if parameters['compressor model'] == 'normal compressor':
                self.model = parameters['compressor model']
                self.n_stages = 1
                self.delta_P = 0      # [-] pressure losses
                
                self.P_target = numpy.zeros(self.n_stages)
                
                for i in range(self.n_stages*2):          # point 0 is the input to the compressor, point 1 is the outlet of the compressor
    
                    self.P_points.append(self.P_in)
                    self.T_points.append(self.T_in)
                    self.h_points.append(PropsSI('H', 'T', self.T_points[i], 'P', self.P_points[i]*100000, self.fluid)/1000)
                    self.s_points.append(PropsSI('S', 'T', self.T_points[i], 'P', self.P_points[i]*100000, self.fluid)/1000)
                    self.rho_points.append(PropsSI('D', 'T', self.T_points[i], 'P', self.P_points[i]*100000, self.fluid))
                
                self.comp_power_list = []
                self.comp_lav_spec = []          
                self.IC_power_list = []
                self.comp_beta_targ = self.beta_max**(1/self.n_stages)
                                               
                self.comp_power_list.append(0.)
                self.comp_lav_spec.append(0.)
                self.IC_power_list.append(0.)
                self.P_target = self.P_in*self.comp_beta_targ**self.n_stages
                P_in = self.P_points[0]
                h_in = self.h_points[0]
                s_in = self.s_points[0]
                T_in = self.T_points[0]
                rho_in = self.rho_points[0]
                
                # compression from 0 to 1
                P_out = self.comp_beta_targ * P_in
                s_out_iso = s_in
                h_out_iso = PropsSI('H', 'P', P_out*100000, 'S', s_out_iso*1000, self.fluid)/1000     # [kJ/kg]
                self.T_is = PropsSI('T', 'P', P_out*100000, 'S', s_out_iso*1000, self.fluid)
                self.eta_is = (self.comp_beta_targ**self.epsilon - 1)/(self.comp_beta_targ**(self.epsilon/self.eta_pol)-1)
                h_out = (h_out_iso - h_in) / self.eta_is + h_in
                
                self.T_points[1] = PropsSI('T', 'P', P_out*100000, 'H', h_out*1000, self.fluid)
                self.P_points[1] = P_out
                self.h_points[1] = h_out
                self.s_points[1] = PropsSI('S', 'P', P_out*100000, 'H', h_out*1000, self.fluid)/1000
                self.rho_points[1] = PropsSI('D', 'P', P_out*100000, 'H', h_out*1000, self.fluid)
                
                self.comp_lav_spec[0] = h_out - h_in   # [kJ/kg] specific work of compression given in and out conditions
                
                if 'flow_rate' in parameters:    # if flow rate is defined as input - compressor nominal power is defined as a consequence
                    self.flow_rate = parameters['flow_rate']                    # [kg/h] hourly mass flow rate
                    self.Npower = self.flow_rate/3600 * self.comp_lav_spec[0]   # [kW] compressor nominal power 
                    self.IC_power_list = 0
                    
                elif 'Power' in parameters:      # if Nominal power is defined as input - compressor mass flow rate is defined as a consequence
                    self.Npower = parameters['Power']                       # [kW] compressor nominal power
                    self.flow_rate = self.Npower/self.comp_lav_spec[0]*3600 # [kg/h] hourly mass flow rate
                    self.IC_power_list = 0
                
                ################################################################################################
    
    
            # Single-stage compressor operation with refrigeration 
            
            if parameters['compressor model'] == 'compressor with refrigeration':
                self.model = parameters['compressor model']
                if 'Power' in parameters: self.Npower = parameters['Power']
                self.n_stages = 1
                self.T_IC = parameters['T_IC']
                self.delta_P = parameters['pressure losses IC']  
                
                self.P_target = numpy.zeros(self.n_stages)
                
                for i in range(self.n_stages*2+1):           # point 0 is the input to the compressor, point 1 is the outlet of the compressor
    
                    self.P_points.append(self.P_in)
                    self.T_points.append(self.T_in)
                    self.h_points.append(PropsSI('H', 'T', self.T_points[i], 'P', self.P_points[i]*100000, self.fluid)/1000)
                    self.s_points.append(PropsSI('S', 'T', self.T_points[i], 'P', self.P_points[i]*100000, self.fluid)/1000)
                    self.rho_points.append(PropsSI('D', 'T', self.T_points[i], 'P', self.P_points[i]*100000, self.fluid))
                
                self.comp_power_list = []
                self.comp_lav_spec = []          
                self.IC_power_list = []
                self.comp_beta_targ = self.beta_max**(1/self.n_stages)       
                                        
                self.comp_lav_spec.append(0.)
                self.comp_power_list.append(0.)
                self.IC_power_list.append(0.)
                self.P_target = self.P_in*self.comp_beta_targ**self.n_stages
                P_in = self.P_points[0]
                h_in = self.h_points[0]
                s_in = self.s_points[0]
                T_in = self.T_points[0]
                rho_in = self.rho_points[0]
                
                # compression 0 - 1
                P_out_ic = self.P_target
                P_out = (P_out_ic)/(1-self.delta_P)
                self.beta_new = P_out/P_in
                s_out_iso = s_in
                h_out_iso = PropsSI('H', 'P', P_out*100000, 'S', s_out_iso*1000, self.fluid)/1000
                self.T_is = PropsSI('T', 'P', P_out*100000, 'S', s_out_iso*1000, self.fluid)
                self.eta_is = (self.comp_beta_targ**self.epsilon - 1)/(self.comp_beta_targ**(self.epsilon/self.eta_pol)-1)
                h_out = (h_out_iso - h_in) / self.eta_is + h_in
                
                self.T_points[1] = PropsSI('T', 'P', P_out*100000, 'H', h_out*1000, self.fluid)
                self.P_points[1] = P_out
                self.h_points[1] = h_out
                self.s_points[1] = PropsSI('S', 'P', P_out*100000, 'H', h_out*1000, self.fluid)/1000
                self.rho_points[1] = PropsSI('D', 'P', P_out*100000, 'H', h_out*1000, self.fluid)
                
                # interrefrigeration 1 - 2
                self.P_points[2] = P_out_ic
                self.T_points[2] = self.T_IC
                self.h_points[2] = PropsSI('H', 'P', P_out_ic*100000, 'T', self.T_IC, self.fluid)/1000
                self.s_points[2] = PropsSI('S', 'P', P_out_ic*100000, 'T', self.T_IC, self.fluid)/1000
                self.rho_points[2] = PropsSI('D', 'P', P_out_ic*100000, 'T', self.T_IC, self.fluid)
                
                h_out_ic = PropsSI('H', 'P', P_out_ic*100000, 'T', self.T_points[2], self.fluid)/1000
                
                self.comp_lav_spec[0] = h_out - h_in    # [kJ/kg] specific work of compression given in and out conditions
                self.delta_H = h_out - h_out_ic
                
                if 'flow_rate' in parameters:       # if flow rate is defined as input - compressor nominal power is defined as a consequence
                    self.maxflowrate = parameters['flow_rate']                            # [kg/h] hourly mass flow rate
                elif self.maxflowrate:  
                    self.Npower = self.maxflowrate/3600 * self.comp_lav_spec[0]           # [kW] compressor nominal power 
                    self.IC_power_list = self.maxflowrate/3600 * self.delta_H             # [kW] Total heat to be removed by the cooling system
                    
                elif 'Power' in parameters:     # if nominal power is defined as input - compressor mass flow rate is defined as a consequence
                    self.Npower = parameters['Power']                                   # [kW] compressor nominal power
                    self.maxflowrate = self.Npower/self.comp_lav_spec[0]*3600             # [kg/h] hourly mass flow rate
                    self.IC_power_list = self.maxflowrate/3600 * self.delta_H    # [kW] Total heat to be removed by the cooling system
    
                ########################################################################################à
                
            # Multistage compressor operation with refrigeration
            if parameters['compressor model'] == 'multistage compressor with refrigeration':
                
                self.model = parameters['compressor model']
                self.n_stages = parameters['n_stages']      
                self.T_IC = parameters['T_IC']
                self.delta_P = parameters['pressure losses IC']  
                self.T_is_points = []
                self.P_target = numpy.zeros(self.n_stages)
                self.eta_is = numpy.zeros(self.n_stages)
                self.beta_new = numpy.zeros(self.n_stages)
                self.delta_H = numpy.zeros(self.n_stages)
                
                for i in range(self.n_stages*2 + 1):            # point 0 is the input to the compressor, point 1 is the outlet of the compressor
    
                    self.P_points.append(self.P_in)
                    self.T_points.append(self.T_in)
                    self.h_points.append(PropsSI('H', 'T', self.T_points[i], 'P', self.P_points[i]*100000, self.fluid)/1000)
                    self.s_points.append(PropsSI('S', 'T', self.T_points[i], 'P', self.P_points[i]*100000, self.fluid)/1000)
                    self.rho_points.append(PropsSI('D', 'T', self.T_points[i], 'P', self.P_points[i]*100000, self.fluid))
                
                self.comp_lav_spec = []     
                self.comp_power_list = []
                self.IC_power_list = []
                self.comp_beta_targ = numpy.power(self.beta_max, 1 / self.n_stages)
                
                P_ref = 2
                T_ref = 298.15
                deltaT_app = 5
                                            
                for i in range(self.n_stages):
                                               
                    self.comp_lav_spec.append(0.)
                    self.comp_power_list.append(0.)
                    self.IC_power_list.append(0.)
                    self.P_target[i] = self.P_in*self.comp_beta_targ**(i+1)
                    self.delta_H = numpy.zeros(self.n_stages)
                    P_in = self.P_points[i*2]
                    h_in = self.h_points[i*2]
                    s_in = self.s_points[i*2]
                    T_in = self.T_points[i*2]
                    rho_in = self.rho_points[i*2]  
                    
                    # compression stage i
                    P_out_ic = self.P_target[i]
                    P_out = (P_out_ic)/(1-self.delta_P)
                    self.beta_new[i] = P_out/P_in
                    s_out_iso = s_in
                    h_out_iso = PropsSI('H', 'P', P_out*100000, 'S', s_out_iso*1000, self.fluid)/1000
                    self.T_is_points.append(PropsSI('T', 'P', P_out*100000, 'S', s_out_iso*1000, self.fluid))
                    self.eta_is[i] = (self.beta_new[i]**self.epsilon - 1)/(self.beta_new[i]**(self.epsilon/self.eta_pol)-1)
                    h_out = (h_out_iso - h_in) / self.eta_is[i] + h_in
                    
                    self.T_points[i*2 + 1] = PropsSI('T', 'P', P_out*100000, 'H', h_out*1000, self.fluid)
                    self.P_points[i*2 + 1] = P_out
                    self.h_points[i*2 + 1] = h_out
                    self.s_points[i*2 + 1] = PropsSI('S', 'P', P_out*100000, 'H', h_out*1000, self.fluid)/1000
                    self.rho_points[i*2 + 1] = PropsSI('D', 'P', P_out*100000, 'H', h_out*1000, self.fluid)
                    
                    # interrefrigeration stage i
                    self.P_points[i*2 + 2] = P_out_ic
                    self.T_points[i*2 + 2] = self.T_IC
                    self.h_points[i*2 + 2] = PropsSI('H', 'P', P_out_ic*100000, 'T', self.T_IC, self.fluid)/1000
                    self.s_points[i*2 + 2] = PropsSI('S', 'P', P_out_ic*100000, 'T', self.T_IC, self.fluid)/1000
                    self.rho_points[i*2 + 2] = PropsSI('D', 'P', P_out_ic*100000, 'T', self.T_IC, self.fluid)
                    
                    h_out_ic = PropsSI('H', 'P', P_out_ic*100000, 'T', self.T_IC, self.fluid)/1000
                    self.delta_H[i] = h_out - h_out_ic
                    self.comp_lav_spec[i] = h_out - h_in
                    
                if 'flow_rate' in parameters:       # if flow rate is defined as input - compressor nominal power is defined as a consequence
                    self.maxflowrate = parameters['flow_rate']                            # [kg/h] hourly mass flow rate
                    self.Npower = round(self.maxflowrate/3600*sum(self.comp_lav_spec))    # [kW] compressor nominal power 
                    self.IC_power_list = self.maxflowrate/3600 * sum(self.delta_H)             # [kW] Total heat to be removed by the cooling system
                elif self.maxflowrate:   
                    self.Npower = round(self.maxflowrate/3600*sum(self.comp_lav_spec))    # [kW] compressor nominal power 
                    self.IC_power_list = self.maxflowrate/3600 * sum(self.delta_H)             # [kW] Total heat to be removed by the cooling system
                elif 'Power' in parameters:         # if nominal power is defined as input - compressor mass flow rate is defined as a consequence
                    self.Npower = parameters['Power']
                    self.maxflowrate = self.Npower/sum(self.comp_lav_spec)     # [kg/s] Nominal mas flow rate
                    self.IC_power_list = self.maxflowrate * sum(self.delta_H)      # [kJ] Total heat to be removed by the cooling system

    ##########################################################################################################################################################
                    
                # # isoentalpic transformation, control from T to dispenser
                # Pout_dispenser = 350  # [bar]
                # self.P_points.append(Pout_dispenser) #Pout_dispenser bar, serbatoio trenio
                # self.h_points.append(h_out_ic)
                # self.T_points.append(PropsSI('T', 'P', Pout_dispenser*100000, 'H', h_out_ic*1000, self.fluid))
                # self.s_points.append(PropsSI('S', 'P', Pout_dispenser*100000, 'H', h_out_ic*1000, self.fluid)/1000)
                # self.rho_points.append(PropsSI('D', 'P', Pout_dispenser*100000, 'H', h_out_ic*1000, self.fluid))
                
    ################################################################################################################################################################
                # # Alternative calculation for compressor power (without IC at last compressor) Ref:The TechnoEconomics of Hydrogen Compression TECHNICAL BRIEF
                
                # if 'flow_rate' in parameters:       # if flow rate is defined as input - compressor nominal power is defined as a consequence
                #     m_day       = self.flow_rate*3600*24 # [kgH2/day]
                #     T_in        = self.T_in         # [K]
                #     P_in        = self.P_in         # [bar]
                #     P_out       = self.P_max        # [bar]
                #     beta_stage  = self.beta_new[0]  # compressore a diaframma
                #     eta_is      = self.eta_is[0]
                #     eta_motor   = 0.99              # efficiency of electric engine
                #     R_uni       = self.R_univ
                    
                #     N = round(m.log(P_out/P_in)/(m.log(beta_stage)))   # Number of stages of the compressor
                #     T_out = T_in*(1+((P_out/P_in)**(self.epsilon/N)-1)/eta_is)
                #     P_avg = (P_out+P_in)/2
                #     T_avg = (T_out+T_in)/2
                #     Z = PropsSI('Z', 'T', T_avg, 'P', P_avg*100000, 'Hydrogen')
                #     MM = c.H2MOLMASS
                #     q_m = (m_day/MM)/(3600*24)  # [mol/s]
                #     Power = N*1/self.epsilon*q_m/eta_is*T_in*Z*R_uni*((P_out/P_in)**(self.epsilon/N)-1)/1000 #kW
                    
                #     Power_shaft = Power/eta_motor
                #     print(Power_shaft)  
                    
        def thermodynamics_points(self):
            
            """
            Costruzione dei punti termodinamici
            
            """
    
            f = []
            sim_hour = 1
    
            if self.model == 'normal compressor':
                print('è solo una trasformazione adiabatica')
            else:
                points = self.n_stages*2+2
    
            
                #Parameters:
                # fluid_name (string or AbstractState) – The name of the fluid to be plotted or a state instance
                # graph_type (string) – The graph type to be plotted, like “PH” or “TS”
                # axis (matplotlib.pyplot.gca(), Optional) – The current axis system to be plotted to. Default: create a new axis system
                # fig (matplotlib.pyplot.figure(), Optional) – The current figure to be plotted to. Default: create a new figure
                # unit_system (string, ['EUR','KSI','SI']) – Select the units used for the plotting. ‘EUR’ is bar, kJ, C; ‘KSI’ is kPa, kJ, K; ‘SI’ is Pa, J, K
                # tp_limits (string, ['NONE','DEF','ACHP','ORC']) – Select the limits in T and p.
     
                plt.rcParams['figure.dpi'] = 1000
                plot = PropertyPlot('Hydrogen', 'Ts', unit_system='EUR')
                plot.set_axis_limits([0, 50, -250, 275])
                plot.calc_isolines(CoolProp.iQ, num=11)
                plot.calc_isolines(CoolProp.iP, iso_range=[35,350,525], num=3, rounding=True)
                plot.draw()
    
                for i in range(points):
                    f.append(Compressor(inp_test, sim_hour))
                    # plt.f(prova[0].P_points[i], prova[0].h_points[i], marker="o", color="red")
                    plt.plot(f[0].s_points[i], f[0].T_points[i]-273.15, marker="o", color="red")
                    plt.annotate(i+1, (f[0].s_points[i]+0.5, f[0].T_points[i]-273.15+10))
                plt.show()
            
                plot = PropertyPlot('Hydrogen', 'Ts', unit_system='EUR')
                plot.set_axis_limits([25, 35, 25, 200])
                plot.calc_isolines(CoolProp.iQ, num=11)
                plot.calc_isolines(CoolProp.iP, iso_range=[525,self.P_points[5]], num=2, rounding=True)
                plot.draw()
                plt.plot(f[0].s_points[5], f[0].T_points[5]-273.15, marker="o", color="red")
                plt.plot(f[0].s_points[6], f[0].T_points[6]-273.15, marker="o", color="red")
                plt.annotate(6, (f[0].s_points[5]+0.25, f[0].T_points[5]-283.15), fontsize = 14)
                plt.annotate(7, (f[0].s_points[6]+0.25, f[0].T_points[6]-263.15), fontsize = 14)
            
                plt.show()
                                   
        def parametric_stage(self):
                
            P_max = numpy.linspace(70,875,24)
            n_stages = numpy.linspace(1,4,4)
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
            
            # flow_rate = 0.05
            flow_rate_H2 = self.flow_rate
            flow_rate = flow_rate_H2*self.M_air/self.M_H2
            # flow_rate_H2 = self.M_H2/self.M_air*flow_rate
            
            beta = numpy.linspace(2,25,24)
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
                        'Power': numpy.zeros(len(beta))
                        }
            
            fluid = 'Hydrogen'
            
            prop_fluid = {'P': [],
                        'T': [],
                        'h': [],
                        's': [],
                        'rho': [],
                        'eta_is': [],
                        'Delta H': [],
                        'Power': numpy.zeros(len(beta))
                        }
            
            s_in_air = self.CP_air*m.log(T_in/T_0)+self.R_air*m.log(P_in/P_0)
            s_in_fluid = PropsSI('S', 'T', T_in, 'P', P_in*100000, fluid)/1000
            h_in_air = self.CP_air*(T_in-T_0)
            h_in_fluid = PropsSI('H', 'T', T_in, 'P', P_in*100000, fluid)/1000
            
            eta_pol = 0.7 
            
            eta_is_air = numpy.zeros(len(beta))
            eta_is_fluid = numpy.zeros(len(beta))
            
            self.power_fluid_corr = numpy.zeros(len(beta))
            self.fract_power = numpy.zeros(len(beta))
            
            for i in range(len(beta)):
                prop_air['P'].append(P_in*beta[i])
                prop_air['eta_is'].append((beta[i]**self.epsilon_air - 1)/(beta[i]**(self.epsilon_air/eta_pol)-1))
                T_is_a = T_in*beta[i]**self.epsilon_air
                h_is_a = self.CP_air*(T_is_a-T_0)
                prop_air['h'].append((h_is_a - h_in_air) / prop_air['eta_is'][i] + h_in_air)
                prop_air['T'].append(T_in+(prop_air['h'][i]-h_0_air)/self.CP_air)
                prop_air['s'].append(self.CP_air*m.log(prop_air['T'][i]/T_0)+self.R_air*m.log(prop_air['P'][i]/P_0))
                prop_air['Delta H'].append(prop_air['h'][i] - h_in_air)
                prop_air['Power'][i]=((prop_air['h'][i] - h_in_air)*flow_rate)
                
                
                prop_fluid['P'].append(P_in*beta[i])
                prop_fluid['eta_is'].append((beta[i]**self.epsilon - 1)/(beta[i]**(self.epsilon/eta_pol)-1))
                h_is_f = PropsSI('H', 'P', prop_fluid['P'][i]*100000, 'S', s_in_fluid*1000, fluid)/1000
                prop_fluid['h'].append((h_is_f - h_in_fluid) / prop_fluid['eta_is'][i] + h_in_fluid)
                prop_fluid['Delta H'].append(prop_fluid['h'][i] - h_in_fluid)
                prop_fluid['Power'][i] = ((prop_fluid['h'][i] - h_in_fluid)*flow_rate)
                
                self.fract_power[i] =  prop_fluid['Power'][i]/ prop_air['Power'][i]
                
                self.power_fluid_corr[i] = flow_rate_H2*(prop_fluid['h'][i] - h_in_fluid)
    
            
            # plt.figure(dpi=1000)    
            # plt.plot(beta, prop_air['Delta H'], ls = (0,(5,1)), label = "Air")
            # plt.plot(beta, prop_fluid['Delta H'], ls = (0,(5,1)), label = fluid)
            # plt.plot(beta[13], prop_air['Delta H'][13] , marker="o", color="purple")
            # plt.plot(beta[13], prop_fluid['Delta H'][13] , marker="o", color="purple")
            # plt.text(beta[13],prop_air['Delta H'][13]+300,str(round(prop_air['Delta H'][13])))
            # plt.text(beta[13],prop_fluid['Delta H'][13]-500,str(round(prop_fluid['Delta H'][13])))
    
            
            # plt.legend(fontsize=8)
            # plt.xlabel('$\\beta$ [-]')
            # plt.ylabel('$\\Delta$h [kJ/kg]')
    
            # plt.title('$\\Delta$h vs $\\beta$')
            # plt.grid()       
            # plt.show()
            
            # plt.figure(dpi=1000)    
            # plt.plot(beta, prop_air['Power'], ls = (0,(5,1)), label = "Air")
            # plt.plot(beta, prop_fluid['Power'], ls = (0,(5,1)), label = fluid)
            # plt.plot(beta[13], prop_air['Power'][13] , marker="o", color="purple")
            # plt.plot(beta[13], prop_fluid['Power'][13] , marker="o", color="purple")
            # plt.text(beta[13],prop_air['Power'][13]+15,str(round(prop_air['Power'][13])))
            # plt.text(beta[13],prop_fluid['Power'][13]-30,str(round(prop_fluid['Power'][13])))
    
            
            # plt.legend(fontsize=8)
            # plt.xlabel('$\\beta$ [-]')
            # plt.ylabel('W$_{c}$ [kW]')
    
            # plt.title('W$_{c}$ vs $\\beta$')
            # plt.grid()
            # plt.show()
            
            # plt.figure(dpi=1000)    
            # plt.plot(beta, self.fract_power)
            # plt.xlabel('$\\beta$ [-]')
            # plt.ylabel('$\\Delta$h$_{H2}$ / $\\Delta$h$_{air}$ [-]')
    
            # plt.title('$\\Delta$h$_{H2}$ / $\\Delta$h$_{air}$ vs $\\beta$')
            # plt.grid()
    
            # plt.show()
            
            
            
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
                          'Power': numpy.zeros(len(beta))
                          }
    
    
            beta_air = numpy.linspace(2,25,24)    
            cost = numpy.zeros(len(beta_air))
    
            beta_fluid = numpy.zeros(len(beta_air))
            self.delta_beta = numpy.zeros(len(beta_air))
            beta_corr = numpy.zeros(len(beta_air))
            n_stages = numpy.zeros(len(beta_air))
            
            costi_corr = {'Z_1': numpy.zeros(len(beta_air)),
                          'Z_2': numpy.zeros(len(beta_air))
                         }
    
            costi = {'Z_1': numpy.zeros(len(beta_air)),
                     'Z_2': numpy.zeros(len(beta_air))
                    }    
    
            for n in range(len(beta_air)):
                cost[n] = ((beta_air[n]**(self.epsilon_air)-1)*self.CP_air*T_in/0.85)
                beta_fluid[n] = (1+0.85/(self.CP_H2*T_in)*(cost[n]))**(1/self.epsilon)
                self.delta_beta[n] = (beta_air[n])/(beta_fluid[n])
                beta_corr[n] = self.delta_beta[n]*beta_air[n]
                n_stages[n] = m.log(beta_air[n])/m.log(beta_fluid[n])
            
            
            #Punti termodinamici dell'aria con il beta corretto#
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
                
            
                
            # plt.figure(dpi=1000)    
            
            # # plt.scatter(cost, beta_air,s=20,color='tab:blue',edgecolors='k',zorder=3, label = 'cost correlation 1')
            # # plt.scatter(cost, beta_fluid,s=20,color='tab:orange',edgecolors='k',zorder=3, label = 'cost correlation 2')
            # plt.scatter(cost, self.delta_beta,s=20,color='tab:pink',edgecolors='k',zorder=3, label = '$\\eta$$_{is,c}$ = 0.85')
    
            # plt.legend(fontsize=10)
            # plt.xlabel('$\\psi$u$^{2}$')
            # plt.ylabel('$\\beta$$_{air}$/$\\beta$$_{H2}$ [-]')
            # plt.title('$\\beta$$_{air}$/$\\beta$$_{H2}$ vs $\\psi$u$^{2}$')
            # plt.grid()
            # plt.show()
            
            # plt.figure(dpi=1000)
            # plt.scatter(beta_air, n_stages, s=20, color='tab:pink', edgecolor='k')
            # plt.legend(fontsize=10)
            # plt.xlabel('$\\beta$$_{air}$ [-]')
            # plt.ylabel('n$_{stages, equi}$')
    
            # # plt.title('$\\beta$$_{air}$/$\\beta$$_{H2}$ vs $\\psi$u$^{2}$')
            # plt.grid()
    
            # plt.show()
            
    
            for n in range(len(beta_corr)):
                # m air con beta corr #
                costi_corr['Z_1'][n] = (39.5*flow_rate*beta_corr[n])/(0.9-0.85)*(m.log(beta_corr[n]))
                costi_corr['Z_2'][n] = 44.71*flow_rate/(0.95-0.85)*beta_corr[n]*m.log(beta_corr[n])
                # m air con beta air (2,3,4....) #
                costi['Z_1'][n] = (39.5*flow_rate*beta_air[n])/(0.9-0.85)*(m.log(beta_air[n]))
                costi['Z_2'][n] = 44.71*flow_rate/(0.95-0.85)*beta_air[n]*m.log(beta_air[n])
            
                
            # plt.figure(dpi = 1000)
            # plt.plot(beta_air, costi_corr['Z_1'], color='gold', marker = 's', zorder=3, label = 'Saghafifar - Gadalla')
            # plt.plot(beta_air, costi_corr['Z_2'], color='darkkhaki', marker = 's', zorder=3, label = 'Roosen et al')
            # plt.legend(fontsize=8, bbox_to_anchor=(1.04, 1), borderaxespad=0)
            # plt.xlabel('$\\beta$ [-]')
            # plt.ylabel('Costi [€]')
            # plt.title('$\\beta$ vs Costi')
            # plt.grid()
    
            # plt.show()
            
            # plt.figure(dpi = 1000)
            # plt.plot(beta_corr, costi_corr['Z_1'], color='gold', marker = 's', zorder=3, label = 'Saghafifar - Gadalla')
            # plt.plot(beta_corr, costi_corr['Z_2'], color='darkkhaki', marker = 's', zorder=3, label = 'Roosen et al')
            # plt.legend(fontsize=8, bbox_to_anchor=(1.04, 1), borderaxespad=0)
            # plt.xlabel('$\\beta$$_{corr}$ [-]')
            # plt.ylabel('Costi [€]')
            # plt.title('$\\beta$$_{corr}$ vs Costi')
            # plt.grid()
    
            # plt.show()
            
            # plt.figure(dpi = 1000)
            # plt.plot(beta_air, costi['Z_1'], color='gold', marker = 's', zorder=3, label = 'Saghafifar - Gadalla')
            # plt.plot(beta_air, costi['Z_2'], color='darkkhaki', marker = 's', zorder=3, label = 'Roosen et al')
            # plt.legend(fontsize=8, bbox_to_anchor=(1.04, 1), borderaxespad=0)
            # plt.xlabel('$\\beta$ [-]')
            # plt.ylabel('Costi [€]')
            # plt.title('$\\beta$ vs Costi')
            # plt.grid()
    
            # plt.show()
            
    ###########################################################################################################################################################
            
            # costs expressed as a function of mass flow rate and integration of compression ratio (beta)
            
            costi_air = {'Z_1': numpy.zeros(len(beta)), 'Z_2': numpy.zeros(len(beta)), 'Z_3': numpy.zeros(len(beta)), 'Z_4': numpy.zeros(len(beta)),
                          'Z_5': numpy.zeros(len(beta)),'Z_6': numpy.zeros(len(beta)), 'Z_7': numpy.zeros(len(beta)), 'Z_8': numpy.zeros(len(beta)),
                          'Z_9': numpy.zeros(len(beta))
                          }
            
            costi_fluid = {'Z_1': numpy.zeros(len(beta)),'Z_2': numpy.zeros(len(beta)),'Z_3': numpy.zeros(len(beta)), 'Z_4': numpy.zeros(len(beta)),
                          'Z_5': numpy.zeros(len(beta)), 'Z_6': numpy.zeros(len(beta)),'Z_7': numpy.zeros(len(beta)), 'Z_8': numpy.zeros(len(beta)),
                          'Z_9': numpy.zeros(len(beta))
                          }
            
            #correzione della portata#
            costi_fluid_corr = {'Z_1': numpy.zeros(len(beta)), 'Z_2': numpy.zeros(len(beta)), 'Z_3': numpy.zeros(len(beta)), 'Z_4': numpy.zeros(len(beta)),
                                'Z_5': numpy.zeros(len(beta)), 'Z_6': numpy.zeros(len(beta)), 'Z_7': numpy.zeros(len(beta)), 'Z_8': numpy.zeros(len(beta)),
                                'Z_9': numpy.zeros(len(beta))
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
            # # Grafici parità di beta - costi#
            
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
            
    ###############################################################################################################################################################################
            #tolgo la relazione di costi Z_9 e Z_1 perchè vanno fuori scala e faccio fitting -------- > Alla fine non l'ho utilizzato
            
            # y = []
            
            # y.extend(costi_fluid_corr['Z_2'])
            # y.extend(costi_fluid_corr['Z_3'])
            # y.extend(costi_fluid_corr['Z_5'])
            # y.extend(costi_fluid_corr['Z_6'])
            # y.extend(costi_fluid_corr['Z_7'])
            # y.extend(costi_fluid_corr['Z_8'])
            # y.extend(costi_corr['Z_1'])
            # y.extend(costi_corr['Z_2'])
            
            # x = []
            # for i in range(int(len(y)/24)):
            #     x.extend(beta_air)
                    
            # def mapping_function(x,a,b,c,d):
            #     fit = a*x + b*x**2 + c*x**3 + d
            #     return fit
            
            # opt_param = curve_fit(mapping_function,x,y)
            # a,b,c,d = opt_param[0][0], opt_param[0][1], opt_param[0][2], opt_param[0][3]
            
            # plt.scatter(x,y,color='green', cmap='viridis', s=10, alpha=0.5, marker='X', label='Data')
            # x_line = numpy.arange(min(x),max(x),1)
            # y_line = mapping_function(x_line, a, b, c, d)
            # plt.plot(x_line, y_line, '--', color = 'red')
            # plt.xlabel('prova 1')
            # plt.ylabel('prova 2')
            # plt.title('fitt')
            # plt.legend(loc = 'lower right', prop = {'size': 8})
            # plt.show
            
            ###################################
            # curve = numpy.polyfit(x,y,3)
            # poly = numpy.poly1d(curve)  #equazione di curve#
            # print(curve)
            # new_x = []
            # new_y = []
            # for i in range(len(beta_air)):
            #     new_x.append(i+2)
            #     calc = poly(i+2)
            #     new_y.append(calc)
                
            # R = r2_score(y, poly(x))
            # plt.plot(new_x, new_y)
            # plt.show
                   
            
    ####################################################################################################################################################
            # Grafici parità di Potenza - costi#
            
            # for n in sorted(costi_air):
            #     if a <= 2:
            #         plt.plot(prop_air['Power'], costi_air[n], marker = 'o', color = color[a], zorder=3, label = label[a])
            #         a +=1 
            # #     elif a>=4 and a < 9:
            # #         plt.plot(prop_air['Power'], costi_air[n], marker = '^', color = color[a], zorder=3, label = label[a])
            # #         a +=1 
            # #     else:
            # #         a +=1
            
            # # plt.plot(prop_air['Power'], costi['Z_1'], marker = 's', color = color[9], zorder=3, label = label[9])
            # # plt.plot(prop_air['Power'], costi['Z_2'], marker = 's', color = color[10], zorder=3, label = label[10])
         
        
            # plt.legend(fontsize=8, bbox_to_anchor=(1.04, 1), borderaxespad=0)
            # plt.xlabel('W$_{c}$ [kW]')
            # plt.ylabel('Costi [€]')
            # plt.xlim([0,30])
            # plt.ylim([0,50000])
            # plt.title('W$_{c}$ vs Costi, m'+'\u0307'+'$_{air}$= '+str(round(flow_rate,4)))
            # plt.grid()
    
            # plt.show()
    
            # plt.figure(dpi=1000)    
            # for n in sorted(costi_fluid):
            #     if b <= 2:
            #         plt.plot(prop_fluid['Power'], costi_fluid[n], marker = 'o', color = color[b], zorder=3, label = label[b])
            #         b +=1 
            #     elif b>=4:
            #         plt.plot(prop_fluid['Power'], costi_fluid[n], marker = '^', color = color[b], zorder=3, label = label[b])
            #         b +=1
            #     else:
            #         b +=1
                    
            # plt.plot(prop_fluid['Power'], costi_corr['Z_1']*self.M_air/self.M_H2, marker = 's', color = color[9], zorder=3, label = label[9])
            # plt.plot(prop_fluid['Power'], costi_corr['Z_2']*self.M_air/self.M_H2, marker = 's', color = color[10], zorder=3, label = label[10])  
    
            # plt.legend(fontsize=8, bbox_to_anchor=(1.04, 1), borderaxespad=0)
            # plt.xlabel('W$_{c}$ [kW]')
            # plt.ylabel('Costi [€]')
            # plt.title('W$_{c}$ vs Costi, m'+'\u0307'+'$_{H2}$= '+str(flow_rate))
            # plt.grid()
    
            # plt.show()
            
    
            
            # plt.figure(dpi=1000)    
            # for n in sorted(costi_fluid_corr):
            #     if c <= 2:
            #         plt.plot(self.power_fluid_corr, costi_fluid_corr[n], marker = 'o', color = color[c], zorder=3, label = label[c])
            #         c +=1  
            #     elif c>=4:
            #         plt.plot(self.power_fluid_corr, costi_fluid_corr[n], marker = '^', color = color[c], zorder=3, label = label[c])
            #         c +=1
            #     else:
            #         c +=1
         
                    
            # # plt.plot(self.power_fluid_corr, costi_corr['Z_1'], marker = 's', color = color[9], zorder=3, label = label[9])
            # # plt.plot(self.power_fluid_corr, costi_corr['Z_2'], marker = 's', color = color[10], zorder=3, label = label[10])
      
            # # plt.plot(prop_air_2['Power'], costi_corr['Z_1'], marker = 's', color = color[9], zorder=3, label = label[9])
            # # plt.plot(prop_air_2['Power'], costi_corr['Z_2'], marker = 's', color = color[10], zorder=3, label = label[10])
      
    
            # plt.legend(fontsize=8, bbox_to_anchor=(1.04, 1), borderaxespad=0)
            # # plt.xlim([0,30])
            # # plt.ylim([0,50000])
            # plt.xlabel('W$_{c}$ [kW]')
            # plt.ylabel('Costi [€]')
            # plt.title('W$_{c}$ vs Costi, m'+'\u0307'+'$_{H2}$= '+str(round(flow_rate_H2,5)))
            # plt.grid()
    
            # plt.show()     
            
            
            
    ##########################################################################################################################################################################################################
            
            #Grafico costi aria + costi idrogeno equivalente + costi aria equivalente e beta corretto
                    
            # plt.figure(dpi=600)    
    
            # plt.plot(prop_air['Power'], costi_air['Z_1'], marker = 'o', color = color[0], zorder=3, label = label[0])        
            # plt.plot(prop_air['Power'], costi_air['Z_2'], marker = 'o', color = color[1], zorder=3, label = label[1])   
            # plt.plot(prop_air['Power'], costi_air['Z_3'], marker = 'o', color = color[2], zorder=3, label = label[2])   
            # # plt.plot(self.power_fluid_corr, costi_fluid_corr['Z_5'], marker = '^', color = color[4], zorder=3, label = label[4])
            # # plt.plot(self.power_fluid_corr, costi_fluid_corr['Z_6'], marker = '^', color = color[5], zorder=3, label = label[5])
            # # plt.plot(self.power_fluid_corr, costi_fluid_corr['Z_7'], marker = '^', color = color[6], zorder=3, label = label[6])
            # # plt.plot(self.power_fluid_corr, costi_fluid_corr['Z_8'], marker = '^', color = color[7], zorder=3, label = label[7])
            # # plt.plot(self.power_fluid_corr, costi_fluid_corr['Z_9'], marker = '^', color = color[8], zorder=3, label = label[8])
            # # plt.plot(prop_air_2['Power'], costi_corr['Z_1'], marker = 's', color = color[9], zorder=3, label = label[9])        
            # # plt.plot(prop_air_2['Power'], costi_corr['Z_2'], marker = 's', color = color[10] ,zorder=3, label = label[10])   
    
            # plt.legend(fontsize=8, bbox_to_anchor=(1.04, 1), borderaxespad=0)
            # plt.xlabel('W$_{c}$ [kW]')
            # plt.xlim([0,30])
            # plt.ylim([0,50000])
            # plt.ylabel('Costi [€]')
            # plt.title('W$_{c}$ vs Costi, m'+'\u0307'+'$_{H2}$= '+str(round(flow_rate_H2,5)))
            # plt.grid()
            # plt.show() 
            
            # plt.figure(dpi=1000)    
    
            # plt.plot(beta_air, costi_air['Z_1'], marker = 'o', color = color[0], zorder=3, label = label[0])        
            # plt.plot(beta_air, costi_air['Z_2'], marker = 'o', color = color[1], zorder=3, label = label[1])   
            # plt.plot(beta_air, costi_air['Z_3'], marker = 'o', color = color[2], zorder=3, label = label[2])   
            # plt.plot(beta_air, costi_fluid_corr['Z_5'], marker = '^', color = color[4], zorder=3, label = label[4])
            # plt.plot(beta_air, costi_fluid_corr['Z_6'], marker = '^', color = color[5], zorder=3, label = label[5])
            # plt.plot(beta_air, costi_fluid_corr['Z_7'], marker = '^', color = color[6], zorder=3, label = label[6])
            # plt.plot(beta_air, costi_fluid_corr['Z_8'], marker = '^', color = color[7], zorder=3, label = label[7])
            # plt.plot(beta_air, costi_fluid_corr['Z_9'], marker = '^', color = color[8], zorder=3, label = label[8])
            # plt.plot(beta_air, costi_corr['Z_1'], marker = 's', color = color[9], zorder=3, label = label[9])        
            # plt.plot(beta_air, costi_corr['Z_2'], marker = 's', color = color[10], zorder=3, label = label[10])   
           
            # plt.legend(fontsize=8, bbox_to_anchor=(1.04, 1), borderaxespad=0)
            # plt.xlabel('$\\beta$ [-]')
            # plt.ylabel('Costi [€]')
            # plt.title('$\\beta$ vs Costi, m'+'\u0307'+'$_{H2}$= '+str(round(flow_rate_H2,5)))
            # plt.grid()
    
            # plt.show() 
            
############################################################################################################################################################################################################    


    def use(self, h, available_hyd_lp=False, storable_hydrogen_hp=False, massflowrate=False):
        
        """
        compressed hydrogen
    
        storable_hydrogen_hp: float storable hydrogen H tank max_capacity - SOC[h-1] [kg] or maximum absorbable production if there is no tank
        available_hyd_lp: float available hydrogen H tank SOC[h-1] [kg]
        massflowrate: float hydrogen flow rate to be processed in the timestep [kg/h]
        h : timestep float timestep in hours [h]

        output : 
        float hydrogen compressed in the timestep [kg]    
        float electricity absorbed that hour [kWh]
        float cooling needs [kWh]
        """
        if self.model == 'simple compressor':
            
            self.hyd[h] = massflowrate
            e_absorbed  = self.en_cons*massflowrate
            t_absorbed  = 0                             # cooling not considered in this model
            
            return(self.hyd[h],-e_absorbed,t_absorbed)
        
        if self.model == "normal compressor":
            
            self.hyd[h] = self.flow_rate*self.timestep*3600 
            e_absorbed = (self.Npower)*self.timestep           # [kWh]energia del compressore in un'ora di utilizzo
            t_absorbed = 0
            
            if self.hyd[h] > available_hyd_lp:        # if the is not enough space in the H tank to store the hydrogen (H tank is nearly full)
                self.hyd[h] = 0
                e_absorbed = 0
                t_absorbed = 0
            
            return(self.hyd[h],-e_absorbed,-t_absorbed)
                
        elif self.model == 'multistage compressor with refrigeration' or self.model == 'compressor with refrigeration':
            
            if massflowrate == False:       # compressor working in between a buffer and an HPH storage. Only full load operations allowed. 
                self.hyd[h] = self.maxflowrate*self.timestep        # [kg/h] hydrogen mass flow rate
                t_absorbed = (self.IC_power_list)*self.timestep     # [kWh] of required cooling
                e_absorbed = (self.Npower)*self.timestep            # [kWh] absorbed energy 
                    
                if self.hyd[h] > available_hyd_lp:                  # if the is not enough space in the H tank to store the hydrogen (H tank is nearly full)
                    self.hyd[h] = 0
                    e_absorbed = 0
                    t_absorbed = 0
                    
                elif self.hyd[h] > storable_hydrogen_hp:            # if mass flow rate is higher than the storable amount of hydrogen
                    self.hyd[h] = 0
                    e_absorbed = 0
                    t_absorbed = 0
                            
                return(self.hyd[h],-e_absorbed,-t_absorbed)
            
            else:   # working with a single hydrogen tank. Simplification of partial load functioning - Can be upgraded
                self.hyd[h] = massflowrate
                e_absorbed  = massflowrate/3600*sum(self.comp_lav_spec)
                t_absorbed  = massflowrate/3600*sum(self.delta_H)
                
                return(self.hyd[h],-e_absorbed,-t_absorbed)
                
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
            'refud': dict
                'rate': float, percentage of initial investment which will be rimbursed [%]
                'years': int, years for reimbursment
            'replacement': dict
                'rate': float, replacement cost as a percentage of the initial investment [%]
                'years': int, after how many years it will be replaced
        """
        compressor_costs = {                                                                        # Ref: https://hsweb.hs.uni-hamburg.de/projects/star-formation/hydrogen/P2H_Full_Study_FCHJU.pdf
                            'Gas Flow [kg/h]': [1,10,100,1000,10000],                               # dataset to be updated 
                            'Compressed Tanks [€/(kg/h)]': [63460,29007,13259,6016,2770]
                            }
        
        costdf = pd.DataFrame(compressor_costs)
        
        x = costdf['Gas Flow [kg/h]']
        y = costdf['Compressed Tanks [€/(kg/h)]']
        
        f = interp1d(x, y)
        
        specific_cost = f(self.maxflowrate)  # [€/(kg/h)]
        
        tech_cost = {key: value for key, value in tech_cost.items()}

        if self.model == 'simple compressor':
            exchange_rate   = 0.91                                           # [2015USD/2023€]  exchange rate between USD and €
            correlation1    = (51901*(self.maxflowrate**0.65))*exchange_rate # Ref: https://www.sciencedirect.com/science/article/pii/S0360319919330022?via%3Dihub
            correlation = correlation1
        elif self.model != 'simple compressor':
            size = self.Npower
            IF              = 1.3                                            # [-] Installation Factor
            correlation2    = 63684.6*(self.Npower**0.4603)*IF               # Ref: https://transitionaccelerator.ca/wp-content/uploads/2021/10/TA-Briefs-1.2-The-Techno-Economics-of-Hydrogen-Compression-FINALPDF.pdf
            correlation = correlation2
        
        if tech_cost['cost per unit'] == 'default price correlation':
            # C = specific_cost*self.maxflowrate
            C = correlation                        # [€] 
        else:
            C = size * tech_cost['cost per unit']   # [€]

        tech_cost['total cost'] = tech_cost.pop('cost per unit')
        tech_cost['total cost'] = C
        tech_cost['OeM'] = tech_cost['OeM'] *C /100 # €

        self.cost = tech_cost

###########################################################################################################################################################################

if __name__ == "__main__":
    
    from matplotlib.patches import Patch

    inp_test = {    'P_out' : 525,
                    'P_in'  : 35,
                    'compressor model' : "multistage compressor with refrigeration",
                    # 'compressor model' : "normal compressor",
                    # 'flow_rate': 0.0034799,
                    # 'flow_rate': 0.0277777778,
                    'flow_rate': 170,
                    'fluid': 'Hydrogen',
                    'T_in': 343.15,
                    'pressure losses IC': 0.00,
                    'T_IC': 308.15,
                    'n_stages': 3
                    }

    sim_hour = 8760
    a = []
    jojo = []
    
    fluid  = 'Hydrogen'
    jojo = Compressor(inp_test, sim_hour)
    # a = jojo.thermal_power()
    # b = jojo.comp_power
    # c = jojo.thermodynamics_points()
    f = jojo.use(1, 1500, 10000)

    
    # plo = jojo.parametric_stage()
    # pp = jojo.fluid_vs_air()
    
    # plot m = 100 kg/h, tank 350 kg e 3500 kg####
    LCOH = [10.269,10,621,10.148,10.334,0,0]
    c = numpy.zeros(10)
    beta = inp_test['P_out']/inp_test['P_in']
    
    CP_H2 = 14.3060971
    CP_air = 1.0063081425141251
    epsilon = 0.28843752220862523
    epsilon_air = 0.28661432606258475
    M_H2 = 2.01588
    M_air = 28.96547
    flow_rate_H2 = inp_test['flow_rate']/3600
    flow_rate_air = flow_rate_H2*M_air/M_H2
    # flow_rate_air = inp_test['flow_rate']

    
    cost = ((beta**(epsilon_air)-1)*CP_air*inp_test['T_in']/0.85)
    beta_fluid = (1+0.85/(CP_H2*inp_test['T_in'])*(cost))**(1/epsilon)
    delta_beta = beta/beta_fluid
    beta_corr = delta_beta*beta
    
    
    # c[0] = 91562*(13.163/455)**0.67
    # c[1] = 9642.2*13.163**0.46
    # c[2] = 10167.5*13.163**0.46
    # c[3] = (39.5*flow_rate_air*beta)/(0.9-0.85)*(m.log(beta))
    # c[4] = 44.71*flow_rate_air/(0.95-0.85)*beta*m.log(beta)
    # c[5] = 0
    # c[6] = 0
    # c[7] = 0
    # c[8] = 0
    # c[9] = 0
    
    c[0] = 91562*(192.87/455)**0.67
    c[1] = 9642.2*192.87**0.46
    c[2] = 10167.5*192.87**0.46
    c[3] = (39.5*flow_rate_air*beta_corr)/(0.9-0.85)*(m.log(beta_corr))
    c[4] = 44.71*flow_rate_air/(0.95-0.85)*beta_corr*m.log(beta_corr)
    c[5] = 25421*192.87**0.61
    c[6] = 40038*192.87**0.61
    c[7] = 82790*192.87**0.4603
    c[8] = 80873*192.87**0.6038
    c[9] = 6186.6*192.87**0.8355


    
    
    color = ['orange', 'lightgreen', 'blue', 'gold', 'darkkhaki'
               , 'brown', 'pink', 'indigo', 'turquoise', 'salmon'
             ]
    label = ['Sadeghi', 'Parikhani', 'Rezayanet', 'Saghafifar -\nGadalla','Roosen'
               , 'Syed', 'Blazquez-Diaz', 'HDSAM -\nReciprocating -\n350 bar', 'HDSAM -\nReciprocating- \n700 bar', 'HDSAM -\nCentifugal'    
             ]
    p =  ('Sadeghi', 'Parikhani', 'Rezayanet', 'Saghafifar -Gadalla','Roosen'
            , 'Syed', 'Blazquez-Diaz', 'HDSAM Reciprocating 350 bar', 'HDSAM Reciprocating 700 bar',  'HDSAM Centifugal'
          )
    
    plt.figure(dpi = 1000)
    plt.rc('xtick', labelsize=6)    # fontsize of the tick labels
    plt.bar(label, c, color = color)
    cmap = dict(zip(c, color))
    patches = [Patch(color=v, label=k) for k, v in cmap.items()]
    plt.legend(title = 'Correlazioni di costo', labels = p, handles=patches, bbox_to_anchor=(1.04, 0.5), loc='center left', borderaxespad=0)
    plt.grid(color = 'grey', linestyle = '-', linewidth = 0.5, alpha = 0.4)
    plt.ylabel('Z$_{comp}$ [€]', fontsize = 8)
    plt.ylim([0,2000000])
    plt.show()
    
