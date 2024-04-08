import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
from sklearn.linear_model import LinearRegression
from numpy import log as ln
import math
import os
import sys 
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),os.path.pardir)))   # temorarily adding constants module path 
from core import constants as c

class electrolyzer:
    
    def __init__(self,parameters,timestep_number,timestep=False):
        """
        Create an electrolyzer object
    
        parameters : dictionary
            'Npower': float nominal power [kW]
            'number of modules': str  number of modules in the stack [-]
            'stack model': str 'Enapter 2.1','McLyzer 800' are aviable or 'PEM General'
            'minimum_load': 0-1 float [%], if specified, minimum load the electrolyser must operate at
            'strategy': str - 'full-time'. Electrolyzers operational 24/7, grid connection must be present. 
                            - 'hydrogen-first'. Electrolyzers working only when renewable power is available, 
                               prioritizing production of hydrogen over electricity production from RES
                      
        output : electrolyzer object able to:
            abrosrb electricity and water and produce hydrogen and oxygen .use(p)
        """

        self.model      = parameters['stack model']         # [-] selected electorlyzer model
        self.Npower     = parameters['Npower']              # [kW] float nominal power of electrolyzers installed capacity for the location
        if self.model in ['PEM General','Alkaline'] and self.Npower> 1000:
            raise ValueError(f"Warning: {self.Npower} kW of rated power has been selected for the single PEM or Alkaline electrolyser module (Detailed Models). \n\
            The maximum capacity is 1000 kW.\n\
            Options to fix the problem: \n\
                (a) -  Global electrolysis capacity can be increased by adding more modules in electrolyser parameters in studycase.json")
        if self.model == 'Alkaline' and self.Npower < 10:
            raise ValueError(f"Warning: {self.Npower} kW of rated power has been selected for the Alkaline electrolyzer model.\n\
            The minimum capacity for this model is 10 kW according to the manufacturer's specifications.\nPlease select a nominal power of 10 kW or higher for the Alkaline model.")
        self.n_modules          = parameters['number of modules']       # [-] number of modules in the stack
        self.strategy           = parameters['strategy']                # definig operational strategy for electrolyzers
        self.only_renewables    = parameters['only_renewables']
        self.min_load           = parameters.get('minimum_load', 0) # if 'minimum load' is not specified as model input, the default value of 0 is selected by default
        self.ageing             = parameters.get('ageing', False)   # if 'ageing' is not specified as model input, the default value is set to False 
        self.min_year   = c.MINUTES_YEAR                    # [min/year]    number of minutes in one year
        self.min_week   = c.MINUTES_WEEK                    # [min/week]    number of minutes in one week
        self.min_month  = c.MINUTES_MONTH                   # [min/month]   number of minutes in one month
        self.rhoNrh2    = c.H2NDENSITY                      # [kg/m^3]      hydrogen density under normal conditions
        self.rhoStdh2   = c.H2SDENSITY                      # [kg/m^3]      hydrogen density under standard conditions
        self.H2MolMass  = c.H2MOLMASS                       # [kg/mol]      Hydrogen molar mass
        self.H2VolMol   = c.H2NMOLVOL                       # [Nm^3/mol]    Hydrogen molar volume at Normal conditions (T = 0°C, P = 101325 Pa)
        self.O2MolMass  = c.O2MOLMASS                       # [kg/mol]      Hydrogen molar mass
        self.h2oMolMass = c.H2OMOLMASS                      # [kg/mol]      Water molar mass
        self.rhoStdh2o  = c.H2OSDENSITY                     # [kg/m3]       H2O density @ T = 15°C p = 101325 Pa
        self.H2_lhv     = c.LHVH2                           # [MJ/kg]       Hydrogen Lower Heating Value
        self.lhv_nvol   = c.LHV_H2NVOL                      # [kWh/Nm3]     Hydrogen volumetric LHV under normal conditions
        self.watercons  = 0.015                             # [m^3/kgH2]    cubic meters of water consumed per kg of produced H2. Fixed value of 15 l of H2O per kg of H2. https://doi.org/10.1016/j.rset.2021.100005
        self.cp_water   = c.CP_WATER*1000                   # [J/kgK]      Water Mass specific constant pressure specific heat
        # self.CF         = np.zeros(timestep_number)      # [%] electrolyer stack Capacity Factor
        self.cost = False # will be updated with tec_cost()
        
        if timestep == False: 
            self.timestep           = c.timestep            # [min] simulation timestep if code is launched from main
            self.timestep_number    = c.timestep_number     # [-]   number of timestep if code is launched from main  
        else:
            self.timestep           = timestep              # [min] simulation timestep if code is launched from electrolyzer.py
            self.timestep_number    = timestep_number       # [-]   number of timestep if code is launched from electrolyzer.py
        
        self.timesteps_year =  self.min_year/self.timestep  # [-] number of simulation steps in one year for the given simulaiton inputs
        self.timesteps_week =  self.min_week/self.timestep  # [-] number of simulation steps in one year for the given simulaiton inputs
        
        "2H2O --> 2H2 + O2"                             # Electorlysis chemical reaction general
        
        self.oxy = self.O2MolMass/(2*self.H2MolMass)    # [gO2/gH2] g O2 produced every g H2 

        if parameters['stack model'] == 'Enapter 2.1':      # https://www.enapter.com/it/newsroom/enapter-reveals-new-electrolyser-el-2-1
            stack_operative_power_consumption = 2.4         # [kW] Single module nominal power
            stack_production_rate_L = 500/3600              # [Nl/s]
            stack_production_rate = stack_production_rate_L / 1000 * self.rhoNrh2 # [kg/s]
            self.spec_cons = stack_operative_power_consumption/ stack_production_rate_L*1000                 # [kJ/Nm3] electrolyzer specific consumption at nominal mass flow rate production
            self.production_rate = stack_production_rate * (self.Npower / stack_operative_power_consumption) # [kg/s] actual production defined in a simplified way 
                                                                                                             # by dividing the selected electorlyzer size in the location
                                                                                                             # for the single unit nominal power
                                                                                                             
            self.MaxPowerStack = self.Npower*self.n_modules     # [kW] electrolyzer stack total power
            self.eff = self.lhv_nvol/self.spec_cons             # [-] LHV efficiency 
            self.n_modules_used=np.zeros(timestep_number)      # array containing modules used at each timestep
            self.EFF = np.zeros(timestep_number)               # keeping track of the elecrolyzer efficiency over the simulation
            self.EFF_last_module = np.zeros(timestep_number)       # last module efficiency array initialization
        
        if parameters['stack model'] == 'McLyzer 800':   # https://mcphy.com/it/apparecchiature-e-servizi/elettrolizzatori/large/
            stack_operative_power_consumption = 4000     # [kW] Single module nominal power
            stack_production_rate_m3 = 800/3600               # [Nm3/s]
            stack_production_rate = stack_production_rate_m3 * self.rhoNrh2   # [kg/s]
            self.spec_cons = stack_operative_power_consumption/ stack_production_rate_L                      # [kJ/Nm3] 
            self.production_rate = stack_production_rate * (self.Npower / stack_operative_power_consumption) # [kg/s] actual production defined in a simplified way 
                                                                                                             # by dividing the selected electorlyzer size in the location                                                                                                     # for the single unit nominal power      
                                                                                                             # for the single unit nominal power
            self.MaxPowerStack = self.Npower*self.n_modules     # [kW] electrolyzer stack total power            self.eff = self.lhv_nvol/self.spec_cons             # [-] LHV efficiency 
            self.eff = self.lhv_nvol/self.spec_cons             # [-] LHV efficiency
            self.n_modules_used=np.zeros(timestep_number)      # array containing modules used at each timestep
            self.EFF = np.zeros(timestep_number)               # keeping track of the elecrolyzer efficiency over the simulation
            self.EFF_last_module = np.zeros(timestep_number)       # last module efficiency array initialization
        
        if parameters['stack model'] == 'Hylizer 1000':  # https://www.ie-net.be/sites/default/files/Presentatie%205_Baudouin%20de%20Lannoy.pdf
            stack_operative_power_consumption = 5000     # [kW] Single module nominal power
            stack_production_rate_m3 = 1000/3600              # [Nm3/s]
            stack_production_rate = stack_production_rate_m3 * self.rhoNrh2   # [kg/s]
            self.spec_cons = stack_operative_power_consumption/ stack_production_rate_L                      # [kJ/Nm3] 
            self.production_rate = stack_production_rate * (self.Npower / stack_operative_power_consumption) # [kg/s] actual production defined in a simplified way 
                                                                                                             # by dividing the selected electorlyzer size in the location
                                                                                                             # for the single unit nominal powe
            self.MaxPowerStack = self.Npower*self.n_modules     # [kW] electrolyzer stack total power                                                                             
            self.eff = self.lhv_nvol/self.spec_cons             # [-] LHV efficiency 
            self.n_modules_used=np.zeros(timestep_number)      # array containing modules used at each timestep
            self.EFF = np.zeros(timestep_number)               # keeping track of the elecrolyzer efficiency over the simulation
            self.EFF_last_module = np.zeros(timestep_number)       # last module efficiency array initialization
        
        if parameters['stack model'] == 'PEM General':
            
            self.EFF                    = np.zeros(timestep_number)     # keeping track of the elecrolyzer efficiency over the simulation
            self.wat_cons               = np.zeros(timestep_number)     # water consumption array initialization
            self.EFF_last_module        = np.zeros(timestep_number)     # last module efficiency array initialization
            self.wat_cons_last_module   = np.zeros(timestep_number)     # last module water consumption initialization
            self.n_modules_used         = np.zeros(timestep_number)     # array containing modules used at each timestep
            self.cell_currdens          = np.zeros(timestep_number)     # cell current density at every hour
            Runiv                       = c.R_UNIVERSAL                 # [J/(mol*K)] Molar ideal gas constant
            self.FaradayConst           = c.FARADAY                     # [C/mol]     Faraday constant
            self.LHVh2                  = c.LHVH2                       # [MJ/kg]     H2 LHV
            self.HHVh2Mol               = c.HHVH2MOL                    # [kJ/mol]    H2 HHV molar
            
            # Math costants
            self.eNepero      = c.NEPERO         # [-]         Euler's number
            # Ambient conditions 
            self.AmbTemp      = c.AMBTEMP        # [K]         Standard ambient temperature - 15 °C

            """
            Anode  : 2H2O --> O2 + 4H+ + 4e-\n

            Cathode: 4H+ + 4e- --> 2H2"

            Overall reaction: 2H2O(l) + electrical energy --> 2H2(g) + O2(g)"
            
            # https://doi.org/10.1016/j.ijhydene.2008.11.083    # Electrolyzer Parameters - depending on the different types of electrolyzers chosen - model to be specified also in class - name
            """
            # self.MembThickness       = 250           # [micron]   
            self.MembThickness       = 158.1842         # [micron]
            self.Lambda              = 20               # [-] Cell mositure content
            self.AnodeCurrDensity    = 0.00013          # [A/cm^2]
            self.AnodePressure       = 101325           # [Pa]
            self.CathodeCurrDensity  = 0.001            # [A/cm^2]
            self.CTCanode            = 0.6              # [-] Charge transfer coefficient - Anode   https://www.sciencedirect.com/science/article/pii/S0360319918309017
            self.CTCcathode          = 0.45             # [-] Charge transfer coefficient - Cathode //      //      //
            self.CurrDensityMin      = 0.005            # [A/cm^2]
            self.OperatingTemp       = 273.15 + 70      # [K]
            self.OperatingPress      = 3000000          # [Pa]
            self.MinInputPower       = 0.1*self.Npower  # [kW] minimum input power chosen as 10% of module nominal power 
            self.MaxPowerStack       = self.n_modules*self.Npower           # [kW] electrolyzer stack total power
            self.min_partial_load    = self.min_load*self.MaxPowerStack     # [kW] minimum operational load for the electrolyzer during simulation
            self.nc                  = 10+int((self.Npower/1000)*(35-10))   # For a power range between 0kW and 1000kW the number of cells in the stack varies between 10 and 35 
            self.CurrDensityMax      = 2.1+(self.Npower/1000)*(3-2.1)       # For a power range between 0kW and 1000kW the maximum current density varies between 2.1 and 3 A/cm2 
            
            #    self.CurrDensityMax = 2            # [A/cm^2] https://www.sciencedirect.com/science/article/pii/S266638642030151X#:~:text=In%20contrast%2C%20PEM%20electrolyzers%20experience,at%20high%20current%20density%20operations.&text=While%20commercial%20electrolyzers%20typically%20operate,reported%20by%20Lewinski%20et%20al.
          
            # Computing the electrolyzer polarization curve V-i
            
            'POLARIZATION CURVE'
            
            Ndatapoints = 1000                     # Number of points used to compute the polarization curve 
    
            self.CurrDensityMax_id = self.CurrDensityMax+1    # Calculations done with a higher CurrDensMax because the cell mass transport isn't valid in correspondence of CurrDensMax
            self.CellCurrDensity   = np.linspace(self.CurrDensityMin,self.CurrDensityMax_id,Ndatapoints)  # [A/cm^2]
            self.CellVoltage       = np.zeros(Ndatapoints)    # [V] 

            'Polarization (V-i) curve calculation'
            
            # V is obtained by summing 4 different contributions 
                
            '1- Cell open curcuit voltage'
            
            pH2O = self.eNepero**(11.676-(3816.44 /(self.OperatingTemp-46.13)))  # [atm] 
            pO2  = self.AnodePressure/101325 - pH2O                              # [atm]
            pH2  = self.OperatingPress/101325 - pH2O                             # [atm]
            
            Ecell = 1.229-0.85e-3*(self.OperatingTemp-298.15)+4.3085e-5*self.OperatingTemp*ln(pH2*(pO2**0.5)/pH2O)   # [V] Cell open curcuit voltage
            
            for i in range(0,Ndatapoints):
                
                '2- Cell activation overpotential'
                
                Vact_an  = (Runiv*self.OperatingTemp*ln(self.CellCurrDensity[i]/self.AnodeCurrDensity))/(2*self.CTCanode*self.FaradayConst)     # [V]
                Vact_cat = Runiv*self.OperatingTemp*ln(self.CellCurrDensity[i]/self.CathodeCurrDensity)/(4*self.CTCcathode*self.FaradayConst)   # [V]
                  
                Vact = Vact_an + Vact_cat                                      # [V] Cell activation overpotential
            
                '3- Cell mass transport overpotential'
                
                MaxCurrDensity_2 = self.CurrDensityMax_id + 0.00001            # [A/cm^2]
                Vdiff = -Runiv*self.OperatingTemp*ln(1-self.CellCurrDensity[i]/MaxCurrDensity_2)/(2*self.FaradayConst)  # [V] Cell mass transport overpotential
            
                '4- Cell Ohmic losses'
                  
                MembConductivity = (0.005139*self.Lambda-0.00326)*self.eNepero**(1268*(1/303 - 1/self.OperatingTemp))   # [S/cm]
                Rcell = (self.MembThickness/10000)/MembConductivity                                                     # [cm^2/S]
                  
                Vohmic = self.CellCurrDensity[i]*Rcell                                                                  # [V] Cell Ohmic losses
                
                'Resulting stack voltage'
                
                self.CellVoltage[i]=Ecell+Vact+Vdiff+Vohmic                    # [V] - Cell voltage
            
            self.Voltage = self.nc*self.CellVoltage                            # [V] - Module voltage
            
            Ndata                  = int(Ndatapoints*(self.CurrDensityMax)/(self.CurrDensityMax_id))  
            self.CurrentDensityMax = self.CellCurrDensity[Ndata]
            self.CellCurrDensity   = self.CellCurrDensity[:Ndata]
            self.Voltage           = self.Voltage[:Ndata]
            
            self.CellArea          = (self.Npower/(self.CurrDensityMax*1e-3*self.CellVoltage[Ndata-1]))/self.nc       # [cm^2] cell active area  file:///C:/Users/Andrea/Downloads/1-s2.0-S0360319913002607-main.pdf  up to 5000cm2 (Fig.15 and Table A)            
            self.Current           = self.CellCurrDensity*self.CellArea # [A]

            'Interpolation of calculated functioning points to detect the best fit-function for i-V curve'

            self.num    = Ndatapoints                                            # Number of intervals to be considered for the interpolation
            self.x2     = np.linspace(0.05,max(self.CellCurrDensity),self.num)    # Setting xlim for range of validity of LinRegression Calculation - Only for plot-related reasons 

            # Interpolation
            self.iV1    = interp1d(self.CellCurrDensity,self.Voltage,bounds_error=False)                # Linear spline 1-D interpolation                                                                                                                                                                                                                                                
           
            # Defining Electrolyzer Max Power Consumption
            Power_inp = []                                      # [kW] Initializing power input series
            for i in range(len(self.Current)):
                pot = (self.Current[i]*self.Voltage[i])/1000    # [kW]
                Power_inp.append(pot)
                    
            # Interpolation
            self.PI = interp1d(Power_inp,self.Current,bounds_error=False,fill_value='extrapolate')        # Linear spline 1-D interpolation
            
            'Single module H2 production'
            hydrogen = []
            etaFar   = []
            etaEle   = []
            
            for i in range(len(Power_inp)):
               
                'Electrolyzer efficiency'    
                      
                etaFaraday  = 9.95e-1*self.eNepero**((-9.5788-0.0555*self.OperatingTemp)/(self.Current[i]/(self.CellArea*self.nc/10000))+ \
                              (1502.7083-70.8005*self.OperatingTemp)/((self.Current[i]/(self.CellArea*self.nc/10000))**2))                   # [-] Faraday efficiency
                    
                Vstack      = (Power_inp[i]/self.Current[i])*1000        # [V] Stack operating voltage
                etaElectr   = self.nc*self.LHVh2*1e6*self.H2MolMass*etaFaraday/(2*Vstack*self.FaradayConst)                                   # [-] Electric efficiency
                                                                                                                                            # etaElectr - ref. to justify 'etaElectr' term presence in the equation --- https://www.taylorfrancis.com/books/edit/10.1201/b19096/pem-electrolysis-hydrogen-production-hui-li-haijiang-wang-dmitri-bessarabov-nana-zhao
                etaFar.append(etaFaraday)
                etaEle.append(etaElectr)       
                 
                'Hydrogen Production'
                      
                HydroProdMol  = (etaFaraday*self.nc*self.Current[i])/(2*self.FaradayConst)   # [mol/s] (Guilbert 2020)                                                                                                                                                                                                                                                                                                                                                                                                              
                hyd           = HydroProdMol*self.H2MolMass                                  # [kg/s] hydrogen produced                     
                hyd_vol       = hyd/self.rhoNrh2                                             # [Nm^3/s] hydrogen produced        
                deltaHydrogen = hyd_vol*self.LHVh2*self.rhoNrh2*1000                         # [kW] power produced, in the form of hydrogen
                
                hydrogen.append(hyd)     
            
            self.maxh2prod          = round(max(hydrogen),8)            # [kg/s] maximum amount of produced hydrogen for the considered module
            self.maxh2prod_stack    = self.maxh2prod*self.n_modules     # [kg/s] maximum amount of produced hydrogen for the considered stack 
          
            self.etaF       = interp1d(hydrogen,etaFar,bounds_error=False,fill_value='extrapolate')       # Linear spline 1-D interpolation -> produced H2 - Faraday efficiency
            self.etaEle     = interp1d(hydrogen,etaEle,bounds_error=False,fill_value='extrapolate')       # Linear spline 1-D interpolation -> produced H2 - Electric efficiency
            self.h2P        = interp1d(hydrogen,Power_inp,bounds_error=False,fill_value='extrapolate')    # Linear spline 1-D interpolation -> produced H2 - Power consumption 
            self.P2h        = interp1d(Power_inp,hydrogen,bounds_error=False,fill_value='extrapolate')    # Linear spline 1-D interpolation -> Power consumption - Produced H2
            self.PetaEle    = interp1d(Power_inp,etaEle,bounds_error=False,fill_value='extrapolate')      # Linear spline 1-D interpolation -> Power consumption - Electric efficiency
            
        if parameters['stack model'] == 'Alkaline':
            '''
            Alkaline Electorlyzer - McPhy model
            '''
            self.EFF                    = np.zeros(timestep_number)     # keeping track of the elecrolyzer efficiency over the simulation
            self.wat_cons               = np.zeros(timestep_number)     # water consumption array initialization
            self.EFF_last_module        = np.zeros(timestep_number)     # last module efficiency array initialization
            self.wat_cons_last_module   = np.zeros(timestep_number)     # last module water consumption initialization
            self.n_modules_used         = np.zeros(timestep_number)     # array containing modules used at each timestep
            self.cell_currdens          = np.zeros(timestep_number)     # cell current density at every hour
            self.AmbTemp                = c.AMBTEMP                     # [K]         Standard ambient temperature - 15 °C
            if self.ageing:             # if ageing effects are being considered
                self.stack = {
                                'T[°C]': np.zeros(timestep_number),                            # Initialize empty array to track module temperature for each timestep 
                                'Activation[-]': np.zeros(timestep_number),                    # Initialize an array to track module activation (1 for on, 0 for off) for each timestep
                                'Pol_curve_history': [],                                       # Initialize an empty list to keep track of polarization curve shifts during utilization
                                'Module_efficiency[-]': [],                                    # Initialize an empty list to keep track of module efficiency over time
                                'Conversion_factor_op[kg/MWh]': np.zeros(timestep_number),        # Initialize an array to keep track of performance evolution
                                'Conversion_factor_rated[kg/MWh]': np.zeros(timestep_number),        # Initialize an array to keep track of performance evolution
                                'hydrogen_production[kg/s]': np.zeros(timestep_number),        # Initialize an array to keep track of hydrogen production
                                'I_op[A]': np.zeros(timestep_number),                          # Initialize an array to keep track of operating current
                                'V_op[V]': np.zeros(timestep_number)                           # Initialize an array to keep track of operating voltage
                                }
                
                # for i in range(self.n_modules):                         # iterate over the number of modules specified in the parameters
                    # self.stack[f"module_{i+1}"] = {'T[°C]'      : np.zeros(timestep_number),    # for each module, an empty array is created to track parameter behaviour for every timestep of the simulation
                    #                                'Activation' : np.zeros(timestep_number)}  
                
            """
            Anode: 2OH- --> 1/2 O2 + H2O + 2e-

            Cathode: 2H2O + 2e- --> H2 + 2OH-

            Overall reaction: 2H2O(l) + electrical energy --> 2H2(g) + O2(g)
            
            # refs for modelling approach:  1 - https://pdf.sciencedirectassets.com/282073/1-s2.0-S2212017314X00058/1-s2.0-S2212017314001686/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjELz%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJHMEUCIQCLqAxPMnlJl7XlWV23NW2woNNkz2eP9nGUbx1mBdFnRwIgH9hmF%2BcYiOjxMrKrA6Jl12fl0y98n9%2BVWlCbUedvprYqvAUI1f%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FARAFGgwwNTkwMDM1NDY4NjUiDPz1fO4DDWcPwUXNrCqQBc%2FiTeOEkirkQzgE9kWPPP2Qb59lZrWOK%2BW9MZUrZIUuFuR2X48VX44QDFZamkCW05TNpk7KoxjxWC8ybCBZxhWvWybSnwXPd9B%2B2Nbix5jH4HCfTl80uYz6Hjl8OKLSt6QpwearvSn90Jm0d6KekO2o48dM5xoCBqSmeo0xLER3hW6c6gs4Kp8ag0BsG8QvVVc86%2FxYNPawAZTcKVyQacqQe0GM%2BgJBSpMg%2Bc6WLSbds8O7qxzayNvY3NBkoMr9bojZsHIpY3Yhpj5v1ccAONRV8aOoR9uy40hGfG4GNo8oAD%2BZ3Tks%2B0LkPc5fVVyHLQ1z91QzSkwCUgJ31C5s10PkKXJurbfSnPqskUb8ixt%2F1INEfuuYDsfeLRQWmZ9BJNdwX33qYcQ7L6Zr0M3y5DOFNwJrMGUIwNaaduqtSgDZ4CQLVUGZvnTaPCxgAEE88l7Vv3BBQiA9%2B6ErvH3neXcA2fNmiOINcd0stXq%2FMS5zzos%2FH5%2Fhl%2FGyR8%2BqWUpX7CHBb2PfcjKnzsbKM11atR7P7q2xd1Tamf6fu3Bsf%2BYrQU%2FCI9GJhcblyJPIGOt4nwO2mqdB3VFpafQPns1Pg2RPLhMzGRwzIKjx9Mlr0CPNsbrUApCrNI8wyacvBsefyWQAcAc%2F259ssAVC8o9XikBoLEHta2QAkqWlkx5VYWDfA%2FS9y8VsxtJwlb%2F90rjRHfA%2FlPUsyjC%2FFo%2BJ7pMKaq%2FC9DJwJPlgCwoAz23QeaG37RN3n4UgrEaGuoAQ6Kr%2Fr8C7ifyBhkM1BgkB09Z3%2Fjddw0s0dVYLeZBSDqGMjDh3KxHlZDLTaKOLaS%2BvDxKScbDhgPkXaJbrp1elwQZKZdMbEyI3vJ1TSthLS9XjcvcsMMGnlbAGOrEBYmG1DzkwA4ZFTjrC8Cgo8w1gcXXKVPaaBWvWrAa4coycEqDHMDRhzBDYJ%2FJVIwzh%2FB0X6wtvA7UY9Xxqnq2Q7BONGrdUGSRZHGNOKguA86TNog62EyXRyCAlRJnoj3OVyyny2m36kpgtfd2AlbT98LSUC4Z71oH7feAIXIytn0JGJKIkH8VnrpLeFhpSjOeMolkGIv6gEtae%2FgzYi0IudUNjk1nrKA4X%2BSiyaH%2Fw%2BAhc&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20240328T120208Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTY6G6REGFR%2F20240328%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=ded9d214c0be69f013fc02c186d231db3195e81d891e14886df1b4d2e786d285&hash=d927908279169f3621d66403d1b2dba982156ff44938cf009dbfbf53bdf8788a&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S2212017314001686&tid=spdf-a036668c-85e6-48a0-8791-c8a6ef8434d2&sid=692f72be52846740520893a042dac0589099gxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=13125d50055550565753&rr=86b77833ae1083b8&cc=it
                                            2 - https://iopscience.iop.org/article/10.1088/1757-899X/1279/1/012005/pdf
                                            3 - https://www.sciencedirect.com/science/article/pii/S266620272100063X?ref=pdf_download&fr=RR-2&rr=86bf99327e193757
                                            4 - https://www.sciencedirect.com/science/article/pii/S0360319902000332
                                            5 - https://www.sv-jme.eu/?ns_articles_pdf=/ns_articles/files/ojs/858/public/858-6610-1-PB.pdf&id=3045
                                            6 - https://assets.siemens-energy.com/siemens/assets/api/uuid:a33a8c39-b694-4d91-a0b5-4d8c9464e96c/efficiency-white-paper.pdf
            """
            Runiv                       = c.R_UNIVERSAL                 # [J/(mol*K)] Molar ideal gas constant
            self.FaradayConst           = c.FARADAY                     # [C/mol]     Faraday's constant
            self.LHVh2                  = c.LHVH2                       # [MJ/kg]     H2 LHV
            self.HHVh2Mol               = c.HHVH2MOL                    # [kJ/mol]    H2 HHV molar
            
            # Math costants
            self.eNepero    = c.NEPERO          # [-]       Euler's number
            # Ambient conditions 
            self.T_ref      = 298               # [K]       Standard ambient temperature - 15 °C
            
            self.design_T       = 273.15 + 71   # [K]       Design Operating Temperature
            self.design_P       = 3000000       # [Pa]      Design Operating Pressure
            self.design_Pbar    = 30            # [bar]     Design Operating Pressure
            self.KOH    = 30        # [-] potassium hydroxide wt% in the electorlyte
            self.n      = 2         # [-] number of electrons
            self.ΔS0    = 162.86   # [J/(mol*K)] Standard state entropy change at 25° C
            self.ΔH0    = 286.02   # [J/(mol*K)] Standard state enthalpy change at 25° C
            self.ΔG0    = self.ΔH0*1000 - self.T_ref*self.ΔS0   # [J/(mol*K)] Standard state enthalpy change at 25° C
            
            self.Erev       = self.ΔG0/(2*self.FaradayConst)       # [V] minimum reversible voltage required for electrolysis reaction to occur at STP
            self.m          = self.KOH*(183.1221 - 0.56845*self.design_T + 984.5679*np.exp(self.KOH/115.96277))/5610.5          # [mol/kg] molarity of the electrolyte
            self.P0_H2O_bar = np.exp(37.043 - 6275.7 / self.design_T - 3.4159 * np.log(self.design_T))                          # [bar] vapor pressure of pure water
            self.PH2O_bar   = np.exp(0.016214 - 0.13802*self.m + 0.19330*np.sqrt(self.m) + 1.0239*np.log(self.P0_H2O_bar))          # [bar] water vapor pressure over the electrolyte
            self.P0_H2O     = self.P0_H2O_bar*100000    # [Pa] conversion bar ---> Pascal
            self.PH2O       = self.PH2O_bar*100000      # [Pa] conversion bar ---> Pascal
            
            # Coefficients for overvoltage on electrodes 
            self.overvolt_c = {'t1': 0.501,     # [m^2/A]
                               't2': 8.424,     # [m^2°C/A]
                               't3': 247.3,     # [m^2°C/A]
                               's' : 0.024}    # [V]
            
            # Coefficients for ohmic resistance of electrolyte
            self.ohmic_r = {'r1': 7.80e-5,      # [Ωm^2]
                            'r2': -0.98e-7}     # [Ωm^2/°C]
            
            self.CellArea           = 0.5       # [m^2] single cell area 
            self.cell_power         = 10        # [kW] single cell rated power (ref. https://linkinghub.elsevier.com/retrieve/pii/S0960148123003725)
            self.rated_h2prod       = 1.9       # [Nm^3/h] rated hydrogen production of the single cell as reported in ref. https://linkinghub.elsevier.com/retrieve/pii/S0960148123003725
            self.nc                 = round(self.Npower/self.cell_power)    # number of cells present in a single module (max rated power= 1000 kW)  
            self.CurrDensityMax_id  = 1000      # [mA/cm^2] 
            self.CurrDensityMax_Am2 = self.CurrDensityMax_id*10             # [A/m^2]
            self.CurrDensityMax     = self.CurrDensityMax_id/1000           # [A/cm^2]
            self.CellVoltage_limit  = 2.3       # [V] imposed limit for maximum cell voltage over degradation
            
            # Computing the electrolyzer polarization curve V-i
            'POLARIZATION CURVE'
            
            Ndatapoints = 3000      # [-] number of points used to compute the polarization curve 
    
            self.CellCurrDensity        = np.linspace(0,self.CurrDensityMax_Am2,Ndatapoints)    # [A/m^2]
            self.CellCurrDensity_Acm2   = self.CellCurrDensity/1e4                              # [A/cm^2]
            self.CellVoltage            = np.zeros(Ndatapoints)                                 # [V]
        
            'Polarization (V-i) curve calculation'
            
            # V is obtained by summing 3 different contributions 
            '1- Reversible/Open Circuit voltage'
            Vrev = self.Erev + (self.design_T-self.T_ref)*self.ΔS0/(self.n*self.FaradayConst) + \
                    (Runiv*self.design_T)/(2*self.FaradayConst)*(np.log((self.design_Pbar - self.PH2O_bar)**(3/2) / (self.PH2O_bar/self.P0_H2O_bar)))    
            
            for i in range(0,Ndatapoints):
                '2- Cell activation overpotential'
                Vact  = self.overvolt_c['s']*np.log((self.overvolt_c['t1']+(self.overvolt_c['t2']/self.design_T)+(self.overvolt_c['t3']/(self.design_T**(2))))*self.CellCurrDensity[i]+1)   # [V]
            
                '3- Cell Ohmic losses'
                Vohm = (self.ohmic_r['r1'] + self.ohmic_r['r2']*self.design_T)*self.CellCurrDensity[i]
                
                self.CellVoltage[i] = Vrev + Vact +Vohm         # [V] - cell voltage
            
            self.CellVoltage = np.round(self.CellVoltage,4)     # [V] - cell voltage rounded      
            self.Voltage = self.nc*self.CellVoltage             # [V] - module voltage
            self.Current = self.CellCurrDensity*self.CellArea   # [A] - module current
            self.Current = np.round(self.Current)               # [A] - module current rounded
            # Interpolation of calculated functioning points to detect the best fit-function for i-V curve'
            self.num    = Ndatapoints                           # [-] number of intervals to be considered for the interpolation
            self.x2     = np.linspace(0,max(self.CellCurrDensity),self.num)# setting xlim for range of validity of LinRegression Calculation - Only for plot-related reasons 

            # Interpolation
            self.iV1    = interp1d(self.CellCurrDensity,self.Voltage,bounds_error=False)    # Linear spline 1-D interpolation CurrDensity-Voltage                                                                                                                                                                                                                                        
            self.IV     = interp1d(self.Current,self.Voltage,bounds_error=False)            # Linear spline 1-D interpolation Current-Voltage                                                                                                                                                                                                                                                
           
            # Defining Electrolyzer Max Power Consumption
            Power_inp = []                                      # [kW] Initializing power input series
            for i in range(len(self.Current)):
                pot = (self.Current[i]*self.Voltage[i])/1000    # [kW]
                Power_inp.append(pot)
                        
            # Interpolation
            self.PI = interp1d(Power_inp,self.Current,bounds_error=False,fill_value='extrapolate')        # Linear spline 1-D interpolation - Power-Current
            self.PV = interp1d(Power_inp,self.Voltage,bounds_error=False,fill_value='extrapolate')        # Linear spline 1-D interpolation - Power-Voltage
            
            self.MaxPowerModule = round(max(Power_inp),2)                       # [kW] max power input for the considered module
            self.MinInputPower  = self.MaxPowerModule*self.min_load             # [kW] min power input for the considered module based on specified constraints
            self.CurrMin        = round(self.PI(self.MinInputPower).item(),2)   # [A] min module current value based on operational minimum imposed
            self.CurrDensityMin = round(self.CurrMin/self.CellArea,2)           # [A/m^2] min module current density derived from previous calculations            
            self.MaxPowerStack  = self.n_modules*self.MaxPowerModule            # [kW] electrolyzer stack total power

            # Faraday Efficiency Parameters --->>> ref. https://pdf.sciencedirectassets.com/282073/1-s2.0-S2212017314X00058/1-s2.0-S2212017314001686/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEJD%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJIMEYCIQDgqxwVjJxiHLLBODMw03%2FHIYvT4hvP8iWsuOphCOxhygIhAOmdhVlRWzCgYTkqJbbkQlEojL%2FyZ65wc9fJQ8QEVi%2FPKrsFCKn%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEQBRoMMDU5MDAzNTQ2ODY1Igy%2FakbWeXhGmUJ2l%2B4qjwVmU0QmbIFx2ql02k%2FCIRG56%2BJINPVd5yg4O7TW1w6W4OVHbtsFp36ElAow2tb6XVsl3xawjf2kDlCNhflC%2BuA6KZGrxejl2jzxM3f%2FIU9IlTb%2BltwVqoG6s6dNSKtGwexeW%2F5%2BZu%2FzpjNkMgrwoWUH0bPmoK1rXYIrQ8vzePI%2BLL7lPiLgdYZ%2B1jiikIyG74%2F6aewIreubveyqcLLpICb1X435YEXvEXzy62AyZszRsQV8Ev%2B1V84IFG4TClJgGw7s4qHyRbUjoOqRvGwSKJ1lSrjw%2B17%2FZJy0oa5%2BgkEjVYjoLhV82kpQsYmlzn31blKXxy546Phq2c4G5t1M8O4vlxvXybiMx2NSaLpheiCHBX8KHaMmaYAQqBkQfsklcpfp25oMDF5h3iTIvTg0gHbODI1StZOTXA1ySJ1NGiDy9Z%2Brel1VGgVHaWty7vvdFCpZKZk%2Bv%2BjuUmuLZlwkOATPGzt%2FM07QBGpJkjsL40gL3Xgj7khUK6brJ2Ko25GnxEuX3sQr01KKX8ZmsVI%2Bilm7Hz9FwA6nw1uO6RSBxTCjf9r5WHQ6GkcWnF0MJaahGwxWIWC2aI%2FguZS9qckNeCyFJePfXgNa3HTVVeUjNch45trpuBs4hRbtdQKZd9ttwA8s%2BRrDilDBzBOcH01iPNAOrDcwm6H62VhneAei6sz5TUwbIyLiAoNYpqiZSBOYYOkjiO4UrwSto5BkUAT4dRnSt4DFLafF4rvIkWfS%2BFByKHj6SbUsrMzVZ55JUEvhsoF50lEO2eNIhYBr1tRfIbGYOKfBFi%2BesQ1pxqok%2Fa0XUx0lKAX%2BMMeARg3tZBpVxbOTrlh%2BtvvUdyrPyLVBtkEubeYo7SBFJ5rTFX8nWrVrMLTUi7AGOrABgfQNQtzpMPjdQVvcJ9Q95AZsp3PIJhcg5gZq3quesElVagQuqcFOdgyNny4iwvO4ajRmKB7RYlRiEmCxQnPnxML4MQRuo5Cb0h09Ekj5EPuTUc%2Fqpxn1uQAABA7P%2FfiUprSXcMPZSnc7gu2xRNI9u6nQqaE6mD2JqAJI9RT077olvcXjpDBO1yQfmyvwuNcaQWMPjCd9osQ9KFBGG3s5ps3xX92uW0OoghSueZ3iJVM%3D&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20240326T162812Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYXCFU7G7J%2F20240326%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=c6fa76575d41faee4c908f59cb708d7dce934ea53d388fa2ca2a52a1450ea0b3&hash=c52a55ec47589bc504e35e9490ba46e4685d9e596b258c1b380345487e8d6af2&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S2212017314001686&tid=spdf-d132cc46-e88d-42de-a9b3-2401b0568a91&sid=0dbca1a288a0a1443d7a9560727a5bb103aegxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=13105c53025f5e5c5157&rr=86a88330df6fba8c&cc=it
            Temp    = [40,60,80]            # [°C]
            f1      = [150,200,250]         # [mA^2/cm^4]
            f2      = [0.990,0.985,0.980]   # [mA^2/cm^4]

            t_op = self.design_T -273.15    # [°C] operational temperature
            
            self.effe1 = interp1d(Temp,f1)  # interpolating function for coefficient 1
            self.effe2 = interp1d(Temp,f2)  # interpolating function for coefficient 2
            
            f1_= self.effe1(t_op)           # [mA^2/cm^4] coeff 1 based on operating temp
            f2_= self.effe2(t_op)           # [mA^2/cm^4] coeff 2 based on operating temp
            
            'Faraday efficiency'
            self.eta_F = []                 # [-] Faraday efficiency empty list
            # Faraday efficiency computation
            for i in range(len(self.CellCurrDensity)):
                
                eta = ((self.CellCurrDensity[i]/10)**2/(f1_+(self.CellCurrDensity[i]/10)**2))*f2_
                self.eta_F.append(eta)
            
            'Hydrogen produciton'
            # hydrogen production rate - single cell and module
            hydrogen    = []
            for i in range(len(self.Current)):
                hyd = (self.eta_F[i]*self.Current[i])/(2*self.FaradayConst)  # [mol/s] hydrogen production
                hydrogen.append(hyd)
            
            h2_prodcellvol      = [h * self.H2VolMol for h in hydrogen]         # [Nm^3/s] single cell volumetric hydrogen productoin
            h2_prodcellmass     = [h * self.rhoNrh2 for h in h2_prodcellvol]    # [kg/s] single cell hydrogen mass production
            
            h2_prodmodulemass   = [h * self.nc for h in h2_prodcellmass]        # [kg/s] single module hydrogen mass production
            
            'Electrolyzer efficiency'
            # efficiency - single cell and module
            self.eta_cell   = []                 # [-] single cell efficiency empty list
            self.eta_module = []                 # [-] module efficiency empty list
            for i in range(len(self.Current)):
                if self.Current[i]/1000*self.CellVoltage[i] == 0: # denominator check
                    self.eta_cell.append(0)
                    self.eta_module.append(0)
                else: 
                    eff = (self.eta_F[i]*h2_prodcellmass[i]*(self.H2_lhv*1000))/((self.Current[i]/1000)*self.CellVoltage[i])    # [-] electrolytic cell efficiency (Eq. 10) - ref. https://www.sv-jme.eu/?ns_articles_pdf=/ns_articles/files/ojs/858/public/858-6610-1-PB.pdf&id=3045
                    self.eta_cell.append(eff)
                    eff1 = (self.eta_F[i]*h2_prodmodulemass[i]*(self.H2_lhv*1000))/((self.Current[i]/1000)*self.Voltage[i])     # [-] module efficiency
                    self.eta_module.append(eff1)
            
            self.maxh2prod          = round(max(h2_prodmodulemass),8)   # [kg/s] maximum amount of produced hydrogen for the considered module
            self.maxh2prod_stack    = self.maxh2prod*self.n_modules     # [kg/s] maximum amount of produced hydrogen for the considered stack 
            self.Σ                  = (self.maxh2prod*3600)/(self.Npower/1000)   # [kg/MWh] ideal converison factor
            self.eta_module         = np.array(self.eta_module)

            if self.ageing:             # if ageing effects are being considered
                self.stack['Pol_curve_history'].append(self.Voltage) # saving ideal polarization curve as first element to keep track og ageing effects
                self.stack['Module_efficiency[-]'].append(self.eta_module) # saving ideal efficiency curve
                self.stack['T[°C]'][0] = self.design_T -273.15  # [°C] initialising electorlyser temperature
            
            'Functions for predicting the operating behaviour'
            # interpolating functions
            self.etaF       = interp1d(h2_prodmodulemass,self.eta_F,bounds_error=False,fill_value='extrapolate')        # Linear spline 1-D interpolation -> produced H2 - Faraday efficiency
            self.etaEle     = interp1d(h2_prodmodulemass,self.eta_module,bounds_error=False,fill_value='extrapolate')   # Linear spline 1-D interpolation -> produced H2 - Electric efficiency
            self.h2P        = interp1d(h2_prodmodulemass,Power_inp,bounds_error=False,fill_value='extrapolate')         # Linear spline 1-D interpolation -> produced H2 - Power consumption 
            self.P2h        = interp1d(Power_inp,h2_prodmodulemass,bounds_error=False,fill_value='extrapolate')         # Linear spline 1-D interpolation -> Power consumption - Produced H2
            self.PetaEle    = interp1d(Power_inp,self.eta_module,bounds_error=False,fill_value='extrapolate')           # Linear spline 1-D interpolation -> Power consumption - Electric efficiency

            
    def h2power(self,h2):
        """
        Inverse function that computes the power consumption and Faraday efficiency
        corresponding to a required amount of produced hydrogen by means of interpolating functions.

        Parameters
        ----------
        h2 : float hydrogen to be produced [kg/s]

        Returns
        -------
        hyd: float hydrogen produced [kg/s] (same as input 'h2')
        P_absorbed : float electricity absorbed [kW]
        etaElectr : module efficiency [-]
        watCons : module water consumption [m^3/s]
        CellCurrDensity1 : single cell current density [A/cm^2]
        etaFaraday : Faraday efficiency [-]

        """
        if 0 <= h2 <= self.maxh2prod:
            
            hyd        = h2                                                                             # [kg/s]
            P_absorbed = self.h2P(hyd)                                                                  # [kW] required power
            CellCurrDensity1 = self.PI(P_absorbed)/self.CellArea                                        # [A/cm^2] Cell working current density
            etaElectr = self.etaEle(hyd)                                                                # [-] Electrolyzer efficiency   
            # watCons   = hyd_vol*self.rhoStdh2*self.h2oMolMass/self.H2MolMass/etaElectr/self.rhoStdh2o   # [m^3/s] water used by the electrolyzer - volume calculated @ 15°C & Pamb  
            watCons   = hyd*self.watercons                                                              # [m^3/s] water used by the electrolyzer - volume calculated @ 15°C & Pamb  
        
            if P_absorbed >=  self.MinInputPower :   
                
                pass
            
            else:      # If absorbed power is too low, produced hydrogen is zero
                
                hyd               = 0
                P_absorbed        = 0
                etaElectr         = 0
                watCons           = 0
                CellCurrDensity1  = 0
        
        return(hyd,P_absorbed,etaElectr,watCons,CellCurrDensity1)
    
             
    def plot_polarizationpts(self):
        
        if self.model in ['PEM General','Alkaline']: 
            # i-V plot
            x = np.linspace(min(self.CellCurrDensity),max(self.CellCurrDensity),self.num) 
            plt.figure(dpi=300)
            plt.plot(self.CellCurrDensity,self.Voltage,label='data', color='b',marker='.', linestyle='None', mec='r', markersize=7, markerfacecolor='white', zorder=0)
            plt.axvline(x=self.CurrDensityMin,color='b',linestyle=':', label= 'Lower Functioning Boundary', zorder=3, linewidth=1.2)   # Ref. doi: 10.1016/j.ijhydene.2008.11.083
            plt.plot(x,self.iV1(x),label='linear', linestyle='--') 
            plt.grid()
            plt.legend(fontsize=8)
            if self.model == 'Alkline':
                plt.xlabel('Cell Current Density [A/m$^{2}$]')
            else:
                plt.xlabel('Cell Current Density [A/cm$^{2}$]')
            plt.ylabel('Stak Voltage [V]')
            plt.title(f'{self.model}EL Polarization Curve - SplineInterp' )
            plt.show()
            
            # i-V plot
            x = np.linspace(min(self.CellCurrDensity),max(self.CellCurrDensity),self.num) 
            plt.figure(dpi=300, figsize = (5,3.5))
            # plt.plot(self.CellCurrDensity,self.CellVoltage,label='data', color='b',marker='.', linestyle='None', mec='r', markersize=7, markerfacecolor='white', zorder=0)
            plt.plot(x,self.iV1(x)/self.nc, color = 'steelblue') 
            # plt.plot(x,self.iV1(x),label='Model', color = 'steelblue') 
            plt.axvline(x=self.CurrDensityMin,color='k',linestyle=':', label= 'CurDens$_\mathregular{min}$', zorder=3, linewidth=1.2)   # Ref. doi: 10.1016/j.ijhydene.2008.11.083
            plt.grid(alpha = 0.5)
            plt.legend(fontsize=8)
            if self.model == 'Alkline':
                plt.xlabel('Cell Current Density [A/m$^{2}$]')
            else:
                plt.xlabel('Cell Current Density [A/cm$^{2}$]')
            plt.ylabel('Cell Voltage [V]')
            # plt.title('PEMEL Polarization Curve - SplineInterp' )
            plt.show()
        
        else: 
            print('Polarization curve not available')
        
                 
    def plot_linregression(self): 
        
        if self.model == 'PEM General':
    
            'Linear Regression'
            x1 = self.CellCurrDensity.reshape((-1,1))
            y1 = self.Voltage
            
            model = LinearRegression().fit(x1,y1)
            r_sq_linreg = model.score(x1,y1)
            print('Coeff. of determination:', r_sq_linreg)
            
            self.coeff_A = model.intercept_      # Obtaining the calculated intercept for the linear fit - returns a scalar
            print('intercept:', self.coeff_A)
            self.coeff_B = model.coef_           # Obtaining the calculated slope for the linear fit - returns an array with only 1 value
            print('slope:', self.coeff_B)
            
            Volt_LinReg = self.coeff_A + self.coeff_B*self.CellCurrDensity
               
            # i-V plot Linear Regression 
            plt.figure(dpi=1000)
            plt.plot(self.CellCurrDensity,self.Voltage,label='data', color='b',marker='.', linestyle='None', mec='r', markersize=7, markerfacecolor='white', zorder=0)
            plt.axvline(x=0.05,color='b',linestyle=':', label= 'Lower Functioning Boundary', zorder=3, linewidth=1.2)   # Ref. doi: 10.1016/j.ijhydene.2008.11.083
            plt.plot(self.x2,Volt_LinReg,label='Least Squares LinReg') 
            plt.grid()
            plt.legend(fontsize=8)
            plt.xlabel('Cell Current Density [A cm$^{-2}$]')
            plt.ylabel('Stak Voltage [V]')
            plt.title('PEMEL Polarization Curve - LinReg' )
            plt.show()  

        else: 
            print('Polarization curve not available')            
    
    
    def ageing(self,step,power,Text):
        """
        Calculates the impact of ageing on the electrolyzer, adjusting its performance over time.
        Ageing effects are modeled as increases in operational voltage due to both time and temperature,
        which in turn affect hydrogen production efficiency.
    
        Parameters:
        - step (int): Current simulation step indicating the operational time.
        - power (float): Electrical power input to the electrolyzer [kW].
        - Text (float | None): External temperature [°C]; defaults to a standard value if None.
    
        Returns:
        - hyd_produced (float): Adjusted hydrogen production [kg/s], accounting for ageing effects.
    
        This function updates the electrolyzer's polarization curve to reflect degradation, then uses
        this updated curve to determine the new operational parameters, including the hydrogen production rate.
        """
        
        V_inctime   = (3*1e-6)/60    # [V/min]   voltage increase - time degradation in μV
        V_incTemp   = 5*1e-3         # [5mV/°C]  voltage increase - thermal degradation
        
        'Ideal behaviour'
        Iop_id    = self.PI(power)      # [A] module operating ideal current based on system power input
        Vop_id    = self.PV(power)      # [V] module operating ideal voltage based on system power input 
        H2op_id   = self.P2h(power)     # [kg/s] module oprating ideal hydrogen production based on system power input 
        
        if step == 0:
            hydprod_prev = 0                        # [kg/s] effect on thermal degradation of hydrogen produced at step 0 are neglected
            Iop_prev     = 0                        # [A]                  
            Vop_prev     = 0                        # [V]
            Tel_prev     = self.design_T-273.15     # [°C]            
        else:
            hydprod_prev = self.stack['hydrogen_production[kg/s]'][step-1] # [kg/s] hydrogen prod. Effects of production at step-1 manifesting on polarization curve at current step
            Iop_prev     = self.stack['I_op[A]'][step-1]    # [A]                  
            Vop_prev     = self.stack['V_op[V]'][step-1]    # [V]
            Tel_prev     = self.stack['T[°C]'][step-1]      # [°C]
            
        temp = el.thermal_effects(Tel_prev,hydprod_prev,Iop_prev,Vop_prev,Text)     # [°C] module temperature at current timestep
        
        'Computing ageing phenomena' 
        V_time          = (V_inctime*self.timestep)*self.nc     # [V] voltge time degradation tha may occur for the considered step in the simulation if the electorlyzer is turned on
        
        # operating time counter
        operation_time  = sum(self.stack['Activation[-]'])
        
        V_thermal       = (V_incTemp*((self.design_T-273.15)-temp))*self.nc    # [V] voltage thermal degradation
        
        # updating polarization curve
        polarization_curve_new = self.Voltage + V_time*operation_time + V_thermal # [V] self.Voltage represents the design polarization curve
        
        # limit on time degradation for single cell voltage
        if max(polarization_curve_new-V_thermal)/self.nc > self.CellVoltage_limit:
            print('Electorlyzer module voltage exceeds safe limits due to ageing. Module must be replaced')
        
        #link between cell current (= stack current) and cell voltage: polarization curve
        IV_new   = interp1d(self.Current,polarization_curve_new) # Linear spline 1-D interpolation - updating I-V function for ageing effect
        V_op     = IV_new(Iop_id)       # [V] operational voltage accounting for ageing effect 
        
        ageing_factor_op    = Vop_id/V_op     # [-] ageing factor expressed as the ratio between operational and ideal voltage for the considered current. Denominator increases over time
        ageing_factor_rated = max(self.Voltage)/max(polarization_curve_new) # [-] ageing factor for functioning at rated power
        hyd_produced        = H2op_id*ageing_factor_op # [kg/s] hydrogen produced in operative conditions accounting for ageing effect 

        if hyd_produced > 0: # if the elctrolyzer has been activated at current step
            self.stack['Activation[-]'][step] = 1
        
        self.stack['Conversion_factor_op[kg/MWh]'][step]        = self.Σ*ageing_factor_op       # [kg/MWh]  ideal converison factor
        self.stack['Conversion_factor_rated[kg/MWh]'][step]     = self.Σ*ageing_factor_rated    # [kg/MWh]  ideal converison factor
        self.stack['hydrogen_production[kg/s]'][step]           = hyd_produced                  # [kg/s]    hydrogen produced in the timestep
        self.stack['I_op[A]'][step]                             = Iop_id                        # [A]       operating current
        self.stack['V_op[V]'][step]                             = V_op                          # [V]       operating voltage
        if step == 0:
            pass
        else:
            self.stack['T[°C]'][step]                           = temp                          # [°C]       operating voltage
        
        if step % self.timesteps_year == 0:
            self.stack['Pol_curve_history'].append(polarization_curve_new)
            self.stack['Module_efficiency[-]'].append(self.eta_module*(self.Voltage/polarization_curve_new))         # [kg/MWh] ideal converison factor
            print(f'Year {int(step/self.timesteps_year)}')
        
        return hyd_produced
    
    def thermal_effects (self,T_el,hydrogen,Iop,Vop,Text):
        """
        Calculates the updated temperature of the electrolyzer module considering thermal effects.
          
        This function models the thermal behavior of the electrolyzer by accounting for both heat generation
        from the exothermic reaction (considering thermal losses) and heat loss to the environment. It uses a simplified
        thermal model that considers the geometry of the electrolyzer, material properties, and operational conditions
        to estimate the heat transfer and subsequent temperature change.
          
        Parameters:
        - T_el [°C]: Current temperature of the electrolyzer.
        - hydrogen [kg/s]: Amount of hydrogen being produced, which influences the thermal power generated within the electrolyzer.
        - Iop [A]: Operating current, used to calculate the thermal power generated from electrical losses above the thermoneutral voltage.
        - Vop [V]: Operating voltage, used in conjunction with the operating current to calculate the thermal power generated.
        - Text [°C] (optional): External ambient temperature. If not provided, a default value is used.
          
        Returns:
        - Updated temperature [°C] of the electrolyzer after accounting for the thermal effects during the operation.
          
        Note: The function scales the electrolyzer's geometry based on the number of cells and uses material properties
        and heat transfer coefficients to estimate heat losses. It also considers the impact of insulation and adjusts
        the temperature based on the heat generated and lost during the operation.
        """
        
        # 1 MW module - design dimensions
        design_lenght   = 3                 # [m] design length of the gas-liquid separator
        design_radius   = 0.3               # [m] internal radius
        ncells_des      = 100               # [-] number of cells for 1 MW module
        
        scale_factor    = ncells_des/self.nc  # [-] scale factor of the configuration. As of current version, it can only been smaller of 1 MW module
        scaled_lenght   = design_lenght*(scale_factor**(1/3))       # [m] scale of the geometry accoring to the SF
        
        π = math.pi                     # [-]
        op_time = self.timestep*60      # [s]  simulation timestep in seconds
        T_op    = self.design_T-273.15  # [°C] design operating temperature
        if Text:
            T_ext = Text                # [°C] external temperature
        else:
            T_ext = 25                  # [°C] external temperature 
        V_tn = 1.48*self.nc             # [V]  thermoneutral voltage https://www.scopus.com/record/display.uri?eid=2-s2.0-77958033092&origin=inward
        
        'module geometry'
        r1 = design_radius*(scale_factor**(1/3)) # [m] module radius       
        s1 = 0.004         # [m] thickness of the electrolyzer container
        r3 = 1             # [m]  container internal radius
        s2 = 0.005         # [m] container thickness
        
        'heat transfeer coefficeints'
        h1 = 100          # [W/ m^2K]   internal convection between water (H2O + 30% KOH) - tank
        h2 = 10           # [W/ m^2K]   convection tank-container
        h3 = 20           # [W/ m^2K]   external convection container-air
        k1 = 52           # [W/ mK]     steel tank conduction 
        k2 = 52           # [W/ mK]     steel container conduction     
        
        'container insulation'
        s2 = 0.2        # [m]     insulation layer thickness
        k2 = 0.05       # [W/ mK] insulation layer conduction
            
        'electrolyte properties'
        m_elect = scaled_lenght*r1*r1*π*1000/2  # [kg] of H2O in gas-liquid separator (half water, half gas)
        c_elect = self.cp_water                 # [J/kg*K]   water specific heat

        a = h1*2*π*r1*scaled_lenght
        b = k1*2*π*scaled_lenght/np.log((r1 + s1)/r1)   
        c = h2*2*π*(r1 + s1)*scaled_lenght
        d = h2*2*π*r3*scaled_lenght
        e = k2*2*π*scaled_lenght/np.log((r3 + s2)/r3)
        f = h3*2*π*(r3 + s2)*scaled_lenght
        
        Qloss = (T_el-T_ext)/(1/a + 1/b + 1/c + 1/d + 1/e + 1/f) # [W] thermal power lost to the environment

        if hydrogen > 0:  
            #thermal power generated form the stack
            Qgain   = (Vop-V_tn)*(Iop)                                  # [V]*[A] = [W] produce thermal power
            temp    = T_el+(op_time/(m_elect*c_elect))*(Qgain-Qloss)    # [°C]
            temp    = min(temp,T_op) # if temp > design temperature cooling system takes it down to design temperature
        else:
            temp = T_el - (op_time/(m_elect*c_elect))*Qloss             # [°C]
        return temp 
    
    
    # Computing Electrolyzers performances via Spline Interpolation  
    def use(self,step,storable_hydrogen=False,p=False,hydrog=False,Text=None):
        """
        Electorlyzers stack and single modules operational parameters
    
        step : timestep float timestep in hours [step]
        p : float > 0 electricity provided to the electrolyzer [kW]
        storable_hydrogen : float storable hydrogen H tank max_capacity - SOC[h-1] [kg] or maximum absorbable production if there is no tank

        output : 
            
        hyd: float hydrogen produced [kg/s]    
        P_absorbed: float electricity absorbed [kW]
        
        """
        if self.strategy == 'hydrogen-first' and hydrog == False:       # defined strategy is either to work only with renewable energy or to prioritize its consumption while interacting also with the electricity grid
            
            if self.model not in ['PEM General','Alkaline']:
                
                P_absorbed = min(self.MaxPowerStack,p)                                      # [kW]
                hyd = self.production_rate*(P_absorbed/self.MaxPowerStack)                  # [kg/s] (P_absorbed / self.Npower) represents the fraction of the maximum possible amount of produced hydrogen
                oxygen = hyd*self.oxy                                                       # [kg/s] Oxygen produced as electorlysis by-product
                # watCons = hyd*self.h2oMolMass/self.H2MolMass/self.eff/self.rhoStdh2o      # [m^3/s] Water consumption for electrolysis process
                watCons   = hyd*self.watercons                                              # [m^3/s] Water consumption for electrolysis process - volume calculated @ 15°C & Pamb
                self.EFF[step] = self.eff
                self.EFF_last_module[step] = self.eff
                if P_absorbed == 0:
                    self.n_modules_used[step] = 0
                else:
                    self.n_modules_used[step] = round(self.MaxPowerStack/P_absorbed)
                
                max_hyd_storable = storable_hydrogen/(self.timestep*60)                     # [kg/s] Maximum storable hydrogen flow rate considering the chosen timestep 
                
                if hyd > max_hyd_storable:     # if there is not enough space in the H tank to store the hydrogen (H tank is nearly full)
                    hyd         = 0                                # turn off the electrolyzer
                    P_absorbed  = 0
                    oxygen      = 0
                    watCons     = 0
                    self.n_modules_used[step] = 0       
                    
                return(hyd,-P_absorbed,oxygen,-watCons) # return hydrogen supplied, electricity absorbed, oxygen produced, water consumed
                        
                    
            elif self.model in ['PEM General','Alkaline']:   
                
                'Defining the working point of the electrolyzer by spline interpolation:'
                
                max_hyd_storable = storable_hydrogen/(self.timestep*60)                           # Maximum storable hydrogen flow rate considering the chosen timestep [kg/s]
                
                # PowerInput [kW] - Electric Power directed to the electrolyzer
                if p <= self.Npower:                      # if available power is lower than single module nominal power
                    
                    P_absorbed = p
                    hyd,P_absorbed,etaElectr,watCons,CellCurrden,hydrogen = electrolyzer.use1(self,step,p,max_hyd_storable,Text=None)
                    oxygen = hyd*self.oxy                 # [kg/s] Oxygen produced as electorlysis by-product 
                    if hyd > 0:
                        self.n_modules_used[step] = 1
                    else:
                        self.n_modules_used[step] = 0
                    self.EFF[step]           = etaElectr     # [-] single module efficiency
                    self.cell_currdens[step] = CellCurrden
                    self.wat_cons[step]      = watCons
                    
                if p > self.Npower:      # if available power is higher than nominal one, i.e., more modules can be used
      
                    n_modules_used = min(self.n_modules,int(p/self.Npower))
                    P_absorbed     = self.Npower                 # power absorbed by the single module          
                    hyd,P_absorbed,etaElectr,watCons,CellCurrden,hydrogen = electrolyzer.use1(self,step,P_absorbed,max_hyd_storable,Text) 
                    oxygen = hyd*self.oxy                        # [kg/s] Oxygen produced as electorlysis by-product 
                    if hydrogen == hyd:                          # if produced hydrogen is equal to the maximum producible one, i.e., if the produced hydrogen is lower than the storable one. Otherwise it means that only one module can be used
                        hyd_11 = np.zeros(n_modules_used+1)      # creating the array where index represents nr of modules and value is the produced hydrogen summing all the modules production
                        for i in range (n_modules_used+1):
                            hyd_11[i] = hydrogen*i                  # Saving the producible hydrogen for each number of n_modules_used
                            if hyd_11[i] > max_hyd_storable and hyd_11[i-1] < max_hyd_storable: # if, using i modules, the total amount of producible hydrogen is higher than storable one 
                                n_modules_used = i-1                
                                hyd_1 = hyd_11[n_modules_used]                 # Hydrogen produced using i-1 modules
                                P_absorbed_1 = P_absorbed*(n_modules_used)     # Total power absorbed using i-1 modules 
                                watCons_1 = watCons*(n_modules_used) 
                                self.EFF[step] = etaElectr                        # work efficiency of modules working at nominal power 
                                self.cell_currdens[step] = CellCurrden
                                
                                hyd_remained = max_hyd_storable-hyd_1
                                hyd,P_absorbed,etaElectr,watCons,CellCurrDensity1 = electrolyzer.h2power(self,hyd_remained)
                                if hyd > 0:
                                    n_modules_used = n_modules_used+1       # considering the module working at partial load
                                    self.EFF_last_module[step] = etaElectr     # work efficiency of the last module working with the remaining power
                                    self.wat_cons_last_module[step] = watCons  # water consumption // // // // //
                                hyd = hyd+hyd_1
                                oxygen = hyd*self.oxy                       # [kg/s] Oxygen produced as electorlysis by-product 
                                P_absorbed = P_absorbed_1+ P_absorbed       # abs value
                                watCons = watCons_1+watCons
                                self.wat_cons[step] = watCons
                                self.n_modules_used[step] = n_modules_used
        
                        if hyd_11[-1] <= max_hyd_storable:                  # if, using n_modules, the total amount of producible hydrogen is lower than storable one  
                            hyd_1 = hyd*n_modules_used                       # total amount of H2 produced by modules working at full load
                            P_absorbed_1 = P_absorbed*n_modules_used         # total power absorbed    // // // // // 
                            watCons_1 = watCons*n_modules_used                 # total water consumption // // // // // 
                            self.EFF[step] = etaElectr                          # work efficiency of modules working at nominal power 
                            self.cell_currdens[step] = CellCurrden
    
                            if p <= self.MaxPowerStack:                      # if available power is lower than total installed power
                                P_remained = p-self.Npower*n_modules_used    # remaining power after considering modules at full load
                                remained_storable_hydrogen = max_hyd_storable-hyd_1 #[kg/s]
                                hyd,P_absorbed,etaElectr,watCons,CellCurrden,hydrogen = electrolyzer.use1(self,step,P_remained,remained_storable_hydrogen,Text=None)
                                if P_absorbed > 0:                           # if remaining power is higher than 10% of Npower, last module working at partial load is activated
                                    n_modules_used = n_modules_used+1        # considering the module working at partial load
                                    self.EFF_last_module[step] = etaElectr      # work efficiency of the last module working with the remaining power
                                    self.wat_cons_last_module[step] = watCons   # water consumption // // // // // 
                                self.n_modules_used[step] = n_modules_used
                                watCons = watCons_1+watCons
                                self.wat_cons[step] = watCons
                                hyd = hyd+hyd_1
                                oxygen = hyd*self.oxy                       # [kg/s] Oxygen produced as electorlysis by-product 
                                P_absorbed = P_absorbed_1+P_absorbed
                                    
                            else:                                   # if available power is higher than total installed power it means that there are no more exploitable modules, thus no more hydrogen can be produced 
                                P_absorbed = P_absorbed_1
                                hyd = hyd_1
                                watCons = watCons_1
                                oxygen = hyd*self.oxy                       # [kg/s] Oxygen produced as electorlysis by-product 
                                self.wat_cons[step] = watCons
                                self.n_modules_used[step] = self.n_modules
                 
                    else:                                       # Equal or less than one modules used
                        if hyd > 0:
                            self.n_modules_used[step] = 1
                        else:
                            self.n_modules_used[step] = 0
                        self.EFF[step] = etaElectr
                        self.cell_currdens[step] = CellCurrden
                        self.wat_cons[step]      = watCons
            
            # elif self.model == 'Alkaline':
        
        
        elif self.strategy == 'hydrogen-first' and p == False:                   # if defined strategy is to work with electricity from renewables and from grid
            if self.model != 'PEM General':
                u=3  #TO BE DONE
                
            elif self.model in ['PEM General','Alkaline']:
                if hydrog <= self.maxh2prod:      # if requested hydrogen is lower than the maximum amount producible by a single module
                    hyd,P_absorbed,etaElectr,watCons,CellCurrden = electrolyzer.h2power(self,hydrog)
                    oxygen = hyd*self.oxy                 # [kg/s] Oxygen produced as electorlysis by-product 
                    if hyd > 0:
                        self.n_modules_used[step] = 1
                    else:
                        self.n_modules_used[step] = 0
                    self.EFF[step]           = etaElectr     # [-] single module efficiency
                    self.cell_currdens[step] = CellCurrden
                    self.wat_cons[step]      = watCons
                                        
                elif self.maxh2prod < hydrog <= self.maxh2prod_stack:      # if requested hydrogen is higher than the maximum one generable by a single module, i.e., more modules can be used      
                    hyd,P_absorbed,etaElectr,watCons,CellCurrden = electrolyzer.h2power(self,self.maxh2prod) # hydrogen to be produced by the single module  
                    hyd_11 = np.zeros(self.n_modules+1)      # creating the array where index represents nr of modules and value is the produced hydrogen summing all the modules production
                    for i in range (self.n_modules+1):
                        hyd_11[i] = hyd*i                    # Saving the producible hydrogen for each number of n_modules_used
                        if hyd_11[i] > hydrog and hyd_11[i-1] < hydrog: # if, using i modules, the total amount of producible hydrogen is higher than the target one 
                            n_modules_used = i-1                
                            hyd_1 = hyd_11[n_modules_used]                 # Hydrogen produced using i-1 modules
                            P_absorbed_1 = P_absorbed*(n_modules_used)     # Total power absorbed using i-1 modules 
                            watCons_1 = watCons*(n_modules_used) 
                            self.EFF[step] = etaElectr                        # work efficiency of modules working at nominal power 
                            self.cell_currdens[step] = CellCurrden
                            
                            hyd_remained = hydrog-hyd_1
                            hyd,P_absorbed,etaElectr,watCons,CellCurrden = electrolyzer.h2power(self,hyd_remained) # hydrogen left to be produced by the last module
                            if hyd > 0:
                                n_modules_used = n_modules_used+1       # considering the module working at partial load
                                self.EFF_last_module[step] = etaElectr     # work efficiency of the last module working with the remaining power
                                self.wat_cons_last_module[step] = watCons  # water consumption // // // // //
                            hyd = hyd+hyd_1
                            oxygen = hyd*self.oxy                       # [kg/s] Oxygen produced as electorlysis by-product 
                            P_absorbed = P_absorbed_1+ P_absorbed       # abs value
                            watCons = watCons_1+watCons
                            self.wat_cons[step] = watCons
                            self.n_modules_used[step] = n_modules_used
                            break
                        
                elif hydrog >= self.maxh2prod_stack:         # if, using n_modules, the total amount of producible hydrogen is lower than the target one  
                    hyd1,P_absorbed1,etaElectr,watCons1,CellCurrden1 = electrolyzer.h2power(self,self.maxh2prod) # hydrogen to be produced by the single module  
                    hyd = hyd1*self.n_modules                      # total amount of H2 produced by modules working at full load
                    P_absorbed = P_absorbed1*self.n_modules        # total power absorbed    // // // // // 
                    watCons = watCons1*self.n_modules                # total water consumption // // // // // 
                    oxygen = hyd*self.oxy                           # [kg/s] Oxygen produced as electorlysis by-product   
                    self.n_modules_used[step] = self.n_modules
                    self.EFF[step] = etaElectr                         # work efficiency of modules working at nominal power 
                    self.cell_currdens[step] = CellCurrden1
                    self.wat_cons[step] = watCons
                   
                    
        elif self.strategy == 'full-time':   # if electrolyzers are supposed to work full-time for their operation
            if self.model != 'PEM General':
                P_absorbed = self.Npower                                 # [kW]
                hyd = self.production_rate*(P_absorbed / self.Npower)    # [kg/s] (P_absorbed / self.Npower) represents the fraction of the maximum possible 
                oxygen = hyd*self.oxy                                    # [kg/s] Oxygen produced as electorlysis by-product 
                # watCons   = hyd*self.h2oMolMass/self.H2MolMass/etaElectr/self.rhoStdh2o   # [m^3/s] water used by the electrolyzer - volume calculated @ 15°C & Pamb  
                watCons   = hyd*self.watercons                           # [m^3/s] water used by the electrolyzer - volume calculated @ 15°C & Pamb                                        
                self.n_modules_used[step] = self.n_modules
                self.EFF[step] = self.eff
                self.EFF_last_module[step] = self.eff
            
            elif self.model == 'PEM General':
                max_hyd_storable = storable_hydrogen/(self.timestep*60)                           # Maximum storable hydrogen flow rate considering the chosen timestep [kg/s]
                P_absorbed = self.Npower
                hyd,P_absorbed,etaElectr,watCons,CellCurrden,hydrogen = electrolyzer.use1(self,step,P_absorbed,max_hyd_storable,Text=None)
                P_absorbed = P_absorbed*self.n_modules
                hyd = hyd*self.n_modules
                watCons = watCons*self.n_modules
                oxygen = hyd*self.oxy                                    # [kg/s] Oxygen produced as electorlysis by-product                                                        
                self.n_modules_used[step] = self.n_modules
                self.EFF[step]           = etaElectr     # [-] single module efficiency
                self.cell_currdens[step] = CellCurrden
                self.wat_cons[step]      = watCons
                
     
        return (hyd,-P_absorbed,oxygen,-watCons)
        
               
    def use1(self,step,p,max_hyd_storable,Text):
        """
        This function calculates the performances of a single electrolyzer module.
        
        Parameters
        ----------
        p : float > 0 electricity provided to the electrolyzer in the timestep [kW]
        max_hyd_storable : float maximum hydrogen flow rate (H tank max_capacity - SOC[h-1])/(timestep*60) [kg/s] or maximum absorbable production if there is no tank

        Returns
        -------
        hyd: float hydrogen produced [kg/s]    
        P_absorbed: float electricity absorbed in the timestep [kW]
        etaElectr : module efficiency [-]
        watCons : module water consumption [m^3/s]
        CellCurrDensity1 : single cell current density [A/cm^2]

        """
        
        P_absorbed        = p                                   # [kW]
        CellCurrDensity1  = self.PI(P_absorbed)/self.CellArea   # [A/cm^2]

        # Checking if resulting current density is high enough for the electrolyzer to start, otherwise hydrogen prod = 0
        if CellCurrDensity1 < self.CurrDensityMin or P_absorbed < self.MinInputPower:      # condition for operability set for current density and input power 
            etaElectr        = 0      # [-] electrolyzer efficiency
            hyd              = 0      # [kg/s] hydrogen produced in the considered timestep
            P_absorbed       = 0      # [kW] absorbed power
            watCons          = 0      # [m^3/s] water volumetric consumption
            CellCurrDensity1 = 0
            hydrogen         = hyd
            if self.ageing:
                hyd         = electrolyzer.ageing(self,step,P_absorbed,Text)
                hydrogen         = hyd                                         # [kg/s] hydrogen produced       

        else:
            if self.ageing:
                'Electrolyzer efficiency'                      
                etaElectr   = self.PetaEle(P_absorbed)
                'Hydrogen Production'
                hyd         = electrolyzer.ageing(self,step,P_absorbed,Text)                                         # [kg/s] hydrogen produced       
                hydrogen    = hyd
                'Water consumption' 
                # watCons = hyd_vol*self.rhoStdh2*self.h2oMolMass/self.H2MolMass/etaElectr/self.rhoStdh2o   # [m^3/s] water used by the electrolyzer - volume calculated @ 15°C & Pamb
                                                                                                            # Useful references for water consumption in electrolysis
                                                                                                            # 1) - 15 l of H2O per kg of H2. https://doi.org/10.1016/j.rset.2021.100005
                                                                                                            # 2) - 9 kg of H2O per kg of H2 minimum/18-24 kg of H2O per kg of H2 considering de-mineralization/ up to 32. https://energypost.eu/hydrogen-production-in-2050-how-much-water-will-74ej-need/
                watCons  = hyd*self.watercons                                                               # [m^3/s] water used by the electrolyzer - volume calculated @ 15°C & Pamb                   
            else:
                'Electrolyzer efficiency'                      
                etaElectr   = self.PetaEle(P_absorbed)
                'Hydrogen Production'
                hyd         = self.P2h(P_absorbed)                                         # [kg/s] hydrogen produced       
                hydrogen    = hyd
                'Water consumption' 
                # watCons = hyd_vol*self.rhoStdh2*self.h2oMolMass/self.H2MolMass/etaElectr/self.rhoStdh2o   # [m^3/s] water used by the electrolyzer - volume calculated @ 15°C & Pamb
                                                                                                            # Useful references for water consumption in electrolysis
                                                                                                            # 1) - 15 l of H2O per kg of H2. https://doi.org/10.1016/j.rset.2021.100005
                                                                                                            # 2) - 9 kg of H2O per kg of H2 minimum/18-24 kg of H2O per kg of H2 considering de-mineralization/ up to 32. https://energypost.eu/hydrogen-production-in-2050-how-much-water-will-74ej-need/
                watCons  = hyd*self.watercons                                                               # [m^3/s] water used by the electrolyzer - volume calculated @ 15°C & Pamb                                        
        
        if hyd > max_hyd_storable:           # if there is not enough space in the H tank to store the hydrogen (H tank is nearly full)
        
            hyd,P_absorbed,etaElectr,watCons,CellCurrDensity1 = electrolyzer.h2power(self,max_hyd_storable) # defining the power requirements to produce storable hydrogen amount
       
        return (hyd,P_absorbed,etaElectr,watCons,CellCurrDensity1,hydrogen)           


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

        size = self.Npower * self.n_modules 
        
        if tech_cost['cost per unit'] == 'default price correlation': # ref: https://doi.org/10.1016/j.ijhydene.2023.04.100 PEMEL cost, 2030 year of reference
            C0 = 1185.69 # €/kW
            scale_factor = 0.925 # 0:1
            C =  C0*(size**scale_factor)
        else:
            C = size*tech_cost['cost per unit']
            
        if tech_cost['OeM'] == 'default price correlation': # ref: https://doi.org/10.1016/j.ijhydene.2023.04.100, 2030 year of reference
            C0 = 349.8 # €/kW
            scale_factor = -0.305
            OeM = (C0*(size**scale_factor))*size
        else:
            OeM = tech_cost['OeM'] *C /100 # €

        tech_cost['total cost'] = tech_cost.pop('cost per unit')
        tech_cost['total cost'] = C
        tech_cost['OeM'] = OeM

        self.cost = tech_cost    

##########################################################################################

if __name__ == "__main__":
    
    """
    Functional test
    """
    inp_test = {  
                  "Npower": 1000,
                  "number of modules": 1,
                  "stack model": 'Alkaline',
                  "minimum_load": 0.2,
                  "ageing": False,
                  "strategy": 'hydrogen-first',
                  "only_renewables":True
                }
    
    
    sim_steps = 1000           # [step] simulated period of time - usually it's 1 year minimum
    timestep = 60              # [min] selected timestep
    storable_hydrogen = 10000                      # [kg] Available capacity in tank for H2 storage at timestep 'step'

    el = electrolyzer(inp_test,sim_steps,timestep=timestep)        # creating electrolyzer object
    el.plot_polarizationpts()                    # plot example
    
    'Test 1 - Tailored ascending power input'

    flow  = np.linspace(0,el.Npower*6,sim_steps)  # [kW] power input - ascending series
    flow1 = np.linspace(0,el.Npower,sim_steps)  # [kW] power input - ascending series

    hyd         = np.zeros(len(flow))           # [kg] produced hydrogen
    p_absorbed  = np.zeros(len(flow))           # [kW] absorbed power
    oxygen      = np.zeros(len(flow))           # [kg] produced oxygen
    water       = np.zeros(len(flow))           # [m3] consumed water
    
    for i in range(len(flow)):
        hyd[i],p_absorbed[i],oxygen[i],water[i] = el.use(i,storable_hydrogen=storable_hydrogen,p=flow[i])

    fig, ax = plt.subplots(dpi=600)
    ax.plot(-p_absorbed,el.n_modules_used,color='tab:green',zorder=3)
    ax.set_xlabel('Absorbed Power [kW]')
    ax.set_ylabel('Active electrolyzer modules [-]')
    ax.grid()
    ax.set_title('Electrolyzer Stack - Absorbed Power')
    
    fig, ax = plt.subplots(dpi=600)
    ax.plot(flow,el.n_modules_used,color='indianred',zorder=3)
    ax.set_xlabel('Input Power [kW]')
    ax.set_ylabel('Active electrolyzer modules [-]')
    ax.grid()
    ax.set_title('Electrolyzer Stack - Input Power')
        
    # print('Produced Hydrogen in the timestep [kg]: \n\n',hyd)
    # print('\nAbsorbed power [kWh]: \n\n',p_absorbed)     
    # print('\nElectrolyzer Efficiency [-]: \n\n',el.EFF)                # electrolyzer efficiency at every hour     

    cellarea = el.CellArea
    if el.model == 'Alkaline':
        cellarea = el.CellArea*10000
    nompower = el.Npower
    
    plt.figure(dpi=600)
    plt.scatter(el.cell_currdens,el.EFF,s=20,color='tab:orange',edgecolors='k',zorder=3)
    plt.title("Efficiency vs Power Input")
    # plt.ylim([0,0.8]) 
    textstr = '\n'.join((
        r'$CellArea=%.1f$ $cm^{2}$' % (cellarea,),
        r'$P_{nom}= %.1f$ kW' % (nompower,),
        r'$i_{max}= %.1f$ A $cm^{-2}$' % (el.CurrDensityMax,),
        r'$n_{cell}= %.0f$' % (el.nc,)))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    plt.text(max(el.cell_currdens)/2,0.2,textstr,fontsize=10,va='bottom',backgroundcolor='none', bbox=props)
    plt.grid()
    plt.xlabel('Curr dens [A/cm2]')
    plt.ylabel('$\\eta$')
    
    for i in range(len(flow1)):
        
        hyd[i],p_absorbed[i],oxygen[i],water[i] = el.use(i,storable_hydrogen=storable_hydrogen,p=flow1[i])
        
    plt.figure(dpi=600)
    plt.plot(flow1,el.EFF)
    plt.title("Electrolyzer Module Efficiency")
    # plt.ylim([0,0.8]) 
    textstr = '\n'.join((
        r'$CellArea=%.1f$ $cm^{2}$' % (cellarea,),
        r'$P_{nom}= %.1f$ kW' % (nompower,),
        r'$i_{max}= %.1f$ A $cm^{-2}$' % (el.CurrDensityMax,),
        r'$n_{cell}= %.0f$' % (el.nc,)))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    plt.text(el.Npower/2,0.2,textstr,fontsize=10,va='bottom',backgroundcolor='none', bbox=props)
    plt.grid()
    plt.xlabel('Input Power [kW]')
    plt.ylabel('$\\eta$')    

    min_load = el.min_load
    first_component = int(min_load*sim_steps)
    fig, ax = plt.subplots(dpi =300, figsize = (5,3.5))
    ax2 = ax.twinx()
    ax.plot(flow1[first_component:-1],el.EFF[first_component:-1],label='Efficiency', color = '#eb4034')
    ax2.plot(flow1,hyd,label='H2$_\mathregular{prod}$', color ='#4ea3f2')
    plt.axvline(x=min_load*inp_test['Npower'],color='black',linestyle='--',zorder = 3)
    ax.set_xlabel('Power input [kW]')
    ax.set_ylabel('$\\eta$ [-]')
    ax.set_ylim(0,None)
    ax2.set_ylim(0,None)    
    ax2.set_ylabel('Produced hydrogen [kg/s]')
    ax.grid(alpha = 0.5)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax2.legend(h1+h2, l1+l2, loc='best', fontsize = 'small')
        
    'Test 2 - Random power input'

    flow = np.random.uniform(0.08*el.Npower,5.2*el.Npower,sim_steps)   # [kWh] randomic power input as example
    
    hyd          = np.zeros(len(flow))          # [kg] produced hydrogen
    p_absorbed    = np.zeros(len(flow))          # [kW] absorbed hydrogen
    oxygen       = np.zeros(len(flow))          # [kg] produced oxygen
    water        = np.zeros(len(flow))          # [m3] consumed water
    
    for i in range(len(flow)):
        
        hyd[i],p_absorbed[i],oxygen[i],water[i] = el.use(i,storable_hydrogen=storable_hydrogen,p=flow[i])
        
    fig, ax = plt.subplots(dpi=1000)
    ax2 = ax.twinx() 
    ax.bar(np.arange(sim_steps)-0.2,el.EFF,width=0.35,zorder=3,edgecolor='k',label='$1^{st}$ module efficiency', alpha =0.8)
    ax.bar(np.arange(sim_steps)+0.,el.EFF_last_module,width=0.35,zorder=3, edgecolor = 'k',align='edge',label='Last module efficiency',alpha =0.8)
    ax2.scatter(np.arange(sim_steps),flow,color ='limegreen',s=25,edgecolors='k',label='Available Power')
    ax.set_ylim(None,0.8)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1+h2, l1+l2, loc='lower center',bbox_to_anchor=(0.5, 1.08), ncol =3, fontsize ='small')
    ax.set_xlabel('Time [step]')
    ax.set_ylabel('Efficiency [-]')
    ax2.set_ylabel('Power Input [kW]')
    ax.grid()
    ax.set_title('Electrolyzer Stack functioning behaviour')


#%%
    """
    Evaluation of ageing effects - test
    """
    if el.model == 'Alkaline':
        
        inp_test["ageing"] = True
        
        sim_steps = 8760*15              # [step] simulated period of time - usually it's 1 year minimum
        timestep = 60                   # [min] selected timestep
        
        el = electrolyzer(inp_test,sim_steps,timestep=timestep)        # creating electrolyzer object
        
        # Temperature array to test temperature variation effect 
        timestep_day        = (60*24)/timestep     # [-] number of timesteps in one day
        daily_fluctuations  = np.sin(np.linspace(0, 2*np.pi, int(timestep_day)))
        temp_min            = -5 # [°C]
        temp_max            = 30 # [°C]
        daily_fluctuations  = ((daily_fluctuations+1)/2)*(temp_max-temp_min)+temp_min
        ext_temp            = np.tile(daily_fluctuations,int(sim_steps/timestep_day))
        random_variation    = np.random.uniform(-2.5,2.5,size=ext_temp.shape)
        ext_tempv           = ext_temp + random_variation # [°C] ambient temperature array to be fed to the electorlyzer model
        
        flow = np.random.uniform(0.08*el.Npower,1.2*el.Npower,sim_steps)   # [kW] randomic power input as example

        hyd          = np.zeros(len(flow))          # [kg] produced hydrogen
        p_absorbed    = np.zeros(len(flow))          # [kW] absorbed hydrogen
        oxygen       = np.zeros(len(flow))          # [kg] produced oxygen
        water        = np.zeros(len(flow))          # [m3] consumed water
        
        for i in range(len(flow)):
            hyd[i],p_absorbed[i],oxygen[i],water[i] = el.use(i,storable_hydrogen=storable_hydrogen,p=flow[i],Text=ext_tempv[i])
             
#%%
        'Polarization Curve History'
        fig, ax = plt.subplots(dpi=1000) 
        
        for index,item in enumerate(el.stack['Pol_curve_history']):
            ax.plot(el.CellCurrDensity_Acm2,item,label=f'{index}')
        
        ax.legend(title='Year',fontsize='small')
        ax.set_xlabel('Current density [A/cm$^{2}$]')
        ax.set_ylabel('Stack Voltage [V]')
        # ax.set_ylim(135.10,135.3)
        ax.grid(alpha=0.5,zorder=-1)    
        plt.show()
        
        'Efficiency History'
        fig, ax = plt.subplots(dpi=1000) 
        
        for index,item in enumerate(el.stack['Module_efficiency[-]']):
            ax.plot(el.CellCurrDensity_Acm2,item,label=f'{index}')
        
        ax.legend(title='Year',fontsize='small')
        ax.set_xlabel('Current density [A/cm$^{2}$]')
        ax.set_ylabel('Module efficiency [-]')
        ax.grid(alpha=0.5,zorder=-1)    
        plt.show()
        
        'Conversion Factor Update'
        fig, ax = plt.subplots(dpi=1000)     
        ax.scatter(np.arange(sim_steps),el.stack['Conversion_factor_op[kg/MWh]'], s=0.1, label='operational')
        ax.scatter(np.arange(sim_steps),el.stack['Conversion_factor_rated[kg/MWh]'], s=0.1, label='rated')
        # ax.plot(np.arange(sim_steps),el.stack['Conversion_factor_rated[kg/MWh]'], linewidth=0.1, label='rated')
        # ax.plot(np.arange(sim_steps),el.stack['Conversion_factor_op[kg/MWh]'], linewidth=0.05, label='operational')
        ax.legend()
        ax.set_ylabel('Σ [kg/MWh]')
        ax.set_xlabel('Time [step]')
        ax.grid(alpha=0.5,zorder=-10)    
        plt.show()
