import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
import math
from numpy import log as ln
from sklearn.linear_model import LinearRegression
import os
import sys 
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),os.path.pardir)))   # temporarily adding constants module path 
import constants as c
import pickle
import scipy.fft
import scipy.optimize

class fuel_cell:
    
    def __init__(self,parameters,simulation_hours):
        """
        Create a Fuel Cell object
    
        parameters : dictionary
            'Npower': float nominal power [kW]
                      
        output : Fuel cell object able to:
            absosrb hydrogen and produce electricity .use(e)
        """
        
        self.model               = parameters['stack model']  # [-]  Fuel cell model
        self.Npower              = parameters['Npower']       # [kW] Single module nominal power   
        self.ageing              = parameters['ageing']       # bool - calucate ageing   
        self.timestep = 1                                     # [h]  simulation timestep
        self.ageing_day = 7     # [days] How often ageing has to bee calculated? 
        self.check_ageing=0     # ageing calculation check. If already computed or not. 
        self.pol_curves = []    # list initialization for corrected curves. Aging effects affecting pol curve behaviour
        
        #########################################
        if self.model == 'FCS-C5000':           # this model is based on FCS-C5000 characteristic curves https://www.horizonfuelcell.com/hseries
            self.nc = 120 * self.Npower/5       # number of cells (120 cells are required to have a nominal power of 5000 W)
        
            # characteristic curves
            I=[0, 10, 20, 30, 40, 50, 60, 70, 80]                        # [A]
            v=[0.96, 0.83, 0.796, 0.762, 0.728, 0.694, 0.66, 0.6, 0.5]   # [V]
            V=[]                                                         # [V] Stack voltage list 
            P=[]                                                         # [W] Stack output power list
            for i in range(len(I)):
                volt=self.nc*v[i]
                V.append(volt)
                pot=I[i]*volt/1000                                       # [W] -> [kW] conversion
                P.append(pot)                                            # [kW] Stack power
            
            self.IP=interp1d(I,P,kind='cubic',bounds_error=False,fill_value='extrapolate')  
            self.Pmax=max(P)
           
            # PI inverse curve that will be used to find the operating point
            I=[]
            P=[]
            for a in range(800):
                i=a/10
                pot=self.IP(i)
                if pot>self.Pmax:
                    break  
                I.append(i)
                P.append(pot)
            self.PI=interp1d(P,I,kind='cubic',bounds_error=False,fill_value='extrapolate')

        ###########################################
        if self.model == 'PEM General':
            
            # Model validation with experimental data for the specifications below:
                
            'Model: NEDSTACK FCS 13-XXL'      # https://nedstack.com/sites/default/files/2022-07/nedstack-fcs-13-xxl-gen-2.9-datasheet-rev01.pdf
            'Rated nominal power:  13.6 kW'
            
            self.EFF       = np.zeros(simulation_hours)              # [-]  Keeping track of fuel cell efficiency
            self.VOLT      = np.zeros(simulation_hours)              # [V]  Keeping track of single cell working voltage - necessary for ageing calculations
            self.CURR_DENS = np.zeros(simulation_hours)              # [A]  Keeping track of single cell working current - necessary for ageing calculations
            self.EFF_last_module     = np.zeros(simulation_hours)    # [-]  Keeping track of the elecrolyzer last module efficiency over the simulation
            self.n_modules_used      = np.zeros(simulation_hours)    # [-]  Number of modules active at each timestep 

            
            self.rhoStdh2        = c.H2SDENSITY           # [kg/Sm3]    PropsSI('D', 'T', 288.15, 'P', 101325, 'H2') H2  density @ T = 15°C p = 101325 Pa
            self.rhoStdh2o       = c.H2OSDENSITY          # [kg/m3]     H2O density @ T = 15°C p = 101325 Pa
            self.Runiv           = c.R_UNIVERSAL          # [J/(mol*K)]
            self.Rh2             = c.R_H2                 # [J/(kg*K)] 
            self.FaradayConst    = c.FARADAY              # [C/mol]     Faraday constant
            self.deltaG0         = c.GIBBS                # [kJ/mol]    Gibbs free energy @ T = 25°C p = 101325 Pa
            self.GammaPerfectGas = c.GAMMA                # [-]         Gamma = cp/cv  
            self.LHVh2           = c.LHVH2                # [MJ/kg]     H2 LHV
            self.HHVh2           = c.HHVH2                # [MJ/kg]     H2 HHV
            self.HHVh2Mol        = c.HHVH2MOL             # [kJ/mol]    H2 HHV molar
            self.cpH2O           = c.CP_WATER             # [kJ/(kgK)]  Water specific heat
            self.h2oMolMass      = c.H2OMOLMASS           # [kg/mol]    Water molar mass
            self.H2MolMass       = c.H2MOLMASS            # [kg/mol]    Hydrogen molar mass

            # Math costants
            self.eNepero      = c.NEPERO                  # [-]         Euler's number

            # Ambient conditions 
            self.AmbTemp      = c.AMBTEMP                 # [K]         Standard ambient temperature - 15 °C
     
            "H2 --> 2H+ + 2e" 
          
            # FuelCell  Parameters 
            
            self.Lambda              = 23                                  # [-]      Cell mositure content
            self.FC_AnodeCurrDens    = 0.000000009                         # [A/cm^2] Anode current density
            self.FC_CathodeCurrDens  = 0.001                               # [A/cm^2] Cathode current density
            self.FC_OperatingTemp    = 273.15 + 62                         # [K]      https://www.google.com/search?q=working+temperature+PEM+fuel+cell&oq=working+temperature+PEM+fuel+cell+&aqs=chrome..69i57j0i512l9.10383j0j7&sourceid=chrome&ie=UTF-8
            self.FC_FuelPress        = 101325 + 25000                      # [Pa]     Fuel supply pressure (Anode) 
            self.FC_AirPress         = 101325                              # [Pa]     Air supply pressure (Cathode)
            self.FC_MinCurrDens      = 0.0001                              # [A/cm^2] FC min. current density      
            self.CTC                 = 0.4                                 # [-]      Charge Transfer Coefficient  - if Nominal Power < 6 kW: CTC = 0.45
                                                                           #                                  - if //   //   //  > 6 kW: CTC = 0.4
            self.MembThickness       = 250                                 # [μm]     Fuel cell membrane thickness - if Nominal Power < 6 kW: MembThickness = 100 
                                                                           #                                    - if //   //   //  > 6 kW: MembThickness = 145 
            self.FC_MaxCurrent       = 230                                 # [A]      Value taken from datasheet. Current value at which max power is delivered
            
            self.n_modules           = parameters['number of modules']     # [-]      Number of modules constituting the stack
            self.MaxPowerTot         = self.n_modules*self.Npower          # [kW]     Stack maximum power output
         
            self.nc                  = 90 + int((self.Npower/1000)*(250-90))   # For a power range between 0kW and 1000kW the number of cells in the stack varies between 90 and 250 
            self.FC_MaxCurrDens      = 1.2 + (self.Npower/1000)*(1.3-1.2)      # For a power range between 0kW and 1000kW the maximum current density varies between 1.2 and 1.3 A/cm2 
            self.FC_NominalPower     = self.Npower                 # [kW] Max Power - Nominal Power at max current
            
            # Varying the number of cells, module efficiency remains unchanged
            # nc = 96 and FC_MaxCurrDens = 1.2 are derived from the above-mentioned datasheet and are specific for the specified size of 13.6 kW 
            
            
            'POLARIZATION CURVE'
            
            Ndatapoints= 1000                    # Number of points used to compute the polarization curve 
            
            # Saving different losses contruibutions to be considered in the polarization curve
            
            self.OCpotential = []                # [V] Open circuit voltage
            self.ActLosses   = []                # [V] Activation losses
            self.OhmLosses   = []                # [V] Ohmic losses
            self.CellVolt = np.zeros(Ndatapoints)
            self.CellCurrDensity = np.linspace(self.FC_MinCurrDens,self.FC_MaxCurrDens,Ndatapoints)   
            
            'Polarization (V-i) curve calculation'
            
            # V is obtained by summing 3 different contributions (concentration losses are considered to be negligible)

            '1- Cell open curcuit voltage'
            
            pO2 = (self.FC_AirPress*0.21)/101325  # [atm] accounting for partial pressure of oxigen in air
            pH2O = 1                              # [atm]
            pH2 = self.FC_FuelPress/101325        # [atm]
            
            Ecell = 1.229 -0.85e-3*(self.FC_OperatingTemp-298.15) + 4.3085e-5*self.FC_OperatingTemp*ln(pH2*(pO2**0.5)/pH2O)  # [V] Open circuit voltage
            self.OCpotential = [Ecell for i in range(Ndatapoints)]
                            
            for i in range(0,Ndatapoints):
                                                        
                '2- Activation losses'
              
                Vact_cat = -self.Runiv*self.FC_OperatingTemp*np.log10(self.FC_CathodeCurrDens)/(self.CTC*4*self.FaradayConst)+self.Runiv*self.FC_OperatingTemp*np.log10(self.CellCurrDensity[i])/(self.CTC*4*self.FaradayConst) #[V]
                Vact_an = -self.Runiv*self.FC_OperatingTemp*np.log10(self.FC_AnodeCurrDens)/(self.CTC*2*self.FaradayConst)+self.Runiv*self.FC_OperatingTemp*np.log10(self.CellCurrDensity[i])/(self.CTC*2.*self.FaradayConst)   #[V]
          
                Vact = Vact_cat +Vact_an               # [V] Activation losses
                
                self.ActLosses.append(Ecell-Vact)
           
                '3- Ohmic losses'
         
                rho_m = (181.6*(1+0.03*(self.CellCurrDensity[i])+0.062*((self.FC_OperatingTemp/303)**2)*(self.CellCurrDensity[i])**2.5))/ \
                    ((self.Lambda-0.634-3*(self.CellCurrDensity[i]))*self.eNepero**(4.18*(self.FC_OperatingTemp-303)/self.FC_OperatingTemp))      #[Ohm*cm] specific membrane recistence
              
                Rm = rho_m*self.MembThickness/10000    # [Ohm/cm2]  Cell resistance depending on temperature and moisture content (Lambda) 
         
                Vohm = self.CellCurrDensity[i]*Rm      # [V] Ohmic losses
                
                self.OhmLosses.append(Ecell-Vact-Vohm)
                
                'SINGLE CELL VOLTAGE'
                
                self.CellVolt[i] = (Ecell-Vact-Vohm)                                  # [V]    Cell Voltage
            
            self.Vmin_FC=self.CellVolt[-1]*self.nc                                    # [V]    Minimum value for working voltage
            self.FC_CellArea=self.Npower*1000/(self.Vmin_FC*self.FC_MaxCurrDens)      # [cm^2] FC cell active area

            'RESULTING STACK VOLTAGE'
         
            self.Voltage = self.CellVolt*self.nc                                      # [V] Stack Voltage

            'Interpolation of  polarization curve: defining the fit-function for i-V curve'
            
            self.num = Ndatapoints                                                    # [-] number of intervals to be considered for the interpolation
            self.x = np.linspace(self.CellCurrDensity[0],self.CellCurrDensity[-1],self.num) 
            
            # Interpolating functions
            
            self.iV1 = interp1d(self.CellCurrDensity,self.Voltage,bounds_error=False,fill_value='extrapolate')    # Linear spline 1-D interpolation
            self.iV2 = interp1d(self.CellCurrDensity,self.CellVolt,bounds_error=False,fill_value='extrapolate')   # Linear spline 1-D interpolation
            
            # Creating the reverse curve IP - necessary to define the exact functioning point        
            
            self.Current = self.CellCurrDensity*self.FC_CellArea    # [A] Defining the current value: same both for the single cell and the full stack!
            
            # Defining Fuel Cell Max Power Generation
            
            'Fuel Cell Max Power Consumption'
            
            FC_power = []
            for i in range(len(self.Current)):
                pot=self.Current[i]*self.Voltage[i]/1000       # [kW] Power
                FC_power.append(pot)
            self.MinPower = min(FC_power)
            
            self.IP = interp1d(self.Current,FC_power,bounds_error=False,fill_value='extrapolate')          
            self.P  = []
            for i in range(len(self.Current)):
                power = (self.iV1(self.CellCurrDensity[i])*self.Current[i])/1000     # [kW] Resolving the equation system via interpolation
                self.P.append(power)                                                 # [kW] Output power values varying current
                
            self.PI=interp1d(self.P,self.Current,bounds_error=False,fill_value='extrapolate')  # Interpolating function returning Current if interrogated with Power 
            
            'Single module electricity production'
            
            # Building lists of values needed for interpolation functions
            
            hydrogen             = []
            electricity_produced = []
            eta_FuelCell         = []
            FC_Heat_produced     = []
            
            for i in range(len(FC_power)):
               
                p_required = FC_power[i]                                   # [kWh] --> equivalent to kW for the considered timestep
                
                FC_CellCurrDensity = self.PI(p_required)/self.FC_CellArea  # [A/cm^2] current density value at which the fuel cell is working 
        
                Current = FC_CellCurrDensity*self.FC_CellArea     # [A] FuelCell Stack operating current 
                
                FC_Vstack = self.iV1(FC_CellCurrDensity)          # [V] Stack operating voltage
                V_cell = FC_Vstack/self.nc                        # [V] Single cell operating voltage
        
                'Computing FC efficiency and hydrogen energy demand'    
                
                pO2 = (self.FC_AirPress*0.21)/101325              # [atm]
                pH2O = 1                                          # [atm]
                pH2 = self.FC_FuelPress/101325                    # [atm]  
                
                Ecell = 1.229-0.85e-3*(self.FC_OperatingTemp-298.15) + 4.3085e-5*self.FC_OperatingTemp*ln(pH2*(pO2**0.5)/pH2O)   # [V] Open circuit voltage 
                deltaG = self.deltaG0 - self.Runiv*self.FC_OperatingTemp*ln(pH2*math.sqrt(pO2)/pH2O)/(2*self.FaradayConst)       # [kJ/mol] Gibbs free energy at actual conditions
              
                eta_voltage = FC_Vstack/(Ecell*self.nc)           # [-] Voltage efficiency 
                eta_th = - deltaG/self.HHVh2Mol                   # [-] Thermodynamic efficiency
             
                etaFC = eta_th*eta_voltage                        # [-] FC efficiency

                'Hydrogen demand'
             
                FC_HydroCons = Current*self.nc*3600/(95719.25*1000)/self.rhoStdh2*self.timestep   # [Sm3] (Chavan 2017)
                hyd = FC_HydroCons*self.rhoStdh2/self.timestep                                    # [kg/h]
                FC_deltaHydrogen = - FC_HydroCons*self.rhoStdh2*self.HHVh2*1000/3600              # [kWh]
             
             
                'Process heat, that can be recovered'
              
                FC_Heat = ((1.481*self.nc)/FC_Vstack-1)*p_required*self.timestep     # [kWh] --> equivalent to kW for the considered timestep of 1h
                    
                hydrogen.append(hyd)                       # [kg]  produced hydrogen
                electricity_produced.append(p_required)    # [kWh] output power
                eta_FuelCell.append(etaFC)                 # [-]   fc efficiency
                FC_Heat_produced.append(FC_Heat)           # [kWh] co-product heat
                
            self.maxh2used = max(hydrogen)                 # [kg] maximum amount of exploitable hydrogen for the considered module
          
            self.etaFuelCell = interp1d(hydrogen,eta_FuelCell,bounds_error=False,fill_value='extrapolate')          # Linear spline 1-D interpolation -> H2 consumption - FC efficiency
            self.h2P         = interp1d(hydrogen,electricity_produced,bounds_error=False,fill_value='extrapolate')  # Linear spline 1-D interpolation -> H2 consumption - produced electricity
            self.FC_Heat     = interp1d(hydrogen,FC_Heat_produced,bounds_error=False,fill_value='extrapolate')      # Linear spline 1-D interpolation -> H2 consumption - produced heat

        ####################################   
        if self.model == 'SOFC':
            
            self.EFF=np.zeros(simulation_hours)                    # [-]      keeping track of fuel cell efficiency
            
            self.EFF_last_module     = np.zeros(simulation_hours)  # [-]      keeping track of the elecrolyzer last active module efficiency over the simulation
            self.n_modules_used      = np.zeros(simulation_hours)  # [-]      Number of modules active at each timestep 
            
            self.e                            = c.NEPERO           # [-]      Euler's number

            self.FC_OperatingTemp             = 273.15+800         # [K]      Operating temperature
            self.FC_RefTemp                   = 273.15+750         # [K]      Operating temperature reference
            self.FC_FuelPress                 = 116000             # [Pa]     Fuel supply pressure (Anode)  
            self.FC_AirPress                  = 101325/101325      # [atm]    Air supply pressure (Cathode)
            # self.FC_MinCurrDens               = 0.0001              # [A/cm^2] FC min. current density - arbitrary
            self.FC_MinCurrDens               = 0.04194            # [A/cm^2] FC min. current density - derived from experimental plot

            self.FC_NominalPower              = self.Npower        # [W]      FC nominal power
            self.FC_AlfaAnode                 = 0.55               # [-]      Charge transfer coefficient for anode
            self.FC_AlfaCathode               = 0.33               # [-]      Charge transfer coefficient for cathode 
            self.FC_ActivationEnergyAnode     = 110                # [kJ/mol] Anonde activation energy  -> ref.:http://dx.doi.org/10.1016/j.desal.2017.02.013
            self.FC_ActivationEnergyCathode   = 160                # [kJ/mol] Cathode activation energy -> ref.: //        //        //        //        //
            self.FC_ActivationCoeff           = 10                 # [-]      Activation coefficient 
            self.FC_ExchangeCurrDensChAn      = 0.530              # [A/cm^2] Exchange current density channel anode 
            self.FC_ExchangeCurrDensChCat     = 0.200              # [A/cm^2] Exchange current density channel cathode 
            self.FC_ExchangeCurrDensCathode   = self.FC_ExchangeCurrDensChCat*self.e**(((self.FC_ActivationCoeff*self.FC_ActivationEnergyCathode)\
                                                /c.R_UNIVERSAL)*((1/self.FC_RefTemp)-(1/self.FC_OperatingTemp)))                                   # [mA/cm^2]
            self.FC_ModFactor                 = 2                  # [-]      Ohmic losses modification factor 
            self.ThicknessElectrolyte         = 0.00003            # [m]      Electrolyte thickness
            self.ElectrolyteCostant           = 50                 # [K/ohm*m]
            self.ActivationEnergyElectrolyte  = 9*10**(7)          # [kJ/mol] Electrolyte activation energy
            self.FC_LimitCurrDens             = 6                  # [A/cm^2]
            self.DiffusionVolumeH20           = 13.1               # [-] 
            self.DiffusionVolumeH2            = 6.1                # [-]
            self.HHVh2Mol                     = c.HHVH2MOL         # [kJ/mol] Hydrogen Higher Heating Value - Molar
            self.rhoStdh2                     = c.H2SDENSITY       # [kg/Sm3] Hydrogen density at Standard conditions 
            self.HHVh2                        = c.HHVH2            # [MJ/kg]  Hydrogen  higher heating value      
            
            h2oMolMass       = c.H2OMOLMASS   # [kg/mol]     Water molar mass
            H2MolMass        = c.H2MOLMASS    # [kg/mol]     Hydrogen molar mass
            H2MolStdEntropy  = c.H2MOL_S_E    # [J/K*mol]    Specific molar entropy (gaseous phase)
            O2MolStdEntropy  = c.O2MOL_S_E    # [J/K*mol]    Specific molar entropy
            H20MolStdEntropy = c.H2OMOL_S_E   # [J/K*mol]    Specific molar entropy
            
            self.n_modules           = parameters['number of modules']  # [-]      Number of modules constituting the stack
            self.MaxPowerTot         = self.n_modules*self.Npower       # [kW]     Max power output 
            
            self.nc                  = 90 + int((self.Npower/1000)*(250-90))       # For a power range between 0kW and 1000kW the number of cells in the stack varies between 90 and 250 
            self.FC_MaxCurrDens      = 1.2 + (self.Npower/1000)*(1.3-1.2)          # For a power range between 0kW and 1000kW the maximum current density varies between 1.2 and 1.3 A/cm2 
            
            # Varying the number of cells, module efficiency remains unchanged
            # nc = 77 and FC_MaxCurrDens = 1.23457 are derived from the above-mentioned datasheet
            
            'POLARIZATION CURVE'
            
            Ndatapoints = 1000                        # Number of points used to compute the polarization curve 
            
            self.DeltaV_ohm = np.zeros(Ndatapoints)   # Ohmic losses
            self.DeltaV_con = np.zeros(Ndatapoints)   # Concentration losses
            self.Ecell = []                           # Open circuit voltage
            self.CellVoltage = np.zeros(Ndatapoints)
            self.CellCurrDensity = np.linspace(self.FC_MinCurrDens,self.FC_MaxCurrDens,Ndatapoints)  
            
            pH2 = self.FC_FuelPress/101325                  # [atm]
            pH2O = 1                                        # [atm]
            pO2 = (self.FC_AirPress*0.21*101325)/101325     # [atm]
            
            # V is obtained by summing 4 different contributions 
            
            'Polarization (V-i) curve calculation'
        
            '1 - Open circuit potential'                
            
            Ecell = 1.19 + ((c.R_UNIVERSAL*self.FC_OperatingTemp)/(2*c.FARADAY))*ln(pH2*(pO2**0.5)/pH2O)  # [V] Open circuit voltage
            self.Ecell  = [Ecell for i in range(Ndatapoints)]

            for i in range(0,Ndatapoints):
                                
                '2 - Ohmic losses'
                
                self.DeltaV_ohm[i] = (self.CellCurrDensity[i]*self.ThicknessElectrolyte*self.FC_OperatingTemp)/\
                                     (9000*self.e**(-100000/(c.R_UNIVERSAL*self.FC_OperatingTemp)))                   # [V] (100000 activation energy in [kJ/mol], 9000 electrolyte constant)
               
                '3 - Concentration losses'   # Significant losses only for cathode, activation losses (minimal contribution for SOFC) are also present in the following formula
                
                self.DeltaV_con[i] = ((c.R_UNIVERSAL*self.FC_OperatingTemp)/(2*c.FARADAY))*\
                                     ln((self.CellCurrDensity[i]/(self.FC_ExchangeCurrDensCathode*self.FC_AirPress*(0.21-0.0008*10000*self.CellCurrDensity[i]*c.R_UNIVERSAL*self.FC_OperatingTemp/(4*c.FARADAY*101325*0.00002)))))        # [V] 
                
                '4 - Cell voltage'
                
                self.CellVoltage[i] = self.Ecell[i]-self.DeltaV_ohm[i]-self.DeltaV_con[i]
                
            self.Vmin_FC_stack = self.CellVoltage[-1]*self.nc                                 # [V]    Minimum value for working voltage
            self.FC_CellArea = self.Npower*1000/(self.Vmin_FC_stack*self.FC_MaxCurrDens)      # [cm^2] FC cell active area
            self.FC_Pmax = self.Npower                                                        # [kW] Max output power

            '5- Stack voltage'
            
            self.Voltage = self.nc*self.CellVoltage

            'Interpolation of  polarization curve: defining the fit-function for i-V curve'    
          
            self.num = Ndatapoints                                          # [-] number of intervals to be considered for the interpolation
            self.x = np.linspace(self.CellCurrDensity[0],self.CellCurrDensity[-1],self.num)
                    
            # Defining different interpolation methods
            self.iV1 = interp1d(self.CellCurrDensity,self.Voltage,bounds_error=False,fill_value='extrapolate')          # Linear spline 1-D interpolation
            
            # Creating the reverse curve IP - necessary to define the exact functioning point
            self.Current = self.CellCurrDensity*self.FC_CellArea

            'Fuel Cell Max Power Consumption'
            
            FC_power = []
            for i in range(len(self.Current)):
                 pot = self.Current[i]*self.Voltage[i]/1000     # [kW] Power
                 FC_power.append(pot)
            self.MinPower = min(FC_power)

            
            self.IP=interp1d(self.Current,FC_power,kind='cubic',bounds_error=False,fill_value='extrapolate') 
            self.P = np.zeros(Ndatapoints)
            
            for i in range (len(self.Current)):
                  self.P[i] = (self.iV1(self.CellCurrDensity[i])*self.Current[i])/1000   # [kW] power output 
            
            self.PI=interp1d(self.P,self.Current,bounds_error=False,fill_value='extrapolate')   # Interpolating function returning Current if interrogated with Power 

            'Single module electricity production'
   
            # Building lists of values needed for interpolation functions    
    
            hydrogen = []
            electricity_produced = []
            eta_FuelCell   = []
            FC_Heat_produced = []
            
            for i in range(len(FC_power)):
                
                p_required = FC_power[i]                                   # [kWh] --> equivalent to kW for the considered timestep                
                FC_CellCurrDensity = self.PI(p_required)/self.FC_CellArea  # [A/cm^2] current density value at which the fuel cell is working 
                
                Current = FC_CellCurrDensity*self.FC_CellArea     # [A] FuelCell Stack operating current
                
                FC_Vstack= self.iV1(FC_CellCurrDensity)           # [V] Stack operating voltage
    
                'Computing FC efficiency and hydrogen energy demand'    
                
                pO2 = (self.FC_AirPress*0.21*101325)/101325      # [atm]
                pH2O = 1                                         # [atm]
                pH2 = self.FC_FuelPress/101325                   # [atm]  
                
                Ecell = 1.19 + ((c.R_UNIVERSAL*self.FC_OperatingTemp)/(2*c.FARADAY))*ln(pH2*(pO2**0.5)/pH2O)
                DeltaV_OHM = (FC_CellCurrDensity*self.ThicknessElectrolyte*self.FC_OperatingTemp)/(9000*self.e**(-100000/(c.R_UNIVERSAL*self.FC_OperatingTemp)))
                DeltaV_CON = ((c.R_UNIVERSAL*self.FC_OperatingTemp)/(2*c.FARADAY))*ln((FC_CellCurrDensity/(self.FC_ExchangeCurrDensCathode*self.FC_AirPress*(0.21-0.0008*10000*FC_CellCurrDensity*c.R_UNIVERSAL*self.FC_OperatingTemp/(4*c.FARADAY*101325*0.00002)))))
    
                TotalLoss = DeltaV_OHM + DeltaV_CON   # [V]
                
                deltaG = c.GIBBS - c.R_UNIVERSAL*self.FC_OperatingTemp*ln(pH2*math.sqrt(pO2)/pH2O)/(2*c.FARADAY)    # [kJ/mol] Gibbs free energy at actual conditions
              
                eta_voltage = FC_Vstack/(Ecell*self.nc)           # [-] Voltage efficiency 
                eta_th = c.GIBBS/self.HHVh2Mol                    # [-] Thermodynamic efficiency
                etaFC = -eta_th*eta_voltage                       # [-] FC efficiency
                
                'Hydrogen demand'
                
                FC_HydroCons = ((Current*self.nc*3600)/(c.FARADAY*1000))/(self.rhoStdh2*self.timestep)     # [Sm^3]
                hyd=FC_HydroCons*self.rhoStdh2/self.timestep                                               # [kg/h]
                FC_deltaHydrogen = - FC_HydroCons*self.rhoStdh2*self.HHVh2*1000/3600                       # [kWh]
                
                z = Current/(2*c.FARADAY)   # [mol/s]
                
                DeltaS = -((H20MolStdEntropy-(O2MolStdEntropy/2)-H2MolStdEntropy)+(c.R_UNIVERSAL/2)*ln((pH2**2)*pO2/(pH2O**2))) # [J/mol*K]
                
                FC_Heat = (((z*self.FC_OperatingTemp*DeltaS + Current*TotalLoss)*self.timestep)*self.nc)/1000   # [kWh] --> equivalent to kW for the considered timestep of 1h

                hydrogen.append(hyd)                       # [kg]  produced hydrogen
                electricity_produced.append(p_required)    # [kWh] output power
                eta_FuelCell.append(etaFC)                 # [-]   fc efficiency
                FC_Heat_produced.append(FC_Heat)           # [kWh] co-product heat
                
            self.maxh2used = max(hydrogen)                                                          # [kg] maximum amount of exploitable hydrogen for the considered module
          
            self.etaFuelCell = interp1d(hydrogen,eta_FuelCell,bounds_error=False,fill_value='extrapolate')       # Linear spline 1-D interpolation -> H2 consumption - FC efficiency
            self.h2P  = interp1d(hydrogen,electricity_produced,bounds_error=False,fill_value='extrapolate')      # Linear spline 1-D interpolation -> H2 consumption - produced electricity
            self.FC_Heat = interp1d(hydrogen,FC_Heat_produced,bounds_error=False,fill_value='extrapolate')       # Linear spline 1-D interpolation -> H2 consumption - produced heat  
            
    def h2power(self,h2):
        """
        Inverse function that computes the h2 consumption and Fuel cell efficiency
        corresponding to a required amount of electricity by means of interpolating functions.

        Parameters
        ----------
        h2 : float exploitable hydrogen to produce electricity in the timestep [kg]

        Returns
        -------
        hyd: float exploited hydrogen in the timestep [kg] (same as input 'h2')
        p_required : float electricity produced in the timestep [kWh]
        FC_Heat : Heat produced by Fuel cell operation in the time step [kWh]
        etaFC : module efficiency [-]

        """
        if 0 <= h2 <= self.maxh2used:
            
            hyd                = h2                                     # [kg] hydrogen consumption
            p_required         = self.h2P(hyd)                          # [kW] coverable electric power
            FC_CellCurrDensity = self.PI(p_required)/self.FC_CellArea   # [A/cm^2] current density value at which the fuel cell is working 
            
            # Checking if resulting current density is high enough for the fuel cell to start, otherwise hydrogen used = 0
    
            if FC_CellCurrDensity < self.FC_MinCurrDens:      # condition for operability set for current density 
                                                                                          
                etaFC         = 0       # [-]      fuel cell efficiency
                hyd           = 0       # [kg]     hydrogen used in the considered timestep
                Current       = 0       # [A]      Operational Current
                p_required    = 0       # [kWh]    required energy - when timestep is kept at 1 h kWh = kW
                FC_Heat       = 0       # [kWh]    thermal energy used
                FC_CellCurrDensity = 0  # [A/cm^2] current density
            
            else: 
            
                etaFC      = self.etaFuelCell(hyd)  # [-]   FC efficiency
                FC_Heat    = self.FC_Heat(hyd)      # [kWh] FC produced heat
                
        # elif h2 < 0:   #Andrea--->se non metto questo else e h<0 dà errore sul return
            
        #     etaFC         = 0       # [-]      fuel cell efficiency
        #     hyd           = 0       # [kg]     hydrogen used in the considered timestep
        #     Current       = 0       # [A]      Operational Current
        #     p_required    = 0       # [kWh]    required energy - when timestep is kept at 1 h kWh = kW
        #     FC_Heat       = 0       # [kWh]    thermal energy used
        #     FC_CellCurrDensity = 0  # [A/cm^2] current density
        
        # elif h2 > self.maxh2used:   #Andrea--->se non metto questo else e h>hmax dà errore sul return
            
        #     hyd                = self.maxh2used                                     # [kg] hydrogen consumption
        #     p_required         = self.h2P(hyd)                          # [kW] coverable electric power
        #     FC_CellCurrDensity = self.PI(p_required)/self.FC_CellArea   # [A/cm^2] current density value at which the fuel cell is working 
        #     etaFC      = self.etaFuelCell(hyd)  # [-]   FC efficiency
        #     FC_Heat    = self.FC_Heat(hyd)      # [kWh] FC produced heat
            
            return(hyd,p_required,FC_Heat,etaFC)

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
        
        if tech_cost['cost per unit'] == 'default price correlation':
            C0 = 1500 # €/kW
            scale_factor = 0.8 # 0:1
            C = size * C0 **  scale_factor
        else:
            C = size * tech_cost['cost per unit']

        tech_cost['total cost'] = tech_cost.pop('cost per unit')
        tech_cost['total cost'] = C
        tech_cost['OeM'] = tech_cost['OeM'] *C /100 # €


        self.cost = tech_cost    
#%%                     
    def plot_polarizationpts(self):
         
          if self.model == 'PEM General': 
              
              fig=plt.figure(figsize=(10,8),dpi=1000)
              fig.suptitle("PEMFC Polarization Curve STACK - P ={}".format(round(self.Npower,1)) +" kW")
               
              plt.plot(self.CellCurrDensity,self.Voltage,label='data', color='b',marker='.', linestyle='None', mec='r', markersize=7, markerfacecolor='white', zorder=0)
              plt.plot(self.x,self.iV1(self.x),label='linear interp',linestyle='--') 
              plt.grid()
              plt.legend(fontsize=8)
              plt.xlabel('Cell Current Density [A cm$^{-2}$]')
              plt.ylabel('Stak Voltage [V]')
              plt.title('PEMFC Polarization Curve (V-i)' )
              
              fig=plt.figure(figsize=(8,8),dpi=1000)
              fig.suptitle("PEMFC Polarization Curve STACK - P ={}".format(round(self.Npower,1)) +" kW")
            
              ax_1=fig.add_subplot(221)
              ax_1.plot(self.CellCurrDensity,self.Voltage,label='data', color='b',marker='.', linestyle='None', mec='r', markersize=7, markerfacecolor='white', zorder=0)
              ax_1.plot(self.x,self.iV1(self.x),label='linear interp',linestyle='--') 
              ax_1.grid()
              ax_1.legend(fontsize=8)
              ax_1.set_xlabel('Cell Current Density [A cm$^{-2}$]')
              ax_1.set_ylabel('Stak Voltage [V]')
              ax_1.set_title('PEMFC Polarization Curve (V-i)' )
            
              ax_2=fig.add_subplot(222)
              ax_2.plot(self.Current,self.Voltage,label='data', color='b',marker='.', linestyle='None', mec='r', markersize=7, markerfacecolor='white', zorder=0)
              #ax_2.plot(self.x,self.iV1(self.x),label='linear',linestyle='--') 
              ax_2.grid()
              ax_2.legend(fontsize=8)
              ax_2.set_xlabel('Cell Current [A]')
              ax_2.set_ylabel('Stak Voltage [V]')
              ax_2.set_title('PEMFC Polarization Curve (V-I)' )
         
              plt.tight_layout()
              plt.show()
             
             
              'Polarization Curve'
             
              plt.figure(dpi=600, figsize=(9,5))
              plt.plot(self.CellCurrDensity,self.CellVolt,label='data', color='b',marker='.', linestyle='None', mec='r', markersize=7, markerfacecolor='white', zorder=0)
              plt.plot(self.CellCurrDensity,self.OCpotential,label='Open Circuit Potential')
              plt.plot(self.CellCurrDensity,self.ActLosses,label='Activation Losses')
              plt.plot(self.CellCurrDensity,self.OhmLosses,label='Ohmic Losses')
              plt.plot(self.x,self.iV2(self.x),label='linear interp',linestyle='--') 
              plt.grid()
              plt.legend(fontsize=8)
              plt.xlabel('Cell Current Density [A cm$^{-2}$]')
              plt.ylabel('Stak Voltage [V]')
              plt.title('PEMFC Polarization Curve (V-i)' )
                
          elif self.model == 'SOFC':
              
              plt.figure(dpi=1000)
              plt.plot(self.CellCurrDensity,self.Voltage,label='data', color='b',marker='.', linestyle='None', mec='r', markersize=7, markerfacecolor='white', zorder=0)
              plt.plot(self.x,self.iV1(self.x),label='linear', linestyle='--') 
              plt.grid()
              plt.legend(fontsize=8)
              plt.xlabel('Cell Current Density [A cm$^{-2}$]')
              plt.ylabel('Stack Voltage [V]')
              plt.title('SOFC Polarization Curve' )
              plt.show()  
                  
              'Polarization Curve'
                 
              plt.figure(dpi=600, figsize=(9,5))
              plt.plot(self.CellCurrDensity,self.Voltage/self.nc,label='data', color='b',marker='.', linestyle='None', mec='r', markersize=7, markerfacecolor='white', zorder=0)
              plt.plot(self.CellCurrDensity,self.Ecell,label='Open Circuit Potential')
              plt.plot(self.CellCurrDensity,self.DeltaV_ohm,label='Ohmic Losses')
              plt.plot(self.CellCurrDensity,self.DeltaV_con,label='Concentration Losses')
              plt.grid()
              plt.legend(fontsize=8)
              plt.xlabel('Cell Current Density [A cm$^{-2}$]')
              plt.ylabel('Cell Voltage [V]')
              plt.title('SOFC Polarization Curve (i-V)' )
          else: 
              print('Polarization curve not available') 
              
                  
    def plot_stackperformanceSOFC(self):
        
        if self.model == 'SOFC':

            'Experimental values'
            Corrente = [3.4,20.0,25.5,30.6,40.1,55.2,65.3,80.4,90.4,100.0]                   # [A]
            VoltaggioStack = [89.0,78.1,75.7,73.7,70.1,64.8,61.3,56.3,53.0,49.9]             # [V]
            Potenza = [302.4,1563.0,1934.4,2253.5,2813.4,3576.9,4005.2,4527.0,4793.8,4990.0] # [W]
            x=np.linspace(3.4,100,100)                                                       # [A]
                
            interpolazione1=interp1d(Corrente,VoltaggioStack,bounds_error=None,kind='cubic',fill_value='extrapolate')
            interpolazione2=interp1d(Corrente,Potenza, bounds_error=None, kind='cubic', fill_value='extrapolate')
                
            fig,ax = plt.subplots(dpi=600,figsize=(9,5))
            ax2 =ax.twinx()
            ax.plot(self.Current,self.Voltage,label='V$_\mathregular{stack}$ Model')
            ax.plot(Corrente,VoltaggioStack,linestyle='None',marker='.',mec='r',markersize=10)
            ax.plot(x,interpolazione1(x),label='V data',marker='.', markersize=2)
            ax.axvline(x=100,linestyle='--',color='k')
            ax2.axhline(y=4953,linestyle='--',color='k')
            ax.grid()
            ax2.plot(self.Current,self.P*1000,label='P$_\mathregular{stack}$ Model',color='tab:red')
            ax2.plot(Corrente,Potenza,linestyle='None', marker='.',mec='b',markersize=10)
            ax2.plot(x,interpolazione2(x),label='P data',color='tab:orange',marker='.', markersize=2)
            h1, l1 = ax.get_legend_handles_labels()
            h2, l2 = ax2.get_legend_handles_labels()
            ax.legend(h1+h2, l1+l2, loc=4)
            ax.set_ylim(0,None)
            ax.set_xlabel('Current [A]')
            ax.set_ylabel('Stack Tension [V]')
            ax2.set_ylabel('P$_{el}$ [W]')       
            plt.title('SOFC Stack Performance (I-V)')
         
              
#%%             
    def plot_stackperformancePEM(self):
        
        if self.model == 'PEM General':
        
            'Datasheet FCS 13-XXL Gen 2.9'
            
            Current      = [0,40,80,120,160,200,230,250]
            StackVoltage = [94,78,73,69,66,62,59,57]
            StackPower   = [0,3.1,5.8,8.3,10.5,12.4,13.6,14.1]
            x = np.linspace(0,250,100)
            
            interp = interp1d(Current,StackVoltage,bounds_error=None,kind='cubic',fill_value='extrapolate')
            interp1= interp1d(Current,StackPower,bounds_error=None,kind='cubic',fill_value='extrapolate')
            
            'Polarization Curve'
            
            fig,ax = plt.subplots(dpi=600,figsize=(9,5))
            ax2 =ax.twinx()
            ax.plot(self.Current,self.Voltage,label='V$_\mathregular{stack}$ Model')
            ax.plot(Current,StackVoltage,linestyle='None',marker='.',mec='r',markersize=10)
            ax.plot(x,interp(x),label='V datasheet',marker='.', markersize=2)
            ax.axvline(x=230,linestyle='--',color='k')
            ax2.axhline(y=13.6,linestyle='--',color='k')
            ax.grid()
            ax2.plot(self.Current,self.P,label='P$_\mathregular{stack}$ Model',color='tab:red')
            ax2.plot(Current,StackPower,linestyle='None', marker='.',mec='b',markersize=10)
            ax2.plot(x,interp1(x),label='P datasheet',color='tab:orange',marker='.', markersize=2)
            h1, l1 = ax.get_legend_handles_labels()
            h2, l2 = ax2.get_legend_handles_labels()
            ax.legend(h1+h2, l1+l2, loc=4)
            ax.set_ylim(0,None)
            ax.set_xlabel('Current [A]')
            ax.set_ylabel('Stak Voltage [V]')
            ax2.set_ylabel('Stack Power [kW]')       
            plt.title('PEMFC STACK Performance (V-I)')
        
#%%     
    def plot_linregression(self): 
         
         if self.model == 'PEM General':
             
             'Linear Regression'
             x1 = self.CellCurrDensity.reshape((-1,1))
             y1 = self.Voltage
             
             model = LinearRegression().fit(x1,y1)
             r_sq_linreg = model.score(x1,y1)
             print('Coeff. of Determination:', r_sq_linreg)
             
             self.coeff_A = model.intercept_      # Obtaining the calculated intercept for the linear fit - returns a scalar
             print('intercept:', self.coeff_A)
             self.coeff_B = model.coef_           # Obtaining the calculated slope for the linear fit - returns an array with only 1 value
             print('slope:', self.coeff_B)
             
             Volt_LinReg = self.coeff_A + self.coeff_B*self.CellCurrDensity
             
             # i-V plot Linear Regression 
             
             fig=plt.figure(figsize=(8,8),dpi=1000)
             fig.suptitle("Prestazioni FC da {}".format(round(self.FC_NominalPower,1)) +" kW")
             
             ax_1=fig.add_subplot(121)
             ax_1.plot(self.CellCurrDensity,self.Voltage,label='data', color='b',marker='.', linestyle='None', mec='r', markersize=7, markerfacecolor='white', zorder=0)
             ax_1.plot(self.CellCurrDensity,Volt_LinReg,label='Least Squares LinReg') 
             ax_1.grid()
             ax_1.legend(fontsize=8)
             ax_1.set_xlabel('Cell Current Density [A cm$^{-2}$]')
             ax_1.set_ylabel('Stak Voltage [V]')
             ax_1.set_title('PEMFC Polarization Curve (V-i)' )
             
             ax_2=fig.add_subplot(122)
             ax_2.plot(self.Current,self.Voltage,label='data', color='b',marker='.', linestyle='None', mec='r', markersize=7, markerfacecolor='white', zorder=0)
             #ax_2.plot(self.Current,Volt_LinReg,label='Least Squares LinReg') 
             ax_2.grid()
             ax_2.legend(fontsize=8)
             ax_2.set_xlabel('Cell Current [A]')
             ax_2.set_ylabel('Stak Voltage [V]')
             ax_2.set_title('PEMFC Polarization Curve (V-I)' )
             
             plt.tight_layout()
             plt.show()
             
         else: 
             print('Polarization curve not available') 
             
#################################################################################            
#%%    # Computing FuelCell performance via Spline Interpolation 
       
    def use(self,h,e,available_hyd):
        """
        The Fuel cell can absorb hydrogen and produce electricity: H2 --> 2H+ + 2e
        
        h: int hour under examination [h]
        e: float < 0 energy required  [kWh]
        available_hyd: float available hydrogen H tank SOC[h-1] [kg]
      
        output : hydrogen absorbed and electricity supplied that hour [kg]
        """
        ##########################
        if self.model=='FCS-C5000':
           
            p_required = -e                    # [kWh]
            
            p = min(p_required,self.Npower)    # [kW] how much electricity can be absorbed by the fuel cell absorb
            
            # Calculating the operating point on the characteristic curves
            I=self.PI(p)
           
            # calculate the hydrogen consumed by each single cell
            qe=I/96485                         # [mol_e/s] Faraday
            qH=qe*0.5*3600                     # [mol_H2/h] the moles of H2 are 1/2 the moles of electrons. 3600s in one hour.
            QH=qH*2.058                        # [g_H2/h] H2 molar mass = 2.058
            
            hyd = QH*self.nc/1000 # total stack hydrogen [kg/hr]
            
            if hyd > available_hyd: # if not enough hydrogen is available to meet demand (H tank is nearly empty)
                hyd = 0
                p = 0
                # turn off the fuel cell
                
            return (-hyd,p,0,0) # return hydrogen absorbed [kg] and electricity required [kWh]

        ###############################
        if self.model=='PEM General' or self.model=='SOFC':
            
            # PowerOutput [kW] - Electric Power required from the fuel cell
            if abs(e) <= self.Npower:
                
                e_required = e
                hyd,p,FC_Heat,etaFC,hydrogen = fuel_cell.use1(self,h,e_required,available_hyd)
                if abs(hyd) > 0:
                    self.n_modules_used[h] = 1
                else:
                    self.n_modules_used[h] = 0
                self.EFF[h]   = etaFC
    
                            
            if abs(e) > self.Npower:      # if Electric Power required is higher than nominal one, i.e., more modules can be used

                n_modules_used = min(self.n_modules,int(abs(e)/self.Npower))
                e_required =  -self.Npower                     # power required by the single module   
                hyd,p,FC_Heat,etaFC,hydrogen = fuel_cell.use1(self,h,e_required,available_hyd)
                if hydrogen == hyd:                            # if required hydrogen is equal to the maximum requestable one, i.e., if the required hydrogen is lower than the available one. Otherwise it means that only one module can be used
                    hyd_11 = np.zeros(n_modules_used+1)
                    for i in range (n_modules_used+1):
                        hyd_11[i] = hydrogen*i                   # Saving the requested hydrogen for each number of n_modules_used
                        if abs(hyd_11[i]) > available_hyd and abs(hyd_11[i-1]) < available_hyd: #if, using i modules, the total amount of requested hydrogen is higher than the available one 
                            hyd_1 = hyd_11[i-1]                 # Hydrogen requested using i-1 modules
                            p_1 = p*(i-1)                       # Total power requested using i-1 modules 
                            FC_Heat_1 = FC_Heat*(i-1)
                            self.EFF[h] = etaFC                 # work efficiency of modules working at nominal power 
                        
                            remained_available_hyd = available_hyd-abs(hyd_1)
                            hyd,p,FC_Heat,etaFC = fuel_cell.h2power(self,remained_available_hyd)
                            if abs(hyd) > 0:
                                self.n_modules_used[h] = i
                                self.EFF_last_module[h] = etaFC
                            else:
                                self.n_modules_used[h] = i-1
                            hyd = -hyd+hyd_1
                            p = p+p_1
                            FC_Heat = FC_Heat+FC_Heat_1
          
                    if abs(hyd_11[-1]) <= available_hyd:       # if, using n_modules, the total amount of requested hydrogen is lower than available one  
                        hyd_1 = hyd*n_modules_used             # total amount of H2 requested by modules working at full load
                        p_1 = p*n_modules_used                 # total power produced    // // // // //  
                        FC_Heat_1 = FC_Heat*n_modules_used
                        self.EFF[h] = etaFC                    # work efficiency of modules working at nominal power 
                        
                        if abs(e) <= self.MaxPowerTot:                             # if available energy is lower than total installed power
                            e_excess_req = -(abs(e)-self.Npower*n_modules_used)    # remaining requested power after considering modules at full load
                            remained_available_hyd = available_hyd-abs(hyd_1)
                            hyd,p,FC_Heat,etaFC,hydrogen = fuel_cell.use1(self,h,e_excess_req,remained_available_hyd)
                            if abs(hyd) > 0:
                                self.n_modules_used[h] = n_modules_used+1
                                self.EFF_last_module[h] = etaFC
                            else:
                                self.n_modules_used[h] = n_modules_used
                            hyd = hyd+hyd_1
                            p = p+p_1
                            FC_Heat = FC_Heat+FC_Heat_1   
                            
                        else:                     # if available energy is higher than total installed power it means that there are no more exploitable modules, thus no more hydrogen can be produced 
                            p = p_1
                            hyd = hyd_1
                            FC_Heat = FC_Heat_1
                            self.n_modules_used[h] = self.n_modules
                            
                else:                             # only one module used because requested hydrogen from one module is higher than the available one 
                    if abs(hyd) > 0:
                        self.n_modules_used[h] = 1
                    else:
                        self.n_modules_used[h] = 0
                    self.EFF[h] = etaFC

        return (hyd,p,FC_Heat)  # return hydrogen absorbed [kg] electricity required [kWh] and heat as a co-product [kWh]
        
            
    def use1(self,h,e,available_hyd):
         
        'Finding the working point of the FuelCell by explicitly solving the system:'
        
        if self.model == 'PEM General':
            
            p_required = -e                                             # [kWh] --> equivalent to kW for the considered timestep
            FC_CellCurrDensity = self.PI(p_required)/self.FC_CellArea   # [A/cm^2] current density value at which the fuel cell is working 
            
            # Checking if resulting current density is high enough for the fuel cell to start, otherwise hydrogen used = 0
    
            if FC_CellCurrDensity < self.FC_MinCurrDens:      # condition for operability set for current density 
                                                                                          
                etaFC              = 0      # [-] fuel cell efficiency
                hyd                = 0      # [kg] hydrogen used in the considered timestep
                Current            = 0      # [A] Operational Current
                p_required         = 0      # [kWh] required energy - when timestep is kept at 1 h kWh = kW
                FC_Heat            = 0      # [kWh] thermal energy used
                FC_CellCurrDensity = 0
                hydrogen           = hyd
            
            else: 
                
                Current = FC_CellCurrDensity*self.FC_CellArea     # [A] FuelCell Stack operating current 
                
                FC_Vstack= self.iV1(FC_CellCurrDensity)           # [V] Stack operating voltage
                V_cell=FC_Vstack/self.nc
                # FC_Vstack_ = (p/Current)*1000                     # [V] voltage check
        
                'Computing FC efficiency and  hydrogen energy demand'    
                
                pO2 = (self.FC_AirPress*0.21)/101325              # [atm]
                pH2O = 1                                          # [atm]
                pH2 = self.FC_FuelPress/101325                    # [atm]  
                
                Ecell = 1.229-0.85e-3*(self.FC_OperatingTemp-298.15) + 4.3085e-5*self.FC_OperatingTemp*ln(pH2*(pO2**0.5)/pH2O)   # [V] Open circuit voltage 
                deltaG = self.deltaG0 - self.Runiv*self.FC_OperatingTemp*ln(pH2*math.sqrt(pO2)/pH2O)/(2*self.FaradayConst)       # [kJ/mol] Gibbs free energy at actual conditions
              
                eta_voltage = FC_Vstack/(Ecell*self.nc)           # [-] Voltage efficiency 
                eta_th = - deltaG/self.HHVh2Mol                   # [-] Thermodynamic efficiency
             
                etaFC = eta_th*eta_voltage                        # [-] FC efficiency
                
                self.VOLT[h]      = V_cell                        # [V]      Cell voltage history
                self.CURR_DENS[h] = FC_CellCurrDensity            # [A/cm2]  Cell current density history
                
                'Hydrogen demand'
             
                FC_HydroCons = Current*self.nc*3600/(95719.25*1000)/self.rhoStdh2*self.timestep   # [Sm3] (Chavan 2017)
                hyd = FC_HydroCons*self.rhoStdh2/self.timestep                                    # [kg/h]
                FC_deltaHydrogen = - FC_HydroCons*self.rhoStdh2*self.HHVh2*1000/3600              # [kWh]
             
                'Process heat, that can be recovered'
              
                FC_Heat = ((1.481*self.nc)/FC_Vstack-1)*p_required*self.timestep               # [kWh] --> equivalent to kW for the considered timestep of 1h
                
                hydrogen = -hyd   # [kg] check value to be implemented in use function
                
                if hyd > available_hyd: # if not enough hydrogen is available to meet demand (H tank is nearly empty)
                    
                    # defining the electric load that can be covered with the available hydrogen
                    hyd,p_required,FC_Heat,etaFC = fuel_cell.h2power(self,available_hyd)

            return (-hyd,p_required,FC_Heat,etaFC,hydrogen)  # return hydrogen absorbed [kg] electricity required [kWh] and heat as a co-product [kWh]

        if self.model == 'SOFC':
            
            h2oMolMass       = c.H2OMOLMASS   # [kg/mol]     Water molar mass
            H2MolMass        = c.H2MOLMASS    # [kg/mol]     Hydrogen molar mass
            H2MolStdEntropy  = c.H2MOL_S_E    # [J/K*mol]    Specific molar entropy (gaseous phase)
            O2MolStdEntropy  = c.O2MOL_S_E    # [J/K*mol]    Specific molar entropy
            H20MolStdEntropy = c.H2OMOL_S_E   # [J/K*mol]    Specific molar entropy
            
            
            p_required=-e                                               # [kWh] --> equivalent to kW for the considered timestep
            
            FC_CellCurrDensity = self.PI(p_required)/self.FC_CellArea   # [A/cm^2] current density value at which the fuel cell is working 
            
            # Checking if resulting current density is high enough for the fuel cell to start, otherwise hydrogen used = 0
    
            if FC_CellCurrDensity < self.FC_MinCurrDens:      # condition for operability set for current density 
                                                                                          
                etaFC         = 0      # [-] fuel cell efficiency
                hyd           = 0      # [kg] hydrogen used in the considered timestep
                Current       = 0      # [A] Operational Current
                p_required    = 0      # [kWh] required energy - when timestep is kept at 1 h kWh = kW
                FC_Heat       = 0      # [kWh] thermal energy used
                FC_CellCurrDensity = 0
                hydrogen = hyd
            
            else: 
 
                Current = FC_CellCurrDensity*self.FC_CellArea      # [A] FuelCell Stack operating current
             
                FC_Vstack= self.iV1(FC_CellCurrDensity)            # [V] Stack operating voltage
    
                'Computing FC efficiency and hydrogen energy demand'    
                
                pO2 = (self.FC_AirPress*0.21*101325)/101325       # [atm]
                pH2O = 1                                          # [atm]
                pH2 = self.FC_FuelPress/101325                    # [atm]  
                
                Ecell = 1.19 + ((c.R_UNIVERSAL*self.FC_OperatingTemp)/(2*c.FARADAY))*ln(pH2*(pO2**0.5)/pH2O)
                DeltaV_OHM = (FC_CellCurrDensity*self.ThicknessElectrolyte*self.FC_OperatingTemp)/(9000*self.e**(-100000/(c.R_UNIVERSAL*self.FC_OperatingTemp)))
                DeltaV_CON = ((c.R_UNIVERSAL*self.FC_OperatingTemp)/(2*c.FARADAY))*ln((FC_CellCurrDensity/(self.FC_ExchangeCurrDensCathode*self.FC_AirPress*(0.21-0.0008*10000*FC_CellCurrDensity*c.R_UNIVERSAL*self.FC_OperatingTemp/(4*c.FARADAY*101325*0.00002)))))
    
                TotalLoss = DeltaV_OHM + DeltaV_CON   # [V]
                
                deltaG = c.GIBBS - c.R_UNIVERSAL*self.FC_OperatingTemp*ln(pH2*math.sqrt(pO2)/pH2O)/(2*c.FARADAY)       # [kJ/mol] Gibbs free energy at actual conditions
              
                eta_voltage = FC_Vstack/(Ecell*self.nc)           # [-] Voltage efficiency 
                eta_th = c.GIBBS/self.HHVh2Mol                    # [-] Thermodynamic efficiency
             
                etaFC = -eta_th*eta_voltage                       # [-] FC efficiency
                
                'Hydrogen demand'
                
                FC_HydroCons     = ((Current*self.nc*3600)/(c.FARADAY*1000))/(self.rhoStdh2*self.timestep)     # [Sm^3] volumetric hydrogen consumption in the timestep
                hyd              = FC_HydroCons*self.rhoStdh2/self.timestep                                    # [kg/h] hourly hydrogen consumption
                FC_deltaHydrogen = - FC_HydroCons*self.rhoStdh2*self.HHVh2*1000/3600                           # [kWh]
                
                z = Current/(2*c.FARADAY)   # [mol/s]
                
                DeltaS = -((H20MolStdEntropy-(O2MolStdEntropy/2)-H2MolStdEntropy)+(c.R_UNIVERSAL/2)*ln((pH2**2)*pO2/(pH2O**2))) # [J/mol*K]
                
                FC_Heat = (((z*self.FC_OperatingTemp*DeltaS + Current*TotalLoss)*self.timestep)*self.nc)/1000   # [kWh] --> equivalent to kW for the considered timestep of 1h
            
                hydrogen = -hyd          # [kg] check value to be implemented in use function
                
                if hyd > available_hyd: # if not enough hydrogen is available to meet demand (H tank is nearly empty)
                    
                    # defining the electric load that can be covered with the available hydrogen
                    hyd,p_required,FC_Heat,etaFC = fuel_cell.h2power(self,available_hyd)

            return (-hyd,p_required,FC_Heat,etaFC,hydrogen)  # return hydrogen absorbed [kg] electricity required [kWh] and heat as a co-product [kWh]


    def calculate_ageing(self,economic_data):
        
        # ref. study  https://doi.org/10.1016/j.ijhydene.2017.02.146
        
        load_i = self.CURR_DENS[:8760]          # defining fuel cell load profile - current for the 1st year of operation 
        load_V = self.VOLT[:8760]               # defining fuel cell load profile - voltage for the 1st year of operation 
        self.i_round=np.round(load_i,3)
        self.i_pos=np.where(self.i_round!=0)    # non-null values
        self.v_round=np.around(load_V,3)
        self.v_pos=np.where(self.v_round>0)[0]
        h_utilizzo_year=len(self.i_pos[0])      # working hours over 1 year of operations
        self.replacements=[]

        # Parameters definition

        #current

        k_1=25      # current costant (paper)

        #voltage

        k_2=5       # voltage constant (paper)
        self.v_0=self.Voltage[-1]/self.nc       #Voltaggio minimo cella
        self.vol_max=self.Voltage[0]/self.nc    #Voltaggio massimo cella
        self.v_L=0.6*self.vol_max               #valore minimo range ottimale 
        self.v_U=0.8*self.vol_max               #valore massimo range ottimale 

        #life time

        self.v_rated=self.vol_max

        v_dot_data=np.array([4,1,2,11,6,50,50,75,200,300,260,400])    #dataset V_dot (paper) "migliore"
        phi_data=np.array([1,1,1,1,1,3,4.3,5,7,7.5,9.2,9])              #dataset phi (paper) "migliore"
        polyfit_forced= [4.1353, 0.0925,0]                                # funzione interpolante "migliore"--->mi serve phi=1.121 per avere 20000 h (v_deg=5.3)
        
        # v_dot_data=np.array([4,1,2,11,25,30,6,50,50,75,200,300,260,400])    #dataset V_dot (paper)
        # phi_data=np.array([1,1,1,1,1,1,1,3,4.3,5,7,7.5,9.2,9])              #dataset phi (paper)
        
        paper_fit=[40,-46]                                                  #funzione interpolante lineare v_dot(phi) dal paper--->mi serve phi=0.93 (impossibile) per avere 20000 h (v_deg=5.3)

        # polyfit_forced= [3.8913, 2.0875,0]
        
        exp_fit=scipy.optimize.curve_fit(lambda t,a,b: a*np.exp(b*t),  phi_data,  v_dot_data)   #funzione interpolante esponenziale v_dot(phi)
        a=exp_fit[0][0]
        b=exp_fit[0][1]

        def v_phi_fitting_plot():
            x=np.arange(11)
            h=np.poly1d(paper_fit)
            t=h(x)
            f=np.poly1d(polyfit_forced)
            y=f(x)
            g=a*np.exp(b*x)
            fig=plt.figure(dpi=1000)
            ax=fig.add_subplot(111)
            ax.plot(phi_data,v_dot_data,marker="o",linestyle="",label='data')
            ax.plot(x,t,label='linear fitting [],R^2=0,8738')
            ax.plot(x,y,label='polynomial fitting,R^2=0,9132')
            ax.plot(x,g,label='exponential fitting,R^2=0,7784')
            ax.set_xlabel('phi [-]')
            ax.set_ylabel('v_deg_rate [10e-6V/h]')
            ax.set_title('Fitting of degradation rate')
            ax.legend(fontsize=9)
            ax.grid()
            plt.show()

        v_phi_fitting_plot()
        
        ## Load current

        #DFT current

        DFT_i=scipy.fft.fft(load_i)     # Fast Fourier Transorm (Discretized)

        N = len(DFT_i)
        n = np.arange(N)

        Fs = 1 / (60*60)    # sampling frequency
        T = N/Fs
        freq = n/T 
        
        if Fs<freq[-1]:
            f_max=Fs        #paper
        else:
            f_max=freq[-1]

        DFT_i=scipy.fft.fft(load_i)     # Fast Fourier Transorm (Discretized)
        DFT_i_mag=np.abs(DFT_i)/N

        def DFT_current_plot():
            plt.figure(dpi=1000)      
            plt.plot(freq[0:int(N/2+1)],2*DFT_i_mag[0:int(N/2+1)])
            plt.grid()
            plt.xlim(0,f_max)
            plt.ylabel('|FFT_i(J)| [A/cm^2]')
            plt.xlabel('Frequency (Hz)')
            plt.title('Single-Sided Amplitude Spectrum of J(t)')
            plt.show() 
        
        DFT_current_plot()
        
        #weight current

        w_curr=[k_1,1]

        #current value

        grad=3
        p_Cln=np.polyfit(freq,np.abs(DFT_i),grad)  #costanti della funzione F(omega) di terzo grado trovata col fitting di DFT sulla frequenza
        prod=np.polyint(np.poly1d(p_Cln)*np.poly1d(w_curr))           #integro il prodotto tra F(omega) e w_curr indefinito
        I=np.polyval(prod,f_max)-np.polyval(prod,0)    #rendo l'integrale definito tra f_max e 0
        self.m_curr=1/T*I+1
        print(self.m_curr)

        ## Load voltage

        # Voltage histogram

        self.bins = np.arange(self.v_0-0.01,self.vol_max+0.01,0.01) # define some self.bins that cover the range of interest
        self.delta_v=self.bins[1]-self.bins[0]
        self.v_counts=np.histogram(self.v_round,self.bins)[0]        
        self.H_v=self.v_counts/len(self.v_pos)
        self.H_v_tot=self.H_v.sum()
        
        def v_count_plot():
            fig=plt.figure(dpi=1000)
            #self.v_counts=plt.hist(self.v_round, self.bins, density=False, facecolor='b', alpha=0.75, edgecolor='black',rwidth=0.8)[0]
            plt.stairs(self.v_counts,self.bins,fill=True,facecolor='b',edgecolor='black')
            plt.vlines(x=self.v_L,ymin=0,ymax=max(self.v_counts),color="red",linestyle="dashed",label="v_L")
            plt.vlines(x=self.v_U,ymin=0,ymax=max(self.v_counts),color="red",linestyle="dashed",label="v_U")
            plt.xlim([self.v_0-0.1 ,self.vol_max+0.1])
            plt.xlabel('Operation Voltage Range [V]')
            plt.ylabel('Voltage counts')
            #plt.xticks([self.v_L,self.v_U],['v_L','v_U'],weight='bold')
            plt.text(self.v_L,-100,"v_L",weight='bold',horizontalalignment='center')
            plt.text(self.v_U,-100,"v_U",weight='bold',horizontalalignment='center')
            plt.grid()
            plt.show()
            
        v_count_plot()
        
        def H_v_plot():
            fig1=plt.figure(dpi=1000)
            plt.stairs(self.H_v*100,self.bins,fill=True,facecolor="b",edgecolor='black')
            plt.vlines(x=self.v_L,ymin=0,ymax=max(self.H_v)*100,color="red",linestyle="dashed",label="v_L")
            plt.vlines(x=self.v_U,ymin=0,ymax=max(self.H_v)*100,color="red",linestyle="dashed",label="v_U")
            plt.xlim([self.v_0-0.1 ,self.vol_max+0.1])
            plt.xlabel('Operation Voltage Range [V]')
            plt.ylabel('Voltage distribution [%]')
            plt.text(self.v_L,-1,"v_L",weight='bold',horizontalalignment='center')
            plt.text(self.v_U,-1,"v_U",weight='bold',horizontalalignment='center')
            plt.grid()
            plt.show()
        
        H_v_plot()
        
        self.bins=self.bins[1:]
        
        # definition of w_vol function

        w_vol_1= lambda v: (k_2*(v-self.v_L))**2+1   #if v<self.v_L

        w_vol_2=1                               #if self.v_L<v<self.v_U

        w_vol_3= lambda v: (k_2*(v-self.v_U))**2+1   #if v>self.v_U
        
        # voltage value (summation method)
        
        
        m=0
        for indice in range(len(self.bins)):
            
            if self.bins[indice]<self.v_L:
                val=self.H_v[indice]*w_vol_1(self.bins[indice]-self.delta_v/2)
                
            elif self.bins[indice]>=self.v_L and self.bins[indice]<=self.v_U:
                val=self.H_v[indice]*w_vol_2
                
            elif self.bins[indice]>self.v_U:
                val=self.H_v[indice]*w_vol_3(self.bins[indice]-self.delta_v/2)
                
            m+=val
            
        self.m_vol=m
        print("m_vol="+str(self.m_vol))
        
        ## Characteristic value of load

        self.phi=self.m_curr*self.m_vol
        print("phi="+str(self.phi))

        ## Life time definition

        v_deg_rate_paper=paper_fit[0]*self.phi+paper_fit[1]
        self.v_deg_rate_pol=(polyfit_forced[0]*(self.phi**2))+(polyfit_forced[1]*self.phi)+polyfit_forced[2] #[microV/h]
        v_deg_rate_exp=a*np.exp(self.phi*b)
        print("v_deg="+str(self.v_deg_rate_pol))
        
        T_cell_h_paper=(0.1*self.v_rated)/(v_deg_rate_paper/1e6)   #[h]
        T_cell_y_paper=T_cell_h_paper/h_utilizzo_year

        self.T_cell_h_pol=(0.1*self.v_rated)/(self.v_deg_rate_pol/1e6)   #[h]
        self.T_cell_y_pol=self.T_cell_h_pol/h_utilizzo_year
        print("T_cell="+str(self.T_cell_h_pol))

        T_cell_h_exp=(0.1*self.v_rated)/(v_deg_rate_exp/1e6)   #[h]
        T_cell_y_exp=T_cell_h_exp/h_utilizzo_year

        n_replacements=int(economic_data['investment years']*(8760/self.T_cell_h_pol))
        for n in np.arange(n_replacements):
            self.replacements.append(self.T_cell_h_pol*(1+n))
        
        
        return self.replacements
    
    def plot_corrected_polcurves(self,h):
        
        fig=plt.figure(figsize=(12,8),dpi=1000)
        fig.suptitle("Prestazioni FC da {}".format(round(self.FC_NominalPower,1)) +" kW con ageing",fontsize=25)
        
        ax_1=fig.add_subplot(111)
        ax_1.plot(self.CellCurrDensity,self.Voltage,label='BoL', color='b',linewidth=3.0)
        ax_1.plot(self.x,self.iV1(self.x),label='h='+str(h),color='r',linewidth=3.0)
        #ax_1.plot(self.x,self.iV1(self.x),label="EoL",color='r',linewidth=3.0) 
        ax_1.grid()
        ax_1.legend(fontsize=20)
        ax_1.set_xlabel('Cell Current Density [A cm$^{-2}$]',fontsize=20)
        ax_1.set_ylabel('Stak Voltage [V]',fontsize=20)
        ax_1.set_title('PEMFC Polarization Curve (V-i)',fontsize=20 )
        
        plt.tight_layout()
        plt.show()
    
#%%##########################################################################################

if __name__ == "__main__":
    
    """
    Functional test
    """
    
    inp_test = {'Npower': 1000,
                "number of modules": 4,
                'stack model':'SOFC',
                'ageing': True
                }
    
    sim_hours=60                               # [h] simulated period of time - usually it's 1 year minimum
    time=np.arange(sim_hours)
    
    fc = fuel_cell(inp_test,sim_hours)         # creating fuel cell object
    # fc.plot_polarizationpts()                  # cell polarization curve
    
    available_hydrogen = 1000                  # [kg] hydrogen available in the storage system

    if fc.model == 'PEM General' and fc.Npower == 13.6:
        fc.plot_stackperformancePEM()                   # PEMFC stack performance curve
    if fc.model == 'SOFC' and fc.Npower == 5:
        fc.plot_stackperformanceSOFC()                  # SOFC stack performance curve

    'Test 1 - Tailored ascending power input'

    flow  = - np.linspace(0,fc.Npower*6,sim_hours)  # [kW] power input - ascending series
    flow1 = - np.linspace(0,fc.Npower,sim_hours)    # [kW] power input - ascending series

    hyd_used = np.zeros(len(flow))      # [kg] hydrogen used by fuel cell
    P_el     = np.zeros(len(flow))      # [kW] electricity produced
    P_th     = np.zeros(len(flow))      # [kW] produced heat
    
    for h in range(len(flow1)):
        hyd_used[h],P_el[h],P_th[h] = fc.use(h,flow1[h],available_hydrogen)
        
    # fc.EFF[fc.EFF == 0] = math.nan      # activate to avoid representation o '0' values when fuel cell is turned off
    
    fig=plt.figure(figsize=(8,8),dpi=1000)
    fig.suptitle("{} ({} kW) performance".format(inp_test['stack model'],round(fc.FC_NominalPower,1)))
    
    PI=fig.add_subplot(211)
    PI.plot(-flow1,P_th,label="Thermal Power") 
    PI.axvline(x=fc.MinPower,linestyle=':',label= 'Lower Functioning Boundary', zorder=3, linewidth = 2)   
    PI.set_title("Heat vs Electric Power")
    PI.grid()
    PI.set_xlabel("Power Output [kW]")
    PI.set_ylabel("Heat Output [kW]")
    PI.legend(fontsize=15)
    
    ETA=fig.add_subplot(212)
    ETA.scatter(-flow1,fc.EFF,label="Efficiency",color="green",edgecolors='k')
    ETA.axvline(x=fc.MinPower,color='tab:blue',linestyle=':',label= 'Lower Functioning Boundary', zorder=3, linewidth = 2)   
    ETA.set_title("Efficiency vs Power")
    ETA.grid()
    ETA.set_xlabel("Power Output [kW]")
    ETA.set_ylabel("Efficiency [-]")
    ETA.legend(fontsize=15)
    
    plt.tight_layout()
    plt.show()
    
    
    fig=plt.figure(figsize=(8,8),dpi=1000)
    fig.suptitle("{} ({} kW) performance".format(inp_test['stack model'],round(fc.FC_NominalPower,1)))
    
    PI=fig.add_subplot(211)
    PI.plot(-flow1,-hyd_used)
    PI.set_title("H$_{2}$ Consumption vs Power")
    PI.grid()
    PI.set_xlabel("Power Output [kW]")
    PI.set_ylabel("Hydrogen consumption [kg/h]")
    # PI.legend(fontsize=15)
    
    ETA=fig.add_subplot(212)
    ETA.scatter(-flow1,fc.EFF,label="Efficiency",color="green",edgecolors='k',zorder =3)
    ETA.axvline(x=fc.MinPower,color='tab:blue',linestyle=':',label= 'Lower Functioning Boundary', zorder=3, linewidth = 2)   
    ETA.set_title("Efficiency vs Power")
    ETA.grid()
    ETA.set_xlabel("Power Output [kW]")
    ETA.set_ylabel("Efficiency [-]")
    ETA.legend(fontsize=15)
    
    plt.tight_layout()
    plt.show()
    
    fig, ax = plt.subplots(dpi=600)
    ax.scatter(-flow1,fc.EFF, edgecolors='k', zorder = 3)
    ax.set_title("Fuel Cell Module Efficiency")
    textstr = '\n'.join((
        r'$CellArea=%.1f$ $cm^{2}$' % (fc.FC_CellArea,),
        r'$P_{nom}= %.1f$ kW' % (fc.Npower,),
        r'$i_{max}= %.1f$ A $cm^{-2}$' % (fc.FC_MaxCurrDens,),
        r'$n_{cell}= %.0f$' % (fc.nc,)))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax.text(fc.Npower/2,0.2,textstr,fontsize=10,va='bottom',backgroundcolor='none', bbox=props)
    ax.grid()
    ax.set_xlabel('Power Output [kW]')
    ax.set_ylabel('$\\eta$') 

   
    for h in range(len(flow)):
        hyd_used[h],P_el[h],P_th[h] = fc.use(h,flow[h],available_hydrogen)
        
    fig, ax = plt.subplots(dpi=600)
    ax.bar(-flow,fc.n_modules_used,color='tab:green',width=60,zorder=3)
    ax.set_xlabel('Required Power [kW]')
    ax.set_ylabel('Active modules [-]')
    ax.grid()
    ax.set_title('Fuel Cell Stack - Nr of working modules')
  
    
    'Test 2 - Random power input'

    fd   = -np.random.uniform(0.08*fc.Npower,5.2*fc.Npower,sim_hours)   # [kW] power required from the Fuel Cell - random values
    
    for h in range(len(fd)):
        hyd_used[h],P_el[h],P_th[h] = fc.use(h,fd[h],available_hydrogen)
        
    fig, ax = plt.subplots(dpi=1000)
    ax2 = ax.twinx() 
    ax.bar(np.arange(sim_hours)-0.2,fc.EFF,width=0.35,zorder=3,edgecolor='k',label='$1^{st}$ module efficiency', alpha =0.8)
    ax.bar(np.arange(sim_hours)+0.,fc.EFF_last_module,width=0.35,zorder=3, edgecolor = 'k',align='edge',label='Last module efficiency',alpha =0.8)
    ax2.scatter(np.arange(sim_hours),-fd,color ='limegreen',s=25,edgecolors='k',label='Required Power')
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1+h2, l1+l2, loc='lower center',bbox_to_anchor=(0.5, 1.08), ncol =3, fontsize ='small')
    ax.set_xlabel('Time [h]')
    ax.set_ylabel('Efficiency [-]')
    ax2.set_ylabel('Power Output [kW]')
    ax.grid()
    ax.set_title('Fuel Cell Stack functioning behaviour')
    
    num = 24   # number of hours to be represented in the plot below
    
    fig, ax = plt.subplots(dpi=1000)
    ax.bar(np.arange(num)-0.2,P_el[:num],width=0.35,zorder=3,color='lightseagreen',edgecolor='k',label='Electric output', alpha =0.8)
    ax.bar(np.arange(num)+0.,P_th[:num],width=0.35,zorder=3,color='indianred',edgecolor='k',align='edge',label='Thermal output',alpha =0.8)
    ax.legend(loc='lower center',bbox_to_anchor=(0.5, 1.08), ncol =3, fontsize ='small')
    ax.set_xlabel('Time [h]')
    ax.set_ylabel('Power [kW]')
    ax2.set_ylabel('Power Output [kW]')
    ax.grid()
    ax.set_title('Fuel Cell Stack functioning behaviour')