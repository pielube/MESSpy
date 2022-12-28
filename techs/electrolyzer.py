import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
from sklearn.linear_model import LinearRegression
from numpy import log as ln
import os
import sys 
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),os.path.pardir)))   # temorarily adding constants module path 
import constants as c

class electrolyzer:
    
    def __init__(self,parameters,simulation_hours):
        """
        Create an electrolyzer object
    
        parameters : dictionary
            'Npower': float nominal power [kW] - optional
            'stack model': str 'Enapter 2.1','McLyzer 800' are aviable or 'PEM General'
                      
        output : electrolyzer object able to:
            abrosrb electricity and produce hydrogen .use(e)
        """
        
        self.rhoNrh2 = c.H2NDENSITY          # [kg/m^3] hydrogen density under normal condition
        self.Npower = parameters['Npower']   # [kW] float nominal power of electrolyzers installed capacity for the location
    

        if parameters['stack model'] == 'Enapter 2.1':   # https://www.enapter.com/it/newsroom/enapter-reveals-new-electrolyser-el-2-1
            self.model= parameters['stack model'] 
            stack_operative_power_consumption = 2.4      # [kW] Single module nominal power
            stack_production_rate_L = 500                # [Nl/h]
            stack_production_rate = stack_production_rate_L / 1000 * self.rhoNrh2 # [kg/h]
            self.production_rate = stack_production_rate * (self.Npower / stack_operative_power_consumption) # [kg/h] actual production defined in a simplified way 
                                                                                                             # by dividing the selected electorlyzer size in the location
                                                                                                             # for the single unit nominal power
       
        if parameters['stack model'] == 'McLyzer 800':   # https://mcphy.com/it/apparecchiature-e-servizi/elettrolizzatori/large/
            self.model= parameters['stack model']      
            stack_operative_power_consumption = 4000     # [kW] Single module nominal power
            stack_production_rate_m3 = 800               # [Nm3/hr]
            stack_production_rate = stack_production_rate_m3 * self.rhoNrh2   # [kg/hr]
            self.production_rate = stack_production_rate * (self.Npower / stack_operative_power_consumption) # [kg/h] actual production defined in a simplified way 
                                                                                                             # by dividing the selected electorlyzer size in the location                                                                                                     # for the single unit nominal power      
                                                                                                             # for the single unit nominal power
    
        if parameters['stack model'] == 'Hylizer 1000':  # https://www.ie-net.be/sites/default/files/Presentatie%205_Baudouin%20de%20Lannoy.pdf
            self.model= parameters['stack model']      
            stack_operative_power_consumption = 5000     # [kW] Single module nominal power
            stack_production_rate_m3 = 1000              # [Nm3/hr]
            stack_production_rate = stack_production_rate_m3 * self.rhoNrh2   # [kg/hr]
            self.production_rate = stack_production_rate * (self.Npower / stack_operative_power_consumption) # [kg/h] actual production defined in a simplified way 
                                                                                                             # by dividing the selected electorlyzer size in the location
                                                                                                             # for the single unit nominal powe
        if parameters['stack model'] == 'PEM General':
            
            self.EFF = np.zeros(simulation_hours)   # keeping track of the elecrolyzer efficiency over the simulation
            self.wat_cons = np.zeros(simulation_hours)
            self.EFF_last_module = np.zeros(simulation_hours)   
            self.wat_cons_last_module = np.zeros(simulation_hours)
            self.n_modules_used=np.zeros(simulation_hours)
            self.cell_currdens = np.zeros(simulation_hours)   # cell current density at every hour
            
            self.model= parameters['stack model']
            self.rhoStdh2o    = c.H2OSDENSITY    # [kg/m3]     H2O density @ T = 15°C p = 101325 Pa
            Runiv             = c.R_UNIVERSAL    # [J/(mol*K)] Molar ideal gas constant
            self.FaradayConst = c.FARADAY        # [C/mol]     Faraday constant
            self.LHVh2        = c.LHVH2          # [MJ/kg]     H2 LHV
            self.HHVh2Mol     = c.HHVH2MOL       # [kJ/mol]    H2 HHV molar
            self.h2oMolMass   = c.H2OMOLMASS     # [kg/mol]    Water molar mass
            self.H2MolMass    = c.H2MOLMASS      # [kg/mol]    Hydrogen molar mass

            # Math costants
            self.eNepero      = c.NEPERO         # [-]         Euler's number
            # Ambient conditions 
            self.AmbTemp      = c.AMBTEMP        # [K]         Standard ambient temperature - 15 °C
           
            # At current development stage it is taken for granted we are working with hourly balances 
            self.timestep = 1              # [h]

                
            "2H2O --> 2H2 + O2"
            # https://doi.org/10.1016/j.ijhydene.2008.11.083    # Electrolyzer Parameters - depending on the different types of electrolyzers chosen - model to be specified also in class - name
        
            self.nc                  = 10            # number of cells in the stack
#            self.MembThickness       = 250           # [micron]   
            self.MembThickness       = 158.1842      # [micron]
            self.Lambda              = 20            # [-] Cell mositure content
            self.AnodeCurrDensity    = 0.00013       # [A/cm^2]
            self.AnodePressure       = 101325        # [Pa]
            self.CathodeCurrDensity  = 0.001         # [A/cm^2]
            self.CTCanode            = 0.6           # [-] Charge transfer coefficient - Anode   https://www.sciencedirect.com/science/article/pii/S0360319918309017
            self.CTCcathode          = 0.45          # [-] Charge transfer coefficient - Cathode //      //      //
            self.CurrDensityMin      = 0.005         # [A/cm^2]
            self.OperatingTemp       = 273.15 + 70   # [K]
            self.OperatingPress      = 4000000       # [Pa]
            self.n_modules           = parameters['number of modules']
            self.MaxPowerStack         = self.n_modules*self.Npower
            
            
            
            self.nc=10+int((self.Npower/1000)*(35-10))
            self.CurrDensityMax=2.1+(self.Npower/1000)*(3-2.1)
            # if 0<=self.Npower<=10:
            #    self.nc = 10
            #    self.CurrDensityMax = 2            # [A/cm^2] https://www.sciencedirect.com/science/article/pii/S266638642030151X#:~:text=In%20contrast%2C%20PEM%20electrolyzers%20experience,at%20high%20current%20density%20operations.&text=While%20commercial%20electrolyzers%20typically%20operate,reported%20by%20Lewinski%20et%20al.
                
            # if 10<self.Npower<=50:
            #    self.nc = 20
            #    self.CurrDensityMax = 2.375        # [A/cm^2]
                  
            # if 50<self.Npower<=200:
            #    self.nc = 30
            #    self.CurrDensityMax = 2.75         # [A/cm^2] 
                  
            # if 200<self.Npower<=500:
            #    self.nc = 40
            #    self.CurrDensityMax = 3.125        # [A/cm^2] 
             
            # if 500<self.Npower<=1000:
            #    self.nc = 50
            #    self.CurrDensityMax = 3.5          # [A/cm^2] 
           
            # Computing the electrolyzer polarization curve V-i
            
            'POLARIZATION CURVE'
            
            Ndatapoints = 3000                     # Number of points used to compute the polarization curve 
    
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
            self.Current           = self.CellCurrDensity*self.CellArea

            'Interpolation of calculated functioning points to detect the best fit-function for i-V curve'

            self.num = Ndatapoints                                            # Number of intervals to be considered for the interpolation
            self.x2 = np.linspace(0.05,max(self.CellCurrDensity),self.num)    # Setting xlim for range of validity of LinRegression Calculation - Only for plot-related reasons 

            # Interpolation
           
            self.iV1 = interp1d(self.CellCurrDensity,self.Voltage,bounds_error=False)#,fill_value='extrapolate')                # Linear spline 1-D interpolation                                                                                                                                                                                                                                                

            # Defining Electrolyzer Max Power Consumption
            
            Power_inp = []                                  # [kW] Initializing power input series
            for i in range(len(self.Current)):
                
                pot=(self.Current[i]*self.Voltage[i])/1000  # [kW]
                Power_inp.append(pot) 

            # PowerNominal = max(Power_inp)              # [kW] Max Power Consumption of Electrolyzer computed
            # print(PowerNominal)                        # declared nominal power check
                    
            # Interpolation
    
            self.PI=interp1d(Power_inp,self.Current,bounds_error=False,fill_value='extrapolate')        # Linear spline 1-D interpolation
            
            'Single module H2 production'

            hydrogen = []
            etaFar   = []
            
            for i in range(len(Power_inp)):
               
                'Electrolyzer efficiency and hydrogen output'    
                      
                etaFaraday = 9.95e-1*self.eNepero**((-9.5788-0.0555*self.OperatingTemp)/(self.Current[i]/(self.CellArea*self.nc/10000))+ \
                            (1502.7083-70.8005*self.OperatingTemp)/((self.Current[i]/(self.CellArea*self.nc/10000))**2))                   # [-] Faraday efficiency
                
                etaFar.append(etaFaraday)       
                 
                'Hydrogen Production'
                      
                HydroProdMol  = (etaFaraday*self.nc*self.Current[i]*3600)/(2*self.FaradayConst)     # [mol/h] (Guilbert 2020)  
                HydroMol      = HydroProdMol*self.timestep                                          # [mol] necessary for implementing compressor model, to check if its fundamental or can be removed
                hyd           = HydroMol*self.H2MolMass                                             # [kg] hydrogen produced in the considered timestep
                hydrogen.append(hyd)     
            
            self.maxh2prod = max(hydrogen)                                                          # [kg] maximum amount of produced hydrogen for the considered module
          
            
            self.etaF = interp1d(hydrogen,etaFar,bounds_error=False,fill_value='extrapolate')       # Linear spline 1-D interpolation -> produced H2 - Faraday efficiency
            self.h2P  = interp1d(hydrogen,Power_inp,bounds_error=False,fill_value='extrapolate')    # Linear spline 1-D interpolation -> produced H2 - Power consumption 
                        
            
    def h2power(self,h2):
        """
        Inverse function that computes the power consumption and Faraday efficiency
        corresponding to a required amount of produced hydrogen by means of interpolating functions.

        Parameters
        ----------
        h2 : float hydrogen to be produced in the timestep [kg]

        Returns
        -------
        hyd: float hydrogen produced in the timestep [kg] (same as input 'h2')
        e_absorbed : float electricity absorbed in the timestep [kWh]
        etaElectr : module efficiency [-]
        watCons : module water consumption [m^3]
        CellCurrDensity1 : single cell current density [A/cm^2]
        etaFaraday : Faraday efficiency [-]

        """
        if 0 < h2 <= self.maxh2prod:
            
            hyd        = h2
            e_absorbed = self.h2P(hyd)                                                                # [kW] required power
            etaF       = self.etaF(hyd)                                                               # [-]  Faraday efficiency
            CellCurrDensity1 = self.PI(e_absorbed)/self.CellArea                                      # [A/cm^2] Cell working current density
            Current   = CellCurrDensity1*self.CellArea                                                # [A] Stack operating current  
            Vstack    = (e_absorbed/Current)*1000                                                     # [V] Stack operating voltage
            etaElectr = self.nc*self.LHVh2*1e6*self.H2MolMass*etaF/(2*Vstack*self.FaradayConst)       # [-] Electrolyzer efficiency
            HydroMol  = (etaF*self.nc*Current*3600)/(2*self.FaradayConst)*self.timestep               # [mol/h] (Guilbert 2020)  
            hyd_vol   = HydroMol*self.H2MolMass/self.rhoNrh2                                          # [Nm^3] hydrogen produced in the considered timestep   
            watCons   = hyd_vol*self.rhoNrh2*self.h2oMolMass/self.H2MolMass/etaElectr/self.rhoStdh2o  # [m^3] water used by the electrolyzer - volume calculated @ 15°C & Pamb  
        
        return(hyd,e_absorbed,etaElectr,watCons,CellCurrDensity1,etaF)
        
             
    def plot_polarizationpts(self):
        
        if self.model == 'PEM General': 
            # i-V plot
            x = np.linspace(min(self.CellCurrDensity),max(self.CellCurrDensity),self.num) 
            plt.figure(dpi=1000)
            plt.plot(self.CellCurrDensity,self.Voltage,label='data', color='b',marker='.', linestyle='None', mec='r', markersize=7, markerfacecolor='white', zorder=0)
            plt.axvline(x=0.05,color='b',linestyle=':', label= 'Lower Functioning Boundary', zorder=3, linewidth=1.2)   # Ref. doi: 10.1016/j.ijhydene.2008.11.083
            plt.plot(x,self.iV1(x),label='linear', linestyle='--') 
            plt.grid()
            plt.legend(fontsize=8)
            plt.xlabel('Cell Current Density [A cm$^{-2}$]')
            plt.ylabel('Stak Voltage [V]')
            plt.title('PEMEL Polarization Curve - SplineInterp' )
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
 
    
    # Computing Electrolyzers performances via Spline Interpolation  
    def use(self,h,e,storable_hydrogen,TankMaxCapacity):
        """
        Electorlyzers stack and single modules operational parameters
    
        h : timestep float timestep in hours [h]
        e : float > 0 electricity provided to the electrolyzer in one hour [kWh]
        storable_hydrogen : float storable hydrogen H tank max_capacity - SOC[h-1] [kg] or maximum absorbable production if there is no tank
        TankMaxCapacity: tank size [kg]

        output : 
            
        hyd: float hydrogen produced in the timestep [kg]    
        e_absorbed: float electricity absorbed in the timestep [kWh]
        
        """
       
        if self.model != 'PEM General':
            
            e_absorbed = min(self.Npower,e)                          # [kWh] when timestep is kept at 1 h kWh = kW
            hyd = self.production_rate*(e_absorbed / self.Npower)    # [kg/h] (e_absorbed / self.Npower) represents the fraction of the maximum possible 
                                                                     # amount of produced hydrogen
            
            if hyd > storable_hydrogen:        # if the is not enough space in the H tank to store the hydrogen (H tank is nearly full)
                hyd = 0
                e_absorbed = 0
                # turn off the electrolyzer
                # this behavior could be solved with more advanced models, necessary inverse production functions.
                
            return(hyd,-e_absorbed) # return hydrogen supplied and electricity absorbed
                    
                
        elif self.model == 'PEM General':   
            
            'Defining the working point of the electrolyzer by spline interpolation:'
            
            # PowerInput [kW] - Electric Power from renewables directed to the electrolyzer
            if e <= self.Npower:                           # if available power is lower than single module Nominal Power
                
                e_absorbed=e
                hyd,e_absorbed,etaElectr,watCons,CellCurrden = electrolyzer.use1(self,e,storable_hydrogen,TankMaxCapacity)
                if hyd != 0:
                    self.n_modules_used[h]=1
                else:
                    self.n_modules_used[h]=0
                self.EFF[h]           = etaElectr
                self.cell_currdens[h] = CellCurrden
                self.wat_cons[h]      = watCons
                
            if self.Npower < e <= self.MaxPowerStack:      # if available power is between single module power and stack power
                
                n_modules_used=int(e/self.Npower)
                e_absorbed=self.Npower                     # power absorbed by the single module          
                hyd,e_absorbed,etaElectr,watCons,CellCurrden = electrolyzer.use1(self,e_absorbed,storable_hydrogen,TankMaxCapacity) 
                hyd_1=hyd*n_modules_used                   # total amount of H2 produced by modules working at full load
                e_absorbed_1=e_absorbed*n_modules_used     # total power absorbed    // // // // // 
                watCons_1=watCons*n_modules_used           # total water consumption // // // // // 
                self.EFF[h] = etaElectr                    # work efficiency of modules working at nominal power 
                self.cell_currdens[h] = CellCurrden
             
                e_remained=e-self.Npower*n_modules_used    # remaining power after considering modules at full load
                hyd,e_absorbed,etaElectr,watCons,CellCurrden = electrolyzer.use1(self,e_remained,storable_hydrogen,TankMaxCapacity)
               
                if -e_absorbed>0:                          # if remaining power is higher than 10% of Npower, last module working at partial load is activated
                    n_modules_used=n_modules_used+1        # considering the module working at partial load
                    self.EFF_last_module[h]=etaElectr      # work efficiency of the last module working with the remaining power
                    self.wat_cons_last_module[h]=watCons   # water consumption // // // // // 
                self.n_modules_used[h]=n_modules_used
                self.wat_cons[h]=watCons_1+watCons
                hyd=hyd+hyd_1
                e_absorbed=e_absorbed_1+e_absorbed
                
            if e > self.MaxPowerStack:                     # if available power greater than stack power
               
                self.n_modules_used[h]=self.n_modules
                e_absorbed=self.Npower
                hyd,e_absorbed,etaElectr,watCons,CellCurrden = electrolyzer.use1(self,e_absorbed,storable_hydrogen,TankMaxCapacity)
                hyd=hyd*self.n_modules
                e_absorbed=e_absorbed*self.n_modules
                self.EFF[h] = etaElectr                
                self.wat_cons[h]=watCons*self.n_modules
                self.cell_currdens[h] = CellCurrden
                
        return (hyd,-e_absorbed)
        
               
    def use1(self,e,storable_hydrogen,TankMaxCapacity):
        """
        This function calculates the performances of a single electrolyzer module.
        

        Parameters
        ----------
        e : float > 0 electricity provided to the electrolyzer in one hour [kWh]
        storable_hydrogen : float storable hydrogen H tank max_capacity - SOC[h-1] [kg] or maximum absorbable production if there is no tank
        TankMaxCapacity: tank size [kg] 

        Returns
        -------
        hyd: float hydrogen produced in the timestep [kg]    
        e_absorbed: float electricity absorbed in the timestep [kWh]
        etaElectr : module efficiency [-]
        watCons : module water consumption [m^3]
        CellCurrDensity1 : single cell current density [A/cm^2]

        """
        
        e_absorbed        = e                                   # [kWh] when timestep is kept at 1 h kWh = kW
        CellCurrDensity1  = self.PI(e_absorbed)/self.CellArea   # [A/cm^2]
        MinInputPower     = 0.1*self.Npower                     # [kW] minimum input power chosen as 10% of module nominal power 

        # Checking if resulting current density is high enough for the electrolyzer to start, otherwise hydrogen prod = 0
    
        if CellCurrDensity1 < self.CurrDensityMin or e_absorbed < MinInputPower:      # condition for operability set for current density and input power 
                                                                                      
            deltaHydrogen = 0      # [kWh] absorbed energy to produce hydrogen
            etaElectr     = 0      # [-] electrolyzer efficiency
            Vstack        = 0      # [V] stack voltage 
            HydroMol      = 0      # [mol] produced moles of hydrogen
            hyd_vol       = 0      # [Nm^3] hydrogen produced in the considered timestep  
            hyd           = 0      # [kg] hydrogen produced in the considered timestep
            etaFaraday    = 0      # [-] Faraday efficiency
            Current       = 0      # [A] Operational Current
            e_absorbed    = 0      # [kWh] absorbed energy - when timestep is kept at 1 h kWh = kW
            watCons       = 0      # [m^3] water volumetric consumption
        
        else:     
            
            Current = self.PI(e_absorbed)              # [A] Stack operating current                     
            Vstack  = (e_absorbed/Current)*1000        # [V] Stack operating voltage
                    
            'Calculating electrolyzer efficiency and hydrogen energy output for the given inputs'    
                  
            etaFaraday = 9.95e-1*self.eNepero**((-9.5788-0.0555*self.OperatingTemp)/(Current/(self.CellArea*self.nc/10000))+ \
                        (1502.7083-70.8005*self.OperatingTemp)/((Current/(self.CellArea*self.nc/10000))**2))                   # [-] Faraday efficiency
                
            etaElectr = self.nc*self.LHVh2*1e6*self.H2MolMass*etaFaraday/(2*Vstack*self.FaradayConst)                          # [-] Electrolyzer efficiency
                
            'Hydrogen Production'
                  
            HydroProdMol  = (etaFaraday*self.nc*Current*3600)/(2*self.FaradayConst)      # [mol/h] (Guilbert 2020)  
            HydroMol      = HydroProdMol*self.timestep                                   # [mol] necessary for implementing compressor model, to check if its fundamental or can be removed
            hyd           = HydroMol*self.H2MolMass                                      # [kg] hydrogen produced in the considered timestep
            hyd_vol       = HydroMol*self.H2MolMass/self.rhoNrh2                         # [Nm^3] hydrogen produced in the considered timestep          
            deltaHydrogen = hyd_vol*self.LHVh2*self.rhoNrh2*(1000/3600)                  # [kWh] Energy produced, in the form of hydrogen 
                  
            'Water consumption' 
            
            watCons = hyd_vol*self.rhoNrh2*self.h2oMolMass/self.H2MolMass/etaElectr/self.rhoStdh2o     # [m^3] water used by the electrolyzer - volume calculated @ 15°C & Pamb           
        
        if hyd > storable_hydrogen:           # if there is not enough space in the H tank to store the hydrogen (H tank is nearly full)
        
            # defining the power requirements to produce storable hydrogen amount 
            
            hyd,e_absorbed,etaElectr,watCons,CellCurrDensity1,etaFaraday = electrolyzer.h2power(self,storable_hydrogen)
     
            if e_absorbed >=  MinInputPower : 
                
                pass
            
            else:
                
                hyd               = 0
                e_absorbed        = 0
                etaElectr         = 0
                watCons           = 0
                CellCurrDensity1  = 0
       
        return (hyd,-e_absorbed,etaElectr,watCons,CellCurrDensity1)            
           
##########################################################################################

if __name__ == "__main__":
    
    """
    Functional test
    """
    
    inp_test = {  
                  "Npower": 1000,
                  "number of modules": 80,
                  "stack model": "PEM General"
                }
    
    sim_hours = 24                              # [h] simulated period of time - usually it's 1 year minimum

    el = electrolyzer(inp_test,sim_hours)       # creating electrolyzer object
    el.plot_polarizationpts()                   # plot example

    storable_hydrogen = 200                     # [kg] Available capacity in tank for H2 storage at hour 'h'
    tank_maxcapacity  = 2000                    # [kg] Tank max capacity
    
    'Test 1 - Tailored ascending power input'

    flow = np.linspace(0,el.Npower,sim_hours)   # [kW] power input - ascending series

    hyd          = np.zeros(len(flow))          # [kg] produced hydrogen
    eabsorbed    = np.zeros(len(flow))          # [kW] absorbed hydrogen
    
    for i in range(len(flow)):
        
        hyd[i],eabsorbed[i] = el.use(i,flow[i],storable_hydrogen,tank_maxcapacity)
            
    print('Produced Hydrogen in the timestep [kg]: \n\n',hyd)
    print('\nAbsorbed energy [kWh]: \n\n',eabsorbed)     
    print('\nElectrolyzer Efficiency [-]: \n\n',el.EFF)                # electrolyzer efficiency at every hour     

    cellarea = el.CellArea
    nompower = el.Npower

    plt.figure(dpi=600)
    plt.plot(flow,el.EFF)
    plt.title("Electrolyzer Module Efficiency")
    # plt.ylim([0,0.8]) 
    textstr = '\n'.join((
        r'$CellArea=%.1f$ $cm^{2}$' % (cellarea,),
        r'$P_{nom}= %.1f$ kW' % (nompower,),
        r'$i_{max}= %.1f$ A $cm^{-2}$' % (el.CurrDensityMax,),
        r'$n_{cell}= %.0f$' % (el.nc,)))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    plt.text(max(flow)/2,0.2,textstr,fontsize=10,va='bottom',backgroundcolor='none', bbox=props)
    plt.grid()
    plt.xlabel('Input Power [kW]')
    plt.ylabel('$\\eta$')    
    
    plt.figure(dpi=600)
    plt.plot(el.cell_currdens,el.EFF,color='tab:orange')
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
    
    
    'Test 2 - Random power input'

    flow = np.random.uniform(0.08*el.Npower,1.2*el.Npower,sim_hours)   # [kWh] randomic power input as example
    
    hyd          = np.zeros(len(flow))          # [kg] produced hydrogen
    eabsorbed    = np.zeros(len(flow))          # [kW] absorbed hydrogen
    
    for i in range(len(flow)):
        
        hyd[i],eabsorbed[i] = el.use(i,flow[i],storable_hydrogen,tank_maxcapacity)
        
    fig, ax = plt.subplots(dpi=600)
    ax2 = ax.twinx() 
    ax.bar(np.arange(sim_hours),el.EFF, width = 0.4, zorder = 3, align = 'edge')
    ax2.bar(np.arange(sim_hours),flow, width = 0.4, zorder = 3, align = 'center', color ='tab:orange')
    ax.set_xlabel('Time [h]')
    ax.set_ylabel('Efficiency [h]')
    ax2.set_ylabel('Power Input [kW]')
    ax.grid()
    ax.set_title('Electrolyzer functioning behaviour')
