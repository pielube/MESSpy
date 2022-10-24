import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
# from CoolProp.CoolProp import PropsSI
import numpy as np
from sklearn.linear_model import LinearRegression
from numpy import log as ln


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
        
        if parameters['Npower']:
            self.Npower = parameters['Npower'] # float nominal power of electrolyzers installed capacity for the location [kW]
        
        H_N_density = 0.08988237638480538  # PropsSI('D', 'T', 273.15, 'P', 101325, 'H2') # hydrogen density under normal condition

        if parameters['stack model'] == 'Enapter 2.1': # https://www.enapter.com/it/newsroom/enapter-reveals-new-electrolyser-el-2-1
            self.model= parameters['stack model'] 
            stack_operative_power_consumption = 2.4    # [kW] Single module nominal power
            stack_production_rate_L = 500              # [NL/h]
            stack_production_rate = stack_production_rate_L / 1000 * H_N_density # [kg/h]
            self.production_rate = stack_production_rate * (self.Npower / stack_operative_power_consumption) # [kg/h] actual production defined in a simplified way 
                                                                                                             # by dividing the selected electorlyzer size in the location  
                                                                                                             # for the single unit nominal power      

        if parameters['stack model'] == 'McLyzer 800': # https://mcphy.com/it/apparecchiature-e-servizi/elettrolizzatori/large/
            self.model= parameters['stack model']      
            stack_operative_power_consumption = 4000   # [kW] Single module nominal power
            stack_production_rate_m3 = 800             # [Nm3/hr]
            stack_production_rate = stack_production_rate_m3 * H_N_density   # [kg/hr]
            self.production_rate = stack_production_rate * (self.Npower / stack_operative_power_consumption) # [kg/h] actual production defined in a simplified way 
                                                                                                             # by dividing the selected electorlyzer size in the location                                                                                                     # for the single unit nominal power      
    
        if parameters['stack model'] == 'PEM General':
            
            self.EFF = np.zeros(simulation_hours)   # so far - easiest way of keeping track of the elecrolyzer efficiency over the simulation
            
            self.model= parameters['stack model']
            self.rhoNrh2      = 0.088707      # [kg/Nm3]    H2  density @ T = 0°C p = 101325 Pa  PropsSI('D', 'T', 273.15, 'P', 1*1e5, 'H2')
            self.rhoStdh2o    = 999.06        # [kg/m3]     H2O density @ T = 15°C p = 101325 Pa
            Runiv             = 8.3144621     # [J/(mol*K)] Molar ideal gas constant
            self.FaradayConst = 96485         # [C/mol]     Faraday constant
            self.LHVh2        = 119.96        # [MJ/kg]     H2 LHV
            self.HHVh2Mol     = 285.83        # [kJ/mol]    H2 HHV molar
            self.h2oMolMass   = 0.01801528    # [kg/mol]    Water molar mass
            self.H2MolMass    = 2.01588e-3    # [kg/mol]    Hydrogen molar mass

            # Math costants
            self.eNepero    = 2.71828182845904523536

            # Ambient conditions 
            self.AmbTemp = 288             # [K] temperature ofthe surrounding environment

            # At current development stage it is taken for granted we are working with hourly balances 
            self.timestep = 1              # [h]

                
            "2H2O --> 2H2 + O2"
            # https://doi.org/10.1016/j.ijhydene.2008.11.083    # Electrolyzer Parameters - depending on the different types of electrolyzers chosen - model to be specified also in class - name
        
            # self.nc                  = parameters['CellNumber']           # number of cells in the stack
            # self.MembThickness       = parameters['MembThickness']        # [micron]
            # self.Lambda              = parameters['Lambda']               # [-] Cell mositure content
            # self.AnodeCurrDensity    = parameters['AnodeCurrDensity']     # [A/cm^2]
            # self.AnodePressure       = parameters['AnodePressure']        # [Pa]
            # self.CathodeCurrDensity  = parameters['CathodeCurrDensity']   # [A/cm^2]
            # self.CTCanode            = parameters['CTCanode']             # [-] Charge transfer coefficient - Anode   https://www.sciencedirect.com/science/article/pii/S0360319918309017
            # self.CTCcathode          = parameters['CTCcathode']           # [-] Charge transfer coefficient - Cathode //      //      //
            # self.OperatingTemp       = parameters['OperatingTemp']        # [K]
            # self.OperatingPress      = parameters['OperatingPress']       # [Pa]
            # self.CurrDensityMax      = parameters['CurrDensityMax']       # [A/cm^2]
            # self.CurrDensityMin      = parameters['CurrDensityMin']       # [A/cm^2]
            # self.CellArea            = parameters['CellArea']             # [cm^2] cell active area 
            # self.PowerNominal        = 5           # [kW] Max power absorbable - this value is now defined through interpolation 

            self.nc                  = 10            # number of cells in the stack
            self.MembThickness       = 250           # [micron]
            self.Lambda              = 20            # [-] Cell mositure content
            self.AnodeCurrDensity    = 0.00013       # [A/cm^2]
            self.AnodePressure       = 101325        # [Pa]
            self.CathodeCurrDensity  = 0.001         # [A/cm^2]
            self.CTCanode            = 0.6           # [-] Charge transfer coefficient - Anode   https://www.sciencedirect.com/science/article/pii/S0360319918309017
            self.CTCcathode          = 0.45          # [-] Charge transfer coefficient - Cathode //      //      //
            self.OperatingTemp       = 273.15 + 60   # [K]
            self.OperatingPress      = 5000000       # [Pa]
            self.CurrDensityMax      = 2             # [A/cm^2]
            self.CurrDensityMin      = 0.005         # [A/cm^2]
            self.CellArea            = 100           # [cm^2] cell active area 

            # Computing the electrolyzer polarization curve V-i
            
            'POLARIZATION CURVE'
            
            Ndatapoints = 101                      # Number of points used to compute the polarization curve 
    
            self.CellCurrDensity= np.linspace(self.CurrDensityMin,self.CurrDensityMax,Ndatapoints)
            self.Voltage=np.zeros(Ndatapoints)
        
    
            for i in range(0,Ndatapoints):
                
                'Polarization (V-i) curve calculation'
                # V is obtained by summing 4 different contributions 
                
                '1- Cell open curcuit voltage'
                
                pH2O = self.eNepero**(11.676-(3816.44 /(self.OperatingTemp-46.13)))  # [atm]
                pO2  = self.AnodePressure/101325 - pH2O                              # [atm]
                pH2  = self.OperatingPress/101325-pH2O                               # [atm]
                
                Ecell = 1.229-0.85e-3*(self.OperatingTemp-298.15)+4.3085e-5*self.OperatingTemp*ln(pH2*(pO2**0.5)/pH2O) #[V] Cell open curcuit voltage
                
                '2- Cell activation overpotential'
                
                Vact_an  = (Runiv*self.OperatingTemp*ln(self.CellCurrDensity[i]/self.AnodeCurrDensity))/(2*self.CTCanode*self.FaradayConst)     # [V]
                Vact_cat = Runiv*self.OperatingTemp*ln(self.CellCurrDensity[i]/self.CathodeCurrDensity)/(4*self.CTCcathode*self.FaradayConst)   # [V]
                  
                Vact = Vact_an + Vact_cat                                      # [V] Cell activation overpotential
            
                '3- Cell mass transport overpotential'
                
                MaxCurrDensity_2 = self.CurrDensityMax + 0.0001
                  
                Vdiff = -Runiv*self.OperatingTemp*ln(1-self.CellCurrDensity[i]/MaxCurrDensity_2)/(2*self.FaradayConst)  # [V] Cell mass transport overpotential
            
                '4- Cell Ohmic losses'
                  
                MembConductivity = (0.005139*self.Lambda-0.00326)*self.eNepero**(1268*(1/303 - 1/self.OperatingTemp))   # [S/cm]
                Rcell = (self.MembThickness/10000)/MembConductivity                                                # [cm^2/S]
                  
                Vohmic = self.CellCurrDensity[i]*Rcell                 # [V] Cell Ohmic losses
                
                'Resulting stack voltage'
                      
                self.Voltage[i] = self.nc*(Ecell+Vact+Vdiff+Vohmic) # [V] - Stack Voltage
            
                
            'Interpolation of calculated functioning points to detect the best fit-function for i-V curve'
             
            self.num = Ndatapoints                                                   # Number of intervals to be considered for the interpolation
            self.x2 = np.linspace(0.05,max(self.CellCurrDensity),self.num)        # Setting xlim for range of validity of LinRegression Calculation - Only for plot-related reasons 
            
            # Interpolation
            self.iV1 = interp1d(self.CellCurrDensity,self.Voltage,bounds_error=False)#,fill_value='extrapolate')                # Linear spline 1-D interpolation
           
            # Creating the reverse curve IP which will be necessary to define the exact functioning point
            self.Current = self.CellCurrDensity*self.CellArea          # Defining the current value: same both for the single cell and the full stack!
           
            # Defining Electrolyzer Max Power Consumption
            Power_inp = []
            for i in range(len(self.Current)):
                pot=(self.Current[i]*self.Voltage[i])/1000  # [kW]
                Power_inp.append(pot)
            
            self.PowerNominal = max(Power_inp)  # [kW] Max Power Consumption of Electrolyzer - Nominal Power   
            print('Nominal Power Electorlyzer = {}'.format(round(self.PowerNominal,1)))
            self.IP=interp1d(self.Current,Power_inp,bounds_error=False,fill_value='extrapolate')  
                
            P = []
            for i in range (len(self.Current)):
                power = (self.iV1(self.CellCurrDensity[i])*self.Current[i])/1000   # [kW]
                if power > self.PowerNominal:
                    break
                P.append(power) 
            
            # Interpolation
          
            #self.PI=interp1d(P,self.Current,kind='cubic',bounds_error=False,fill_value='extrapolate')  # Cubic spline 1-D interpolation
            self.PI=interp1d(P,self.Current,bounds_error=False,fill_value='extrapolate')                # Linear spline 1-D interpolation
            #self.PI=PchipInterpolator(P,self.Current)                                                  # PCHIP 1-D monotonic cubic interpolation
                 
            # interval = np.linspace(min(P),max(P),num)      # Current interval - useful for self.PI interpolation
            # #print('INTERVALLO---------------------------------------------', interval)
            # plt.figure(dpi=1000)
            # plt.plot(interval,self.PI(interval))
            # plt.title('Curva P-I prova')
    
    # Step eventuale successivo definire tali funzioni secondo il principio di ereditarietà delle classi 
    # (https://www.programmareinpython.it/video-corso-python-programmazione-a-oggetti/03-ereditarieta/)
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
    def use(self,h,e,storable_hydrogen):
        """
        Produce hydrogen
    
        e : float > 0 electricity provided to the electrolyzer in one hour [kWh]
        storable_hydrogen : float storable hydrogen H tank max_capacity - SOC[h-1] [kg] or maximum absorbable production if there is no tank
        h : timestep float timestep in hours [h]

        output : 
        float hydrogen produced that hour [kg]    
        float electricity absorbed that hour [kWh]
        """
       
        if self.model == 'Enapter 2.1' or self.model == 'McLyzer 800':
            
            e_absorbed = min(self.Npower,e)    # [kWh] when timestep is kept at 1 h kWh = kW
            hyd = self.production_rate * (e_absorbed / self.Npower)  # [kg/h] (e_absorbed / self.Npower) represents the fraction of the maximum possible 
                                                                     # amount of produced hydrogen
            
            if hyd > storable_hydrogen:        # if the is not enough space in the H tank to store the hydrogen (H tank is nearly full)
                hyd = 0
                e_absorbed = 0
                # turn off the electrolyzer
                # this behavior could be solved with more advanced models, necessary inverse production functions.
                
            return(hyd,-e_absorbed) # return hydrogen supplied and electricity absorbed
                    
                
        elif self.model == 'PEM General':    
            
            'Defining the working point of the electrolyzer by spline interpolation:'
               # PowerInput=CellCurrDensity*CellActiveArea*Vstack
               # Vstack=PowerInput/CellCurrDensity
            
            # PowerInput [kW] - Electric Power from renewables directed to the electrolyzer
            e_absorbed = min(self.PowerNominal,e)
  
            
            CellCurrDensity1 = self.PI(e_absorbed)/self.CellArea # [A/cm^2]
            
            # Checking if resulting current density is high enough for the electrolyzer to start, otherwise hydrogen prod = 0
            
            if CellCurrDensity1 < self.CurrDensityMin:
                
                deltaHydrogen = 0      # [kWh] absorbed energy to produce hydrogen
                etaElectr     = 0      # [-] electrolyzer efficiency
                WatCons       = 0      # [m^3] water volumetric consumption
                Vstack        = 0      # [V] stack voltage 
                HydroMol      = 0      # [mol] produced moles of hydrogen
                hyd_vol       = 0      # [Nm^3] hydrogen produced in the considered timestep  
                hyd           = 0      # [kg] hydrogen produced in the considered timestep
                etaFaraday    = 0      # [-] Faraday efficiency
                Current       = 0      # [A] Operational Current
                
            else:     
                
                Current = CellCurrDensity1*self.CellArea  # [A] Stack operating current  
    
                Vstack = (e_absorbed/Current)*1000        # [V] Stack operating voltage
                
                'Calculating electrolyzer efficiency and hydrogen energy output for the given inputs'    
              
                etaFaraday = 9.95e-1*self.eNepero**((-9.5788-0.0555*self.AmbTemp)/(Current/(self.CellArea*self.nc/10000))+ \
                             (1502.7083-70.8005*self.AmbTemp)/((Current/(self.CellArea*self.nc/10000))**2))                   # [-] Faraday efficiency
            
                etaElectr = self.nc*self.LHVh2*1e6*self.H2MolMass*etaFaraday/(2*Vstack*self.FaradayConst)                     # [-] Electrolyzer efficiency
                self.EFF[h] = etaElectr
            
                'Hydrogen Production'
              
                HydroProdMol  = (etaFaraday*self.nc*Current*3600)/(2*self.FaradayConst)      # [mol/h] (Guilbert 2020)  
                HydroMol      = HydroProdMol*self.timestep                                   # [mol] necessary for implementing compressor model, to check if its fundamental or can be removed
                hyd           = HydroMol*self.H2MolMass                                      # [kg] hydrogen produced in the considered timestep
                hyd_vol       = HydroMol*self.H2MolMass/self.rhoNrh2                         # [Nm^3] hydrogen produced in the considered timestep          
                deltaHydrogen = hyd_vol*self.LHVh2*self.rhoNrh2*(1000/3600)                  # [kWh] Energy produced, in the form of hydrogen 
               
                'Water consumption' 
                WatCons = hyd_vol*self.rhoNrh2*self.h2oMolMass/self.H2MolMass/etaElectr/self.rhoStdh2o      # [m^3] water used by the electrolyzer - volume calculated @ 15°C & Pamb           
            if hyd > storable_hydrogen:  # if the is not enough space in the H tank to store the hydrogen (H tank is nearly full)
                hyd = 0 
                e_absorbed = 0
               
            
            return (hyd,-e_absorbed)     

##########################################################################################

if __name__ == "__main__":
    
    """
    Functional test
    """
    
    inp_test = {  "Npower": False,
                  "stack model": "PEM General",
                  }
    
    sim_hours = 24                          # [h] simulated period of time - usually it's 1 year minimum

    el = electrolyzer(inp_test,sim_hours)   # creating electrolyzer object
    el.plot_polarizationpts()               # plot example

    flow = np.random.uniform(1,10,sim_hours)
    
    for h in range(24):
        print(el.use(h,flow[h],99999))
        
        
    print(el.EFF)                           # electrolyzer efficiency at every hour
