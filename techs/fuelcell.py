import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
import math
from numpy import log as ln
from sklearn.linear_model import LinearRegression
import os
import sys 
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),os.path.pardir)))   # temorarily adding constants module path 
import constants as c

class fuel_cell:
    
    def __init__(self,parameters,simulation_hours):
        """
        Create a Fuel Cell object
    
        parameters : dictionary
            'Npower': float nominal power [kW]
                      
        output : Fuel cell object able to:
            absosrb hydrogen and produce electricity .use(e)
        """
        
       
        self.Npower = parameters['Npower'] * 1000 # W
        
        # this model is based on FCS-C5000 characteristic curves https://www.horizonfuelcell.com/hseries
        if parameters['stack model']=='FCS-C5000':
            self.model=parameters['stack model']
            self.nc = 120 * self.Npower / 5000 # # number of cells (120 cells are required to have a nominal power of 5000 W)
        
            # characteristic curves
            I=[0, 10, 20, 30, 40, 50, 60, 70, 80]
            v=[0.96, 0.83, 0.796, 0.762, 0.728, 0.694, 0.66, 0.6, 0.5]
            V=[]
            P=[]
            for i in range(len(I)):
                volt=self.nc*v[i]
                V.append(volt)
                pot=I[i]*volt
                P.append(pot)
            
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
        
        if parameters['stack model']=='PEM General':
            
            'Model: NEDSTACK FCS 13-XXL' # https://nedstack.com/sites/default/files/2022-07/nedstack-fcs-13-xxl-gen-2.9-datasheet-rev01.pdf
            
            self.EFF       = np.zeros(simulation_hours)   # [-]    Keeping track of fuel cell efficiency
            self.VOLT      = np.zeros(simulation_hours)   # [V]    Keeping track of single cell working voltage - necessary for ageing calculations
            self.CURR_DENS = np.zeros(simulation_hours)   # [A]    Keeping track of single cell working current - necessary for ageing calculations
            
            self.model=parameters['stack model']
            
            self.rhoStdh2        = c.H2SDENSITY     # [kg/Sm3]    PropsSI('D', 'T', 288.15, 'P', 101325, 'H2') H2  density @ T = 15°C p = 101325 Pa
            self.rhoStdh2o       = c.H2OSDENSITY    # [kg/m3]     H2O density @ T = 15°C p = 101325 Pa
            self.Runiv           = c.R_UNIVERSAL    # [J/(mol*K)]
            self.Rh2             = c.R_H2           # [J/(kg*K)] 
            self.FaradayConst    = c.FARADAY        # [C/mol]     Faraday constant
            self.deltaG0         = c.GIBBS          # [kJ/mol]    Gibbs free energy @ T = 25°C p = 101325 Pa
            self.GammaPerfectGas = c.GAMMA          # [-]         Gamma = cp/cv  
            self.LHVh2           = c.LHVH2          # [MJ/kg]     H2 LHV
            self.HHVh2           = c.HHVH2          # [MJ/kg]     H2 HHV
            self.HHVh2Mol        = c.HHVH2MOL       # [kJ/mol]    H2 HHV molar
            self.cpH2O           = c.CP_WATER       # [kJ/(kgK)]  Water specific heat
            self.h2oMolMass      = c.H2OMOLMASS     # [kg/mol]    Water molar mass
            self.H2MolMass       = c.H2MOLMASS      # [kg/mol]    Hydrogen molar mass

            # Math costants
            self.eNepero      = c.NEPERO         # [-]         Euler's number

            # Ambient conditions 
            self.AmbTemp      = c.AMBTEMP        # [K]         Standard ambient temperature - 15 °C

            # At current development stage it is taken for granted we are working with hourly balances 

            self.timestep = 1 # [h]

            
            "H2 --> 2H+ + 2e" 
            
            # FuelCellctrolyzer Parameters - depending on different technologies and sizes 
            
            self.nc                  = 96                                  # [-] number of cells in the stack 
            self.Lambda              = 23                                  # [-] Cell mositure content
            self.FC_AnodeCurrDens    = 0.000000009                         # [A/cm^2] Anode current density
            self.FC_CathodeCurrDens  = 0.001                               # [A/cm^2] Cathode current density
            self.FC_OperatingTemp    = 273.15 + 62                         # [K]   https://www.google.com/search?q=working+temperature+PEM+fuel+cell&oq=working+temperature+PEM+fuel+cell+&aqs=chrome..69i57j0i512l9.10383j0j7&sourceid=chrome&ie=UTF-8
            self.FC_FuelPress        = 101325 + 25000                      # [Pa] Fuel supply pressure (Anode) 
            self.FC_AirPress         = 101325                              # [Pa] Air supply pressure (Cathode)
            self.FC_MinCurrDens      = 0.0001                              # [A/cm^2] FC min. current density      
            self.FC_MaxCurrDens      = 1.2                                 # [A/cm^2] FC max. current density      
            self.Vmin_FC             = 59                                  # [V] Minimum value for working voltage - derived from datasheet
            self.CTC                 = 0.4                                 # [-] Charge Transfer Coefficient  - if Nominal Power < 6 kW: CTC = 0.45
                                                                           #                                  - if //   //   //  > 6 kW: CTC = 0.4
            self.MembThickness       = 250                                 # !!! [μm] fuel cell membrane thickness - if Nominal Power < 6 kW: MembThickness = 100 
                                                                           #                                    - if //   //   //  > 6 kW: MembThickness = 145 
            self.NomPower            = 13600                               # [W] datasheet specification                                            
            self.FC_MaxCurrent       = 230                                 # [A] value taken from datasheet. Current value at which max power is delivered
            self.FC_CellArea         = self.NomPower/ \
                                      (self.Vmin_FC*self.FC_MaxCurrDens)   # [cm^2] FC cell active area - Not given by the constructor, has to be calculated as follows -> CellArea= FC_NominalPower/(Vmin*FC_MaxCurrDens)          
            
            
            # Computing the single cell polarization curve V-i
            
            'POLARIZATION CURVE'
            
            Ndatapoints= 101         # Number of points used to compute the polarization curve 
            
            self.OCpotential = []
            self.ActLosses   = []
            self.OhmLosses   = []
            self.CellVolt = np.zeros(Ndatapoints)
            self.Voltage  = np.zeros(Ndatapoints)
            self.CellCurrDensity = np.linspace(self.FC_MinCurrDens,self.FC_MaxCurrDens,Ndatapoints)   
            
            for i in range(0,Ndatapoints):
                
                'Polarization (V-i) curve calculation'
                
                # V is obtained by summing 3 different contributions (concentration losses are considered to be negligible)
                
                '1- Cell open curcuit voltage'
                
                pO2 = (self.FC_AirPress*0.21)/101325  # [atm] accounting for partial pressure of oxigen in air
                pH2O = 1                              # [atm]
                pH2 = self.FC_FuelPress/101325        # [atm]
            
                CellCurr = self.CellCurrDensity[i]*self.FC_CellArea   # [A] Current in the cell
            
                Ecell = 1.229 -0.85e-3*(self.FC_OperatingTemp-298.15) + 4.3085e-5*self.FC_OperatingTemp*ln(pH2*(pO2**0.5)/pH2O)  # [V] Open circuit voltage
                self.OCpotential.append(Ecell)
                
                '2- Activation losses'
              
                Vact_cat = -self.Runiv*self.FC_OperatingTemp*np.log10(self.FC_CathodeCurrDens)/(self.CTC*4*self.FaradayConst)+self.Runiv*self.FC_OperatingTemp*np.log10(self.CellCurrDensity[i])/(self.CTC*4*self.FaradayConst) #[V]
                Vact_an = -self.Runiv*self.FC_OperatingTemp*np.log10(self.FC_AnodeCurrDens)/(self.CTC*2*self.FaradayConst)+self.Runiv*self.FC_OperatingTemp*np.log10(self.CellCurrDensity[i])/(self.CTC*2.*self.FaradayConst)   #[V]
          
                Vact = Vact_cat +Vact_an              # [V] Activation losses
                
                self.ActLosses.append(Ecell-Vact)
           
                '3- Ohmic losses'
         
                rho_m = (181.6*(1+0.03*(CellCurr/self.FC_CellArea)+0.062*((self.FC_OperatingTemp/303)**2)*(CellCurr/self.FC_CellArea)**2.5))/ \
                    ((self.Lambda-0.634-3*(CellCurr/self.FC_CellArea))*self.eNepero**(4.18*(self.FC_OperatingTemp-303)/self.FC_OperatingTemp))                           #[Ohm*cm] specific membrane recistence
            
              
                Rm = rho_m*self.MembThickness/(10000*self.FC_CellArea)   # [Ohm]
         
                Vohm = CellCurr*Rm  # [V] Ohmic losses
                
                self.OhmLosses.append(Ecell-Vact-Vohm)
          
                
                'SINGLE CELL VOLTAGE'
                
                self.CellVolt[i] = (Ecell-Vact-Vohm)        # [V] Cell Voltage
            
                'RESULTING STACK VOLTAGE'
         
                self.Voltage[i] = self.CellVolt[i]*self.nc  # [V] Stack Voltage
                    
                
            'Interpolation of calculated polarization curve points: defining the fit-function for i-V curve'
            
            self.num = Ndatapoints                          # [-] number of intervals to be considered for the interpolation
            
            self.x = np.linspace(self.CellCurrDensity[0],self.CellCurrDensity[-1],self.num) 
            
            # Interpolating functions
            
            self.iV1 = interp1d(self.CellCurrDensity,self.Voltage,bounds_error=False,fill_value='extrapolate')          # Linear spline 1-D interpolation
            self.iV2 = interp1d(self.CellCurrDensity,self.CellVolt,bounds_error=False,fill_value='extrapolate')         # Linear spline 1-D interpolation
            
             # Creating the reverse curve IP - necessary to define the exact functioning point        
            
            self.Current = self.CellCurrDensity*self.FC_CellArea          # Defining the current value: same both for the single cell and the full stack!

            # Defining Fuel Cell Max Power Generation
            
            FC_power = []
            for i in range(len(self.Current)):
                pot=self.Current[i]*self.Voltage[i]/1000     # [kW] Power
                FC_power.append(pot)
                
            self.FC_NominalPower = max(FC_power)             # [kW] Max Power Consumption of Electrolyzer - Nominal Power at max current  
            # print('Nominal Power FuelCell = {}'.format(round(self.FC_NominalPower,1)))
            
            self.IP=interp1d(self.Current,FC_power,bounds_error=False,fill_value='extrapolate')  
        
            self.P = []
            for i in range(len(self.Current)):
                power = (self.iV1(self.CellCurrDensity[i])*self.Current[i])/1000          # [kW] Resolving the equation system via interpolation
                self.P.append(power)                                                      # [kW] Output power values varying current
                
            self.PI=interp1d(self.P,self.Current,bounds_error=False,fill_value='extrapolate')  # Interpolating function returning Current if interrogated with Power 
                                  
    def plot_polarizationpts(self):
         
          if self.model == 'PEM General': 
              
              fig=plt.figure(figsize=(10,8),dpi=1000)
              fig.suptitle("PEMFC Polarization Curve STACK - P ={}".format(round(self.FC_NominalPower,1)) +" kW")
               
              plt.plot(self.CellCurrDensity,self.Voltage,label='data', color='b',marker='.', linestyle='None', mec='r', markersize=7, markerfacecolor='white', zorder=0)
              plt.plot(self.x,self.iV1(self.x),label='linear interp',linestyle='--') 
              plt.grid()
              plt.legend(fontsize=8)
              plt.xlabel('Cell Current Density [A cm$^{-2}$]')
              plt.ylabel('Stak Voltage [V]')
              plt.title('PEMFC Polarization Curve (V-i)' )
              
              # fig=plt.figure(figsize=(8,8),dpi=1000)
              # fig.suptitle("PEMFC Polarization Curve STACK - P ={}".format(round(self.FC_NominalPower,1)) +" kW")
            
              # ax_1=fig.add_subplot(221)
              # ax_1.plot(self.CellCurrDensity,self.Voltage,label='data', color='b',marker='.', linestyle='None', mec='r', markersize=7, markerfacecolor='white', zorder=0)
              # ax_1.plot(self.x,self.iV1(self.x),label='linear interp',linestyle='--') 
              # ax_1.grid()
              # ax_1.legend(fontsize=8)
              # ax_1.set_xlabel('Cell Current Density [A cm$^{-2}$]')
              # ax_1.set_ylabel('Stak Voltage [V]')
              # ax_1.set_title('PEMFC Polarization Curve (V-i)' )
            
              # ax_2=fig.add_subplot(222)
              # ax_2.plot(self.Current,self.Voltage,label='data', color='b',marker='.', linestyle='None', mec='r', markersize=7, markerfacecolor='white', zorder=0)
              # #ax_2.plot(self.x,self.iV1(self.x),label='linear',linestyle='--') 
              # ax_2.grid()
              # ax_2.legend(fontsize=8)
              # ax_2.set_xlabel('Cell Current [A]')
              # ax_2.set_ylabel('Stak Voltage [V]')
              # ax_2.set_title('PEMFC Polarization Curve (V-I)' )
             
              # plt.tight_layout()
              # plt.show()
             
             
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
                
             
          else:
              
              print('Polarization curve not available')
             
    def plot_stackperformance(self):
        
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
       
    # Computing FuelCell performances via Spline Interpolation 
       
    def use(self,h,e,available_hyd):
        """
        The Fuel cell can absorb hydrogen and produce electricity: H2 --> 2H+ + 2e
    
        e: float < 0 energy required  [kWh]
        available_hyd: float available hydrogen H tank SOC[h-1] [kg]
      
        output : hydrogen absorbed and electricity supplied that hour [kg]
        """
        if self.model=='FCS-C5000':
            
            p_required = -e*1000  # kWh --> W
            
            p = min(p_required,self.Npower) # how much electricity can Fuel Cell absorb? [W]
            
            # find the operating point on the characteristic curves
            I=self.PI(p)
           
            # calculate the hydrogen consumed by each single cell
            qe=I/96485 #[mol_e/s] Faraday
            qH=qe*0.5*3600 #[mol_H2/h] the moles of H2 are 1/2 the moles of electrons. 3600s in one hour.
            QH=qH*2.058 #[g_H2/h] H2 molar mass = 2.058
            
            hyd = QH*self.nc/1000 # total stack hydrogen [kg/hr]
            
            if hyd > available_hyd: # if not enough hydrogen is available to meet demand (H tank is nearly empty)
                hyd = 0
                p = 0
                # turn off the fuel cell
                # this behavior could be solved with more advanced models, necessary inverse production functions.
                
            return (-hyd,p/1000) # return hydrogen absorbed [kg] and electricity required [kWh]
        
        elif self.model=='PEM General':
            
            'Finding the working point of the FuelCell by explicitly solving the system:'
            
            p_required = -e  # [kWh] --> equivalent to kW for the considered timestep
            
            p = min(p_required,self.FC_NominalPower)          # [kW] limiting maximum power that the Fuel Cell can provide 
            
            FC_CellCurrDensity =self.PI(p)/self.FC_CellArea   # [A/cm^2] current density value at which the fuel cell is working 

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
            
            self.EFF[h]       = etaFC                         # [-]      FC efficiency history
            self.VOLT[h]      = V_cell                        # [V]      Cell voltage history
            self.CURR_DENS[h] = FC_CellCurrDensity            # [A/cm2]  Cell current density history
            
            'Hydrogen demand'
         
            FC_HydroCons = Current*self.nc*3600/(95719.25*1000)/self.rhoStdh2*self.timestep   #[Sm3] (Chavan 2017)
            hyd=FC_HydroCons*self.rhoStdh2/self.timestep                                      #[kg/h]
            FC_deltaHydrogen = - FC_HydroCons*self.rhoStdh2*self.HHVh2*1000/3600              #[kWh]
         
         
            'Process heat, that can be recovered'
          
            FC_Heat = ((1.481*self.nc)/FC_Vstack-1)*p*self.timestep               # [kWh] --> equivalent to kW for the considered timestep
            
            
            if hyd > available_hyd: # if not enough hydrogen is available to meet demand (H tank is nearly empty)
                hyd = 0
                p = 0
                # turn off the fuel cell
                # this behavior could be solved with more advanced models, necessary inverse production functions.
                
            # return (-hyd,p,FC_Heat,Current,etaFC) # return hydrogen absorbed [kg] and electricity required [kWh] and heat as a co-product [kWh]
            return (-hyd,p,FC_Heat,Current) # return hydrogen absorbed [kg] and electricity required [kWh] and heat as a co-product [kWh]
            
            
    
    
##########################################################################################

if __name__ == "__main__":
    
    """
    Functional test
    """
    
    inp_test = {'Npower': False,
                'stack model':'PEM General',
                }
    
    sim_hours=96                               # [h] simulated period of time - usually it's 1 year minimum
    time=np.arange(sim_hours)
    
    fc = fuel_cell(inp_test,sim_hours)         # creating fuel cell object
    fc.plot_polarizationpts()                  # cell polarization curve
    fc.plot_stackperformance()                 # stack performance curve
    
    flow =-np.linspace(1,25,sim_hours)         # [kW] power required from the Fuel Cell 
    fd   =-np.random.uniform(0,15,sim_hours)   # [kW] power required from the Fuel Cell - random values

    P_el=[]
    P_th=[]
    I=[]
    for h in range(len(flow)):
        P_el.append(fc.use(h,flow[h],100)[1])
        P_th.append(fc.use(h,flow[h],100)[2])
        I.append(fc.use(h,flow[h],100)[3]) 
        

    fig=plt.figure(figsize=(8,8),dpi=1000)
    fig.suptitle("PEMFC ({} kW) performances".format(round(fc.FC_NominalPower,1)))
    
    PI=fig.add_subplot(211)
    PI.plot(I,P_el,label="P_el",color="blue")
    PI.plot(I,P_th,label="P_th",color="red")
    PI.set_title("P vs I")
    PI.grid()
    PI.set_xlabel("I [A]")
    PI.set_ylabel("P [kW]")
    PI.legend(fontsize=15)
    
    ETA=fig.add_subplot(212)
    ETA.plot(I,fc.EFF,label="Eta",color="green")
    ETA.set_title("ETA vs I")
    ETA.grid()
    ETA.set_xlabel("I [A]")
    ETA.set_ylabel("Eta")
    ETA.legend(fontsize=15)
    
    plt.tight_layout()
    plt.show()
    
    'Behaviour over time'
    
    P_el=[]
    P_th=[]
    I=[]
    for h in range(len(flow)):
        P_el.append(fc.use(h,fd[h],100)[1])
        P_th.append(fc.use(h,fd[h],100)[2])
        I.append(fc.use(h,fd[h],100)[3])
    
    fig=plt.figure(figsize=(8,8),dpi=1000)
    fig.suptitle("PEMFC ({} kW) performances ".format(round(fc.FC_NominalPower,1)))
    
    PI=fig.add_subplot(211)
    PI.plot(time,P_el,label="P_el",color="blue")
    PI.plot(time,P_th,label="P_th",color="red")
    PI.set_title("P over time")
    PI.grid()
    PI.set_xlabel("Time [h]")
    PI.set_ylabel("P [kW]")
    PI.legend(fontsize=15)
    
    ETA=fig.add_subplot(212)
    ETA.plot(time,fc.EFF,label="Eta",color="green")
    ETA.set_title("ETA over time")
    ETA.grid()
    ETA.set_xlabel("Time [h]")
    ETA.set_ylabel("Eta")
    ETA.legend(fontsize=15)
    
    plt.tight_layout()
    plt.show()
    