import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import os
import sys 
from scipy.interpolate import interp1d
from CoolProp.CoolProp import PropsSI
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),os.path.pardir)))   # temorarily adding constants module path 
import constants as c

class chp_gt:    
    
    def __init__(self,parameters,simulation_hours):
        
        """
        Create a GT-based CHP object
    
        parameters : dictionary
            # 'max capacity': float [kg]
            # 'pressure': float [bar]
        
                      
        output : CHP object able to:
            
            # supply or abrosrb hydrogen .use(h,hyd)
            # record the level of charge .LOC
            # calculate its own volume (pressure) .volume(pressure)
        """
        
        self.tech = parameters["Technology"]
        self.fuel = parameters["Fuel"]
        
        if self.tech != "Gas Turbine" or self.fuel != "Hydrogen":
            raise ValueError(""" 
                             Only H2 fuelled GT available as cogeneration system (5.6 MW rated power).
                             Either use "Gas Turbine" and "Hydrogen" keywords for initialization
                             or remove chp_gt technology from the case study                     
                             """)
        
        self.wel=np.zeros(simulation_hours)            # [W]    produced electricity
        self.mH2CHP=np.zeros(simulation_hours)         # [kg/s] hydrogen mass flow rate required by GT + HRSG
        # mH2SG = np.zeros(simulation_hours)        # [kg/s] hydrogen mass flow rate required by SG
        # mH2=np.zeros(simulation_hours)            # [kg/s] hydrogen mass flow rate required by the whole system GT + HRSG + SG
        self.minprod=np.zeros(simulation_hours)        # [kg/s] minimum steam amount producible given weather conditions (tamb)
        self.maxprod=np.zeros(simulation_hours)        # [kg/s] maximum steam amount producible given weather conditions (tamb)
        # steam_SG=np.zeros(simulation_hours)       # steam required from steam generator units
        self.steam_chp=np.zeros(simulation_hours)      # [kg/s] steam produced by CHP system, GT + HRSG components
        self.steam_miss=np.zeros(simulation_hours)     # [kg/s] amount of steam the CHP system has been unable to provide due to its operational limits
        # pump = np.zeros(simulation_hours)         # [kW] pump power consumption
        
        # Required process steam properties
        self.steam = {
                     'Pressure out'       : 20*1e5,
                     'Vapour fraction out': 1,
                     'Pressure in'        : 1*1e5,
                     'Temperature in'     : 303.15                         
                     }

        # Specific hydrogen gas turbine parameters 
        self.GT_param = {
                        'Power[MW]': 5.649,
                        'Gross Efficiency[-]': 0.31135,
                        'Pressure Ratio[-]': 15.62,
                        'Exhaust Mass Flow[kg/s]': 20.7645,
                        'Exhaust temperature[K]': 821.883,
                        'TIT[K]':1530.07 
                        }


        'Data Extraction - Operation Maps for Real Working Conditions CHP System'
        
        main = False                      # check
        if __name__ == "__main__":        # if code is being executed from chp_gt.py script
            os.chdir(r'./chp_maps')
        else: 
            os.chdir(r'./techs/chp_maps') # if code is being executed from main
            main = True

        with pd.ExcelFile('CHPmaps.xlsx') as xls:
        
            self.Wel_map = pd.read_excel(xls,sheet_name='W_el',header=2,nrows= 7,usecols='A:G',index_col='Tamb [°C]')     # [kW]    Net Electric Power Output of the GT MAP
            self.Eta_map = pd.read_excel(xls,sheet_name='Eta_tag',header=2,usecols='A:G',index_col='Tamb [°C]')           # [-]     GT Efficiency MAP
            self.Beta_map = pd.read_excel(xls,sheet_name='Beta',header=2,usecols='A:G',index_col='Tamb [°C]')             # [-]     Compression Ratio MAP
            self.mH2fuel_map = pd.read_excel(xls,sheet_name='m_fuel',header=2,usecols='A:G',index_col='Tamb [°C]')        # [kg/s]  Fuel mass flow rate MAP
            self.Tmax_map = pd.read_excel(xls,sheet_name='Tmax',header=2,usecols='A:G',index_col='Tamb [°C]')             # [K]     Max Cycle Temperature MAP
            self.mExhGas_map = pd.read_excel(xls,sheet_name='m_fumi',header=2,usecols='A:G',index_col='Tamb [°C]')        # [kg/s]  Exhaust Gases MAP
            self.mCoolant_map = pd.read_excel(xls,sheet_name='m_coolant',header=2,usecols='A:G',index_col='Tamb [°C]')    # [kg/s]  GT coolant mass flow rate MAP
            self.Texh_map = pd.read_excel(xls,sheet_name='T_exh',header=2,usecols='A:G',index_col='Tamb [°C]')            # [K]     Exh gases from GT Temperature MAP
            self.Wth_inGT_map = pd.read_excel(xls,sheet_name='Wth_tag_in',header=2,usecols='A:G',index_col='Tamb [°C]')   # [W]     Wth enetering the GT MAP
            self.Qeva_map = pd.read_excel(xls,sheet_name='Qeva',header=2,usecols='A:G',index_col='Tamb [°C]')             # [kW]    Exchanged Heat in Evaporator MAP
            self.Qeco_map = pd.read_excel(xls,sheet_name='Qeco',header=2,usecols='A:G',index_col='Tamb [°C]')             # [W]     Exchanged Heat in Economizer MAP
            self.Tstack_map = pd.read_excel(xls,sheet_name='Tstack',header=2,nrows=7,usecols='A:G',index_col='Tamb [°C]') # [K]     Exhaust Gases Temperature at Stack MAP
            self.DTpp_map = pd.read_excel(xls,sheet_name='∆Tpp',header=2,usecols='A:G',index_col='Tamb [°C]')             # [K]     Delta T Pinch Point MAP
            self.DTsub_map = pd.read_excel(xls,sheet_name='∆Tsub',header=3,usecols='A:G',index_col='Tamb [°C]')           # [K]     Delta T Pinch Point MAP
           
            self.constr = pd.read_excel(xls,sheet_name='System Boundaries',header=15,usecols='C:F')                       # System constraints imported as a DataFrame
            self.constr = self.constr.drop(columns='Tstack')
        
        if main ==  True:                 # if code is being executed from main, change directory back to main
            os.chdir(r'../..')

        'Operation Constraints'
    
        self.Wel_max= self.GT_param['Power[MW]']*1.2    # [kW] Max Electric Power Output - 120% of Nominal Power
        self.Wel_min= self.GT_param['Power[MW]']*0.3    # [kW] Min Electric Power Output - 30% of Nominal Power
        self.Tmax= self.GT_param['TIT[K]']              # [K]  Turbine Inlet Temperature/ Max cycle T - 1530.07 K
        self.Ts_min= 90 + 273.15                   # [K]  Min Temperature at Stack
    
        'System Boundaries - Lines'
    
        self.t_amb = np.array([-10,0,10,15,20,30,40])   # [°C] Temperature range used for simulations
        self.m_vap = np.array([1.5,2,3,4,5,6])          # [kg/s] Steam Mass Flow Rate range used employed in the simulations 
       
        self.limits={}
        self.constr_func=[]
        self.mylines=[]
        self.a1=[]
        self.indexes=[]
        lines=[]
        
        for label in self.constr.columns:
            a = np.array(self.constr[label])
            idx = np.isfinite(a)  & np.isfinite(self.t_amb)            # cleaning up NaN values given that polyfit function cannot handle NaN values
            self.indexes.append(idx)
            # polreg = np.poly1d(np.polyfit(a[idx],self.t_amb[idx],5))   # polynomial regression
            polreg = interp1d(a[idx],self.t_amb[idx],bounds_error=None,fill_value='extrapolate')    # polynomial interpolation
                                                                       # np.polyfit -> It is a fit polynomial p(x) = p[0] * x**deg + … + p[deg] 
                                                                       #               of degree deg to points (x, y). It fits data to a polynomial function
                                                                       #               https://appdividend.com/2022/01/28/numpy-polyfit-method-in-python/
            # polreg2 = np.poly1d(np.polyfit(self.t_amb[idx],a[idx],5))  # polynomial regression 2 - reversed relation/inverse function 
            polreg2 = interp1d(self.t_amb[idx],a[idx],bounds_error=None,fill_value='extrapolate')   # polynomial interpolation 2 - reversed relation/inverse function 
            line = np.linspace(self.t_amb[idx][0],self.t_amb[idx][-1],100)  # regression interval
            a = a[idx]                                                 # Saving != NaN values in a new array - prolly redundant
            myline=np.linspace(a[0],a[-1],100)                         # creating the interval for the polynomial regression
            self.constr_func.append(polreg(myline))                    # saving polynomial regression values - plot sake
            self.mylines.append(myline)                                # saving polynomial regression intervals
            self.a1.append(a)
            self.limits[label]=polreg,polreg2
       
    def map_plot(self):
        
        markers = list(Line2D.markers.keys())
        colors =['tab:blue','tab:orange','tab:blue']
        plt.figure(dpi=500)
        i=0
        for label in self.constr.columns:
            
            plt.scatter(self.a1[i],self.t_amb[self.indexes[i]],label = label[0]+'$_\mathregular{'+label[1:]+'}$', marker=markers[i], edgecolor='k', c=colors[i],zorder=3)
            plt.plot(self.mylines[i],self.limits[label][0](self.mylines[i]), linestyle='-.',c='k', zorder =2)
            i+=1
            
        plt.ylabel(r'T$_\mathregular{amb}$ [°C]')
        plt.xlabel(r'm$_\mathregular{st}$ [kg/s]')    
        plt.axvline(x=1.5,linestyle='--', c='tab:blue') 
        plt.axhline(y=max(self.t_amb),linestyle='--', c='tab:green') 
        plt.axhline(y=min(self.t_amb),linestyle='--', c='tab:green') 
        plt.fill_between(self.mylines[-1],40,self.constr_func[-1],facecolor='tab:blue')             # Wmin
        plt.fill_between(self.mylines[0],-10,self.constr_func[0],facecolor='tab:blue', zorder=1)    # Wmax
        plt.fill_between(self.mylines[0],40,self.constr_func[0], 
                          where =(self.mylines[0] > 4.4),
                          facecolor='tab:orange',zorder=0)
        plt.fill_between(self.mylines[0],40,self.constr_func[0], 
                          where =(self.mylines[0] >= self.mylines[1][0]-0.01) & (self.mylines[0] < 4.5),
                          facecolor='tab:orange',zorder=0)
        plt.fill_between(self.mylines[1],40,self.constr_func[1], facecolor= 'tab:orange',zorder=-1)  # TIT
        plt.legend(ncol=4, bbox_to_anchor = (0.88,1.14), fontsize = 'large')
        plt.grid(alpha=0.3, zorder=-100)
        plt.show()
                
                
    def bilinear_interpolation(self,mappa,v1,v2):
            """
            bilinear interpolation function. It queries performance maps with required load and Tamb and returns system performance
            
            inputs
                mappa:    performance map
                v1 :      float energy carrier request driving the demand [kWh] (electricity or heat)
                v2 :      float air temperature for the considered timestep [°C]
          
            output 
                y :  exact functioning point for the desired quantity
                    
            """
            x1 = np.array(mappa.columns).astype(float)  # dataset x1 (load)
            x2 = np.array(mappa.index).astype(float)    # dataset x2 (t_amb)
            y_ds = mappa.values.astype(float)           # dataset y (map values) 
            y2 = [] 
            for i in range(len(x2)):
                x_tmp = x1
                y_tmp = y_ds[i,:]
                y2.append(np.interp(v1, x_tmp, y_tmp))
            y2 = np.array(y2)    
            y = np.interp(v2, x2, y2)
            return y
    
    
    def use(self,h,t_air,steamdemand):
        """
        The chp system consumes fuel and produces multiple output energy streams as defined by the specific technology
        
        inputs
            h:      int hour to be simulated
            t_air:  float air temperature for the considered timestep [°C]
            demand: float energy carrier request driving the demand [kWh] (electricity,heat or steam [kg/h])
      
        output 
            carrier1:    produced energy stream driving the system functioning [kWh] 
            carrier2:    2nd energy stream produced as co-product [kWh]
            ...
            nth-carrier: nth energy stream produced as co-product [kWh]
                
        """
        demand = abs(steamdemand)/3600                    # converting steam demand from kg/h to kg/s 
        
        bound = []                                        # operational boundaries of the CHP system for the given air temperature
        for label in self.limits:
            bound.append(self.limits[label][1](t_air))    # evaluating real limits of the system depending on ambient temerature
        bound = np.sort(bound)
        # bound[bound<1.5] = 1.5                       # lower operational boundary set in the operative constraints (PER ORA NON è UTILIZZATO)
        bound = bound[0:2]                                # lower and upper operational boundaries values for the required carrier
        self.minprod[h] = bound[0]                        # min producibility for given conditions
        self.maxprod[h] = bound[1]                        # max producibility for given conditions
       
        if demand in pd.Interval(bound[0],bound[1],closed='both'):  # steam demand within the GT + HRSG range for given temperature
            self.wel[h]= self.bilinear_interpolation(self.Wel_map,demand,t_air)  
            self.steam_chp[h]= demand
            # steam_SG[h]= 0   
            self.steam_miss[h]= 0
            self.mH2CHP[h]= self.bilinear_interpolation(self.mH2fuel_map,demand,t_air)
            # mH2SG[h] = 0
            # mH2[h] = mH2CHP[h]
        elif demand > self.maxprod[h]:
            self.wel[h]= self.bilinear_interpolation(self.Wel_map,self.maxprod[h],t_air)
            self.steam_chp[h]=  self.maxprod[h]
            # steam_SG[h]= demand -  maxprod[h]
            self.steam_miss[h] = (demand-self.maxprod[h]) # - SG.use(steam_SG[h],deltaSG)[0]
            self.mH2CHP[h]= self.bilinear_interpolation(self.mH2fuel_map,self.maxprod[h],t_air)
            # mH2SG[h] = SG.use(steam_SG[h],delta)[2]
            # self.mH2[h] = mH2CHP[h] #+ mH2SG[h]
            # pump[h] = SG.pcons(steam_SG[h], 1.01325, 20)
        elif demand < self.minprod[h]:            # GT running to avoid shutdowns (given it would be turned off only for 62/8760 h/y  - 0.06 % of the time)
            self.wel[h]= self.bilinear_interpolation(self.Wel_map,self.minprod[h],t_air)
            self.steam_chp[h]= demand
            # steam_SG[h]= 0  
            self.steam_miss[h] = 0
            # steam_miss[h] = steam_SG[h] - SG.use(steam_SG[h],deltaSG)[0]
            self.mH2CHP[h]= self.bilinear_interpolation(self.mH2fuel_map,self.minprod[h],t_air) # +SG.use(steam_SG[h],delta)[2] 
            # mH2SG[h] = 0
            # mH2[h] = mH2CHP[h]
        
        return(self.steam_chp[h]*3600,self.wel[h],-self.mH2CHP[h]*3600) # return hydrogen supplied
                    


    def delta_h(self,s_in_P,s_in_Q,s_out_P,s_out_T):
        '''
        Parameters
        ----------
        s_in_P :  [Pa] required pressure of steam
        s_in_Q :  [-]  vapor fraction of the inlet steam
        s_out_P : [Pa] water pressure at the end of condesation process - input conditions of Steam Generator TYPE
        s_out_T : [K]  water temperature after condensation - input conditions of Steam Generator

        Returns
        -------
        h2 :      [J/kg] output stream enthalpy
        h1 :      [J/kg] input stream enthalpy
        delta :   [kJ/kg] Delta h of transformation

        '''
        
        h2= PropsSI('H', 'P', s_in_P, 'Q',s_in_Q, 'Water')       # [J/kg]
        h1= PropsSI('H', 'P', s_out_P, 'T', s_out_T , 'Water')   # [J/kg]
        t1= PropsSI('T','P', s_in_P, 'Q',s_in_Q, 'Water')        # [J/kg]
        delta= (h2-h1)/1000
        
        return (h2,h1,delta)        

    def pes(self):
        
        delta = self.delta_h(self.steam['Pressure out'],
                             self.steam['Vapour fraction out'],
                             self.steam['Pressure in'],
                             self.steam['Temperature in'])[2]
        
        Ele_gen     = sum(self.wel/1000)                # [MWh] Electric Energy produced TOT
        Th_gen      = sum((self.steam_chp*delta)/1000)  # [MWh] Thermal Energy produced TOT
        Fuel_energy = sum(self.mH2CHP)*c.LHVH2          # [MWh] Thermal Energy input as fuel
        
        # Raginonamenti validi se si considera Natural Gas, con H2 non hanno senso di essere calcolati tali valori

        Eta_global = (Ele_gen+Th_gen)/Fuel_energy
        print('The Global CHP plant efficiency is:',round(Eta_global,3))
        CHP_etael = Ele_gen/Fuel_energy
        CHP_etath = Th_gen/Fuel_energy
        
        # Italian scenario 
        RefH_eta = 0.87
        RefE_eta = 0.5033
        
        PES = (1-1/((CHP_etath/RefH_eta)+(CHP_etael/RefE_eta)))*100
        print('PES indicator value:' +str(round(PES,3))+ ' %')
        
        # Certificati Bianchi
        
        Etarif_el = 0.46   # rendimiento medio convenzionale parco di produzione elettrica italiano
        Etarif_th = 0.90   # rendimento termico di riferimento (produzione di vapore)
        
        Risp = Ele_gen/Etarif_el + Th_gen/Etarif_th - Fuel_energy    # [MWh]/anno - Risparmio? 
        k= 1.3
        
        cb= Risp *0.086 * k   # tep risparmiate 
        price = 256.84        # [€/tep] prezzo medio ponderato as of 11/01/2022
        revenews = cb*price   # [€/y]
        
        return(Eta_global,PES)
        
    def tech_cost(self,tech_cost):
        """
        Parameters
        ----------
        tech_cost : dict
            'cost per unit': float [€/kWh]
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

        if tech_cost['cost per unit'] == 'default price correlation':
            
            #Cost references: https://understandingchp.com/chp-applications-guide/6-8-rules-of-thumb-for-chp-engineering-and-installation-costs/
            change = 1.183 # [$/€] average 2021
            
            # GT - Impianti di potenza dispense 
            C = 6380*(self.GT_param['Power[MW]']*1000)**0.73               
            OeM = 0.06*(6380*((self.GT_param['Power[MW]']*1000)**0.73)) 
            
            # HRSG - Impianti di potenza dispense 
            C += (68679.85 + 182811.26 + 211764.99+ 6380*((self.GT_param['Power[MW]']*1000)**0.73))*0.04
            OeM += 0.06*(68679.85 + 182811.26 + 211764.99 + 6380*((self.GT_param['Power[MW]']*1000)**0.73))*0.04
            
            # Alternator
            C += 4028*((180/(1.38))**0.58)              
            OeM += 0.06*4028*((180/(1.38))**0.58)
            
# =============================================================================
#             # Pump
#             C += 68679.85/change
#             OeM += 0.06*(68679.85/change)
#             
#             # Eco
#             C += 182811.26/change
#             OeM += 0.06*(182811.26/change)
#             
#             # Eva
#             C += 211764.99/change
#             OeM += 0.06*(211764.99/change)
#             
#             # Condenser
#             C += (6380*(5649**0.73))*0.04
#             OeM += 0.06*((6380*(5649**0.73))*0.04)
# =============================================================================
                   
        else:
            print( "essendo presente solamente un modello di CHP e di una taglia fissa il costo può essere fatto solo con la default price correlation")

            tech_cost['total cost'] = tech_cost.pop('cost per unit')
            tech_cost['total cost'] = C
            tech_cost['OeM'] = OeM
            tech_cost['refund'] = { "rate": 0, "years": 0}
            tech_cost['replacement'] = {"rate": 80, "years": 30}

            self.cost = tech_cost
        
#%%##########################################################################################

if __name__ == "__main__":
    
    """
    Functional test
    """
    
    inp_test = {"Technology": "Gas Turbine",
                "Fuel"      : "Hydrogen"    }       
    
    simulation_hours = 8760                   # 1 year-long simulation
    chp_gt = chp_gt(inp_test,simulation_hours)   # creating chp object
    
    chp_gt.map_plot()                         # displaying CHP operational constraints

    x = np.arange(0,24)       # hours in a day
    simulation_hours = 24     
    days = ['Winter day', 'Spring day','Summer day','Autumn day']
    daily_temperature = {
                         days[0] : [-2.2,-2.7,-3.1,-3.2,-3.6,-3.9,-4.2,-4.3,-4.5,-1.6,1.9,6.3,7.3,7.8,7.5,7.6,6.9,4,0.7,-0.2,-1.1,-1.9,-1.9,-2.6],        # [°C] typical winter day temperatures 
                         days[1] : [4.5,3.9,3.1,2.7,2.6,2.6,2.6,2.6,8,12,12.8,13.7,13.8,14,14.2,14.3,14.3,13.8,13.8,12.2,9.9,7.2,6.3,5.3],                # [°C] typical spring day temperatures
                         days[2] : [18.4,18.1,17.9,17.4,17.4,17.4,17.4,17.4,24,26.2,26.9,27.6,28,28.4,28.6,27.4,29.2,28.7,29.8,27.6,24.6,23,21.9,20.5],   # [°C] typical summer day temperatures
                         days[3] : [20.6,20.6,20.6,20.3,19.5,18.2,16.7,15.7,19,21.7,23.1,23.9,24,24,24.5,24.7,24.3,22.8,21.1,19.6,18.1,16.9,16.1,15.6]    # [°C] typical autumn day temperatures
                         }
    
    for k in daily_temperature:
        steam_demand = np.random.uniform(0.5,5.2,simulation_hours)*3600    # [kg/h] creating a random steam demand array 
       
        for h in range(24):
            chp_gt.use(h,daily_temperature[k][h],steam_demand[h])
        
        fig, ax = plt.subplots(dpi=600)
        ax.plot(x,chp_gt.minprod[:simulation_hours],label ='CHP$_\mathregular{min}$', alpha=0.9)
        ax.plot(x,chp_gt.maxprod[:simulation_hours],label = 'CHP$_\mathregular{max}$',alpha=0.9)
        ax.scatter(x,steam_demand/3600, label = 'm$_\mathregular{stdem}$', alpha =0.8, c='g',s=70,zorder=10,edgecolors='k', linewidths=0.6)
        ax.scatter(x,chp_gt.steam_chp[:simulation_hours], label = 'm$_\mathregular{stprod}$', alpha =0.8, c='b',s=22,zorder=10,edgecolors='k', linewidths=0.6)
        # plt.plot(x,steam_SG, label = 'm$_{st}$ SG', alpha =0.5, c='k')
        ax.grid(zorder=0, alpha= 0.4)
        ax.set_ylabel('m$_{st}$ [kg/s]')
        ax.set_xlabel('Time [h]')
        ax.set_ylim(0,5.2)
        ax.fill_between(x,max(steam_demand)/3600,chp_gt.maxprod[:simulation_hours], facecolor = 'orangered',zorder=0, alpha=0.3,label = 'SG') 
        ax.fill_between(x,min(steam_demand)/3600,chp_gt.minprod[:simulation_hours], facecolor = 'lightblue',zorder=0, alpha=0.3)
        ax.fill_between(x,chp_gt.minprod[:simulation_hours],chp_gt.maxprod[:simulation_hours], facecolor = 'lightblue',zorder=0, alpha=0.8,label = 'CHP')
        ax.legend(ncol=3,bbox_to_anchor = (0.85,-0.15))
        ax.set_title('CHP system behaviour - ' + k)
        plt.show()
     
    chp_gt.pes()

        
        # 'System Boundaries - Year-long simulation'
            
        # x=np.arange(0,simulation_hours)     # hours in one year
        # steam_demand = np.random.uniform(0.5,5.2,simulation_hours)*3600    # [kg/h] creating a random steam demand array 
        # ambient_t    = np.random.uniform(-5,40,simulation_hours)           # [°C]   creating a random ambient temperature array 
        
        # for h in range(simulation_hours):
            
        #     chp_gt.use(h,ambient_t[h],steam_demand[h])
            
            
        # 'System Boundaries - SCATTER PLOT'
        
        # plt.figure(dpi=500)
        # n=24
        # plt.xticks([n*m for m in np.arange(0,361,30)],['        Jan','        Feb','         Mar','        Apr','        May','         Jun','        Jul','        Aug','        Sep','        Oct','        Nov','        Dec',''])
        # plt.ylabel('m$_{st}$ [kg/s]')
        # plt.plot(x,chp_gt.minprod[:simulation_hours],label ='CHP$_\mathregular{min}$', alpha=0.9)
        # plt.plot(x,chp_gt.maxprod[:simulation_hours],label = 'CHP$_\mathregular{max}$',alpha=0.9)
        # # plt.axhline(y=1.5)      # theoretical minimum steam production of the system 
        # plt.scatter(x,steam_demand/3600, label = 'm$_\mathregular{st}$', alpha =0.8, c='g',s=3,zorder=10,edgecolors='k', linewidths=0.6)
        # # plt.plot(x,steam_SG, label = 'm$_{st}$ SG', alpha =0.5, c='k')
        # plt.grid(zorder=0, alpha= 0.4)
        # plt.fill_between(x,max(steam_demand)/3600,chp_gt.maxprod[:simulation_hours], facecolor = 'orangered',zorder=0, alpha=0.3,label = 'SG') 
        # plt.fill_between(x,min(steam_demand)/3600,chp_gt.minprod[:simulation_hours], facecolor = 'lightblue',zorder=0, alpha=0.3)
        # plt.fill_between(x,chp_gt.minprod[:simulation_hours],chp_gt.maxprod[:simulation_hours], facecolor = 'lightblue',zorder=0, alpha=0.8,label = 'CHP')
        # plt.legend(ncol=5,bbox_to_anchor = (0.94,1.12), fontsize='small')
        # plt.show()