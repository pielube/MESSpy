# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 13:58:34 2023

@author: mati/zini
"""
# Packages
## System and folders
import os
import sys 
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),os.path.pardir)))   # temorarily adding constants module path 

## Data handling
import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator

## Data visualization
import matplotlib.pyplot as plt

## Custom
import constants as c


def bilinear_interp(_map,v1,v2):
        """
        bilinear interpolation function. It queries performance maps with required load and Tamb and returns system performance
        
        inputs
            _map :    performance map
            v1 :      float value representing chp load at the given timestep [-]
            v2 :      float air temperature for the considered timestep [°C]
      
        output 
            y :  exact functioning point for the desired quantity
                
        """
        # x1 = np.array(_map.columns).astype(float)  # dataset x1 (load)
        # x2 = np.array(_map.index).astype(float)    # dataset x2 (t_amb)
        # y_ds = _map.values.astype(float)           # dataset y (map values) 
        # y2 = [] 
        # for i in range(len(x2)):
        #     x_tmp = x1
        #     y_tmp = y_ds[i,:]
        #     y2.append(np.interp(v1, x_tmp, y_tmp))
        # y2 = np.array(y2)    
        # y = np.interp(v2, x2, y2)
        
        interp = RegularGridInterpolator((np.array(_map.index),
                                          np.array(_map.columns)),
                                          np.array(_map.values), 
                                          method='linear')
        
        y = interp((v2,v1))
        return y
    
    
def inverse_bilinear_interp(_map, y, t_amb):
    '''
    Parameters
    ----------
    _map  : pd dataframe of float values - performance map
    y     : float value representing the system "Limit" value - according to the selected "Method"
    t_amb : float value of ambient temperature [°C]

    Returns
    -------
    l2    : float value representing working "Load" of the system correspondent to the defined "Limit"\
            for the given tamb

    '''
    i           = 0         # initializing iteration count
    eps         = 100000000 # absolute error
    epsilon     = 0.001     # imposed error treshold
    cond_while  = True
    
    l1 = _map.columns[-1]   # l1 perturbation
    l2 = _map.columns[0]    # l2 perturbation
    
    y1 = bilinear_interp(_map, l1, t_amb)       # right endpoint for the considered method
    
    while cond_while:
        y2 = bilinear_interp(_map, l2, t_amb)   # first iteration: left endpoint for the considered method
        eps = y - y2                            # absolute error calculation between y (limit) and computed left endpoint value
        l = l2 + ((l1-l2)/(y1-y2))*(y-y2)       # 'inverse' calculation of load - linear approach
        if l < _map.columns[0]:                 # avoid "out of range" interpolation
            l = _map.columns[0]                 # l assigned the minimum load value of the given map           
        elif l > _map.columns[-1]:              # avoid "out of range" interpolation
            l = _map.columns[-1]                # l assigned the maximum load value of the given map
        l1 = l2
        y1 = y2
        l2 = l    
        if  i > 1000 or abs(eps) < (epsilon*y):   # 1000 maxnumber of iterations allowed
            cond_while = False
        i += 1
    return l2    


class Chp:
    
    def __init__(self, parameters, simulation_hours):
        self.fuel           = parameters["Fuel"]            # type of fuel fed to the chp
        self.strategy       = parameters["Strategy"]        # parameter on which the operation of the system is based
        self.coproduct      = parameters["Co-product"]      # co-product energy stream
        self.th_out         = parameters["Thermal Output"]  # type of stream into which heat from cobustion is converted/transferred. Steam or hot water
        self.control_param  = parameters["Control Param"]   # control parameters to define operational boundaries of the system
        self.load           = np.zeros(simulation_hours)    # [-] working load of the system 
        self.q_th           = np.zeros(simulation_hours)    # [kWh] thermal output of the chp
        self.w_el           = np.zeros(simulation_hours)    # [kWh] electricity output of the system
        self.m_fuel         = np.zeros(simulation_hours)    # [kg/h] fuel consumption
        self.l_bound        = np.zeros(simulation_hours)    # [-] minimum load producible given working conditions
        self.u_bound        = np.zeros(simulation_hours)    # [-] maximum load producible given working conditions
        self.steam          = np.zeros(simulation_hours)    # [kg/s] steam produced by CHP system
        self.hot_w          = np.zeros(simulation_hours)    # [kWh] hot water produced by CHP system
        self.shutdown       = np.zeros(simulation_hours)    # [0/1] array to keep track of numbers of system shutdowns
        self.performances   = {                                                 # performance parameters dictionary. Self-consumption % of the considered energy streams
                                self.strategy   : np.zeros(simulation_hours),
                                self.coproduct  : np.zeros(simulation_hours),
                              }    

        if self.fuel == 'gas':
            self.LHVfuel = c.LHVNG
        else:
            self.LHVfuel = c.LHVH2
            
        'Data Extraction - Provided Operation Maps for the considered CHP system'
        main = False                      # check
        if __name__ == "__main__":        # if code is being executed from chp_gt.py script
             os.chdir(r'./chp_maps')
        else: 
            os.chdir(r'./techs/chp_maps') # if code is being executed from main
            main = True
        
        with pd.ExcelFile('CHPmaps.xlsx') as xls:
            self.maps =  {
                          "electricity"         : pd.read_excel(xls,sheet_name='W_el',header=2,nrows= 7,usecols='A:H',index_col='Tamb [°C]'),  # [kW] Net Electric Power Output 
                          "process heat"        : pd.read_excel(xls,sheet_name='Q_th',header=2,nrows= 7,usecols='A:H',index_col='Tamb [°C]'),  # [kW] Thermal Power Output
                          "fuel"                : pd.read_excel(xls,sheet_name='m_fuel',header=2,usecols='A:H',index_col='Tamb [°C]'),         # [kg/s] Fuel Mass Flow Rate Consumption 
                          "process steam"       : pd.read_excel(xls,sheet_name='m_steam',header=2,usecols='A:H',index_col='Tamb [°C]'),        # [kg/s] Steam Mass Flow Rate Production
                          "TIT"                 : pd.read_excel(xls,sheet_name='TIT',header=2,usecols='A:H',index_col='Tamb [°C]'),            # [K] Turbine Inlet Temperature 
                          "Tstack"              : pd.read_excel(xls,sheet_name='Tstack',header=2,usecols='A:H',index_col='Tamb [°C]')          # [K] Exhaust Gases Temperatures at Stack 
                          }
        
        if main ==  True:                 # if code is being executed from main, change directory back to main
            os.chdir(r'../..')
            

    def bound(self, method, lim, t_amb):
        """

        Parameters
        ----------
        method : String, method chosen to control the system
        lim    : float value representing the actual chosen limit, as a consequence of Method
        t_amb  : float value of ambient temperature [°C]

        Raises
        ------
        ValueError
            It is not considered to work at a load lower or higher than the given boundaries (0-1)

        Returns
        -------
        lim    : working load limit of the system for the given conditions. Calculated \
                 via inverse bilinear interpolation if not already given as load value

        """
        if method.__eq__("Load"):
            if lim < 0 or lim > 1:
                raise ValueError("Invalid Load bound values: out of range")                                 
            return lim
        else:            
            return inverse_bilinear_interp(self.maps[method], lim, t_amb)
        
    
    def use(self, h, t_amb, demand, demand2, available_fuel = None):    # None
        """
        The chp system consumes fuel and produces multiple output energy streams as defined by the specific technology
        
        inputs
            h:          int hour to be simulated
            t_air:      float air temperature for the considered timestep [°C]
            demand:     float energy carrier request driving the demand [kWh] (electricity,heat or steam [kg/h])
            demand2:    float energy carrier request as a main CHP co-product  [kWh] (electricity,heat or steam [kg/h])
      
        output 
            carrier1:    produced energy stream driving the system functioning [kWh] 
            carrier2:    2nd energy stream produced as co-product [kWh]
            ...
            nth-carrier: nth energy stream produced as co-product [kWh]
                
        """
        demand = abs(demand)     # demand has to be a positive value     

        # Lower functioning bound
        lb_list =  []
        for lb in self.control_param['Lower'].keys():
            lb_list.append(self.bound(self.control_param['Lower'][lb]['Method'], 
                                      self.control_param['Lower'][lb]['Limit'],
                                      t_amb))
        # Upper functioning bound    
        ub_list =  []
        for ub in self.control_param['Upper'].keys():
            ub_list.append(self.bound(self.control_param['Upper'][ub]['Method'], 
                                      self.control_param['Upper'][ub]['Limit'],
                                      t_amb))
        
        self.l_bound[h]  = max(lb_list)   # saving minimum functioning boundary - expressed as load %
        self.u_bound[h]  = min(ub_list)   # saving maximum functioning boundary - expressed as load %
        
        if available_fuel is not None:  # if available_fuel is a parameter to be considered in the analysis (hydrogen case)
        
            if available_fuel <= 0:       # if there is no fuel available, chp system is turned off
                self.load[h]        = 0
                self.m_fuel[h]      = 0
                self.q_th[h]        = 0
                self.w_el[h]        = 0
                self.steam[h]       = 0
                self.hot_w[h]       = 0
                self.shutdown[h]    = 1     # saving system working history,    1 = 'System has been switched off at this timestep' \
                                            #                                   0 = 'System has been running at this timestep'
                for energy_stream in self.performances:     # no beneficial effect from CHP
                    self.performances[energy_stream][h] = 0
                                                
                return (self.steam[h], self.w_el[h], -self.m_fuel[h], self.q_th[h], self.hot_w[h])  # stop function execution and return values
            
            else:       # computing system performances accounting for fuel availability
        
                minfuel = bilinear_interp(self.maps['fuel'],self.l_bound[h],t_amb)    # control needed when working with hydrogen \
                                                                                      # and no external source of fuel available (i.e. no grid connection)  
                if available_fuel < minfuel:   # if available_fuel is lower than minimum fuel required at the minimum load, system is turned off
                    self.load[h]        = 0
                    self.m_fuel[h]      = 0
                    self.q_th[h]        = 0
                    self.w_el[h]        = 0
                    self.steam[h]       = 0
                    self.hot_w[h]       = 0
                    self.shutdown[h]    = 1
                    for energy_stream in self.performances:     # no beneficial effect from CHP
                        self.performances[energy_stream][h] = 0
                    return (self.steam[h], self.w_el[h], -self.m_fuel[h], self.q_th[h], self.hot_w[h])  # stop function execution and return values
                
                else:
                    load  = inverse_bilinear_interp(self.maps[self.strategy], demand, t_amb)    # operating load corresponding to demand at considered timestep    
                    mfuel = bilinear_interp(self.maps['fuel'], load, t_amb)                     # [kg/h] fuel consumption correspondig to defined operating load
                    
                    if mfuel > available_fuel:   # if too much fuel is required compared to what is available
                        load = inverse_bilinear_interp(self.maps['fuel'], available_fuel, t_amb)   # maximum load based on available fuel is computed and system is operated accordingly
                    else:
                        pass         
                
        else:  
            load = inverse_bilinear_interp(self.maps[self.strategy], demand, t_amb) 
            
            
        if load in pd.Interval(self.l_bound[h],self.u_bound[h],closed='both'):   # if load is within the producible range of chp technology
            pass
        
        elif load > self.u_bound[h]:    # if load value is higher of the producible range of chp technology
            load = self.u_bound[h]
        
        elif load < self.l_bound[h]:    # if load value is lower of the producible range of chp technology
            load = self.l_bound[h]
    
        self.load[h]                    = load 
        self.m_fuel[h]                  = bilinear_interp(self.maps['fuel'],load,t_amb)    
        self.q_th[h]                    = bilinear_interp(self.maps['process heat'],load,t_amb)
        self.w_el[h]                    = bilinear_interp(self.maps['electricity'],load,t_amb)
        self.steam[h]                   = bilinear_interp(self.maps['process steam'],load,t_amb)      # [kg/h] steam produced by CHP system
        self.hot_w[h]                   = 0
        # self.parameters[self.strategy]  =                                                   # [kWh] hot water produced by CHP system - active for specific application
        
        return (self.steam[h], self.w_el[h], -self.m_fuel[h], self.q_th[h], self.hot_w[h])
    
    
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
        # Specific hydrogen gas turbine parameters 
        self.GT_param = {
                        'Power[MW]': 5.649,
                        'Gross Efficiency[-]': 0.31135,
                        'Pressure Ratio[-]': 15.62,
                        'Exhaust Mass Flow[kg/s]': 20.7645,
                        'Exhaust temperature[K]': 821.883,
                        'TIT[K]':1530.07 
                        } 
        
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
            
            tech_cost['total cost'] = C
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
        tech_cost['OeM'] = OeM
        tech_cost['refund'] = { "rate": 0, "years": 0}
        tech_cost['replacement'] = {"rate": 80, "years": 30}

        self.cost = tech_cost
            
            
class Absorber:    
    def __init__(self, parameters, simulation_hours):
        self.COP    = parameters["COP"]             # [-] Coefficient of Performance
        self.Npower = parameters["Npower"]          # [kW] Rated power of the absorber
        self.q_cool = np.zeros(simulation_hours)    # [kWh] initializing cold energy array
        self.q_used = np.zeros(simulation_hours)    # [kWh] initializing used energy array
        
   
    def use(self, h, q_in):
        e_absorbed = min(self.Npower,q_in)                          
        self.q_used[h] = e_absorbed
        self.q_cool[h] = e_absorbed*self.COP        # [kWh] when timestep is kept at 1 h kWh = kW   
        return (self.q_cool[h])                     # [kWh] cold energy produced by the absorber
    
    
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

        size = self.Npower 
        
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



# def performance_idx(self):
    
#     Ele_gen     = sum(self.wel/1000)             # [MWh] Electric Energy produced TOT
#     Th_gen      = sum(self.q_th + self.q_cool)   # [MWh] Thermal Energy produced including absorber 
#     Fuel_energy = sum(self.m_fuel)*self.LHVfuel  # [MWh] Thermal Energy input as fuelof CHP 

#     Eta_global = (Ele_gen+Th_gen)/Fuel_energy
    
    
#     CHP_etael = Ele_gen/Fuel_energy
#     CHP_etath = Th_gen/Fuel_energy
    
#     # Italian scenario 
#     RefH_eta = 0.87
#     RefE_eta = 0.5033
#     PES = (1-1/((CHP_etath/RefH_eta)+(CHP_etael/RefE_eta)))*100
    
#     # Certificati Bianchi
    
#     Etaref_el = 0.46   # the average conventional efficiency of the Italian electricity generation system
#     Etaref_th = 0.90   # average thermal conversion efficiency
    
#     Risp = Ele_gen/Etaref_el + Th_gen/Etaref_th - Fuel_energy    # [MWh]/anno - Risparmio? 
#     k= 1.3
    
#     cb= Risp *0.086 * k   # [toe] saved 
#     price = 256.84        # [€/tep] prezzo medio ponderato as of 11/01/2022
#     revs = cb*price       # [€/y] revenues
    
#     performance_idx =   {
#                         "Global efficiency"         : Eta_global,
#                         "PES"                       : PES,
#                         "Certificati bianchi revs"  : revs                
#                         }
    
#     return(performance_idx)
    
#%%##########################################################################################

if __name__ == "__main__":
    
    """ 
    Functional test
    """
    
    inp_test_chp = {"Fuel"          : "hydrogen",           # "hydrogen", "gas",  hydrogen fuel or natural gas
                    "Strategy"      : "process steam",      # "electricity", "process heat", "process steam", "process hot water"...electric/thermal load follow  
                    "Thermal Output": "process steam",      # "process steam" or "process hot water", type of stream into which heat from cobustion is converted/transferred
                    "Co-product"    : "electricity",        # "electricity", "process heat", "process steam", "process hot water"...co-product of system functioning
                    "Control Param" : {"Lower"  : {"1": {"Method" : "electricity",    # "electricity", "process heat", "Load", "TIT", "Tstack"
                                                        "Limit"  : 3000},              # [kW], [kW], [-], [K], [K]
                                                    },
                                        "Upper"  : {"1": {"Method": "Load",    # "electricity", "process heat", "Load", "TIT", "Tstack"
                                                        "Limit" : 1.},         # [kW], [kW], [-], [K], [K]
                                                    "2": {"Method": "TIT",     
                                                        "Limit" : 1530},                                                 
                                                    "3": {"Method": "Tstack",        
                                                        "Limit" : 363}
                                                    }
                                        }
                    }
    
    
    # inp_test_chp = {"Fuel"          : "hydrogen",           # "hydrogen", "gas",  hydrogen fuel or natural gas
    #                 "Strategy"      : "process steam",      # "electricity", "process heat", "process steam", "process hot water"...electric/thermal load follow  
    #                 "Thermal Output": "process steam",      # "process steam" or "process hot water", type of stream into which heat from cobustion is converted/transferred
    #                 "Co-product"    : "electricity",        # "electricity", "process heat", "process steam", "process hot water"...co-product of system functioning
    #                 "Control Param" : {"Lower"  : {"1": {"Method" : "Load",    # "electricity", "process heat", "Load", "TIT", "Tstack"
    #                                                    "Limit"  : 0.2},              # [kW], [kW], [-], [K], [K]
    #                                                 },
    #                                    "Upper"  : {"1": {"Method": "Load",    # "electricity", "process heat", "Load", "TIT", "Tstack"
    #                                                    "Limit" : 0.8}         # [kW], [kW], [-], [K], [K]                                                  
    #                                                }
    #                                     }
    #                 }
    
    inp_test_abs = {"Npower": 100,
                    "COP"   : 0.72}

    simulation_hours = 24     
    
    Chp = Chp(inp_test_chp,simulation_hours)            # creating Chp object
    absorber = Absorber(inp_test_abs, simulation_hours) # creating Absorber object
    
    x = np.arange(0,24)       # hours in a day
    days = ['Winter day', 'Spring day','Summer day','Autumn day']
    daily_temperature = {
                         days[0] : [-2.2,-2.7,-3.1,-3.2,-3.6,-3.9,-4.2,-4.3,-4.5,-1.6,1.9,6.3,7.3,7.8,7.5,7.6,6.9,4,0.7,-0.2,-1.1,-1.9,-1.9,-2.6],        # [°C] typical winter day temperatures 
                         days[1] : [4.5,3.9,3.1,2.7,2.6,2.6,2.6,2.6,8,12,12.8,13.7,13.8,14,14.2,14.3,14.3,13.8,13.8,12.2,9.9,7.2,6.3,5.3],                # [°C] typical spring day temperatures
                         days[2] : [18.4,18.1,17.9,17.4,17.4,17.4,17.4,17.4,24,26.2,26.9,27.6,28,28.4,28.6,27.4,29.2,28.7,29.8,27.6,24.6,23,21.9,20.5],   # [°C] typical summer day temperatures
                         days[3] : [20.6,20.6,20.6,20.3,19.5,18.2,16.7,15.7,19,21.7,23.1,23.9,24,24,24.5,24.7,24.3,22.8,21.1,19.6,18.1,16.9,16.1,15.6]    # [°C] typical autumn day temperatures
                         }
    
    np.random.seed(42)
    steam_demand = np.random.uniform(0.5,5.2,simulation_hours)*3600         # [kg/h] creating a random steam demand array 
    electricity_demand = np.random.uniform(1000,6000,simulation_hours)      # [kWh] creating a random steam demand array 
    available_fuel = []
    
    for k in daily_temperature:
        available_fuel = [5000]   # [kg] initial value for available fuel
        
        for h in range(24):
            # Chp.use(h, daily_temperature[k][h], steam_demand[h], electricity_demand[h])                     # not considering available fuel
            # absorber.use(h,Chp.use(h, daily_temperature[k][h], steam_demand[h], electricity_demand[h])[0])  # not considering available fuel 
            Chp.use(h, daily_temperature[k][h], steam_demand[h], electricity_demand[h], available_fuel[h])                     # considering available fuel
            absorber.use(h,Chp.use(h, daily_temperature[k][h], steam_demand[h], electricity_demand[h], available_fuel[h])[0])  # considering available fuel
            consumption = available_fuel[h] - Chp.m_fuel[h]
            available_fuel.append(consumption)
        
        fig, ax = plt.subplots(dpi=600)
        ax.plot(x,Chp.l_bound[:simulation_hours]*6*3600,label ='Chp$_\mathregular{min}$', alpha=0.9)
        ax.plot(x,Chp.u_bound[:simulation_hours]*6*3600,label = 'Chp$_\mathregular{max}$',alpha=0.9)
        ax.scatter(x,steam_demand, label = 'm$_\mathregular{stdem}$', alpha =0.8, c='g',s=70,zorder=10,edgecolors='k', linewidths=0.6)
        ax.scatter(x,Chp.steam[:simulation_hours], label = 'm$_\mathregular{stprod}$', alpha =0.8, c='b',s=22,zorder=10,edgecolors='k', linewidths=0.6)
        ax.grid(zorder=0, alpha= 0.4)
        ax.set_ylabel('m$_{st}$ [kg/h]')
        ax.set_xlabel('Time [h]')
        if max(steam_demand) > max(Chp.u_bound[:simulation_hours]*6*3600):
            ax.fill_between(x,Chp.u_bound[:simulation_hours]*6*3600,max(steam_demand), facecolor = 'orangered',zorder=0, alpha=0.3,label = 'Not satisfed') 
        ax.fill_between(x,Chp.l_bound[:simulation_hours]*6*3600,Chp.u_bound[:simulation_hours]*6*3600, facecolor = 'lightblue',zorder=0, alpha=0.3,label = 'Chp')
        if min(steam_demand)< min(Chp.l_bound[:simulation_hours]*6*3600):
            ax.fill_between(x,Chp.l_bound[:simulation_hours],Chp.l_bound[:simulation_hours]*6*3600, facecolor = 'orangered',zorder=0, alpha=0.3)
        ax.legend(ncol=3,bbox_to_anchor = (0.85,-0.15))
        ax.set_title('Chp system behaviour - ' + k)
        plt.show()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    