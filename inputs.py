#### ########################################################################## INPUT definition

simulation_years = 30 # int 

# new-installation study case structure
structure = {'p1': {'demand': {'electricity': 'EnCo_kWh_U0.csv'}, # hourly time series 8760 values [kWh]
                    'PV': {'latitude': 43.7, # float
                           'longitude': 11.2, # float
                           'tilt': 30, # float surface tilt [deg]
                           'azimuth': 0, # float azimuth angle 0 = south, 180 = north [deg]
                           'reference year': 2015, # int year [2005 - 2015] used for output data and to get data from PVGIS if TMY = False 
                           'TMY': True, # bool true if data of TMY is to be used
                           'losses': 0.1, # float losses in cables, power inverters, dirt (sometimes snow), over the years loss of power [%]
                           'peakP': 3}, # float peak DC power [kWp]
                    'battery': {'max capacity': 5, # float [kWh]
                                'ageing': False}, # bool true if aging has to be calculated
                    },
             'p2': {'demand': {'electricity': 'EnCo_kWh_U1.csv'}, # hourly time series 8760 values [kWh]
                    'PV': {'latitude': 43.7, # float
                           'longitude': 11.2, # float
                           'tilt': 30, # float surface tilt [deg]
                           'azimuth': 0, # float azimuth angle 0 = south, 180 = north [deg]
                           'reference year': 2015, # int year [2005 - 2015] used for output data and to get data from PVGIS if TMY = False 
                           'TMY': True, # bool true if data of TMY is to be used
                           'losses': 0.1, # float losses in cables, power inverters, dirt (sometimes snow), over the years loss of power [%]
                           'peakP': 1.5}, # float peak DC power [kWp]
                    'electrolyzer': {'Npower': 1, # float nominal power [kW]
                                     'stack model': 'Enapter 2.1'}, # str 'Enapter 2.1' or 'McLyzer 800' are aviable
                    'fuel cell': {'Npower': 3}, # float nominal power [kW]
                    'H tank': {'max capacity': 6, # float [kg]
                               'pressure': 30} # float [bar]
                    }, # hourly time series 8760 values [kWh]
             'c2': {'demand': {'electricity': 'EnCo_kWh_U2.csv'}}, # hourly time series 8760 values [kWh]
             }
      
# pre-installation reference case structure (must have the same location (same name) of the study_case)
structure0 = {'p1': {'demand': {'electricity': 'EnCo_kWh_U0.csv'}}, # hourly time series 8760 values [kWh]
              'p2': {'demand': {'electricity': 'EnCo_kWh_U1.csv'}}, # hourly time series 8760 values [kWh]
              'c2': {'demand': {'electricity': 'EnCo_kWh_U2.csv'}}, # hourly time series 8760 values [kWh]
              }

# economic parameters
economic_data = {'battery': {'total installation cost': 600, # tech cost + installation costs [€]
                             'OeM': 0, # operations and maintenance costs [€/kW/y]
                             'lifetime': 15}, # time after which the technology must be replaced [years]
                 'PV': {'total installation cost': 1250, # [€]
                        'OeM': 10, # [€/kWh/y]
                        'lifetime': 20}, # [years]
                 'electrolyzer': {'total installation cost': 1500, # [€]
                                  'OeM': 0, # [€/kW/y]
                                  'lifetime': 20}, # [years]
                 'H tank': {'total installation cost': 400, # [€]
                            'OeM': 0, # [€/kg/y]
                            'lifetime': 30}, # [years]
                 'fuel cell': {'total installation cost': 1500, # [€]
                               'OeM': 0, # [€/kW/y]
                               'lifetime': 20}, # [years]
                 'refound': {'rate': 50, # rate of initial investemt that will be refound [%]
                             'time': 10}, # refound time [years]
                 'CER': {'collective self consumption incentives': 0.118, # [€/kWh]
                         'costitution cost': 0, # [€]
                         'OeM': 0}, # [€/y]
                 'electricity': {'purchase': 0.22, # [€/kWh]
                                 'sales': 0.04}, # [€/kWh]
                 'hydrogen': {'purchase': 0, # [€/kg]
                                 'sales': 0}, # [€/kWg   ]   
                 'gas': {'purchase': 0, # [€/kg]
                         'sales': 0}, # [€/kWg]
                 'heat': {'purchase': 0, # [€/kg]
                          'sales': 0}, # [€/kWg] 
                 'interest rate': 0.05,} # [%]

