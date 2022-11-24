"""
COSTANTS MODULE 

    This module contains:
            - Main phisical constants
            - Fluid properties          
            
    For the sake of non-ambiguity, all technologies implement the same values by referring to this module
            
"""
#%%

# from CoolProp.CoolProp import PropsSI

#%%

'PHYSICAL & MATHEMATICAL CONSTANTS'

FARADAY     =  96485                    # [C/mol]     Faraday constant
R_UNIVERSAL =  8.3144621                # [J/(mol*K)] Molar ideal gas constant
GAMMA       =  1.4                      # [-]         Gamma ideal gas = cp/cv 
NEPERO      =  2.71828182845904523536   # [-]         Euler's number
AMBTEMP     =  288                      # [K]         Standard ambient temperature - 15 °C
GIBBS       = -237.17                   # [kJ/mol]    Gibbs free energy @ T = 25°C p = 101325 Pa
R_H2        =  4124.2                   # [J/(kgK)]   H2 characteristic constant

#%%

'FLUID PROPERTIES'

H2NDENSITY   =  0.08988237638480538      # [kg/m^3]    Hydrogen density at Normal conditions (T = 0°C, P = 101325 Pa) -> PropsSI('D', 'T', 273.15, 'P', 101325, 'H2') 
H2SDENSITY   =  0.08520493577522305      # [kg/m^3]    Hydrogen density at Standard conditions (T = 15°C, P = 101325 Pa) -> PropsSI('D', 'T', 288.15, 'P', 101325, 'H2') 
H2ONDENSITY  =  999.8437620819061        # [kg/m^3]    Water density at Normal conditions (T = 0°C, P = 101325 Pa) -> PropsSI('D', 'T', 273.16, 'P', 101325, 'Water') 
H2OSDENSITY  =  999.1026214670995        # [kg/m^3]    Water density at Standard conditions (T = 15°C, P = 101325 Pa) -> PropsSI('D', 'T', 288.15, 'P', 101325, 'Water') 
LHVH2        =  119.96                   # [MJ/kg]     Hydrogen Lower Heating Value              - https://www.eniscuola.net/mediateca/caratteristiche-dellidrogeno/
LHVH2VOL     =  10.05                    # [MJ/m^3]    Hydrogen Lower Heating Value - Volumetric
LHVH2MOL     =  241.80                   # [kJ/mol]    Hydrogen Lower Heating Value - Molar
LHV_H2       =  33.33                    # [kWh/kg]    Hydrogen Lower Heating Value - kWh-mass
LHV_H2VOL    =  2.78                     # [kWh/Nm^3]  Hydrogen Lower Heating Value - kWh-volume
HHVH2        =  141.8                    # [MJ/kg]     Hydrogen Higher Heating Value
HHVH2VOL     =  11.89                    # [MJ/m^3]    Hydrogen Higher Heating Value - Volumetric
HHVH2MOL     =  285.83                   # [kJ/mol]    Hydrogen Higher Heating Value - Molar
HHV_H2       =  39.41                    # [kWh/kg]    Hydrogen Higher Heating Value - kWh-mass
HHV_H2VOL    =  3.28                     # [kWh/Nm^3]  Hydrogen Higher Heating Value - kWh-volume
CP_WATER     =  4188.460622611614        # [J/kgK]     Water Mass specific constant pressure specific heat
H2OMOLMASS   =  0.01801528               # [kg/mol]    Water molar mass
H2MOLMASS    =  2.01588e-3               # [kg/mol]    Hydrogen molar mass
