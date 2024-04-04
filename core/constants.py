"""
COSTANTS MODULE 

    This module contains:
            - Main phisical constants
            - Fluid properties
            
    For the sake of non-ambiguity, all technologies implement the same values by referring to this module
            
"""
#%%

#from CoolProp.CoolProp import PropsSI

#%%

'PHYSICAL & MATHEMATICAL CONSTANTS'

FARADAY     =  96485                       # [C/mol]     Faraday constant
R_UNIVERSAL =  8.3144621                   # [J/(mol*K)] Molar ideal gas constant
GAMMA       =  1.4                         # [-]         Gamma ideal gas = cp/cv 
NEPERO      =  2.71828182845904523536      # [-]         Euler's number
AMBTEMP     =  288                         # [K]         Standard ambient temperature - 15 °C
GIBBS       = -237.17                      # [kJ/mol]    Gibbs free energy @ T = 25°C p = 101325 Pa
R_H2        =  4124.2                      # [J/(kgK)]   H2 characteristic constant
kWh2kJ      =  3600                        # [kJ/kWh]    Conversion factor from kWh to kJ


#%%

'FLUID PROPERTIES'

'Hydrogen'

H2MOLMASS   =  2.01588e-3                 # [kg/mol]     Hydrogen molar mass
H2NDENSITY  =  0.08988237638480538        # [kg/Nm^3]    Hydrogen density at Normal conditions (T = 0°C, P = 101325 Pa) -> PropsSI('D', 'T', 273.15, 'P', 101325, 'H2') 
H2SDENSITY  =  0.08520493577522305        # [kg/Sm^3]    Hydrogen density at Standard conditions (T = 15°C, P = 101325 Pa) -> PropsSI('D', 'T', 288.15, 'P', 101325, 'H2') 
H2NMOLVOL   =  H2MOLMASS/H2NDENSITY       # [Nm^3/mol]   Hydrogen molar volume at Normal conditions (T = 0°C, P = 101325 Pa)
LHVH2       =  119.96                     # [MJ/kg]      Hydrogen Lower Heating Value              - https://www.eniscuola.net/mediateca/caratteristiche-dellidrogeno/
LHVH2VOL    =  LHVH2*H2NDENSITY           # [MJ/Nm^3]    Hydrogen Lower Heating Value - Volumetric
LHVH2MOL    =  H2MOLMASS*LHVH2*1e3        # [kJ/mol]     Hydrogen Lower Heating Value - Molar
LHV_H2      =  LHVH2*(1000/3600)          # [kWh/kg]     Hydrogen Lower Heating Value - kWh-mass
LHV_H2NVOL  =  LHV_H2*H2NDENSITY          # [kWh/Nm^3]   Hydrogen Lower Heating Value - kWh-volume at Normal conditions (T = 0°C, P = 101325 Pa) -> PropsSI('D', 'T', 273.15, 'P', 101325, 'H2') 
HHVH2       =  141.8                      # [MJ/kg]      Hydrogen Higher Heating Value
HHVH2VOL    =  HHVH2*H2NDENSITY           # [MJ/m^3]     Hydrogen Higher Heating Value - Volumetric
HHVH2MOL    =  H2MOLMASS*HHVH2*1e3        # [kJ/mol]     Hydrogen Higher Heating Value - Molar
HHV_H2      =  HHVH2*(1000/3600)          # [kWh/kg]     Hydrogen Higher Heating Value - kWh-mass
HHV_H2VOL   =  HHV_H2*H2NDENSITY          # [kWh/Nm^3]   Hydrogen Higher Heating Value - kWh-volume
H2MOL_S_E   =  130.7                      # [J/K*mol]    Hydrogen Standard Entropy - gaseous phase
CP_H2       =  14.3060                    # [kJ/kgK]     Hydrogen mass specific costant pressure heat capacity (T = 25°C, P = 101325 Pa)
CV_H2       =  10.1796                    # [kJ/kgK]     Hydrogen mass specific costant volume heat capacity (T = 25°C, P = 101325 Pa)

'Water'

H2OMOLMASS  =  0.01801528                 # [kg/mol]     Water molar mass
H2ONDENSITY =  999.8437620819061          # [kg/Nm^3]    Water density at Normal conditions (T = 0°C, P = 101325 Pa) -> PropsSI('D', 'T', 273.16, 'P', 101325, 'Water') 
H2OSDENSITY =  999.1026214670995          # [kg/Sm^3]    Water density at Standard conditions (T = 15°C, P = 101325 Pa) -> PropsSI('D', 'T', 288.15, 'P', 101325, 'Water') 
H2OADENSITY =  998.2071504679284          # [kg/Sm^3]    Water density at 20°C (T = 20°C, P = 101325 Pa) -> PropsSI('D', 'T', 293.15, 'P', 101325, 'Water') 
CP_WATER    =  4.188460622611614          # [kJ/kgK]      Water Mass specific constant pressure specific heat
H2OMOL_S_E  =  188.8                      # [J/K*mol]    Water Standard Entropy - gaseous phase

'Natural Gas'

NGMOLMASS   = 16.04                       # [g/mol]      Gas Methane molar mass
NGNDENSITY  = 0.7174587771429166          # [kg/Nm^3]    Methane density at Normal conditions (T = 0°C, P = 101325 Pa) -> PropsSI('D', 'T', 273.15, 'P', 101325, 'CH4')
NGSDENSITY  = 0.6798343282369096          # [kg/Sm^3]    Methane density at Standard conditions (T = 15°C, P = 101325 Pa) -> PropsSI('D', 'T', 288.15, 'P', 101325, 'CH4')
LHVNG       = 47.451                      # [MJ/kg]      Natural Gas Lower Heating Value
LHVNGVOL    = LHVNG*NGSDENSITY            # [MJ/Sm^3]    Natural Gas Lower Heating Value - Volumetric
LHV_NGSVOL  = LHVNGVOL/3.6                # [kWh/Sm^3]   Natural Gas Lower Heating Value at Standard Conditions

'Diesel'

DIDENSITY   = 837                         # [kg/m^3]     Diesel density at Standard conditions (T = 20°C, P = 101325 Pa) -> PropsSI('D', 'T', 288.15, 'P', 101325, 'H2')
LHVDI       = 43.308                      # [MJ/kg]      Diesel Fuel Lower Heating Value                  
LHVDIVOL    = LHVDI*DIDENSITY             # [MJ/m^3]     Diesel Fuel Lower Heating Value                  

'Coal'

LHVCOAL     = 27.842                      # [MJ/kg]      Coal Fuel Lower Heating Value

'Oxygen'

O2MOLMASS   = 32e-3                       # [kg/mol]     Oxygen molar mass
O2MOL_S_E   = 205.1                       # [J/K*mol]    Oxygen Standard Entropy - gaseous phase

'Air'

AIRMOLMASS   = 28.96547e-3                # [kg/mol]     Air molar mass
AIRSDENSITY  = 1.225                      # [kg/Sm^3]    Air density at Standard conditions (T = 15°C, P = 101325 Pa) -> PropsSI('D', 'T', 273.15, 'P', 101325, 'Air')
CP_AIR       = 1.0063                     # [kJ/kgK]     Air mass specific costant pressure specific heat (T = 25°C, P = 101325 Pa)
CV_AIR       = 0.7178                     # [kJ/kgK]     Air mass specific costant volume specific heat (T = 25°C, P = 101325 Pa) 

'Steam'

H1_STEAM800  = 4159.9                     # [kJ/kg]      Steam mass specific enthalpy @ T = 800°C, P = 116000 Pa


#%%

'TIME UNITS'

MINUTES_HOUR    = 60                      # [min/h]     Number of minutes in one hour 
HOURS_YEAR      = 8760                    # [h/year]    Number of hours in one year
MINUTES_YEAR    = MINUTES_HOUR*HOURS_YEAR # [min/year]  Number of minutes in one year 

