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

Faraday     =  96485                    # [C/mol]     Faraday constant
R_universal =  8.3144621                # [J/(mol*K)] Molar ideal gas constant
Gamma       =  1.4                      # [-]         Gamma ideal gas = cp/cv 
Nepero      =  2.71828182845904523536   # [-]         Euler's number
AmbTemp     =  288                      # [K]         Standard ambient temperature - 15 °C
Gibbs       = -237.17                   # [kJ/mol]    Gibbs free energy @ T = 25°C p = 101325 Pa
R_h2        =  4124.2                   # [J/(kgK)]   H2 characteristic constant

#%%

'FLUID PROPERTIES'

h2Ndensity   =  0.08988237638480538      # [kg/m^3]    Hydrogen density at Normal conditions (T = 0°C, P = 101325 Pa) -> PropsSI('D', 'T', 273.15, 'P', 101325, 'H2') 
h2Sdensity   =  0.08520493577522305      # [kg/m^3]    Hydrogen density at Standard conditions (T = 15°C, P = 101325 Pa) -> PropsSI('D', 'T', 288.15, 'P', 101325, 'H2') 
h2oNdensity  =  999.8437620819061        # [kg/m^3]    Water density at Normal conditions (T = 0°C, P = 101325 Pa) -> PropsSI('D', 'T', 273.16, 'P', 101325, 'Water') 
h2oSdensity  =  999.1026214670995        # [kg/m^3]    Water density at Standard conditions (T = 15°C, P = 101325 Pa) -> PropsSI('D', 'T', 288.15, 'P', 101325, 'Water') 
LHVh2        =  119.96                   # [MJ/kg]     Hydrogen Lower Heating Value              - https://www.eniscuola.net/mediateca/caratteristiche-dellidrogeno/
LHVh2Vol     =  10.05                    # [MJ/m^3]    Hydrogen Lower Heating Value - Volumetric
LHVh2Mol     =  241.80                   # [kJ/mol]    Hydrogen Lower Heating Value - Molar
HHVh2        =  141.8                    # [MJ/kg]     Hydrogen Lower Heating Value
HHVh2Vol     =  11.89                    # [MJ/m^3]    Hydrogen Lower Heating Value - Volumetric
HHVh2Mol     =  285.83                   # [kJ/mol]    Hydrogen Lower Heating Value - Molar
cp_water     =  4188.460622611614        # [J/kgK]     Water Mass specific constant pressure specific heat
h2oMolMass   =  0.01801528               # [kg/mol]    Water molar mass
h2MolMass    =  2.01588e-3               # [kg/mol]    Hydrogen molar mass
