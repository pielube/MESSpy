"""
MESSpy - preprocessing

"""
import numpy as np

def PV_surplus(structure):
    for loc in structure:
        if 'heatpump' in structure[loc]:
            structure[loc]['heatpump']['PV surplus'] = True
    return(structure)

def REC_surplus(structure):
    for loc in structure:
        if 'heatpump' in structure[loc]:
            structure[loc]['heatpump']['PV surplus'] = True
            structure[loc]['heatpump']['REC surplus'] = True
    return(structure)
                        
def set_point(structure,t):
    for loc in structure:
        if 'heatpump' in structure[loc]:
            structure[loc]['heatpump']['set point'] = t
    return(structure)

def iTES_size(structure,L):
    for loc in structure:
        if 'heatpump' in structure[loc]:
            structure[loc]['heatpump']['inertial TES volume'] = L
    return(structure)
            
def regulation(structure,t):
    for loc in structure:
        if 'heatpump' in structure[loc]:
            structure[loc]['heatpump']['regulation'] = True
    return(structure)

def PV_size(structure,location_name,p):
    structure[location_name]['PV']['peakP'] = p
    return(structure)

def Bosconi_10x60(structure):
    new_structure = structure.copy()
    for n in np.arange(5):
        n += 1 
        for loc in structure:            
            if loc not in ['potabilizzatore','uffici']:
                new_name = f"{n} {loc}"
                new_structure[new_name] = {'demand': {'electricity': f"{loc}_rp{n}.csv",
                  'heat': f"{loc}_heat.csv",
                  'dhw': f"{loc}_dhw.csv"},
                 'grid': {'electricity': True, 'gas': True},
                 'boiler_ng': {'Ppeak': structure[loc]['boiler_ng']['Ppeak'], 'efficiency': 0.95}}
    return(new_structure)
                
def Bosconi_10x60_HP(structure):
    new_structure = structure.copy()
    for n in np.arange(5):
        n += 1 
        for loc in structure:            
            if loc not in ['potabilizzatore','uffici']:
                new_name = f"{n} {loc}"
                new_structure[new_name] = structure[loc]
    return(new_structure)              
                
                
def add_PV(structure,location_name,peakP,tilt,azimuth):    
    structure[location_name]['PV'] = {"tilt": tilt,"azimuth": azimuth,"losses": 14,"peakP": peakP,	"serie": 2015, "priority": 3}
    return(structure)
                
def add_battery(structure,location_name,kWh):    
    structure[location_name]['battery'] = {"nominal capacity": kWh,"max E-rate":0.5, "efficiency":0.9, "ageing":False, "life cycles":10000, "end life capacity":0.8, "collective":False, 'priority':4}
    return(structure)

def add_HP(structure,location_name,kW):
    structure[location_name]['heatpump'] = {"type":1,"usage": 1,"nom Pth":kW,"t rad heat": 50,"t rad cool": 10,"inertial TES volume": 200, "inertial TES dispersion": 0.36, "PV surplus": True,"REC surplus": False, "priority": 4	}
    if "boiler_ng" in structure[location_name]:
        structure[location_name].pop("boiler_ng")
    if "gas grid" in structure[location_name]:
        structure[location_name].pop("gas grid")
    return(structure)

def change_peakP(structure,location_name,peakP):
    structure[location_name]['PV']['peakP'] = peakP
    structure[location_name]['PV']['serie'] = "GEC_studycase_pv.csv"
    
    return(structure)

def change_peakW(structure,location_name,peakW):
    structure[location_name]['wind']['Npower'] = peakW
    return(structure)

def change_Elesize(structure,location_name,elesize):
    structure[location_name]['electrolyzer']['number of modules'] = elesize
    return(structure)

def change_Htanksize(structure,location_name,tanksize):
    structure[location_name]['H tank']['max capacity'] = tanksize
    return(structure)

def change_electricityprice(energy_market,electricity_price):
    energy_market['electricity']['purchase'] = electricity_price
    return(energy_market)

def change_windeleprice(energy_market,windele_price):
    energy_market['wind electricity']['purchase'] = windele_price
    return(energy_market)
