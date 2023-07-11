"""
MESSpy - preprocessing

"""

def change_peakP(structure,location_name,peakP):
    structure[location_name]['PV']['peakP'] = peakP
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
