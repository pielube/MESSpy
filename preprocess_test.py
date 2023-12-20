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

def change_O2tanksize(structure,location_name,O2_tanksize):
    structure[location_name]['O2 tank']['max capacity'] = O2_tanksize
    return(structure)

def change_O2demand(structure,location_name,amount):
    structure[location_name]['oxygen demand']['amount'] = amount
    return(structure)

def change_O2tankprice(tech_cost,O2tank_price):
    tech_cost['O2 tank']['cost per unit'] = O2tank_price
    return(tech_cost)

def change_O2sellingprice(energy_market,O2price):
    energy_market['oxygen']['sale'] = O2price
    return(energy_market)