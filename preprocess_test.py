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
