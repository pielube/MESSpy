"""
MESSpy - preprocessing

"""

def change_peakP(structure,location_name,peakP):
    structure[location_name]['PV']['peakP'] = peakP
    structure[location_name]['inverter']['peakP'] = peakP
    return(structure)

