import numpy as np
import pickle

def NPV(structure,structure0,study_case,reference_case,economic_data,simulation_years):
    """
    Economic assesment 
    
    economic_data: dictionary
        'tech_name': dictionary (repeat for each installed technologies: PV, battery, electrolyzer, H tank, fuel cell)
            'total installation cost': tech_name cost + installation costs [€]
            'OeM': operations and maintenance costs [€/kW/y]
            'lifetime': time after which the tech_name must be replaced [years]            
        'refound': dictionary inventives definition
            'rate': rate of initial investemt that will be refound [%]
            'time': refound time [years]
        'CER': dictionary REC economic parameters definition
            'collective self consumption incentives': [€/kWh]
            'costitution cost': [€]
            'OeM': [€/y]
        'carrier_name': dictionary (repeat for reach considered carriers: electricity, hydrogen, gas, heat)
            'purchase': [€/kWh] or [€/kg] 
            'sales': [€/kWh] or [€/kg]   
        'interest rate': [%]
        
    study_case: str name of study case results file.pkl
    reference_case: str name of reference case results file.pkl
        
    structure: dictionary of study case structure (see REC.py)
    structure0: dictionary of reference case structure (see REC.py)
                        
    output: NPV of each location in 'economic_assessment.pkl'
        
    """
    
    # open energy balances of study and reference case
    with open('Results/balances_'+study_case+'.pkl', 'rb') as f:
        rec = pickle.load(f)        
    with open('Results/balances_'+reference_case+'.pkl', 'rb') as f:
        rec0 = pickle.load(f)
        
    results = {} # dictionary initialise economic results of each locations
    
    for location_name in structure: # for reach locations
        
        results[location_name] = {} # dictionary initialise economic results
        system = structure[location_name] # location system (see Location.py)
        system0 = structure0[location_name] # same for reference case
                
        CF = np.zeros(simulation_years) # array initialize Casch Flow
        I0 = 0 # initialise initial investment     
        OeM = 0 # initialise OeM        
        
        # each tech has a different sizing parameter:
        size = {'PV': 'peakP', 'battery': 'max capacity', 'electrolyzer': 'Npower', 'fuel cell': 'Npower', 'H tank': 'max capacity'} 
        
        for tech_name in system: # considering each techonlogies in the location
            if tech_name in economic_data: # to not consider 'electricity demand' as a technology and avoid bugs
                
                # Calculate I0, OeM and replacements
                I0 += system[tech_name][size[tech_name]]*economic_data[tech_name]['total installation cost'] # add technology total installation cost to location I0
                OeM += system[tech_name][size[tech_name]]*economic_data[tech_name]['OeM'] # add technology OeM to location OeM
                if economic_data[tech_name]['lifetime'] < simulation_years: # if tech_name replacement happens before the end of the simulation
                    CF[economic_data[tech_name]['lifetime']] += - system[tech_name][size[tech_name]]*economic_data[tech_name]['total installation cost'] # subtract technology replacement to location Cash Flow
               
        # OeM in the reference case (must be subtracted)
        for tech_name in system0: # considering each techonlogies in the locations (reference_case)
            if tech_name in economic_data: # to not consider 'electricity demand' as a technology and avoid bugs
                OeM += - system[tech_name][size[tech_name]]*economic_data[tech_name]['OeM'] # subtract technology OeM to location OeM
                
        # Replacements in the reference case (must be add to CF if a technology is no longer present)
                if tech_name not in system: # if the tech_name is no longer present in the study case
                    if economic_data[tech_name]['lifetime'] < simulation_years: # if tech_name replacement happens before the end of the simulation
                        CF[economic_data[tech_name]['lifetime']] += + system[tech_name][size[tech_name]]*economic_data[tech_name]['total installation cost'] # add technology replacement to location Cash Flow

        CF[:] += - OeM # OeM every year
        
        # energy sold and purchased in study case 
        for carrier in rec[location_name]: # for each carrier (electricity, hydrogen, gas, heat, cool)
            if 'grid' in rec[location_name][carrier]: 
                sold = rec[location_name][carrier]['grid']*economic_data[carrier]['sales']
                sold = np.reshape(sold,(-1,8760))
                CF = CF - sold.sum(axis=1,where=sold<0)
            if 'grid' in rec[location_name][carrier]:
                purchase = rec[location_name][carrier]['grid']*economic_data[carrier]['purchase']
                purchase = np.reshape(purchase,(-1,8760))
                CF = CF - purchase.sum(axis=1,where=purchase>0)
                
        # energy sold and purchased in reference case 
        for carrier in rec0[location_name]: # for each carrier (electricity, hydrogen, gas, heat)
            if 'grid' in rec0[location_name][carrier]: 
                sold = rec0[location_name][carrier]['grid']*economic_data[carrier]['sales']
                sold = np.reshape(sold,(-1,8760))
                CF = CF + sold.sum(axis=1,where=sold<0)
            if 'grid' in rec0[location_name][carrier]:
                purchase = rec0[location_name][carrier]['grid']*economic_data[carrier]['purchase']
                purchase = np.reshape(purchase,(-1,8760))
                CF = CF + purchase.sum(axis=1,where=purchase>0)
                      
        # refound
        yearly_refound = I0*(economic_data['refound']['rate']/100)/economic_data['refound']['time'] # yearly refound [€]
        refounds = np.zeros(simulation_years) # array initialise
        refounds[:min(simulation_years,economic_data['refound']['time'])] = yearly_refound # array repet yearly refond for 
        CF = CF + refounds # add refound to Cash Flow
        
        # calculate NPV
        results[location_name]['NPV'] = np.zeros(simulation_years) # array initialise Net Present Value
        results[location_name]['NPV'][0] = -I0 # NPV at time 0 is - the initial investment
        i = economic_data['interest rate'] # interest rate [%]
        for y in range(1,simulation_years): # for each year
            results[location_name]['NPV'][y] = results[location_name]['NPV'][y-1] + CF[y-1]/(1+i)**y # calculate NPV 
                   
    # RES incentives redistribution???
    
    # save results in Results/economic_assesment.pkl
    with open('Results/economic_assessment.pkl', 'wb') as f:
        pickle.dump(results,f) 
        
    
            
    
            