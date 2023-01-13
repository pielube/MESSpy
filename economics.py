import numpy as np
import pandas as pd
import pickle

def NPV(structure,structure0,study_case,reference_case,economic_data,simulation_years,path):
    """
    Economic assesment 
    
    economic_data: dictionary
        'tech_name': dictionary (repeat for each installed technologies: PV, battery, electrolyzer, H tank, fuel cell)
            'total installation cost': tech_name cost + installation costs [€]
            'OeM': operations and maintenance costs [€/kW/y]
            'lifetime': time after which the tech_name must be replaced [years]            
            'refund': dictionary incentives definition
                'rate': rate of initial investemt that will be refunded [%]
                'time': refunding time [years]
        'REC': dictionary REC economic parameters definition
            'collective self consumption incentives': [€/kWh]
            'incentives redistribution': 0-100 how the incentives are redistributed between prosumers, consumers and REC manger
        'carrier_name': dictionary (repeat for reach considered carrier: electricity, hydrogen, gas)
            'purchase': [€/kWh] (electricity and gas) or [€/kg] (hydrogen)
            'sales': [€/kWh] (electricity and gas) or [€/kg] (hydrogen)
        'interest rate': 0-1
        'inflation rate': 0-1
        'investment year': time horizon for which to evaluate the economic analysis (must be a multiple of simulation_year in general.json)
        
    study_case: str name of study case results file.pkl
    reference_case: str name of reference case results file.pkl
        
    structure: dictionary of study case structure (see REC.py)
    structure0: dictionary of reference case structure (see REC.py)
                        
    output: NPV of each location in 'economic_assessment.pkl'
        
    """  
    
    years_factor = int(economic_data['investment years'] / simulation_years) # this factor is useful to match the length of the energy simulation with the length of the economic investment
    
    # open energy balances of study and reference case
    with open('Results/balances_'+study_case+'.pkl', 'rb') as f:
        rec = pickle.load(f)        
    with open('Results/balances_'+reference_case+'.pkl', 'rb') as f:
        rec0 = pickle.load(f)
        
    results = {}                              # dictionary initialise economic results of each locations
    
    for location_name in structure:           # for reach location
        
        results[location_name] = {}           # dictionary initialise economic results
        system = structure[location_name]     # location system (see Location.py)
        system0 = structure0[location_name]   # same for reference case
        
        results[location_name]['Annual cash flow'] = {}   # dictionary initialise annual cash flow
                
        CF = np.zeros(economic_data['investment years'])  # array initialize Casch Flow
        I0 = 0                                            # initialise initial investment     
        OeM = 0                                           # initialise OeM        
        
        # each tech has a different sizing parameter:
        size = {'PV': 'peakP', 'battery': 'nominal capacity', 'electrolyzer': 'Npower', 'fuel cell': 'Npower', 'H tank': 'max capacity', 'heatpump': 'nom Pth'} 
        
        for tech_name in system:              # considering each techonlogiy in the location
            if tech_name in economic_data:    # to avoid considering 'electricity demand' as a technology and thus avoiding errors
                
                # Calculate I0 and OeM 
                I0 += system[tech_name][size[tech_name]]*economic_data[tech_name]['total installation cost']  # add technology total installation cost to location I0
                OeM += system[tech_name][size[tech_name]]*economic_data[tech_name]['OeM']                     # add technology OeM to location OeM
               
                # replacements 
                if economic_data[tech_name]['lifetime'] == "ageing": # if replacement time is calculated according to ageing
                    with open('Results/ageing_'+study_case+'.pkl', 'rb') as f:
                        age = pickle.load(f)     
                        age = age[location_name][tech_name][0]
                        for a in age:
                            rep_time = int(a/8760)
                            CF[rep_time] += - system[tech_name][size[tech_name]]*economic_data[tech_name]['total installation cost'] # subtract technology replacement to location Cash Flow
                else: # if replacement time is given
                    rep_time = economic_data[tech_name]['lifetime']
                    while rep_time < economic_data['investment years']: # if tech_name replacement happens before the end of the simulation
                        CF[rep_time] += - system[tech_name][size[tech_name]]*economic_data[tech_name]['total installation cost'] # subtract technology replacement to location Cash Flow
                        rep_time += rep_time
                # NB no refund considered for replacements
                        
                # refund
                if economic_data[tech_name]['refund']['years'] == 0:
                    I0 += - system[tech_name][size[tech_name]]*economic_data[tech_name]['total installation cost']*economic_data[tech_name]['refund']['rate']/100
                else:
                    yearly_refund = system[tech_name][size[tech_name]]*economic_data[tech_name]['total installation cost']*(economic_data[tech_name]['refund']['rate']/100)/economic_data[tech_name]['refund']['years'] # yearly refund [€]
                    refunds = np.zeros(economic_data['investment years']) # array initialise
                    refunds[:min(economic_data['investment years'],economic_data[tech_name]['refund']['years'])] = yearly_refund # array repet yearly refond 
                    CF = CF + refunds # add refund to Cash Flow
                
        # OeM in the reference case (must be subtracted)
        for tech_name in system0: # considering each techonlogies in the locations (reference_case)
            if tech_name in economic_data: # to not consider 'electricity demand' as a technology and avoid bugs
                OeM += - system[tech_name][size[tech_name]]*economic_data[tech_name]['OeM'] # subtract technology OeM to location OeM
                
        # Replacements in the reference case (must be add to CF if a technology is no longer present)
                if tech_name not in system: # if the tech_name is no longer present in the study case
                    if economic_data[tech_name]['lifetime'] < economic_data['investment years']: # if tech_name replacement happens before the end of the simulation
                        CF[economic_data[tech_name]['lifetime']] += + system[tech_name][size[tech_name]]*economic_data[tech_name]['total installation cost'] # add technology replacement to location Cash Flow

        CF[:] += - OeM   # OeM every year
        results[location_name]['Annual cash flow']['OeM'] = -OeM     # OeM every year
        
        # energy sold and purchased in study case 
        results[location_name]['Annual cash flow']['Sale'] = {}      #initialise
        results[location_name]['Annual cash flow']['Purchase'] = {}  #initialise
            
        for carrier in rec[location_name]:                           # for each carrier (electricity, hydrogen, gas, heat)
            if 'grid' in rec[location_name][carrier]:  
                
                results[location_name]['Annual cash flow']['Sale'][carrier] = 0        #initialise
                results[location_name]['Annual cash flow']['Purchase'][carrier] = 0    #initialise                
                
                if type(economic_data[carrier]['sale']) == str: # if there is the price serie
                    sale_serie = np.tile(pd.read_csv(path+'/energy_price/'+economic_data[carrier]['sale'])['0'].to_numpy(),int(simulation_years))  
                    sold = rec[location_name][carrier]['grid'] * sale_serie
                else: # if the price is always the same 
                    sold = rec[location_name][carrier]['grid']*economic_data[carrier]['sale'] 
               
                sold = np.tile(sold,years_factor)
                sold = np.reshape(sold,(-1,8760))
                CF = CF - sold.sum(axis=1,where=sold<0)
                results[location_name]['Annual cash flow']['Sale'][carrier] += - sold.sum(axis=1,where=sold<0)
          
                if type(economic_data[carrier]['purchase']) == str: # if there is the price series
                    purchase_serie = np.tile(pd.read_csv(path+'/energy_price/'+economic_data[carrier]['purchase'])['0'].to_numpy(),int(simulation_years))  
                    purchase = rec[location_name][carrier]['grid'] * purchase_serie
                else: # if the price is always the same 
                    purchase = rec[location_name][carrier]['grid']*economic_data[carrier]['purchase']
               
                purchase = np.tile(purchase,years_factor)
                purchase = np.reshape(purchase,(-1,8760))
                CF = CF - purchase.sum(axis=1,where=purchase>0)
                results[location_name]['Annual cash flow']['Purchase'][carrier] += - purchase.sum(axis=1,where=purchase>0)
            
                
        # energy sold and purchased in reference case 
        for carrier in rec0[location_name]: # for each carrier (electricity, hydrogen, gas, heat)
            if 'grid' in rec0[location_name][carrier]: 
                
                if not carrier in results[location_name]['Annual cash flow']['Sale']:
                    results[location_name]['Annual cash flow']['Sale'][carrier] = 0       #initialise
                    results[location_name]['Annual cash flow']['Purchase'][carrier] = 0   #initialise 
               
                if type(economic_data[carrier]['sale']) == str:                           # if there is the price serie
                    sold = rec0[location_name][carrier]['grid'] * sale_serie
                else: # if the price is always the same 
                    sold = rec0[location_name][carrier]['grid']*economic_data[carrier]['sale'] 
                
                sold = np.tile(sold,years_factor)
                sold = np.reshape(sold,(-1,8760))
                CF = CF + sold.sum(axis=1,where=sold<0)
                results[location_name]['Annual cash flow']['Sale'][carrier] += sold.sum(axis=1,where=sold<0)

                if type(economic_data[carrier]['purchase']) == str: # if there is the price serie
                    purchase = rec0[location_name][carrier]['grid'] * purchase_serie
                else: # if the price is always the same 
                    purchase = rec0[location_name][carrier]['grid']*economic_data[carrier]['purchase']
              
                purchase = np.tile(purchase,years_factor)
                purchase = np.reshape(purchase,(-1,8760))
                CF = CF + purchase.sum(axis=1,where=purchase>0)
                results[location_name]['Annual cash flow']['Purchase'][carrier] += purchase.sum(axis=1,where=purchase>0)
                      
                
        # REC incentives redistribution
        results[location_name]['Annual cash flow']['CSC'] = 0 # initialise
        csc = rec[location_name]['electricity']['collective self consumption']
        inc_pro = - csc * economic_data['REC']['incentives redistribution']['producers']/100 * economic_data['REC']['collective self consumption incentives']
        inc_pro = np.tile(inc_pro,years_factor)
        inc_pro = np.reshape(inc_pro,(-1,8760))
        CF = CF + inc_pro.sum(axis=1,where=inc_pro>0)       
        results[location_name]['Annual cash flow']['CSC'] += inc_pro.sum(axis=1,where=inc_pro>0) 
        
        inc_con = csc * economic_data['REC']['incentives redistribution']['consumers']/100 * economic_data['REC']['collective self consumption incentives']
        inc_con= np.tile(inc_con,years_factor)
        inc_con = np.reshape(inc_con,(-1,8760))
        CF = CF + inc_con.sum(axis=1,where=inc_con>0)   
        results[location_name]['Annual cash flow']['CSC'] += inc_con.sum(axis=1,where=inc_con>0)   
        
        # calculate NPV
        results[location_name]['NPV'] = np.zeros(economic_data['investment years']+1) # array initialise Net Present Value
        results[location_name]['NPV'][0] = -I0 # NPV at time 0 is - the initial investment
        i = economic_data['interest rate'] # interest rate [%]
        f = economic_data['inflation rate'] # inflation rate [%]
        for y in range(1,economic_data['investment years']+1): # for each year
            results[location_name]['NPV'][y] = results[location_name]['NPV'][y-1] + CF[y-1]*(1+f)**y/(1+i)**y # calculate NPV 
                             
            
        results[location_name]['CF_tot'] = CF

    # save results in Results/economic_assesment.pkl
    with open(f"Results/economic_assessment_{study_case}.pkl", 'wb') as f:
        pickle.dump(results,f) 
        
        
    
            
    
            