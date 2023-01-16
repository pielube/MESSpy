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
        'interest rate': 0-1 [rate/year]
        'inflation rate': -1-1 [rate/year] cost evolution of each carrier
        'investment year': time horizon for which to evaluate the economic analysis (must be a multiple of simulation_year in general.json)
        
    study_case: str name of study case results file.pkl
    reference_case: str name of reference case results file.pkl
        
    structure: dictionary of study case structure (see REC.py)
    structure0: dictionary of reference case structure (see REC.py)
                        
    output: NPV of each location in 'economic_assessment.pkl'
        
    """  
    
    years_factor = int(economic_data['investment years'] / simulation_years) # this factor is useful to match the length of the energy simulation with the length of the economic investment
    
    # open energy balances of study and reference case
    with open('Results/balances_'+study_case+'.pkl', 'rb') as f:        balances = pickle.load(f)        
    with open('Results/balances_'+reference_case+'.pkl', 'rb') as f:        balances0 = pickle.load(f)
        
    results = {}                              # dictionary initialise economic results of each locations
    
    for location_name in structure:           # for reach location
        
        system = structure[location_name]     # location system (see Location.py)
        system0 = structure0[location_name]   # same for reference case
        
        results[location_name] = {}           # dictionary initialise economic results
       
        # initialise cash flow:
        results[location_name]['CF_refcase'] = {  'OeM': np.zeros(economic_data['investment years']),
                                                  'Purchase': {},
                                                  'Sale': {},
                                                  'Refund': np.zeros(economic_data['investment years']),
                                                  'CSC': np.zeros(economic_data['investment years']),
                                                  'Tot': np.zeros(economic_data['investment years'])}        
        results[location_name]['CF_studycase'] = {'OeM': np.zeros(economic_data['investment years']),
                                                  'Purchase': {},
                                                  'Sale': {},
                                                  'Refund': np.zeros(economic_data['investment years']),
                                                  'CSC': np.zeros(economic_data['investment years']),
                                                  'Tot': np.zeros(economic_data['investment years'])}       
        results[location_name]['CF'] = {  'OeM': np.zeros(economic_data['investment years']),
                                          'Purchase': {},
                                          'Sale': {},
                                          'Refund': np.zeros(economic_data['investment years']),
                                          'CSC': np.zeros(economic_data['investment years']),
                                          'Tot': np.zeros(economic_data['investment years'])} 

        I0 = 0 # initialise initial investment           
        
        # each tech has a different sizing parameter:
        size = {'PV': 'peakP', 'inverter': 'peakP', 'battery': 'nominal capacity', 'electrolyzer': 'Npower', 'fuel cell': 'Npower', 'H tank': 'max capacity', 'heatpump': 'nom Pth'} 
        
        for tech_name in system:              # considering each techonlogiy in the location
            if tech_name in economic_data:    # to avoid considering 'electricity demand' as a technology and thus avoiding errors
                
                # Calculate I0 and OeM 
                if tech_name in ['electrolyzer','fuel cell']:
                    I0 += system[tech_name][size[tech_name]]*system[tech_name]['number of modules']*economic_data[tech_name]['total installation cost']  # add technology total installation cost to location I0
                elif tech_name == 'inverter':
                    I0 += system[tech_name][size[tech_name]]*system[tech_name]['number']*economic_data[tech_name]['total installation cost']  # add technology total installation cost to location I0
                else:
                    I0 += system[tech_name][size[tech_name]]*economic_data[tech_name]['total installation cost']  # add technology total installation cost to location I0


                results[location_name]['CF_studycase']['OeM'][:] += - system[tech_name][size[tech_name]]*economic_data[tech_name]['OeM']  # add technology OeM to location OeM

                # replacements 
                if economic_data[tech_name]['lifetime'] == "ageing": # if replacement time is calculated according to ageing
                    with open('Results/ageing_'+study_case+'.pkl', 'rb') as f:
                        age = pickle.load(f)     
                        age = age[location_name][tech_name][0]
                        for a in age:
                            rep_time = int(a/8760)
                            results[location_name]['CF_studycase']['OeM'][rep_time] += - system[tech_name][size[tech_name]]*economic_data[tech_name]['total installation cost'] # subtract technology replacement to location Cash Flow
                else: # if replacement time is given
                    rep_time = economic_data[tech_name]['lifetime']
                    while rep_time < economic_data['investment years']: # if tech_name replacement happens before the end of the simulation
                        results[location_name]['CF_studycase']['OeM'][rep_time] += - system[tech_name][size[tech_name]]*economic_data[tech_name]['total installation cost'] # subtract technology replacement to location Cash Flow
                        rep_time += rep_time
                # NB no refund considered for replacements
                        
                # refund
                if economic_data[tech_name]['refund']['years'] == 0:
                    I0 += - system[tech_name][size[tech_name]]*economic_data[tech_name]['total installation cost']*economic_data[tech_name]['refund']['rate']/100
                else:
                    yearly_refund = system[tech_name][size[tech_name]]*economic_data[tech_name]['total installation cost']*(economic_data[tech_name]['refund']['rate']/100)/economic_data[tech_name]['refund']['years'] # yearly refund [€]
                    refunds = np.zeros(economic_data['investment years']) # array initialise
                    refunds[:min(economic_data['investment years'],economic_data[tech_name]['refund']['years'])] = yearly_refund # array repet yearly refond 
                    results[location_name]['CF_studycase']['Refund'] += refunds # add refund to Cash Flow
                
        # OeM in the reference case (must be subtracted)
        for tech_name in system0: # considering each techonlogies in the locations (reference_case)
            if tech_name in economic_data: # to not consider 'electricity demand' as a technology and avoid bugs
                results[location_name]['CF_refcase']['OeM'][:] += - system0[tech_name][size[tech_name]]*economic_data[tech_name]['OeM']  # add technology OeM to location OeM
                
        # Replacements in the reference case (must be add to CF if a technology is no longer present)
                if tech_name not in system: # if the tech_name is no longer present in the study case
                    if economic_data[tech_name]['lifetime'] < economic_data['investment years']: # if tech_name replacement happens before the end of the simulation
                        results[location_name]['CF_refcase']['OeM'][economic_data[tech_name]['lifetime']] += - system[tech_name][size[tech_name]]*economic_data[tech_name]['total installation cost'] # add technology replacement to location Cash Flow
        
        # energy sold and purchased in study case 
        for carrier in balances[location_name]:                           # for each carrier (electricity, hydrogen, gas, heat)
            if 'grid' in balances[location_name][carrier]:  
                
                if type(economic_data[carrier]['sale']) == str: # if there is the price serie
                    sale_serie = np.tile(pd.read_csv(path+'/energy_price/'+economic_data[carrier]['sale'])['0'].to_numpy(),int(simulation_years))  
                    sold = balances[location_name][carrier]['grid'] * sale_serie
                else: # if the price is always the same 
                    sold = balances[location_name][carrier]['grid']*economic_data[carrier]['sale'] 
               
                sold = np.tile(sold,years_factor)
                sold = np.reshape(sold,(-1,8760))
                results[location_name]['CF_studycase']['Sale'][carrier] = - sold.sum(axis=1,where=sold<0)
                results[location_name]['CF']['Sale'][carrier] = np.zeros(economic_data['investment years'])
                
                if type(economic_data[carrier]['purchase']) == str: # if there is the price series
                    purchase_serie = np.tile(pd.read_csv(path+'/energy_price/'+economic_data[carrier]['purchase'])['0'].to_numpy(),int(simulation_years))  
                    purchase = balances[location_name][carrier]['grid'] * purchase_serie
                else: # if the price is always the same 
                    purchase = balances[location_name][carrier]['grid']*economic_data[carrier]['purchase']
               
                purchase = np.tile(purchase,years_factor)
                purchase = np.reshape(purchase,(-1,8760))
                results[location_name]['CF_studycase']['Purchase'][carrier] = - purchase.sum(axis=1,where=purchase>0)
                results[location_name]['CF']['Purchase'][carrier] = np.zeros(economic_data['investment years'])
            
                
        # energy sold and purchased in reference case 
        for carrier in balances0[location_name]: # for each carrier (electricity, hydrogen, gas, heat)
            if 'grid' in balances0[location_name][carrier]: 
                
                if type(economic_data[carrier]['sale']) == str:                           # if there is the price serie
                    sold = balances0[location_name][carrier]['grid'] * sale_serie
                else: # if the price is always the same 
                    sold = balances0[location_name][carrier]['grid']*economic_data[carrier]['sale'] 
                
                sold = np.tile(sold,years_factor)
                sold = np.reshape(sold,(-1,8760))
                results[location_name]['CF_refcase']['Sale'][carrier] = -sold.sum(axis=1,where=sold<0)
                results[location_name]['CF']['Sale'][carrier] = np.zeros(economic_data['investment years'])

                if type(economic_data[carrier]['purchase']) == str: # if there is the price serie
                    purchase = balances0[location_name][carrier]['grid'] * purchase_serie
                else: # if the price is always the same 
                    purchase = balances0[location_name][carrier]['grid']*economic_data[carrier]['purchase']
              
                purchase = np.tile(purchase,years_factor)
                purchase = np.reshape(purchase,(-1,8760))
                results[location_name]['CF_refcase']['Purchase'][carrier] = -purchase.sum(axis=1,where=purchase>0)
                results[location_name]['CF']['Purchase'][carrier] = np.zeros(economic_data['investment years'])
                      
                
        # REC incentives redistribution
        csc = balances[location_name]['electricity']['collective self consumption']
        inc_pro = - csc * economic_data['REC']['incentives redistribution']['producers']/100 * economic_data['REC']['collective self consumption incentives']
        inc_pro = np.tile(inc_pro,years_factor)
        inc_pro = np.reshape(inc_pro,(-1,8760))    
        results[location_name]['CF_studycase']['CSC'] += inc_pro.sum(axis=1,where=inc_pro>0) 
        
        inc_con = csc * economic_data['REC']['incentives redistribution']['consumers']/100 * economic_data['REC']['collective self consumption incentives']
        inc_con= np.tile(inc_con,years_factor)
        inc_con = np.reshape(inc_con,(-1,8760))
        results[location_name]['CF_studycase']['CSC'] += inc_con.sum(axis=1,where=inc_con>0)   
        

        # CF update considering inflation on each carrier
        for carrier in economic_data['inflation rate']:
            f = economic_data['inflation rate'][carrier]
            
            if carrier in results[location_name]['CF_studycase']['Purchase']:
                for y in range(economic_data['investment years']):
                    results[location_name]['CF_studycase']['Purchase'][carrier][y] = results[location_name]['CF_studycase']['Purchase'][carrier][y]*(1+f)**y
                    
            if carrier in results[location_name]['CF_studycase']['Sale']:
                for y in range(economic_data['investment years']):
                    results[location_name]['CF_studycase']['Sale'][carrier][y] = results[location_name]['CF_studycase']['Sale'][carrier][y]*(1+f)**y
                    
            if carrier in results[location_name]['CF_refcase']['Purchase']:
                for y in range(economic_data['investment years']):
                    results[location_name]['CF_refcase']['Purchase'][carrier][y] = results[location_name]['CF_refcase']['Purchase'][carrier][y]*(1+f)**y
                    
            if carrier in results[location_name]['CF_refcase']['Sale']:
                for y in range(economic_data['investment years']):
                    results[location_name]['CF_refcase']['Sale'][carrier][y] = results[location_name]['CF_refcase']['Sale'][carrier][y]*(1+f)**y
            
            
        # Calculate CF comparing CF_studycase and CF_refcase and total cash flow calculation
        results[location_name]['CF']['OeM'] += results[location_name]['CF_studycase']['OeM'] -results[location_name]['CF_refcase']['OeM'] 
        results[location_name]['CF']['Refund'] += results[location_name]['CF_studycase']['Refund'] -results[location_name]['CF_refcase']['Refund']
        results[location_name]['CF']['CSC'] += results[location_name]['CF_studycase']['CSC'] -results[location_name]['CF_refcase']['CSC']
        
        results[location_name]['CF_studycase']['Tot'] += results[location_name]['CF_studycase']['OeM']
        results[location_name]['CF_studycase']['Tot'] += results[location_name]['CF_studycase']['Refund']
        results[location_name]['CF_studycase']['Tot'] += results[location_name]['CF_studycase']['CSC']
        results[location_name]['CF_refcase']['Tot'] += results[location_name]['CF_refcase']['OeM']
        results[location_name]['CF_refcase']['Tot'] += results[location_name]['CF_refcase']['Refund']
        results[location_name]['CF_refcase']['Tot'] += results[location_name]['CF_refcase']['CSC']
      
        for carrier in results[location_name]['CF_studycase']['Purchase']:
            results[location_name]['CF_studycase']['Tot'] += results[location_name]['CF_studycase']['Purchase'][carrier]
            results[location_name]['CF']['Purchase'][carrier] += results[location_name]['CF_studycase']['Purchase'][carrier]
       
        for carrier in results[location_name]['CF_refcase']['Purchase']:
            results[location_name]['CF_refcase']['Tot'] += results[location_name]['CF_refcase']['Purchase'][carrier]
            results[location_name]['CF']['Purchase'][carrier] += - results[location_name]['CF_refcase']['Purchase'][carrier]
        
        for carrier in results[location_name]['CF_studycase']['Sale']:
            results[location_name]['CF_studycase']['Tot'] += results[location_name]['CF_studycase']['Sale'][carrier]
            results[location_name]['CF']['Sale'][carrier] += results[location_name]['CF_studycase']['Sale'][carrier]
        
        for carrier in results[location_name]['CF_refcase']['Sale']:
            results[location_name]['CF_refcase']['Tot'] += results[location_name]['CF_refcase']['Sale'][carrier]
            results[location_name]['CF']['Sale'][carrier] += - results[location_name]['CF_refcase']['Sale'][carrier]
        
        results[location_name]['CF']['Tot'] = results[location_name]['CF_studycase']['Tot'] - results[location_name]['CF_refcase']['Tot']
        
    
        # calculate NPV
        results[location_name]['NPV'] = np.zeros(economic_data['investment years']+1) # array initialise Net Present Value
        results[location_name]['NPV'][0] = -I0 # NPV at time 0 is - the initial investment
        i = economic_data['interest rate'] # interest rate [%]

        for y in range(1,economic_data['investment years']+1): # for each year
            results[location_name]['NPV'][y] = results[location_name]['NPV'][y-1] + results[location_name]['CF']['Tot'][y-1]/(1+i)**y # calculate NPV 

    # save results in Results/economic_assesment.pkl
    with open(f"Results/economic_assessment_{study_case}.pkl", 'wb') as f:  pickle.dump(results,f) 
        
        
    
            
    
            