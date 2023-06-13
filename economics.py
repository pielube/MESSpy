import numpy as np
import pandas as pd
import pickle
import os
import json
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

def NPV(file_studycase,name_studycase,name_refcase,economic_data,simulation_years,path,name_economic):  ### WIP da implementare l'opzione in cui non si possiede la risorsa rinnovabile e di conseguenza vanno rivisti i bilanci di immissione in rete poichè non sono flussi di cassa i cui dispongo io in quanto location
    """
    Economic assesment 
    
    economic_data: dictionary
        'REC': dictionary REC economic parameters definition
            'collective self consumption incentives': [€/kWh]
            'incentives redistribution': 0-100 how the incentives are redistributed between prosumers, consumers and REC manger
        'carrier_name': dictionary (repeat for reach considered carrier: electricity, hydrogen, gas)
            'purchase': [€/kWh] (electricity and gas) or [€/kg] (hydrogen)
            'sales': [€/kWh] (electricity and gas) or [€/kg] (hydrogen)
        'interest rate': 0-1 [rate/year]
        'inflation rate': -1-1 [rate/year] cost evolution of each carrier
        'investment year': time horizon for which to evaluate the economic analysis (must be a multiple of simulation_year in general.json)
        
    name_studycase: str name of study case results file.pkl
    name_refcase: str name of reference case results file.pkl
        
    structure: dictionary of study case structure (see REC.py)
    structure0: dictionary of reference case structure (see REC.py)
                        
    output: NPV of each location in 'economic_assessment.pkl'
        
    """  
    
    # open file study_case
    with open(os.path.join(path,f"{file_studycase}.json"),'r') as f:        studycase = json.load(f)
    years_factor = int(economic_data['investment years'] / simulation_years) # this factor is useful to match the length of the energy simulation with the length of the economic investment
    
    # open energy balances of study and reference case
    with open('Results/balances_'+name_studycase+'.pkl', 'rb') as f:        balances = pickle.load(f)        
    with open('Results/balances_'+name_refcase+'.pkl', 'rb') as f:        balances0 = pickle.load(f)
    
    # open cost of componenets of studycase and refcase
    with open('Results/tech_cost_'+name_studycase+'.pkl', 'rb') as f:        tc = pickle.load(f)        
    with open('Results/tech_cost_'+name_refcase+'.pkl', 'rb') as f:        tc0 = pickle.load(f)
        
    results = {}                              # dictionary initialise economic results of each locations
    
    for location_name in tc:           # for reach location
        
        results[location_name] = {}           # dictionary initialise economic results
       
        # initialise cash flow:
        results[location_name]['CF_refcase'] = {  'OeM': np.zeros(economic_data['investment years']),
                                                  'Initial/Final Tank level': np.zeros(economic_data['investment years']),                                                                                                                                                
                                                  'Purchase': {},
                                                  'Sale': {},
                                                  'Refund': np.zeros(economic_data['investment years']),
                                                  'CSC': np.zeros(economic_data['investment years']),
                                                  'Tot': np.zeros(economic_data['investment years'])}        
        results[location_name]['CF_studycase'] = {'OeM': np.zeros(economic_data['investment years']),
                                                  'Initial/Final Tank level': np.zeros(economic_data['investment years']),                                                                                                                                                
                                                  'Purchase': {},
                                                  'Sale': {},
                                                  'Refund': np.zeros(economic_data['investment years']),
                                                  'CSC': np.zeros(economic_data['investment years']),
                                                  'Tot': np.zeros(economic_data['investment years'])}       
        results[location_name]['CF'] = {  'OeM': np.zeros(economic_data['investment years']),
                                          'Initial/Final Tank level': np.zeros(economic_data['investment years']),                                                                                                                                                
                                          'Purchase': {},
                                          'Sale': {},
                                          'Refund': np.zeros(economic_data['investment years']),
                                          'CSC': np.zeros(economic_data['investment years']),
                                          'Tot': np.zeros(economic_data['investment years'])} 

        results[location_name]['I0'] = {} # initialise initial investment           
        
        
        for tech_name in tc[location_name]:              # considering each techonlogiy in the location
        
            results[location_name]['I0'][tech_name] = tc[location_name][tech_name]['total cost'] # I0
            results[location_name]['CF_studycase']['OeM'][:] += - tc[location_name][tech_name]['OeM'] # OeM

            # replacements 
            if tc[location_name][tech_name]['replacement']['years'] == "ageing": # if replacement year is calculated according to ageing
                with open('Results/ageing_'+name_studycase+'.pkl', 'rb') as f:
                    age = pickle.load(f)     
                    age = age[location_name][tech_name][0]
                    for a in age:
                        rep_time = int(a/8760)
                        results[location_name]['CF_studycase']['OeM'][rep_time] += - results[location_name]['I0'][tech_name] * tc[location_name][tech_name]['replacement']['rate']/100 # subtract technology replacement to location Cash Flow
            else: # if replacement time is given
                rep_time = tc[location_name][tech_name]['replacement']['years']
                while rep_time < economic_data['investment years']: # if tech_name replacement happens before the end of the simulation
                    results[location_name]['CF_studycase']['OeM'][rep_time] += - results[location_name]['I0'][tech_name] * tc[location_name][tech_name]['replacement']['rate']/100 # subtract technology replacement to location Cash Flow
                    rep_time += tc[location_name][tech_name]['replacement']['years']
            # NB no refund considered for replacements
                    
            # refund
            if tc[location_name][tech_name]['refund']['years'] == 0:
                results[location_name]['I0'][tech_name] = results[location_name]['I0'][tech_name]*(100-tc[location_name][tech_name]['refund']['rate'])/100
            else:
                yearly_refund = results[location_name]['I0'][tech_name]*tc[location_name][tech_name]['refund']['rate']/100 / tc[location_name][tech_name]['refund']['years'] # yearly refund [€]
                refunds = np.zeros(economic_data['investment years']) # array initialise
                refunds[:min(economic_data['investment years'],tc[location_name][tech_name]['refund']['years'])] = yearly_refund # array repet yearly refond 
                results[location_name]['CF_studycase']['Refund'] += refunds # add refund to Cash Flow
            
        for tech_name in tc0[location_name]:
            
            results[location_name]['CF_refcase']['OeM'][:] += - tc0[location_name][tech_name]['OeM'] # OeM
            
            # replacements 
            if tc0[location_name][tech_name]['replacement']['years'] == "ageing": # if replacement year is calculated according to ageing
                with open('Results/ageing_'+name_refcase+'.pkl', 'rb') as f:
                    age = pickle.load(f)     
                    age = age[location_name][tech_name][0]
                    for a in age:
                        rep_time = int(a/8760)
                        results[location_name]['CF_refcase']['OeM'][rep_time] += - results[location_name]['I0'][tech_name] * tc0[location_name][tech_name]['replacement']['rate']/100 # subtract technology replacement to location Cash Flow
            else: # if replacement time is given
                rep_time = tc0[location_name][tech_name]['replacement']['years']
                while rep_time < economic_data['investment years']: # if tech_name replacement happens before the end of the simulation
                    results[location_name]['CF_refcase']['OeM'][rep_time] += - tc0[location_name][tech_name]['total cost'] * tc0[location_name][tech_name]['replacement']['rate']/100 # subtract technology replacement to location Cash Flow
                    rep_time += tc0[location_name][tech_name]['replacement']['years']
            
            
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
            
        if 'H tank' in studycase[location_name]:
            with open('results/LOC_'+name_studycase+'.pkl', 'rb') as f:
                    loc = pickle.load(f)
            final_tank_level = loc[location_name]['H tank'][-1]
            initial_tank_level = loc[location_name]['H tank'][0]
            tank_difference_level = initial_tank_level-final_tank_level
            if final_tank_level >= initial_tank_level:
                H2_cost = economic_data['hydrogen']['sale'] #€/kg
            else:
                H2_cost = economic_data['hydrogen']['purchase'] #€/kg
                
            results[location_name]['CF_studycase']['Initial/Final Tank level'][:] = -tank_difference_level*H2_cost    
            results[location_name]['CF_studycase']['Tot'] += results[location_name]['CF_studycase']['Initial/Final Tank level']  
            results[location_name]['CF_refcase']['Tot'] += results[location_name]['CF_refcase']['Initial/Final Tank level']
            results[location_name]['CF']['Initial/Final Tank level'] += results[location_name]['CF_studycase']['Initial/Final Tank level'] -results[location_name]['CF_refcase']['Initial/Final Tank level'] 
                 
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
        
    
        # calculate I0
        results[location_name]['I0']['Tot'] = 0
        for tech_name in results[location_name]['I0']:
            if tech_name != 'Tot':
                results[location_name]['I0']['Tot'] += results[location_name]['I0'][tech_name]
            
        # calculate NPV
        results[location_name]['NPV'] = np.zeros(economic_data['investment years']+1) # array initialise Net Present Value
        results[location_name]['NPV'][0] = -results[location_name]['I0']['Tot'] # NPV at time 0 is - the initial investment
        i = economic_data['interest rate'] # interest rate [%]
        
        PBT = -1
        for y in range(1,economic_data['investment years']+1): # for each year
            aux_var = results[location_name]['NPV'].sum(where=results[location_name]['NPV']>0)
            results[location_name]['NPV'][y] = results[location_name]['NPV'][y-1] + results[location_name]['CF']['Tot'][y-1]/(1+i)**y # calculate NPV
            if aux_var == 0 and results[location_name]['NPV'][y-1] < 0 and results[location_name]['NPV'][y] >= 0:
                PBT = y-1+(-results[location_name]['NPV'][y-1]/(-results[location_name]['NPV'][y-1]+results[location_name]['NPV'][y]))
        
        if PBT > 0:
            results[location_name]['PBT'] = PBT
            results[location_name]['PI'] = results[location_name]['NPV'][-1]/results[location_name]['I0']['Tot']
        else:
            results[location_name]['PBT'] = np.nan
            results[location_name]['PI'] = np.nan
        
    # save results in Results/economic_assesment.pkl
    with open(f"Results/economic_assessment_{name_economic}.pkl", 'wb') as f:  pickle.dump(results,f) 
        
            
            
def LCOH (
            structure,    
            name_studycase, 
            economic_data, 
            simulation_years, 
            path, 
            name_output, 
            revenues    = False, 
            refund      = False,
            plot        = True,
            print_      = True, 
            
            ):
    """
    Levelized Cost Of Hydrogen Calculation
    ----------
    name_studycase: str name of study case results file.pkl

    economic_data: dictionary
        'REC': dictionary REC economic parameters definition
            'collective self consumption incentives': [€/kWh]
            'incentives redistribution': 0-100 how the incentives are redistributed between prosumers, consumers and REC manger
        'carrier_name': dictionary (repeat for reach considered carrier: electricity, hydrogen, gas)
            'purchase': [€/kWh] (electricity and gas) or [€/kg] (hydrogen)
            'sales': [€/kWh] (electricity and gas) or [€/kg] (hydrogen)
        'interest rate': 0-1 [rate/year]
        'inflation rate': -1-1 [rate/year] cost evolution of each carrier
        'decommissioning': 0-1 [-] system dismantling as a % of initial construction costs (Total CAPEX)
        'investment year': time horizon for which to evaluate the economic analysis (must be a multiple of simulation_year in general.json)
     
    simulation_years: int number of years considered for the energy balances       
    
    path: str path of the input data folder 
        
    name_output: str name of the file where to save results of interest
    
    revenues: list/tuple of str/bool defining if generated revenues from excess energy streams have to be included in the calculation. In case the carrier name(s) is/are specified. Default = False
    
    refund: boolean value defining refund cash flows have to be included in the calculation. Default = False
                        
    output:  [€/kgH2] float value of LCOH for the considered configuration

    """       
    years_factor = int(economic_data['investment years'] / simulation_years) # this factor is useful to match the length of the energy simulation with the length of the economic investment
    
    # open energy balances of studycase
    with open('Results/balances_'+name_studycase+'.pkl', 'rb') as f:        balances = pickle.load(f)

    # open cost of componenets of studycase
    with open('Results/tech_cost_'+name_studycase+'.pkl', 'rb') as f:       tc = pickle.load(f)   # ricontrollare tc nel caso dell'analisi di sensitività, probabilmente non lo carica bene. 

    # Check for hydrogen carrier to be included in location balances
    if all(len(balances[location_name]['hydrogen']) != 0 for location_name in tc):    # if hydrogen dictionary has values for at least one location
        pass
    else:
        print("Hydrogen carrier not included in the considered case study - LCOH calculation not available")
        return
    
    results = {}    # dictionary initialising global lcoh results of each location
    lcoh    = {}    # dictionary initialising specific lcoh results of each location
    
    for location_name in tc:           # for each location
        
        results[location_name] = {}           # dictionary initialise economic results
        lcoh[location_name] = {}           # dictionary initialise economic results
        
        # initialise cash flow:     
        results[location_name]['CF'] = {  'OeM': np.zeros(economic_data['investment years']),
                                          'Purchase': {},
                                          'Sale': {},
                                          'Refund': np.zeros(economic_data['investment years']),
                                          'Tot': np.zeros(economic_data['investment years'])} 
    
        results[location_name]['I0'] = {} # initialise initial investment
        
        lcoh[location_name]['Capex'] = {}
        lcoh[location_name]['Opex'] = {}
        
        for tech_name in tc[location_name]:              # considering each technology in the location

            results[location_name]['I0'][tech_name] = tc[location_name][tech_name]['total cost'] # I0   
            results[location_name]['CF']['OeM'][:] += - tc[location_name][tech_name]['OeM'] # OeM
            
            lcoh[location_name]['Capex'][tech_name] = np.concatenate(([tc[location_name][tech_name]['total cost']], np.full(economic_data['investment years'], 0))) # Capex array
            lcoh[location_name]['Opex'][tech_name]  = np.repeat(tc[location_name][tech_name]['OeM'],economic_data['investment years']) # Opex array  
            
            # replacements 
            if tc[location_name][tech_name]['replacement']['years'] == "ageing": # if replacement year is calculated according to ageing
                with open('Results/ageing_'+name_studycase+'.pkl', 'rb') as f:
                    age = pickle.load(f)     
                    age = age[location_name][tech_name][0]
                    for a in age:
                        rep_time = int(a/8760)
                        results[location_name]['CF']['OeM'][rep_time] += - results[location_name]['I0'][tech_name] * tc[location_name][tech_name]['replacement']['rate']/100 # subtract technology replacement to location Cash Flow
                        lcoh[location_name]['Capex'][tech_name] += lcoh[location_name]['Capex'][tech_name] * tc[location_name][tech_name]['replacement']['rate']/100 # subtract technology replacement to location Cash Flow
            else: # if replacement time is given
                rep_time = tc[location_name][tech_name]['replacement']['years']
                while rep_time < economic_data['investment years']: # if tech_name replacement happens before the end of the simulation
                    results[location_name]['CF']['OeM'][rep_time-1] += - np.float64(results[location_name]['I0'][tech_name])*np.float(tc[location_name][tech_name]['replacement']['rate'])/np.float(100) # subtract technology replacement to location Cash Flow. -1 subtracted to index because the array is shifted of 1 position during LCOH calculation -> CF[1:]  = results[location_name]['CF']['Tot'].copy()    # [€] Cash Flows
                    lcoh[location_name]['Capex'][tech_name][rep_time] += np.float64(lcoh[location_name]['Capex'][tech_name][0])*np.float(tc[location_name][tech_name]['replacement']['rate'])/np.float64(100) # subtract technology replacement to location Cash Flow. np.float type added to avoid overflow problems 
                    rep_time += tc[location_name][tech_name]['replacement']['years']
            
            if refund: # if refund has to be included in LCOH calculation
                    
                if tc[location_name][tech_name]['refund']['years'] == 0 and tc[location_name][tech_name]['refund']['rate'] == 0:  # no refunds considered
                    pass
                elif tc[location_name][tech_name]['refund']['years'] == 0 and tc[location_name][tech_name]['refund']['rate'] != 0:
                    results[location_name]['I0'][tech_name] = results[location_name]['I0'][tech_name]*((100-tc[location_name][tech_name]['refund']['rate'])/100)
                    lcoh[location_name]['Capex'][tech_name] = lcoh[location_name]['Capex'][tech_name]*((100-tc[location_name][tech_name]['refund']['rate'])/100)
                else:
                    yearly_refund = results[location_name]['I0'][tech_name]*tc[location_name][tech_name]['refund']['rate']/100 / tc[location_name][tech_name]['refund']['years'] # yearly refund [€]
                    refunds = np.zeros(economic_data['investment years']) # array initialise
                    refunds[:min(economic_data['investment years'],tc[location_name][tech_name]['refund']['years'])] = yearly_refund # array repeat yearly refund 
                    results[location_name]['CF']['Refund'] += refunds # add refund to Cash Flow
                    lcoh[location_name]['Opex'][tech_name] -= refunds # add refund to Opex expenditures
                
                results[location_name]['CF']['Tot'] += - results[location_name]['CF']['Refund']
                
        # energy sold and purchased in study case 
        for carrier in balances[location_name]:                           # for each carrier (electricity, hydrogen, gas, heat)
            
            if 'grid' in balances[location_name][carrier]:  
                
                if type(economic_data[carrier]['sale']) == str: # if the price series is given
                    sale_serie = np.tile(pd.read_csv(path+'/energy_price/'+economic_data[carrier]['sale'])['0'].to_numpy(),int(simulation_years))  
                    sold = balances[location_name][carrier]['grid'] * sale_serie
                else: # if the price is always the same 
                    sold = balances[location_name][carrier]['grid']*economic_data[carrier]['sale'] 
               
                sold = np.tile(sold,years_factor)
                sold = np.reshape(sold,(-1,8760))
                results[location_name]['CF']['Sale'][carrier]   = - sold.sum(axis=1,where=sold<0)
                lcoh[location_name]['Opex'][carrier +' sold']   =   sold.sum(axis=1,where=sold<0) 

                if type(economic_data[carrier]['purchase']) == str: # if the price series is given
                    purchase_serie = np.tile(pd.read_csv(path+'/energy_price/'+economic_data[carrier]['purchase'])['0'].to_numpy(),int(simulation_years))  
                    purchase = balances[location_name][carrier]['grid'] * purchase_serie
                else: # if the price is always the same 
                    purchase = balances[location_name][carrier]['grid']*economic_data[carrier]['purchase']
               
                purchase = np.tile(purchase,years_factor)
                purchase = np.reshape(purchase,(-1,8760))
                results[location_name]['CF']['Purchase'][carrier]   = - purchase.sum(axis=1,where=purchase>0)
                lcoh[location_name]['Opex'][carrier +' purchased']  =   purchase.sum(axis=1,where=purchase>0) 

                
            
            if carrier == 'electricity':
                if 'wind electricity' in balances[location_name]['electricity'] and structure[location_name]['wind']['owned'] == False:   
                
                    # if type(economic_data[carrier]['wind electricity']['sale']) == str: # if the price series is given
                    #     sale_serie = np.tile(pd.read_csv(path+'/energy_price/'+economic_data[carrier]['sale'])['0'].to_numpy(),int(simulation_years))  
                    #     sold = balances[location_name][carrier]['wind electricity'] * sale_serie
                    # else: # if the price is always the same 
                    #     sold = balances[location_name][carrier]['wind electricity']*economic_data['wind electricity']['sale'] 
                   
                    # sold = np.tile(sold,years_factor)
                    # sold = np.reshape(sold,(-1,8760))
                    # results[location_name]['CF']['Sale'][carrier] = - sold.sum(axis=1,where=sold<0)
                    
                    if type(economic_data['wind electricity']['purchase']) == str: # if the price series is given
                        purchase_serie = np.tile(pd.read_csv(path+'/energy_price/'+economic_data[carrier]['purchase'])['0'].to_numpy(),int(simulation_years))  
                        purchase = balances[location_name][carrier]['wind electricity'] * purchase_serie
                    else: # if the price is always the same 
                        purchase = balances[location_name][carrier]['wind electricity']*economic_data['wind electricity']['purchase']
                   
                    purchase = np.tile(purchase,years_factor)
                    purchase = np.reshape(purchase,(-1,8760))
                    results[location_name]['CF']['Purchase'][carrier] += - purchase.sum(axis=1,where=purchase>0)
                    lcoh[location_name]['Opex']['wind electricity purchased']  =  purchase.sum(axis=1,where=purchase>0)
                
                if 'pv electricity' in balances[location_name]['electricity'] and structure[location_name]['PV']['owned'] == False:  
                
                #     if type(economic_data[carrier]['sale']) == str: # if the price series is given
                #         sale_serie = np.tile(pd.read_csv(path+'/energy_price/'+economic_data[carrier]['sale'])['0'].to_numpy(),int(simulation_years))  
                #         sold = balances[location_name][carrier]['pv electricity'] * sale_serie
                #     else: # if the price is always the same 
                #         sold = balances[location_name][carrier]['pv electricity']*economic_data[carrier]['sale'] 
                   
                #     sold = np.tile(sold,years_factor)
                #     sold = np.reshape(sold,(-1,8760))
                #     results[location_name]['CF']['Sale'][carrier] = - sold.sum(axis=1,where=sold<0)
                    
                    if type(economic_data['pv electricity']['purchase']) == str: # if the price series is given
                        purchase_serie = np.tile(pd.read_csv(path+'/energy_price/'+economic_data[carrier]['purchase'])['0'].to_numpy(),int(simulation_years))  
                        purchase = balances[location_name][carrier]['pv electricity'] * purchase_serie
                    else: # if the price is always the same 
                        purchase = balances[location_name][carrier]['pv electricity']*economic_data['pv electricity']['purchase']
                   
                    purchase = np.tile(purchase,years_factor)
                    purchase = np.reshape(purchase,(-1,8760))
                    results[location_name]['CF']['Purchase'][carrier] += - purchase.sum(axis=1,where=purchase>0)
                    lcoh[location_name]['Opex']['pv electricity purchased']  =  purchase.sum(axis=1,where=purchase>0)
        
        # CF update considering inflation on each carrier
        for carrier in economic_data['inflation rate']:
            f = economic_data['inflation rate'][carrier]
            
            if carrier in results[location_name]['CF']['Purchase']:
                for y in range(economic_data['investment years']):
                    results[location_name]['CF']['Purchase'][carrier][y] = results[location_name]['CF']['Purchase'][carrier][y]*(1+f)**y
                    
            if carrier in results[location_name]['CF']['Sale']:
                for y in range(economic_data['investment years']):
                    results[location_name]['CF']['Sale'][carrier][y] = results[location_name]['CF']['Sale'][carrier][y]*(1+f)**y
        
            for k in lcoh[location_name]['Opex']: 
                if carrier in k:
                    for y in range(economic_data['investment years']):
                        lcoh[location_name]['Opex'][k][y] = lcoh[location_name]['Opex'][k][y]*(1+f)**y
                        
        # Building Cash Flow final array while changing sign to revenews and expenditures as needed in LCOH formula
        for carrier in results[location_name]['CF']['Purchase']:
            results[location_name]['CF']['Tot'] -= results[location_name]['CF']['Purchase'][carrier]
            
        if revenues:        # if revenues have to be included in LCOH calculation
            for carrier in results[location_name]['CF']['Sale']:
                if carrier in revenues:
                    results[location_name]['CF']['Tot'] -= results[location_name]['CF']['Sale'][carrier]
            # lcoh dictionary
            keys_to_remove = []
            for carrier in list(lcoh[location_name]['Opex'].keys()):
                if 'sold' in carrier and not any(word in revenues for word in carrier.split()):  # deleting revenues array not included in the analysis
                    keys_to_remove.append(carrier)
            for key in keys_to_remove:
                lcoh[location_name]['Opex'].pop(key)
        else:
            keys_to_remove = []
            for carrier in lcoh[location_name]['Opex']:
                if 'sold' in carrier:
                    keys_to_remove.append(carrier)
            for key in keys_to_remove:
                    del lcoh[location_name]['Opex'][key]   
                                
        results[location_name]['CF']['Tot'] +=  - results[location_name]['CF']['OeM']       
                                            
        # calculate I0
        results[location_name]['I0']['Tot'] = 0
        for tech_name in results[location_name]['I0']:
            if tech_name != 'Tot':
                results[location_name]['I0']['Tot'] += results[location_name]['I0'][tech_name]   
                
        for key in lcoh[location_name]['Opex']:
            lcoh[location_name]['Opex'][key] = np.insert(lcoh[location_name]['Opex'][key], 0, 0) # shifting Opex value to year 1 while only Capex is considered in year 0.
            
        lcoh[location_name]['Capex']['Total'] = sum(lcoh[location_name]['Capex'].values())  # Adding Capex total values to dictionary
        lcoh[location_name]['Opex']['Total'] = sum(lcoh[location_name]['Opex'].values())    # Adding Opex total values to dictionary
       
        # LCOH calculation
        
        # Hydrogen produced each year via electrolysis
        produced_hydrogen = [0] + [sum(balances[location_name]['hydrogen']['electrolyzer'])]*economic_data['investment years']  # [kg/y] - No H2 produced in period 0
        r       = economic_data['interest rate']                # [%] interest rate 
        I0      = results[location_name]['I0']['Tot']           # [€] Initial investment at time = 0
        CF      = np.zeros(economic_data['investment years'] +1)                     # Creating an empty array of year_factor + 1 dimension for Cash Flows in order to insert only I0 as first element
        CF[1:]  = results[location_name]['CF']['Tot'].copy()    # [€] Cash Flows
        CF[0]   = I0
        if economic_data['decommissioning'] > 0:
            CF = np.append(CF, I0*economic_data['decommissioning'])
            produced_hydrogen.append(0)
            lcoh[location_name]['Capex']['Total'] = np.append(lcoh[location_name]['Capex']['Total'], I0*economic_data['decommissioning'])
            
        num=[]  # numerator
        den=[]  # denominator
        
        for i in range(len(CF)):

            num.append((CF[i])*(1/(1+r)**i))
            den.append(produced_hydrogen[i]*(1/(1+r)**i))
        
        LCOH = round(sum(num)/sum(den),3)
        
        results[location_name]['LCOH'] = {'Value [€/kgH2]'          : LCOH,
                                          'Discounted Expenditures' : num,
                                          'Discounted Production'   : den}
                
        for key in lcoh[location_name]['Opex']:
            for i in range(len(lcoh[location_name]['Opex'][key])):
                lcoh[location_name]['Opex'][key][i] = lcoh[location_name]['Opex'][key][i]*(1/(1+r)**i)
        
        for key in lcoh[location_name]['Capex']:
            for i in range(len(lcoh[location_name]['Capex'][key])):
                lcoh[location_name]['Capex'][key][i] = lcoh[location_name]['Capex'][key][i]*(1/(1+r)**i)
        
        # specific costs associated with each technology. Data handling.
        capex   = {}
        opex    = {}
        
        for key, value in lcoh[location_name]['Capex'].items():
            capex[key] = np.sum(value)
        for key, value in lcoh[location_name]['Opex'].items():
            opex[key] = np.sum(value)
            
        hydrogen_prod = sum(den)   # actualized hydrogen production

        # Building DataFrame for data visualization      
        df = pd.DataFrame([capex, opex], index=['Capex', 'Opex'])
        df = df.fillna(0)
        df = df.reindex(columns=[col for col in df.columns if col != 'Total'] + ['Total'])  # index 'Total' must be the last of the column
         
        
        df1 = df.copy()
        df1 = df1/hydrogen_prod
        df1 = df1.round(3)
        df1 = df1.drop(columns=df1.columns[(df1.loc['Capex'] == 0) & (df1.loc['Opex'] == 0)])    # removing elements not contributing to final LCOH value
        
        if plot == True:
          
            if revenues:
                for key in df.loc['Opex'].index.tolist():
                    if key.replace(' sold', '') in revenues and 'sold' in key:
                        df.loc['Opex']['Total'] -= df.loc['Opex'][key]
            
            df = df.drop(columns=df.columns[(df.loc['Capex'] == 0) & (df.loc['Opex'] == 0)])    # removing elements not contributing to final LCOH value
            df = df.drop(columns=[col for col in df.columns if 'sold' in col])
            df = df/hydrogen_prod
            df = df.round(3)
      
            colors = ['#0e4d92', '#2380b2', '#5da5c5', '#9cc2a5', '#c8e0a1', '#ebd279', '#e18b4f', '#ff0000', '#00ff00', '#0000ff']
            colors1 = ['#7fcdbb', '#edf8b1']       
     
            labels_capex    = [col for col in df.columns if col != 'Total' and df.loc['Capex', col] != 0]
            labels_opex     = [col for col in df.columns if col != 'Total' and df.loc['Opex', col] != 0]
            labels_outer    = labels_capex + labels_opex
    
            values_capex = df.loc['Capex', labels_capex].tolist()
            values_opex = df.loc['Opex', labels_opex].tolist()
            values_outer = values_capex + values_opex
    
            repeated_labels = [label for label in labels_outer if labels_outer.count(label) > 1]
    
            label_colors = dict(zip(repeated_labels, colors[:len(repeated_labels)]))
            label_colors = {label: colors.pop(0) if label not in label_colors else label_colors[label] for label in labels_outer}
            
            labels_inner = ['Capex','Opex']
            values_inner = [df.loc['Capex','Total'],df.loc['Opex','Total']]
                    
            # Data to plot
            explode = [0.] * len(labels_outer)
            explode1 = [0.] * len(labels_inner)
                                      
            #Plot
            plt.figure(dpi=300)
            plt.pie(values_outer, 
                    labels=labels_outer, 
                    startangle=180,
                    frame=True,
                    radius=5,
                    colors=[label_colors[label] for label in labels_outer],
                    explode = explode, 
                    wedgeprops={'linewidth': 1, 'edgecolor': 'black'},
                    autopct='%1.1f%%',
                    pctdistance=1.15,
                    labeldistance=1.32,
                    textprops={'fontsize': 8})
    
            prop = fm.FontProperties(weight='bold')
            plt.pie(values_inner,
                    labels=labels_inner,
                    startangle=180, 
                    radius=3.5,
                    colors = colors1,
                    explode=explode1,
                    autopct='%1.1f%%',
                    labeldistance=0.04,
                    pctdistance=0.82,
                    textprops={'fontproperties': prop, 'fontsize':9} )
            
            #Draw circle
            centre_circle = plt.Circle((0,0),1.75,color='black', fc='white',linewidth=0)
            fig = plt.gcf()
            fig.gca().add_artist(centre_circle)
            
            column_labels = [label + " [€/kgH$_\mathregular{2}$]" for label in df.index]
            table = plt.table(cellText=df1.T.values,  
                      colLabels=column_labels,    
                      rowLabels=df1.columns,  
                      cellLoc='center',
                      loc='center',
                      bbox=[0.325, -1.25, 0.75, 0.8])  # table dimensions and position
            table.auto_set_font_size(True)
            table.scale(1, 1.)
            
            plt.suptitle(f'LCOH = {LCOH} €/kgH$_2$', y=0.05, fontsize=12, fontweight='bold', ha='center')
            plt.subplots_adjust(bottom=0.18)    
    
            plt.axis('equal')
            # plt.tight_layout()
            plt.show()        
        
        # save results in Results/economic_assesment.pkl
        with open(f"Results/LCOH_assessment_{name_output}.pkl", 'wb') as f:  pickle.dump(results,f)
        
        if print_ == True:
            print(f"The LCOH calculated for the considered scenario results in {LCOH} €/kgH2")
        
        return(LCOH)                                      


         
    
            
            
    
            