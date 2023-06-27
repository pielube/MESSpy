"""
MESSpy - postprocessing
"""
import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import plotly.io as pio
import plotly.graph_objects as go

        
def total_balances(simulation_name,loc,var=None):
    """
    Total balances figures
    
    simulationa_name : str 
    loc : str 
    var : str -> possibility of specifying a single energy carrier
    
    """
    
    with open('results/balances_'+simulation_name+'.pkl', 'rb') as f: balances = pickle.load(f)
    
    ###### load analysis
    
    ###### total energy balances
    
    carriers = ['electricity','heating water','gas','hydrogen','HP hydrogen']
    units    = {'electricity': 'kWh', 'hydrogen': 'kg', 'HP hydrogen': 'kg', 'gas':'kWh', 'heating water':'kWh'}
    
    if var: 
        
        carriers = [var]
        units    = {var : units[var]}
        balance  = 0                      # initializing the variable to visualize the balance at the end of the simulation period

        if var == 'hydrogen' and 'mechanical compressor' in balances[loc][var]:
            balances[loc][var].pop('mechanical compressor')  # dict.values() to be removed as they alter the total hydrogen balance and ar enot supposed to do so. 
            
    for carrier in carriers:
        print('\nTotal '+carrier+' balances '+loc+':\n') 
        for b in balances[loc][carrier]:
            
            positiv=round(balances[loc][carrier][b][balances[loc][carrier][b]>0].sum(),1)
            negativ=round(balances[loc][carrier][b][balances[loc][carrier][b]<0].sum(),1)
            
            if positiv != 0:
                print(b+' '+str(round(positiv,1))+' '+units[carrier])
                balance += round(positiv,1)
                
            if negativ != 0:
                print(b+' '+str(round(negativ,1))+' '+units[carrier])
                balance += round(negativ,1)
                
                
    print(f"Total {carrier} balance {balance} {units[carrier]}")
    print('\n')

def hydrogen_production(simulation_name,loc,var=None):
    """
    Hydrogen production figures and plots
    
    simulationa_name : str 
    loc : str 
    var : str -> possibility of specifying **kargs
    
    """    
    with open('results/balances_'+simulation_name+'.pkl', 'rb') as f: balances = pickle.load(f)
    
    simulation_hours    = len(balances[loc]['hydrogen']['electrolyzer'])
    constantflow        = sum(balances[loc]['hydrogen']['electrolyzer'])/simulation_hours  # [kg/h] constant flow rate deliverable by the system if operaten in supply-led mode
    demand_constant     = [constantflow * (i + 1) for i in range(simulation_hours)]
    # demand_variable     = -np.cumsum(balances[loc]['hydrogen']['demand'])
    production          = np.cumsum(balances[loc]['hydrogen']['electrolyzer'])
    
    fig, ax = plt.subplots(dpi=1000)
    ax.plot(demand_constant, label = 'constant demand')
    # ax.plot(demand_variable, label = 'variable demand')
    ax.plot(production, label = 'production')
    ax.grid(alpha = 0.3, zorder = 0)
    ax.set_ylabel('Hydrogen [kg]')
    ax.set_xlabel('Time [h]')
    xticks = list(np.linspace(0, simulation_hours - 1, 13).astype(int))
    # xticklabels = [str(value) for value in xticks]
    xticklabels = ['         Jan','         Feb','          Mar','         Apr','         May','          Jun','        Jul','          Aug','           Sep','          Oct','          Nov','           Dec','']
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    ax.legend()
    plt.title('Cumulative H$_\mathregular{2}$ production and demand')
    plt.show()
    

def renewables(simulation_name,simulation_years,loc,first_day,last_day,plot=False):  # to be generalized for all locations
    
    with open('results/balances_'+simulation_name+'.pkl', 'rb') as f: balances = pickle.load(f)
    ###### load analysis
    
    ###### total energy balances
    
    if 'wind' not in balances[loc]['electricity'] and 'PV' not in balances[loc]['electricity']:    # if no renewables are included in the location        
        print("No renewable energy sources are present in the considered case study\nAutoconsumption data calculation not available\n")
        return
    else:
        
        windsc  = np.zeros(simulation_years*8760)  # initializing the array to store wind energy self consumption values
        pvsc    = np.zeros(simulation_years*8760)  # initializing the array to store solar energy self consumption values
        demand  = np.zeros(simulation_years*8760)  # initializing the array to store the overall demand driven also by hydrogen technolgies
        x = np.arange(first_day*24,(last_day+1)*24)   
        
        for tech_name in balances[loc]['electricity']:
        # for tech_name in rec.locations[loc].system: # in order to generalize 
            for i in range(len(balances[loc]['electricity'][tech_name])):
                if tech_name != 'grid' and balances[loc]['electricity'][tech_name][i] < 0: 
                    demand[i] += balances[loc]['electricity'][tech_name][i]
        balances[loc]['electricity']['demand'] = demand
                    
        fromgrid = np.maximum(balances[loc]['electricity']['grid'],0)
        intogrid = np.minimum(balances[loc]['electricity']['grid'],0)
        
        if 'wind' in balances[loc]['electricity'] and not 'PV' in  balances[loc]['electricity']:    # wind ha priority over solar power
            
            windsc = np.minimum(balances[loc]['electricity']['wind'],-demand)

            balances[loc]['electricity']['wind electricity'] = {}
            balances[loc]['electricity']['wind electricity'] = windsc

            k = 1
            if plot == True:
                fig, ax = plt.subplots(dpi=1000)
                ax.plot(x,-demand[first_day*24:last_day*24+24], label = 'Demand', linewidth=k, zorder = 2)
                ax.plot(x,windsc[first_day*24:last_day*24+24],  alpha = 1, zorder = 1, label = 'P$_{Wind SC}$', linewidth=k)
                ax.plot(x,fromgrid[first_day*24:last_day*24+24], label = 'From grid', linewidth=0.8)
                ax.grid(alpha = 0.3, zorder = 0)
                if (last_day-first_day) <= 20:
                    ax.set_xticks(list(range(first_day*24, (last_day+1)*24+1,24)))
                    ax.set_xticklabels([str(x) for x in list(range(first_day*24, (last_day+1)*24+1,24))], rotation=45)
                ax.set_ylabel('Energy [kWh]')
                ax.set_xlabel('Time [h]')
                ax.legend()
                plt.tight_layout()
                plt.show()
            
        elif 'wind' in balances[loc]['electricity'] and 'PV' in  balances[loc]['electricity']:
            windsc      = np.minimum(balances[loc]['electricity']['wind'],-demand)
            pvsc        = np.minimum(-np.minimum(windsc+demand, 0),balances[loc]['electricity']['PV'])                
            
            balances[loc]['electricity']['wind electricity'] = {}
            balances[loc]['electricity']['pv electricity'] = {}
            
            balances[loc]['electricity']['wind electricity']    = windsc
            balances[loc]['electricity']['pv electricity']      = pvsc
            
            k = 1
            if plot == True:
                fig, ax = plt.subplots(dpi=1000)
                ax.plot(x,-demand[first_day*24:last_day*24+24], label = 'Demand', linewidth=k, zorder = 2)
                ax.plot(x,windsc[first_day*24:last_day*24+24],  alpha = 1, zorder = 1, label = 'P$_{Wind SC}$', linewidth=k)
                ax.plot(x,fromgrid[first_day*24:last_day*24+24], label = 'From grid', linewidth=0.8)
                ax.plot(x,pvsc[first_day*24:last_day*24+24], label = 'P$_{PV SC}$', linewidth=k, zorder = 3)
                # ax.plot(pvsc[first_day*24:last_day*24+24]+windsc[first_day*24:last_day*24+24], label = 'P$_{Wind+PV SC}$', linewidth=k, zorder = 4)
                ax.grid(alpha = 0.3, zorder = 0)
                if (last_day-first_day) <= 20:
                    ax.set_xticks(list(range(first_day*24, (last_day+1)*24+1,24)))
                    ax.set_xticklabels([str(x) for x in list(range(first_day*24, (last_day+1)*24+1,24))], rotation=45)
                ax.set_xlabel('Time [h]')
                ax.set_ylabel('Energy [kWh]')
                ax.legend()
                plt.show()
        
        if plot == True:
        
            windcumulative = np.cumsum(windsc)
            pvcumulative = np.cumsum(pvsc)
            gridcumulative = np.cumsum(fromgrid)
            demandcumulative = np.cumsum(abs(demand))
            
            fig, ax = plt.subplots(dpi=1000)
            ax.plot(demandcumulative, label = 'Demand', linewidth=k, zorder = 2)
            ax.plot(windcumulative,  alpha = 1, zorder = 1, label = 'P$_{Wind SC}$', linewidth=k)
            ax.plot(gridcumulative, label = 'From grid', linewidth=0.8)
            ax.plot(pvcumulative, label = 'P$_{PV SC}$', linewidth=k, zorder = 3)
            ax.grid(alpha = 0.3, zorder = 0)
            xticks = list(np.linspace(0, len(demandcumulative) - 1, 13).astype(int))
            # xticklabels = [str(value) for value in xticks]
            xticklabels = ['          Jan','         Feb','          Mar','         Apr','         May','          Jun','         Jul','          Aug','           Sep','          Oct','          Nov','           Dec','']
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticklabels)
            ax.set_ylabel('Energy [kWh]')
            ax.set_xlabel('Time [h]')
            ax.legend()
            plt.show()
        
    
        # saving post process results
        
        with open('results/balances_'+simulation_name+".pkl", 'wb') as f: pickle.dump(balances, f)
    
    
    
def ghg_emissions(simulation_name,loc,emission_intensity,ref_year,print_=False):
    """
    Total balances figures
    
    simulationa_name : str 
    loc : str 
    emission_intensity : dict 
    ref_year: str -> necessary to specify the reference year for grid intensity
    
    """
    with open('results/balances_'+simulation_name+'.pkl', 'rb') as f: balances = pickle.load(f)
    
    carrier = 'electricity'
    
    positiv=round(balances[loc][carrier]['grid'][balances[loc][carrier]['grid']>0].sum(),1)   # amount of energy withdrown from the grid during operations. 
    emission_factor= emission_intensity['reference year'][ref_year]                   # [gCO2/kWh] grid intensity
    co2_tot = positiv*emission_factor/1000                          # [kgCO2] total amount of carbon dioxide due to grid electricity utilization
    produced_hyd = sum(balances[loc]['hydrogen']['electrolyzer'])   # [kgH2] total amount of produced hydrogen via in situ electorlysis
    
    h2_ghg = round(co2_tot/produced_hyd,2)                          # [kgCO2/kgH2] GHG intensity of the produced hydrogen
    
    if print_ == True:
        print(f"The H2 GHG intensity calculated for the considered scenario results in {h2_ghg} kgCO2/kgH2")
    
    return h2_ghg
    
    
    
def NPV_plot(study_case):
    ##### economic
    with open('results/economic_assessment_'+study_case+'.pkl', 'rb') as f: economic = pickle.load(f)
    
    plt.figure(dpi=1000)
    
    for loc in economic:    
        y = economic[loc]['NPV']
        x = np.linspace(0,len(y)-1,len(y))
        plt.plot(x,y,label=loc)
        
    plt.legend()
    plt.title("Investments")
    plt.grid()
    plt.ylabel('Net Present Value [€]')
    plt.xlabel('Time [years]')
    plt.xlim(0,len(y)-1)
    #plt.ylim(-15000,15000)
    plt.show()
    
    
def REC_electricity_balance(simulation_name):
    
    with open('results/balances_'+simulation_name+'.pkl', 'rb') as f:
        balances = pickle.load(f)
        
    df = pd.DataFrame(0.00,columns=["Value [kWh]","Value / production [%]","Value / demand [%]"],
                      index=["SC","CSC","Into grid","From grid","Losses","Production","Demand"])
    
    df.loc['CSC']['Value [kWh]'] = sum(balances['REC']['electricity']['collective self consumption'])
    df.loc['Into grid']['Value [kWh]'] = -sum(balances['REC']['electricity']['into grid'])  -df.loc['CSC']['Value [kWh]']
    df.loc['From grid']['Value [kWh]'] = sum(balances['REC']['electricity']['from grid']) -df.loc['CSC']['Value [kWh]']
    
    for loc in balances:
        if loc != 'REC':   
            balance = balances[loc]['electricity']
            df.loc['Demand']['Value [kWh]'] += -sum(balance['demand'])
            if 'PV' in balance:
                df.loc['Production']['Value [kWh]'] += sum(balance['PV'])
            if 'wind' in balance:
                df.loc['Production']['Value [kWh]'] += sum(balance['wind'])
                
    df.loc['SC']['Value [kWh]'] = df.loc['Demand']['Value [kWh]'] - df.loc['From grid']['Value [kWh]'] - df.loc['CSC']['Value [kWh]']
    df.loc['Losses']['Value [kWh]'] = df.loc['Production']['Value [kWh]'] - df.loc['SC']['Value [kWh]'] - df.loc['Into grid']['Value [kWh]'] - df.loc['CSC']['Value [kWh]']
    
    for b in df.index:
        df.loc[b]['Value / production [%]'] = df.loc[b]['Value [kWh]'] / df.loc['Production']['Value [kWh]'] *100
        df.loc[b]['Value / demand [%]'] = df.loc[b]['Value [kWh]'] / df.loc['Demand']['Value [kWh]'] *100
        
    print('\nREC electricity balances: '+simulation_name+'\n') 
    print(df.round(decimals=2))
    print('\n')
    return(df.round(decimals=2))

    
def LOC_plot(simulation_name):
      
    with open('results/LOC_'+simulation_name+'.pkl', 'rb') as f:
        LOC = pickle.load(f)
           
    unit = {'H tank': '[kg]', 'HPH tank': '[kg]', 'battery': '[kWh]', 'inertial TES': '[°C]'}
    
    for location_name in LOC:
        for tech in LOC[location_name]:
            
            plt.figure(dpi=600)                               
            y = LOC[location_name][tech]
            x = np.linspace(0,len(y)-1,len(y))        
            plt.plot(x,y,label=location_name)
            plt.grid()
            plt.ylabel('LOC '+unit[tech])
            plt.xlabel('Time [hours]')
            xticks = list(np.linspace(0, len(x) - 1, 13).astype(int))
            # xticklabels = [str(value) for value in xticks]
            xticklabels = ['          Jan','         Feb','          Mar','         Apr','         May','          Jun','         Jul','          Aug','           Sep','          Oct','          Nov','           Dec','']
            plt.xticks(xticks,xticklabels,rotation=45)
            plt.title(location_name+' '+tech)
            plt.show()
            
            
def storage_control(simulation_name,e_cost=0.30,H_cost=0.05):
    with open('results/LOC_'+simulation_name+'.pkl', 'rb') as f:
        LOC = pickle.load(f)
        
    for location_name in LOC:
        if "battery" in LOC[location_name]:
            
            print(location_name)
            print('battery initial LOC = '+str(LOC[location_name]['battery'][0])+' kWh')
            print('battery final LOC = '+str(LOC[location_name]['battery'][-1])+' kWh')
            
            print('cost of electricity = '+str(round((LOC[location_name]['battery'][0]-LOC[location_name]['battery'][-1])*e_cost,2))+ ' €')
            print("\n")
            
        if "H tank" in LOC[location_name]:
            
            print(location_name)
            print('H tank initial LOC = '+str(round(LOC[location_name]['H tank'][0],2))+ ' kg')
            print('H tank final LOC = '+str(round(LOC[location_name]['H tank'][-1],2))+' kg')
            print('cost of hydrogen = '+str(round((LOC[location_name]['H tank'][0]-LOC[location_name]['H tank'][-1])*H_cost,2))+ ' €')
            print("\n")
            
    
def hourly_balances_electricity(simulation_name,location_name,first_day,last_day,carrier='electricity',width=0.9,collective=0):
    
        with open('Results/balances_'+simulation_name+'.pkl', 'rb') as f:
            balances = pickle.load(f)
            
        balances = balances[location_name][carrier]
        
        if 'demand' in balances:
            load = -balances['demand'][first_day*24:last_day*24+24]
        else:
            load = np.zeros(last_day*24+24 - first_day*24)
            
        if 'chp_gt' in balances:
            chp_gt = balances['chp_gt'][first_day*24:last_day*24+24]
            
        if 'chp' in balances:
            chp = balances['chp'][first_day*24:last_day*24+24]    
        
        if 'PV' in balances:
            pv = balances['PV'][first_day*24:last_day*24+24]
        
        if 'wind' in balances:
            wind = balances['wind'][first_day*24:last_day*24+24]
            
        if 'grid' in balances:                
            into_grid = np.zeros(24*(last_day-first_day+1))
            from_grid = np.zeros(24*(last_day-first_day+1)) 
            for i,e in enumerate(balances['grid'][first_day*24:last_day*24+24]):
                if e > 0:
                    from_grid[i] = e
                else:
                    into_grid[i] = e 
        
        if 'battery' in balances:
            charge_battery = np.zeros(24*(last_day-first_day+1))
            discharge_battery = np.zeros(24*(last_day-first_day+1))
            for i,e in enumerate(balances['battery'][first_day*24:last_day*24+24]):
                if e > 0:
                    discharge_battery[i] = e
                else:
                    charge_battery[i] = e 
        
        if 'electrolyzer' in balances:
            ele = balances['electrolyzer'][first_day*24:last_day*24+24]
        
        if 'mechanical compressor' in balances:
            compressor = balances['mechanical compressor'][first_day*24:last_day*24+24]
        
        if 'fuel cell' in balances:
            fc = balances['fuel cell'][first_day*24:last_day*24+24]

                     
        to_csc = np.zeros(24*(last_day-first_day+1))
        from_csc = np.zeros(24*(last_day-first_day+1))
        for i,e in enumerate(balances['collective self consumption'][first_day*24:last_day*24+24]):
            if e > 0:
                from_csc[i] = e
            else:
                to_csc[i] = e 
            
        x = np.arange(first_day*24,(last_day+1)*24)
        
        fig = plt.figure(dpi=1000)
        from mpl_toolkits.axisartist.axislines import SubplotZero
        ax = SubplotZero(fig, 1, 1, 1)
        fig.add_subplot(ax)
        ax.axis["xzero"].set_visible(True)
        ax.axis["xzero"].label.set_visible(False)
        ax.axis["xzero"].major_ticklabels.set_visible(False)

        for n in ["bottom","top", "right"]:
            ax.axis[n].set_visible(True)
        ax.grid(axis='y', alpha = 0.5, zorder = -4)
        
        
        if 'PV' in balances and 'wind' in balances:
            # windsc      = np.minimum(balances[loc]['electricity']['wind'],-demand)
            # pvsc        = np.minimum(-np.minimum(windsc+demand, 0),balances[loc]['electricity']['PV'])   
            ax.bar(x, np.minimum(load,wind), width,  label='Wind self consumption', color='yellowgreen')
            ax.bar(x, np.maximum(wind-load,np.zeros(len(pv))), width, bottom=load, label='Wind surplus', color='cornflowerblue')
            ax.bar(x, np.maximum(np.minimum(wind+pv,load)-wind,np.zeros(len(pv))), width, bottom= wind, label='PV self consumption', color='tab:green')
            ax.bar(x, np.maximum(wind-load+pv,np.zeros(len(pv)))-np.maximum(wind-load,np.zeros(len(pv))), width, bottom=load+np.maximum(wind-load,np.zeros(len(pv))),label='PV surplus', color='tab:blue')
            
            if not 'battery' in balances and not 'electrolyzer' in balances:
                ax.bar(x, from_grid-from_csc, width ,bottom=pv+from_csc, label='from grid', color='tomato')
                ax.bar(x, np.array(from_csc), width, bottom=pv,  color='gold')
            
                ax.bar(x, np.array(into_grid), width, label='into grid', color='orange')
                ax.bar(x, np.array(to_csc), width, label='collective self consumption',  color='gold')
                
            if 'battery' in balances:
                ax.bar(x, from_grid-from_csc, width, bottom=pv+discharge_battery+from_csc, label='from grid', color='tomato')
                ax.bar(x, np.array(from_csc), width, bottom=discharge_battery+pv,  color='gold')
                ax.bar(x, discharge_battery, width, bottom = pv, label='discharge battery',  color='purple')
                
                ax.bar(x, np.array(into_grid), width, bottom=charge_battery , label='into grid', color='orange')
                if collective == 0:
                    ax.bar(x, np.array(charge_battery), width,  label='charge battery',  color='violet')
                    ax.bar(x, np.array(to_csc), width, bottom=charge_battery, label='collective self consumption',  color='gold')
                else:
                    ax.bar(x, np.array(charge_battery), width, bottom=np.array(to_csc),  label='charge battery',  color='violet')
                    ax.bar(x, np.array(to_csc), width, label='collective self consumption',  color='gold')
                    
            if 'electrolyzer' in balances and 'fuel cell' in balances:
                ax.bar(x, from_grid-from_csc, width ,bottom=pv+fc+from_csc, label='from grid', color='tomato')
                ax.bar(x, np.array(from_csc), width, bottom=pv+fc,  color='y')
            
                ax.bar(x, np.array(into_grid), width,bottom=ele , label='into grid', color='orange')
                ax.bar(x, np.array(ele), width,  label='to electrolyzer',  color='chocolate')
                ax.bar(x, fc, width, bottom = pv, label='from fuel cell',  color='peru')
                ax.bar(x, np.array(to_csc), width, bottom=ele, label='collective self consumption',  color='gold') 
                
            if 'electrolyzer' in balances:
                ax.bar(x, from_grid-from_csc, width ,bottom=wind+pv+from_csc, label='from grid', color='tomato')
                ax.bar(x, np.array(from_csc), width, bottom=pv,  color='y')
                ax.bar(x, np.array(into_grid), width,bottom=ele , label='into grid', color='orange')
                ax.bar(x, np.array(ele), width,  label='to electrolyzer',  color='chocolate')
                ax.bar(x, np.array(to_csc), width, bottom=ele, label='collective self consumption',  color='gold')
            
            if 'mechanical compressor' in balances:
                ax.bar(x, np.array(compressor), width, bottom = ele, label='to compressor',  color='maroon')
                ax.bar(x, np.array(into_grid), width, bottom=ele+compressor , label='into grid', color='orange')
                
        elif 'PV' in balances and not ('chp_gt' and 'chp') in balances:
            ax.bar(x, np.minimum(load,pv), width,  label='PV self consumption', color='yellowgreen')
            ax.bar(x, np.maximum(pv-load,np.zeros(len(pv))), width, bottom=load, label='PV surplus', color='cornflowerblue')
            
            if not 'battery' in balances and not 'electrolyzer' in balances:
                ax.bar(x, from_grid-from_csc, width ,bottom=pv+from_csc, label='from grid', color='tomato')
                ax.bar(x, np.array(from_csc), width, bottom=pv,  color='gold')
            
                ax.bar(x, np.array(into_grid), width, label='into grid', color='orange')
                ax.bar(x, np.array(to_csc), width, label='collective self consumption',  color='gold')
                
            if 'battery' in balances:
                ax.bar(x, from_grid-from_csc, width, bottom=pv+discharge_battery+from_csc, label='from grid', color='tomato')
                ax.bar(x, np.array(from_csc), width, bottom=discharge_battery+pv,  color='gold')
                ax.bar(x, discharge_battery, width, bottom = pv, label='discharge battery',  color='purple')
                
                ax.bar(x, np.array(into_grid), width, bottom=charge_battery , label='into grid', color='orange')
                if collective == 0:
                    ax.bar(x, np.array(charge_battery), width,  label='charge battery',  color='violet')
                    ax.bar(x, np.array(to_csc), width, bottom=charge_battery, label='collective self consumption',  color='gold')
                else:
                    ax.bar(x, np.array(charge_battery), width, bottom=np.array(to_csc),  label='charge battery',  color='violet')
                    ax.bar(x, np.array(to_csc), width, label='collective self consumption',  color='gold')
                    
            if 'electrolyzer' in balances and 'fuel cell' in balances:
                ax.bar(x, from_grid-from_csc, width ,bottom=pv+fc+from_csc, label='from grid', color='tomato')
                ax.bar(x, np.array(from_csc), width, bottom=pv+fc,  color='y')
            
                ax.bar(x, np.array(into_grid), width,bottom=ele , label='into grid', color='orange')
                ax.bar(x, np.array(ele), width,  label='to electrolyzer',  color='chocolate')
                ax.bar(x, fc, width, bottom = pv, label='from fuel cell',  color='peru')
                ax.bar(x, np.array(to_csc), width, bottom=ele, label='collective self consumption',  color='gold')
            
            if 'mechanical compressor' in balances:
                ax.bar(x, np.array(compressor), width, bottom = ele, label='to compressor',  color='maroon')
                ax.bar(x, np.array(into_grid), width, bottom= ele+compressor , label='into grid', color='orange')
           
        elif 'wind' in balances and not ('chp_gt' and 'chp') in balances:
            ax.bar(x, np.minimum(load,wind), width,  label='wind self consumption', color='yellowgreen')
            ax.bar(x, np.maximum(wind-load,np.zeros(len(wind))), width, bottom=load, label='wind surplus', color='cornflowerblue')
            
            if not 'battery' in balances and not 'electrolyzer' in balances:
                ax.bar(x, from_grid-from_csc, width ,bottom=wind+from_csc, label='from grid', color='tomato')
                ax.bar(x, np.array(from_csc), width, bottom=wind,  color='gold')
            
                ax.bar(x, np.array(into_grid), width, label='into grid', color='orange')
                ax.bar(x, np.array(to_csc), width, label='collective self consumption',  color='gold')
                
            if 'battery' in balances:
                ax.bar(x, from_grid-from_csc, width, bottom=wind+discharge_battery+from_csc, label='from grid', color='tomato')
                ax.bar(x, np.array(from_csc), width, bottom=discharge_battery+wind,  color='gold')
                ax.bar(x, discharge_battery, width, bottom = wind, label='discharge battery',  color='purple')
                
                ax.bar(x, np.array(into_grid), width, bottom=charge_battery , label='into grid', color='orange')
                if collective == 0:
                    ax.bar(x, np.array(charge_battery), width,  label='charge battery',  color='violet')
                    ax.bar(x, np.array(to_csc), width, bottom=charge_battery, label='collective self consumption',  color='gold')
                else:
                    ax.bar(x, np.array(charge_battery), width, bottom=np.array(to_csc),  label='charge battery',  color='violet')
                    ax.bar(x, np.array(to_csc), width, label='collective self consumption',  color='gold')
                
            if 'electrolyzer' in balances and 'fuel cell' in balances:
                ax.bar(x, from_grid-from_csc, width ,bottom=wind+fc+from_csc, label='from grid', color='tomato')
                ax.bar(x, np.array(from_csc), width, bottom=wind+fc,  color='y')
            
                ax.bar(x, np.array(into_grid), width,bottom=ele , label='into grid', color='orange')
                ax.bar(x, np.array(ele), width,  label='to electrolyzer',  color='chocolate')
                ax.bar(x, fc, width, bottom = wind, label='from fuel cell',  color='peru')
                ax.bar(x, np.array(to_csc), width, bottom=ele, label='collective self consumption',  color='gold')
                
            if 'electrolyzer' in balances:
                ax.bar(x, from_grid-from_csc, width ,bottom=wind+from_csc, label='from grid', color='tomato')
                ax.bar(x, np.array(from_csc), width, bottom=wind,  color='y')
            
                ax.bar(x, np.array(into_grid), width,bottom=ele , label='into grid', color='orange')
                ax.bar(x, np.array(ele), width,  label='to electrolyzer',  color='chocolate')
                ax.bar(x, np.array(to_csc), width, bottom=ele, label='collective self consumption',  color='gold')
            
            if 'mechanical compressor' in balances:
                ax.bar(x, np.array(compressor), width, bottom = ele, label='to compressor',  color='maroon')
                ax.bar(x, np.array(into_grid), width, bottom= ele+compressor , label='into grid', color='orange')
                
        elif 'chp_gt' in balances:
            ax.bar(x, np.minimum(load,chp_gt), width,  label='CHP self consumption', color='yellowgreen')
            ax.bar(x, np.maximum(chp_gt-load,np.zeros(len(chp_gt))), width, bottom=load, label='CHP surplus', color='cornflowerblue')
            
            if 'PV' in balances: 
                ax.bar(x, np.maximum(np.minimum(chp_gt+pv,load)-chp_gt,np.zeros(len(pv))), width, bottom= chp_gt, label='PV self consumption', color='tab:green')
                ax.bar(x, np.maximum(chp_gt-load+pv,np.zeros(len(pv)))-np.maximum(chp_gt-load,np.zeros(len(pv))), width, bottom=load+np.maximum(chp_gt-load,np.zeros(len(pv))),label='PV surplus', color='tab:blue')
                
            if not 'battery' in balances and not 'electrolyzer' in balances:
                ax.bar(x, from_grid-from_csc, width ,bottom=chp_gt+pv+from_csc, label='from grid', color='tomato')
                ax.bar(x, np.array(from_csc), width, bottom=chp_gt+pv, color='gold')
            
                ax.bar(x, np.array(into_grid), width, label='into grid', color='orange')
                ax.bar(x, np.array(to_csc), width, label='collective self consumption',  color='gold')
                
            if 'battery' in balances:
                ax.bar(x, from_grid-from_csc, width, bottom=pv+discharge_battery+from_csc, label='from grid', color='tomato')
                ax.bar(x, np.array(from_csc), width, bottom=discharge_battery+pv,  color='gold')
                ax.bar(x, discharge_battery, width, bottom = pv, label='discharge battery',  color='purple')
                
                ax.bar(x, np.array(into_grid), width, bottom=charge_battery , label='into grid', color='orange')
                if collective == 0:
                    ax.bar(x, np.array(charge_battery), width,  label='charge battery',  color='violet')
                    ax.bar(x, np.array(to_csc), width, bottom=charge_battery, label='collective self consumption',  color='gold')
                else:
                    ax.bar(x, np.array(charge_battery), width, bottom=np.array(to_csc),  label='charge battery',  color='violet')
                    ax.bar(x, np.array(to_csc), width, label='collective self consumption',  color='gold')
            
            if 'electrolyzer' in balances:
                ax.bar(x, from_grid-from_csc, width ,bottom=chp_gt+pv+from_csc, label='from grid', color='tomato')
                ax.bar(x, np.array(from_csc), width, bottom=pv,  color='y')
            
                ax.bar(x, np.array(into_grid), width,bottom=ele , label='into grid', color='orange')
                ax.bar(x, np.array(ele), width,  label='to electrolyzer',  color='chocolate')
                # ax.bar(x, fc, width, bottom = pv, label='from fuel cell',  color='peru')
                ax.bar(x, np.array(to_csc), width, bottom=ele, label='collective self consumption',  color='gold')  
                
            if 'electrolyzer' in balances and 'fuel cell' in balances:
                ax.bar(x, from_grid-from_csc, width ,bottom=pv+fc+from_csc, label='from grid', color='tomato')
                ax.bar(x, np.array(from_csc), width, bottom=pv+fc,  color='y')
            
                ax.bar(x, np.array(into_grid), width,bottom=ele , label='into grid', color='orange')
                ax.bar(x, np.array(ele), width,  label='to electrolyzer',  color='chocolate')
                ax.bar(x, fc, width, bottom = pv, label='from fuel cell',  color='peru')
                ax.bar(x, np.array(to_csc), width, bottom=ele, label='collective self consumption',  color='gold')
            
            if 'mechanical compressor' in balances:
                ax.bar(x, np.array(compressor), width, bottom = ele, label='to compressor',  color='maroon')
                ax.bar(x, np.array(into_grid), width, bottom= ele+compressor , label='into grid', color='orange')
        
        elif 'chp' in balances:
            ax.bar(x, np.minimum(load,chp), width,  label='CHP self consumption', color='yellowgreen')
            ax.bar(x, np.maximum(chp-load,np.zeros(len(chp))), width, bottom=load, label='CHP surplus', color='cornflowerblue')
            
            if 'PV' in balances: 
                ax.bar(x, np.maximum(np.minimum(chp+pv,load)-chp,np.zeros(len(pv))), width, bottom= chp, label='PV self consumption', color='tab:green')
                ax.bar(x, np.maximum(chp-load+pv,np.zeros(len(pv)))-np.maximum(chp-load,np.zeros(len(pv))), width, bottom=load+np.maximum(chp-load,np.zeros(len(pv))),label='PV surplus', color='tab:blue')
                
            if not 'battery' in balances and not 'electrolyzer' in balances:
                ax.bar(x, from_grid-from_csc, width ,bottom=chp+pv+from_csc, label='from grid', color='tomato')
                ax.bar(x, np.array(from_csc), width, bottom=chp+pv, color='gold')
            
                ax.bar(x, np.array(into_grid), width, label='into grid', color='orange')
                ax.bar(x, np.array(to_csc), width, label='collective self consumption',  color='gold')
                
            if 'battery' in balances:
                ax.bar(x, from_grid-from_csc, width, bottom=pv+discharge_battery+from_csc, label='from grid', color='tomato')
                ax.bar(x, np.array(from_csc), width, bottom=discharge_battery+pv,  color='gold')
                ax.bar(x, discharge_battery, width, bottom = pv, label='discharge battery',  color='purple')
                
                ax.bar(x, np.array(into_grid), width, bottom=charge_battery , label='into grid', color='orange')
                if collective == 0:
                    ax.bar(x, np.array(charge_battery), width,  label='charge battery',  color='violet')
                    ax.bar(x, np.array(to_csc), width, bottom=charge_battery, label='collective self consumption',  color='gold')
                else:
                    ax.bar(x, np.array(charge_battery), width, bottom=np.array(to_csc),  label='charge battery',  color='violet')
                    ax.bar(x, np.array(to_csc), width, label='collective self consumption',  color='gold')
            
            if 'electrolyzer' in balances:
                ax.bar(x, from_grid-from_csc, width ,bottom=chp+pv+from_csc, label='from grid', color='tomato')
                ax.bar(x, np.array(from_csc), width, bottom=pv,  color='y')
            
                ax.bar(x, np.array(into_grid), width,bottom=ele , label='into grid', color='orange')
                ax.bar(x, np.array(ele), width,  label='to electrolyzer',  color='chocolate')
                # ax.bar(x, fc, width, bottom = pv, label='from fuel cell',  color='peru')
                ax.bar(x, np.array(to_csc), width, bottom=ele, label='collective self consumption',  color='gold')  
                
            if 'electrolyzer' in balances and 'fuel cell' in balances:
                ax.bar(x, from_grid-from_csc, width ,bottom=pv+fc+from_csc, label='from grid', color='tomato')
                ax.bar(x, np.array(from_csc), width, bottom=pv+fc,  color='y')
            
                ax.bar(x, np.array(into_grid), width,bottom=ele , label='into grid', color='orange')
                ax.bar(x, np.array(ele), width,  label='to electrolyzer',  color='chocolate')
                ax.bar(x, fc, width, bottom = pv, label='from fuel cell',  color='peru')
                ax.bar(x, np.array(to_csc), width, bottom=ele, label='collective self consumption',  color='gold')
                
            if 'mechanical compressor' in balances:
                ax.bar(x, np.array(compressor), width, bottom = ele, label='to compressor',  color='maroon')
                ax.bar(x, np.array(into_grid), width,bottom=ele+compressor , label='into grid', color='orange')
        else:
            ax.bar(x, from_grid-from_csc, width ,bottom=from_csc, label='from grid', color='tomato')
            ax.bar(x, np.array(from_csc), width, label='collective self consumption',  color='gold')
           # plt.ylim(0,3.3)
        
        plt.title(location_name+' days '+str(first_day)+'-'+str(last_day))
        plt.plot(x,load,'k',label='load')     
        plt.legend(ncol=2, bbox_to_anchor = (1.01,-0.11))
        plt.ylabel("Hourly energy [kWh/h] ")
        plt.xlabel( "Time  [h] ")
        if (last_day-first_day) <= 10:
            plt.xticks(list(range(first_day*24, (last_day+1)*24+1,24)), [str(x) for x in list(range(first_day*24, (last_day+1)*24+1,24))], rotation=45)
        # ax.xaxis.set_tick_params(bottom=True,labelbottom=True)
        #plt.xticks([0,6,12,18,24],['0','6','12','18','24'],fontsize=10,color='g')
        #plt.xticks([0,6,12,18,24,30,36,42,48],['0','6','12','18','24','30','36','42','48'],fontsize=10,color='g')
        #plt.yticks([-2,-1,0,1,2],['-2','-1','0','1','2'])
        plt.show()
        # print(from_csc)
        
def hourly_balances_heat(simulation_name,location_name,first_day,last_day,carrier='heat',width=0.9,collective=0):
    
    with open('results/balances_'+simulation_name+'.pkl', 'rb') as f:
        balances = pickle.load(f)
        
    balances = balances[location_name][carrier]
    
    if 'demand' in balances:
        load = -balances['demand'][first_day*24:last_day*24+24]
    else:
        return(print('Thermal Load not considered in the case study'))
        
    if 'grid' in balances:
        into_grid=balances['grid'][first_day*24:last_day*24+24]
    else:
        into_grid=np.zeros(last_day*24+24-first_day*24)
        
    if 'boiler_ng' in balances:                
        boiler_ng = balances['boiler_ng'][first_day*24:last_day*24+24]
    
    if 'boiler_h2' in balances:                
        boiler_h2 = balances['boiler_h2'][first_day*24:last_day*24+24]
        
    if 'fuel cell' in balances:
        fuel_cell=balances['fuel cell'][first_day*24:last_day*24+24]
                
                 
    to_csc = np.zeros(24*(last_day-first_day+1))
    from_csc = np.zeros(24*(last_day-first_day+1))
    for i,e in enumerate(balances['collective self consumption'][first_day*24:last_day*24+24]):
        if e > 0:
            from_csc[i] = e
        else:
            to_csc[i] = e 
        
    x = np.arange(first_day*24,(last_day+1)*24)
    
    fig = plt.figure(dpi=1000)
    from mpl_toolkits.axisartist.axislines import SubplotZero
    ax = SubplotZero(fig, 1, 1, 1)
    fig.add_subplot(ax)
    ax.axis["xzero"].set_visible(True)
    ax.axis["xzero"].label.set_visible(False)
    ax.axis["xzero"].major_ticklabels.set_visible(False)
    for n in ["bottom", "top", "right"]:
        ax.axis[n].set_visible(True)
    ax.grid(axis='y', alpha = 0.5) 
    
    
    if 'fuel cell' in balances:
        ax.bar(x, fuel_cell, width,  label='FC self consumption', color='yellowgreen')
        ax.bar(x, into_grid, width, bottom=load, label='FC surplus', color='cornflowerblue')
        
        if 'boiler_h2' in balances:
            ax.bar(x, boiler_h2-from_csc, width ,bottom=fuel_cell+from_csc, label='from boiler$_\mathregular{H2}$', color='yellow')
            ax.bar(x, np.array(from_csc), width, bottom=fuel_cell,  color='gold')
            
        if 'boiler_ng' in balances:
            ax.bar(x, boiler_ng-from_csc, width ,bottom=fuel_cell+from_csc, label='from boiler$_\mathregular{NG}$', color='tomato')
            ax.bar(x, np.array(from_csc), width, bottom=fuel_cell,  color='gold')
                
    else:
        ax.bar(x, boiler_ng-from_csc, width ,bottom=from_csc, label='from boiler$_\mathregular{NG}$', color='tomato')
        ax.bar(x, np.array(from_csc), width, label='collective self consumption',  color='gold')
       # plt.ylim(0,3.3)
    
    plt.title(location_name+' days '+str(first_day)+'-'+str(last_day)+' heat balance')
    plt.plot(x,load,'k',label='load')     
    plt.legend(ncol=2, bbox_to_anchor=(1.2, 0))
    plt.ylabel("Hourly energy [kWh/h] ")
    plt.xlabel( "Time  [h] ")
    #plt.xticks([0,6,12,18,24],['0','6','12','18','24'],fontsize=10,color='g')
    #plt.xticks([0,6,12,18,24,30,36,42,48],['0','6','12','18','24','30','36','42','48'],fontsize=10,color='g')
    #plt.yticks([-2,-1,0,1,2],['-2','-1','0','1','2'])
    plt.show()
        
def hourly_balances_steam(simulation_name,location_name,first_day,last_day,carrier='process steam',width=0.9,collective=0):
    
    with open('results/balances_'+simulation_name+'.pkl', 'rb') as f:
        balances = pickle.load(f)
        
    balances = balances[location_name][carrier]
    
    if 'demand' in balances:
        load = -balances['demand'][first_day*24:last_day*24+24]
    else:
        return(print('Thermal Load not considered in the case study'))
        
    if 'grid' in balances:
        into_grid=balances['grid'][first_day*24:last_day*24+24]
    else:
        into_grid=np.zeros(last_day*24+24-first_day*24)
        
    if 'chp_gt' in balances:                
        chp = balances['chp_gt'][first_day*24:last_day*24+24]
    elif 'chp' in balances: 
        chp = balances['chp'][first_day*24:last_day*24+24]    
        
    if 'boiler_ng' in balances:                
        boiler_ng = balances['boiler_ng'][first_day*24:last_day*24+24]
    
    if 'boiler_h2' in balances:                
        boiler_h2 = balances['boiler_h2'][first_day*24:last_day*24+24]
        
    if 'fuel cell' in balances:
        fuel_cell=balances['fuel cell'][first_day*24:last_day*24+24]
                
                 
    to_csc = np.zeros(24*(last_day-first_day+1))
    from_csc = np.zeros(24*(last_day-first_day+1))
    for i,e in enumerate(balances['collective self consumption'][first_day*24:last_day*24+24]):
        if e > 0:
            from_csc[i] = e
        else:
            to_csc[i] = e 
        
    x = np.arange(first_day*24,(last_day+1)*24)
    
    fig = plt.figure(dpi=1000)
    from mpl_toolkits.axisartist.axislines import SubplotZero
    ax = SubplotZero(fig, 1, 1, 1)
    fig.add_subplot(ax)
    ax.axis["xzero"].set_visible(True)
    ax.axis["xzero"].label.set_visible(False)
    ax.axis["xzero"].major_ticklabels.set_visible(False)
    for n in ["bottom", "top", "right"]:
        ax.axis[n].set_visible(True)
    ax.grid(axis='y', alpha = 0.5) 
    
    if 'chp_gt' or 'chp' in balances: 
        
        ax.bar(x,np.minimum(chp, load), width,  label='CHP self consumption', color='yellowgreen')
        ax.bar(x,np.maximum(chp-load,np.zeros(len(chp))), width, bottom = load, label='CHP surplus', color='cornflowerblue')
    
        if 'fuel cell' in balances:
            ax.bar(x, fuel_cell, width,  label='FC self consumption', color='yellowgreen')
            ax.bar(x, into_grid, width, bottom=load, label='FC surplus', color='cornflowerblue')
        
        if 'boiler_h2' in balances:
            ax.bar(x, boiler_h2-from_csc, width ,bottom=fuel_cell+from_csc, label='from boiler$_\mathregular{H2}$', color='yellow')
            ax.bar(x, np.array(from_csc), width, bottom=fuel_cell,  color='gold')
            
        if 'boiler_ng' in balances:
            ax.bar(x, boiler_ng-from_csc, width ,bottom=fuel_cell+from_csc, label='from boiler$_\mathregular{NG}$', color='tomato')
            ax.bar(x, np.array(from_csc), width, bottom=fuel_cell,  color='gold')
                
    else:
        ax.bar(x, boiler_ng-from_csc, width ,bottom=from_csc, label='from boiler$_\mathregular{NG}$', color='tomato')
        ax.bar(x, np.array(from_csc), width, label='collective self consumption',  color='gold')
       # plt.ylim(0,3.3)
    
    plt.title(location_name+' days '+str(first_day)+'-'+str(last_day)+' heat balance')
    plt.plot(x,load,'k',label='load')     
    plt.legend(ncol=2, bbox_to_anchor = (0.90,-0.11))
    plt.ylabel("Hourly process steam [kg/h] ")
    plt.xlabel( "Time  [h] ")
    #plt.xticks([0,6,12,18,24],['0','6','12','18','24'],fontsize=10,color='g')
    #plt.xticks([0,6,12,18,24,30,36,42,48],['0','6','12','18','24','30','36','42','48'],fontsize=10,color='g')
    #plt.yticks([-2,-1,0,1,2],['-2','-1','0','1','2'])
    plt.show()     

       
def csc_allocation_sum(simulation_name):
    with open('results/balances_'+simulation_name+'.pkl', 'rb') as f:
        balances = pickle.load(f)
        
    for location_name in balances:
        csc = balances[location_name]['electricity']['collective self consumption']
        from_csc = csc.sum(where=csc>0)
        to_csc = csc.sum(where=csc<0)
    
        print(f"{location_name} {int(from_csc)}  {int(to_csc)}")
   

def Flows(simulation_name,carrier='electricity'):

    pio.renderers.default='browser'
    with open('results/balances_'+simulation_name+'.pkl', 'rb') as f:
        balances = pickle.load(f)
    
    ### data prepearing
    node_label = ["LV grid","LV grid","MV grid","MV grid"]
    node_color = ['silver','silver','gray','gray']
    source = []
    target = []
    value = []
    link_color = []
    link_label = []
        
    ### add from grid
    node_number = 4
    for loc in balances:
        if loc != 'REC':
            
            node_color.append('brown')
            node_label.append('load '+loc)
            
            source.append(0)
            target.append(node_number)
            value.append(np.sum(balances[loc][carrier]['grid'], where = balances[loc][carrier]['grid'] > 0))
            link_color.append('chocolate')
            link_label.append('from LV grid')
            node_number += 1
    value.append(sum(value)-sum(balances['REC'][carrier]['collective self consumption']))
    source.append(2)
    target.append(0)
    link_color.append('tomato')
    link_label.append('from MV grid')
    
    ### add to grid and self consumption
    node_number2 = 4
    for loc in balances:
        if loc != 'REC':
            if 'PV' in balances[loc][carrier]:
            
                node_color.append('skyblue')
                node_label.append('PV '+loc)
                
                # into grid
                source.append(node_number)
                target.append(1)
                value.append(- (np.sum(balances[loc][carrier]['grid'], where = balances[loc][carrier]['grid'] < 0)))
                link_color.append('cornflowerblue')
                link_label.append('to LV grid')
                
                # self consumption
                source.append(node_number)
                target.append(node_number2)
                
                if 'battery' in balances[loc][carrier]:
                    value.append(-np.sum(balances[loc][carrier]['demand']) - np.sum(balances[loc][carrier]['grid'], where = balances[loc][carrier]['grid'] > 0) -np.sum(balances[loc][carrier]['battery'], where = balances[loc][carrier]['battery'] > 0))
                    link_color.append('yellowgreen')
                    link_label.append('self consumption (directly from PV)')
                    node_number += 1
                    
                    node_color.append('plum')
                    node_label.append('battery '+loc)
                    
                    source.append(node_number)
                    target.append(node_number2)
                    value.append(np.sum(balances[loc][carrier]['battery'], where = balances[loc][carrier]['battery'] > 0))
                    link_color.append('purple')
                    link_label.append('self consumption (from battery)')
                    
                    source.append(node_number-1)
                    target.append(node_number)
                    value.append(-np.sum(balances[loc][carrier]['battery'], where = balances[loc][carrier]['battery'] < 0))
                    link_color.append('violet')
                    link_label.append('to battery')   
                    
                    node_number += 1
                    
                elif 'electrolyzer' in balances[loc][carrier] and 'fuel cell' in balances[loc][carrier]:
                    value.append(-np.sum(balances[loc][carrier]['demand']) - np.sum(balances[loc][carrier]['grid'], where = balances[loc][carrier]['grid'] > 0) -np.sum(balances[loc][carrier]['fuel cell']))
                    link_color.append('yellowgreen')
                    link_label.append('self consumption (directly from PV)')
                    node_number += 1
                    
                    node_color.append('sandybrown')
                    node_label.append('hydrogen '+loc)
                    
                    source.append(node_number)
                    target.append(node_number2)
                    value.append(np.sum(balances[loc][carrier]['fuel cell']))
                    link_color.append('peru')
                    link_label.append('self consumption (from fuel cell)')
                    
                    source.append(node_number-1)
                    target.append(node_number)
                    value.append(-np.sum(balances[loc][carrier]['electrolyzer']))
                    link_color.append('chocolate')
                    link_label.append('to electrolyzer')   
                    
                    node_number += 1
                                        
                else:
                    value.append(-np.sum(balances[loc][carrier]['demand']) - np.sum(balances[loc][carrier]['grid'], where = balances[loc][carrier]['grid'] > 0))
                
                    link_color.append('yellowgreen')
                    link_label.append('self consumption')
                    node_number += 1
                    
                    
                    
                node_number2 += 1 # load
                
    ### LV to MV
    source.append(1)
    source.append(1)
    target.append(0)
    target.append(3)
    value.append(sum(balances['REC'][carrier]['collective self consumption'])) 
    value.append(- np.sum(balances['REC'][carrier]['into grid'])-sum(balances['REC'][carrier]['collective self consumption']))   
    
    link_color.append('gold')
    link_label.append('collective self consumption')
    
    link_color.append('orange')
    link_label.append('to MV grid')
                
    fig = go.Figure(data=[go.Sankey(
        valueformat = ".1f",
        valuesuffix = " kWh",
        arrangement = "snap",
        node = {
                'thickness': 30,
                #'line': dict(color = "black", width = 0.5),
                'label': node_label,
                'color': node_color,
                #["x": [0.7,0.3,1,0],
                #"y": [0.1,0.1,0.1,0.1], # stranamente è ribaltata rispetto ad annotations
                'pad': 100
                },
        link = {
                'source': source, # indices correspond to labels, eg A1, A2, A1, B1, ...
                'target': target,
                'value': value,
                'color': link_color,
                'label': link_label
                })])
    
    fig.update_layout(title_text=simulation_name, font_size=25,
                      annotations=[dict(
                          x=0.,
                          y=0.,
                          text='',
                          showarrow=False
                          )] )
    fig.write_html(f"Results/energy_flows_{simulation_name}.html")
    fig.show()    
   
        
def ele_param(simulation_name,first_day,last_day,plot=False):              # functioning parameter of the electrolyzer over the simulation period
    
      last_day = last_day+1
      with open('results/tech_params_'+simulation_name+'.pkl', 'rb') as f:
          param = pickle.load(f)
      
      units = {'efficiency': '[-]',
               'hourly capacity factor': '[%]'}  
      for location_name in param:
         if 'electrolyzer' in param[location_name].keys():
             for parameter in param[location_name]['electrolyzer']:
                 
                 if param[location_name]['electrolyzer'][parameter].size > 1:    # plotting only arrays among saved parameters
                     if plot == True:     
                         plt.figure(dpi=600)
                         y = param[location_name]['electrolyzer'][parameter][first_day*24:last_day*24]
                         x = np.arange(first_day*24,last_day*24)
                         plt.plot(x,y,label=location_name)
                         plt.grid(alpha = 0.3)
                         plt.ylabel(parameter+' '+units[parameter])
                         plt.xlabel('Time [hours]')
                         plt.xticks(list(range(first_day*24, (last_day+1)*24+1,24)), [str(x) for x in list(range(first_day*24, (last_day+1)*24+1,24))], rotation=45)
                         plt.title(location_name+' electrolyzer')
                         plt.show()
                 else:
                     print(f"The global electrolyzers {parameter} calculated for the considered scenario results in {round(param[location_name]['electrolyzer'][parameter],2)} %")
                     return round(param[location_name]['electrolyzer'][parameter],2)
    
def fc_param(simulation_name,first_day,last_day):
    last_day = last_day+1

    with open('results/tech_params_'+simulation_name+'.pkl', 'rb') as f:
        param = pickle.load(f)
        
    units = {'current density': '[A/cm2]',
             'cell voltage'   : '[V]'
             } 
        
    for location_name in param:
       if 'fuel cell' in param[location_name].keys():
           for parameter in param[location_name]['fuel cell']:
          
               plt.figure(dpi=300)
               y = param[location_name]['fuel cell'][parameter][first_day*24:last_day*24]
               x = np.arange(first_day*24,last_day*24)
               plt.plot(x,y,label=location_name)
               plt.grid()
               plt.ylabel(parameter+' '+units[parameter])
               plt.xlabel('Time [hours]')
               plt.xticks(list(range(first_day*24, (last_day+1)*24+1,24)), [str(x) for x in list(range(first_day*24, (last_day+1)*24+1,24))], rotation=45)
               plt.title(location_name+' fuel cell')
               plt.show() 
        