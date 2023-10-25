"""
MESSpy - postprocessing
"""
import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from core import constants as c
from matplotlib.patches import Patch

        
def location_balance(simulation_name,loc,var=None):
    """
    Total balances figures
    
    simulationa_name : str 
    loc : str 
    var : str -> possibility of specifying a single energy carrier
    
    """
    
    with open('results/pkl/balances_'+simulation_name+'.pkl', 'rb') as f: balances = pickle.load(f)

    ###### load analysis
    
    ###### total energy balances
    
    carriers = ['electricity','heating water','gas','hydrogen','HP hydrogen']
    units    = {'electricity': 'kJ', 'hydrogen': 'kg', 'HP hydrogen': 'kg', 'gas':'kJ', 'heating water':'kJ'}
    
    if var: 
        
        carriers = [var]
        units    = {var : units[var]}

        if var == 'hydrogen' and 'mechanical compressor' in balances[loc][var]:
            balances[loc][var].pop('mechanical compressor')  # dict.values() to be removed as they alter the total hydrogen balance and ar enot supposed to do so. 
            
    for carrier in carriers:
        balance  = 0   # initializing the variable to visualize the balance at the end of the simulation period
        print('\n\n'+loc+' '+carrier+' balance: '+simulation_name+'\n') 

        for b in balances[loc][carrier]:
            
            positiv=balances[loc][carrier][b][balances[loc][carrier][b]>0].sum()*c.P2E
            negativ=balances[loc][carrier][b][balances[loc][carrier][b]<0].sum()*c.P2E
            
            if positiv != 0:
                print(b+' '+str(round(positiv,1))+' '+units[carrier])
                balance += positiv
                
            if negativ != 0:
                print(b+' '+str(round(negativ,1))+' '+units[carrier])
                balance += negativ
            
        
def REC_electricity_balance(simulation_name,noprint=False,mounth=False):
    
    with open('results/pkl/balances_'+simulation_name+'.pkl', 'rb') as f:
        balances = pickle.load(f)
        
    df = pd.DataFrame(0.00,columns=["Value [kWh]","Value / production [%]","Value / demand [%]"],
                      index=["SC","CSC","Into l. grid","From l. grid","Into n. grid","From n. grid","Battery losses","Production","Demand"])
    
    
    csc = balances['REC']['electricity']['collective self consumption']*c.P2E/c.kWh2kJ
    into_grid = -balances['REC']['electricity']['into grid']*c.P2E/c.kWh2kJ
    from_grid = balances['REC']['electricity']['from grid']*c.P2E/c.kWh2kJ
    
    
    dmc = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365] # duration of months: cumulate [days]             
    if mounth: # 1-12
         csc = csc [ int(dmc[mounth-1]*24*60/c.timestep) : int(dmc[mounth]*24*60/c.timestep)]
         into_grid = into_grid [ int(dmc[mounth-1]*24*60/c.timestep) : int(dmc[mounth]*24*60/c.timestep)]
         from_grid = from_grid [ int(dmc[mounth-1]*24*60/c.timestep) : int(dmc[mounth]*24*60/c.timestep)]
    
    df.loc['CSC']['Value [kWh]'] = sum(csc)
    df.loc['Into l. grid']['Value [kWh]'] = sum(into_grid)
    df.loc['From l. grid']['Value [kWh]'] = sum(from_grid)
    df.loc['Into n. grid']['Value [kWh]'] = sum(into_grid)  -sum(csc)
    df.loc['From n. grid']['Value [kWh]'] = sum(from_grid) -sum(csc)
    
    
    for loc in balances:
        if loc != 'REC':   
            balance = balances[loc]['electricity']
            if 'demand' in balance:
                demand = balance['demand']
                if mounth: # 1-12
                    demand = demand [ int(dmc[mounth-1]*24*60/c.timestep) : int(dmc[mounth]*24*60/c.timestep)]
                df.loc['Demand']['Value [kWh]'] += -sum(demand)*c.P2E/c.kWh2kJ
            if 'PV' in balance:
                pv = balance['PV']
                if mounth: # 1-12
                    pv = pv [ int(dmc[mounth-1]*24*60/c.timestep) : int(dmc[mounth]*24*60/c.timestep)]
                df.loc['Production']['Value [kWh]'] += sum(pv)*c.P2E/c.kWh2kJ
            if 'wind' in balance:
                wind = balance['wind']
                if mounth: # 1-12
                    pv = pv [ int(dmc[mounth-1]*24*60/c.timestep) : int(dmc[mounth]*24*60/c.timestep)]
                df.loc['Production']['Value [kWh]'] += sum(wind)*c.P2E/c.kWh2kJ
          
    df.loc['SC']['Value [kWh]'] = df.loc['Demand']['Value [kWh]'] - df.loc['From n. grid']['Value [kWh]'] 
    df.loc['Battery losses']['Value [kWh]'] = df.loc['Production']['Value [kWh]'] - df.loc['SC']['Value [kWh]'] - df.loc['Into n. grid']['Value [kWh]']
    
    for b in df.index:
        df.loc[b]['Value / production [%]'] = df.loc[b]['Value [kWh]'] / df.loc['Production']['Value [kWh]'] *100
        df.loc[b]['Value / demand [%]'] = df.loc[b]['Value [kWh]'] / df.loc['Demand']['Value [kWh]'] *100
           
    if not noprint:
        if mounth:
            print('\n'+str(mounth))
        print('\n\nRenewable Energy Community electricity balance: '+simulation_name+'\n') 
        print(df.astype(int))
        print('\n')
    return(df.astype(int))


def hist_12_balances_pc(simulation_name,ymax):
    
    sc = []
    csc = []
    into_grid = []
    from_grid = []
    losses = []
    
    for m in np.arange(12):
        balances = REC_electricity_balance(simulation_name, noprint=True, mounth=m+1)
        sc.append(balances['Value [kWh]']['SC'])
        csc.append(balances['Value [kWh]']['CSC'])
        into_grid.append(balances['Value [kWh]']['Into n. grid'])
        from_grid.append(balances['Value [kWh]']['From n. grid'])
        losses.append(balances['Value [kWh]']['Battery losses'])
        
    sc = np.array(sc)
    csc = np.array(csc)
    x = np.arange(12)  # the label locations
    width = 0.8  # the width of the bars
    
    # Create a figure with two subplots arranged horizontally (1 row, 2 columns)
    fig, axes = plt.subplots(1, 2, dpi=1000, figsize=(12, 4))
    
    ax1 = axes[0]
    ax2 = axes[1]
    
    # Plot the data for the first subplot
    ax1.bar(x, sc, width, label='Self-Consumption', color='yellowgreen')
    ax1.bar(x, losses, width, bottom=sc, label='Battery losses', color='tab:purple')
    ax1.bar(x, csc, width, bottom=sc+losses, label='Collective-Self-Consumption', color='gold')
    ax1.bar(x, into_grid, width, bottom=sc+csc+losses, label='Into national grid', color='tab:blue')

    ax1.legend()
    ax1.set_ylabel('Energy produced [kWh/month]')
    ax1.set_xticks(np.arange(12))
    ax1.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
        
    ax1.set_ylim(0, ymax)
    ax1.set_title('Production '+simulation_name)
    ax1.grid(axis='y',zorder=-10)
    
    # Plot the data for the second subplot
    ax2.bar(x, sc, width, label='Self-Sufficiency', color='yellowgreen')
    ax2.bar(x, csc, width, bottom=sc, label='Collective-Self-Sufficiency', color='gold')
    ax2.bar(x, from_grid, width, bottom=sc+csc, label='From national grid', color='tomato')
    
    ax2.legend()
    ax2.set_ylabel('Energy consumed [kWh/month]')
    ax2.set_xticks(np.arange(12))
    ax2.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])

    ax2.set_ylim(0, ymax)
    ax2.set_title('Demand '+simulation_name)
    ax2.grid(axis='y',zorder=-10)
    
    # Adjust the layout to prevent overlapping of titles and labels
    plt.tight_layout()
    
    plt.show()
        
def NPV_plot(study_case):
    ##### economic
    with open('results/pkl/economic_assessment_'+study_case+'.pkl', 'rb') as f: economic = pickle.load(f)
    
    plt.figure(dpi=1000)
    
    for loc in economic:    
        y = economic[loc]['NPV']
        x = np.linspace(0,len(y)-1,len(y))
        plt.plot(x,y,label=loc)
    plt.plot(x,np.zeros(len(x)),color='k')
        
    plt.legend()
    plt.title("Investments")
    plt.grid()
    plt.ylabel('Net Present Value [€]')
    plt.xlabel('Time [years]')
    plt.xlim(0,len(y)-1)
    plt.show()
                
def hourly_balances_electricity(simulation_name,location_name,first_day,last_day,carrier='electricity',width=0.9,collective=0):
    
        with open('results/pkl/balances_'+simulation_name+'.pkl', 'rb') as f:
            balances = pickle.load(f)
            
        balances = balances[location_name][carrier]
        
        if 'demand' in balances:
            load = -balances['demand'][first_day*24*60//c.timestep:last_day*24*60//c.timestep+24*60//c.timestep]
        else:
            load = np.zeros(last_day*24*60/c.timestep+24*60/c.timestep - first_day*24*60/c.timestep)
            
        if 'chp_gt' in balances:
            chp_gt = balances['chp_gt'][first_day*24*60//c.timestep:last_day*24*60//c.timestep+24*60//c.timestep]
            
        if 'chp' in balances:
            chp = balances['chp'][first_day*24*60//c.timestep:last_day*24*60//c.timestep+24*60//c.timestep]
        
        if 'PV' in balances:
            pv = balances['PV'][first_day*24*60//c.timestep:last_day*24*60//c.timestep+24*60//c.timestep]
        
        if 'wind' in balances:
            wind = balances['wind'][first_day*24*60//c.timestep:last_day*24*60//c.timestep+24*60//c.timestep]
            
        if 'grid' in balances:                
            into_grid = np.zeros(24*(last_day-first_day+1)*60//c.timestep)
            from_grid = np.zeros(24*(last_day-first_day+1)*60//c.timestep)
            for i,e in enumerate(balances['grid'][first_day*24*60//c.timestep:last_day*24*60//c.timestep+24*60//c.timestep]):
                if e > 0:
                    from_grid[i] = e
                else:
                    into_grid[i] = e 
        
        if 'battery' in balances:
            charge_battery = np.zeros(24*(last_day-first_day+1)*60//c.timestep)
            discharge_battery = np.zeros(24*(last_day-first_day+1)*60//c.timestep)
            for i,e in enumerate(balances['battery'][first_day*24*60//c.timestep:last_day*24*60//c.timestep+24*60//c.timestep]):
                if e > 0:
                    discharge_battery[i] = e
                else:
                    charge_battery[i] = e 
        
        if 'electrolyzer' in balances:
            ele = balances['electrolyzer'][first_day*24*60//c.timestep:last_day*24*60//c.timestep+24*60//c.timestep]
        
        if 'mechanical compressor' in balances:
            compressor = balances['mechanical compressor'][first_day*24*60//c.timestep:last_day*24*60//c.timestep+24*60//c.timestep]
        
        if 'fuel cell' in balances:
            fc = balances['fuel cell'][first_day*24*60//c.timestep:last_day*24*60//c.timestep+24*60//c.timestep]

                     
        to_csc = np.zeros(24*(last_day-first_day+1)*60//c.timestep)
        from_csc = np.zeros(24*(last_day-first_day+1)*60//c.timestep)
        for i,e in enumerate(balances['collective self consumption'][first_day*24*60//c.timestep:last_day*24*60//c.timestep+24*60//c.timestep]):
            if e > 0:
                from_csc[i] = e
            else:
                to_csc[i] = e 
            
        x = np.arange(first_day*24*60/c.timestep,(last_day+1)*24*60/c.timestep)
        
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
        plt.ylabel("Power [kW] ")
        plt.xlabel( f"Time  [hour] ")
        if (last_day-first_day) <= 10:
             plt.xticks(list(range(first_day*24*60//c.timestep, (last_day+1)*24*60//c.timestep+1,24*60//c.timestep)), [str(x) for x in list(range(first_day*24, (last_day+1)*24+1,24))], rotation=45)
        # ax.xaxis.set_tick_params(bottom=True,labelbottom=True)
        #plt.xticks([0,6,12,18,24],['0','6','12','18','24'],fontsize=10,color='g')
        #plt.xticks([0,6,12,18,24,30,36,42,48],['0','6','12','18','24','30','36','42','48'],fontsize=10,color='g')
        #plt.yticks([-2,-1,0,1,2],['-2','-1','0','1','2'])
        plt.show()

def LOC_plot(simulation_name):
      
    with open('results/pkl/LOC_'+simulation_name+'.pkl', 'rb') as f:
        LOC = pickle.load(f)
           
    unit = {'H tank': '[kg]', 'HPH tank': '[kg]', 'battery': '[kJ]', 'inertial TES': '[°C]'}
    
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

def csc_allocation_sum(simulation_name):
    
    print('\n'+'Collective-Self-Consumption proportional contribution') 
    
    with open('results/pkl/balances_'+simulation_name+'.pkl', 'rb') as f:        balances = pickle.load(f)
        
    for location_name in balances:
        csc = balances[location_name]['electricity']['collective self consumption']*c.P2E/c.kWh2kJ
        from_csc = csc.sum(where=csc>0)
        to_csc = csc.sum(where=csc<0)
    
        print(f"{location_name} {int(from_csc)}  {int(to_csc)}")