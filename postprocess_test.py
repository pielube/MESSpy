"""
MESSpy - postprocessing
"""
import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from core import constants as c
import matplotlib.patches as mpatches

        
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
    
    df.loc['CSC', 'Value [kWh]'] = sum(csc)
    df.loc['Into l. grid', 'Value [kWh]'] = sum(into_grid)
    df.loc['From l. grid', 'Value [kWh]'] = sum(from_grid)
    df.loc['Into n. grid', 'Value [kWh]'] = sum(into_grid)  -sum(csc)
    df.loc['From n. grid', 'Value [kWh]'] = sum(from_grid) -sum(csc)
    
    
    for loc in balances:
        if loc != 'REC':   
            balance = balances[loc]['electricity']
            if 'demand' in balance:
                demand = balance['demand']
                if mounth: # 1-12
                    demand = demand [ int(dmc[mounth-1]*24*60/c.timestep) : int(dmc[mounth]*24*60/c.timestep)]
                df.loc['Demand', 'Value [kWh]'] += -sum(demand)*c.P2E/c.kWh2kJ
            if 'PV' in balance:
                pv = balance['PV']
                if mounth: # 1-12
                    pv = pv [ int(dmc[mounth-1]*24*60/c.timestep) : int(dmc[mounth]*24*60/c.timestep)]
                df.loc['Production', 'Value [kWh]'] += sum(pv)*c.P2E/c.kWh2kJ
            if 'wind' in balance:
                wind = balance['wind']
                if mounth: # 1-12
                    pv = pv [ int(dmc[mounth-1]*24*60/c.timestep) : int(dmc[mounth]*24*60/c.timestep)]
                df.loc['Production', 'Value [kWh]'] += sum(wind)*c.P2E/c.kWh2kJ
          
    df.loc['SC', 'Value [kWh]'] = df.loc['Demand', 'Value [kWh]'] - df.loc['From n. grid', 'Value [kWh]'] 
    df.loc['Battery losses', 'Value [kWh]'] = df.loc['Production', 'Value [kWh]'] - df.loc['SC', 'Value [kWh]'] - df.loc['Into n. grid', 'Value [kWh]']
    
    for b in df.index:
        df.loc[b, 'Value / production [%]'] = df.loc[b, 'Value [kWh]'] / df.loc['Production', 'Value [kWh]'] * 100
        df.loc[b, 'Value / demand [%]'] = df.loc[b, 'Value [kWh]'] / df.loc['Demand', 'Value [kWh]'] *100
           
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
            xticklabels = ['          Jan','         Feb','          Mar','         Apr','         May','          Jun','         Jul','          Aug','           Sep','          Oct','          Nov','           Dec','']
            plt.xticks(xticks,xticklabels,rotation=45)
            plt.title(location_name+' '+tech)
            plt.xlim(0,x[-1])
            plt.show()

def csc_allocation_sum(simulation_name):
    
    print('\n'+'Collective-Self-Consumption proportional contribution') 
    
    with open('results/pkl/balances_'+simulation_name+'.pkl', 'rb') as f:        balances = pickle.load(f)
        
    for location_name in balances:
        csc = balances[location_name]['electricity']['collective self consumption']*c.P2E/c.kWh2kJ
        from_csc = csc.sum(where=csc>0)
        to_csc = csc.sum(where=csc<0)
    
        print(f"{location_name} {int(from_csc)}  {int(to_csc)}")
        
def satisfaction_story(sati):
    mode = {0: 'no demand', 
            1: 'demand satisfied by iTES without switching on the HP ',
            2: 'demand satisfied by switching on the HP',
            3: 'demand satisfied by switching on the HP and iTES overheated',
            -1: 'unsatisfied demand because there is not enough power',
            -2: 'unsatisfied demand because i TES cant reach the minimum required temperature in one step',
            -3: 'unsatisfied demand becasue t amb is too cold to heat water to the desired temperature'       
            }
    
    print('\n history of heating demand satisfaction:')
    for m in mode:
        print(f"{mode[m]}: {(sati == m).sum()} steps")
  
def cop(cop):
    plt.figure(dpi=1000)
    plt.scatter(np.arange(len(cop)),cop,s=1)
    plt.ylabel('COP [-]')
    xticks = list(np.linspace(0, len(cop) - 1, 13).astype(int))
    xticklabels = ['          Jan','         Feb','          Mar','         Apr','         May','          Jun','         Jul','          Aug','           Sep','          Oct','          Nov','           Dec','']
    plt.xticks(xticks,xticklabels,rotation=45)
    plt.grid()
    plt.title('HP Coefficient of Performance')
    plt.xlim(0,len(cop))
    plt.show()
    
def heating_demand(demand):
    plt.figure(dpi=1000)
    plt.plot(np.arange(len(demand)),demand,color='red')
    xticks = list(np.linspace(0, len(demand) - 1, 13).astype(int))
    xticklabels = ['          Jan','         Feb','          Mar','         Apr','         May','          Jun','         Jul','          Aug','           Sep','          Oct','          Nov','           Dec','']
    plt.xticks(xticks,xticklabels,rotation=45)
    plt.ylabel('Demand [kW]')
    plt.grid()
    plt.xlim(0,len(demand))
    plt.title('Heating water demand')
    plt.show()
    
def energy_balance_results(studycase,simulation_name,loc,print_=False,plot=False):
    
    """
    Total balances figures and hydrogen-related ghg emissions calculation
    ----------
    studycase : dictionary (all the inputs are optional)
        'location_1_name': inputs required to create a location object (see Location.py)
        'location_2_name': 
            ...
        'location_n_name':
            
    simulation_name : str
    
    simulation_years : str - years to be simulated
    
    loc : str - location_name

    emission_intensity : dict
    
    ref_year: str -> necessary to specify the reference year for grid intensity
                        
    output:  balances_pp: dictionary containing total balances calculation useful for LCOH calculation, NPV calculation and post process plots

    """  

    with open('results/pkl/balances_'+simulation_name+'.pkl', 'rb') as f: balances = pickle.load(f)

    if 'wind' not in balances[loc]['electricity'] and 'PV' not in balances[loc]['electricity']:    # if no renewables are included in the location        
        print("No renewable energy sources are present in the considered case study\nAutoconsumption data calculation not available\n")
        return

    else: 
        windsc      = {}
        pvsc        = {}
        from_grid   = {}
        wind_autoconsumption    = np.zeros(c.timestep_number)
        pv_autoconsumption      = np.zeros(c.timestep_number)
        
        if 'PV' in balances[loc]['electricity']:
            pv_still_available = balances[loc]['electricity']['PV'].copy()
        else:
            pv_still_available = np.zeros(c.timestep_number)
        
        if 'wind' in balances[loc]['electricity']:
            wind_still_available = balances[loc]['electricity']['wind'].copy()
        else:
            wind_still_available = np.zeros(c.timestep_number)
        
        system = dict(sorted(studycase[loc].items(), key=lambda item: item[1]['priority'])) # ordered by priority
                
            
        if 'wind' in balances[loc]['electricity'] and 'PV' not in balances[loc]['electricity']:    
            for tech_name in system: # (which is ordered py priority)
                if tech_name == 'electricity demand':
                    tech_name = 'demand'
                if tech_name in balances[loc]['electricity'] and all(num <= 0 for num in balances[loc]['electricity'][tech_name]) and tech_name != 'collective self consumption':
                    windsc[tech_name]       = np.zeros(c.timestep_number)  # initializing the array to store wind energy self consumption values for tech_name technology
                    from_grid[tech_name]    = np.zeros(c.timestep_number)
                    for i in range(len(balances[loc]['electricity'][tech_name])):
                        if balances[loc]['electricity'][tech_name][i] < 0:
                            windsc[tech_name][i]    = np.minimum(wind_still_available[i],-balances[loc]['electricity'][tech_name][i])
                            from_grid[tech_name][i] = -balances[loc]['electricity'][tech_name][i] - windsc[tech_name][i]
                            wind_still_available[i] = wind_still_available[i] - windsc[tech_name][i]
                            wind_autoconsumption[i] += windsc[tech_name][i]
                    if print_ == True:
                        print('\nwindsc_'+str(tech_name)+'    = '+str(round((sum(windsc[tech_name]*(c.timestep/60))/1e3),2))+' MWh') 
                        print('from_grid_'+str(tech_name)+' = '+str(round((sum(from_grid[tech_name]*(c.timestep/60))/1e3),2))+' MWh')
                    
            wind_surplus = wind_still_available
            if print_ == True:
                print('\nwind_autoconsumption     = '+str(round((sum(wind_autoconsumption*(c.timestep/60))/1e3),2))+' MWh')
                print('wind_surplus             = '+str(round((sum(wind_surplus*(c.timestep/60))/1e3),2))+' MWh')

                    
        elif 'PV' in balances[loc]['electricity'] and 'wind' not in balances[loc]['electricity']:    
            for tech_name in system: # (which is ordered py priority)
                if tech_name == 'electricity demand':
                    tech_name = 'demand'
                if tech_name in balances[loc]['electricity'] and all(num <= 0 for num in balances[loc]['electricity'][tech_name]) and tech_name != 'collective self consumption':
                    pvsc[tech_name] = np.zeros(c.timestep_number)  # initializing the array to store wind energy self consumption values for tech_name technology
                    from_grid[tech_name] = np.zeros(c.timestep_number)
                    for i in range(len(balances[loc]['electricity'][tech_name])):
                        if balances[loc]['electricity'][tech_name][i] < 0:
                            pvsc[tech_name][i] = np.minimum(pv_still_available[i],-balances[loc]['electricity'][tech_name][i])
                            from_grid[tech_name][i] = -balances[loc]['electricity'][tech_name][i] - pvsc[tech_name][i]
                            pv_still_available[i] = pv_still_available[i] - pvsc[tech_name][i]
                            pv_autoconsumption[i] += pvsc[tech_name][i]
                    if print_ == True:
                        print('\npvsc_'+str(tech_name)+'      = '+str(round((sum(pvsc[tech_name]*(c.timestep/60))/1e3),2))+' MWh') 
                        print('from_grid_'+str(tech_name)+' = '+str(round((sum(from_grid[tech_name]*(c.timestep/60))/1e3),2))+' MWh')
                        
            pv_surplus = pv_still_available
            if print_ == True:
                print('\npv_autoconsumption    = '+str(round((sum(pv_autoconsumption*(c.timestep/60))/1e3),2))+' MWh')
                print('pv_surplus           = '+str(round((sum(pv_surplus)*(c.timestep/60)/1e3),2))+' MWh')


        elif 'PV' in balances[loc]['electricity'] and 'wind' in balances[loc]['electricity']:
            
            if list(system.keys()).index('PV') < list(system.keys()).index('wind'): # PV comes before wind
                
                for tech_name in system: # (which is ordered py priority)
                    if tech_name == 'electricity demand':
                        tech_name = 'demand'
                    if tech_name in balances[loc]['electricity'] and all(num <= 0 for num in balances[loc]['electricity'][tech_name]) and tech_name != 'collective self consumption':
                        pvsc[tech_name] = np.zeros(c.timestep_number)  # initializing the array to store pv energy self consumption values for tech_name technology
                        windsc[tech_name] = np.zeros(c.timestep_number)  # initializing the array to store wind energy self consumption values for tech_name technology
                        from_grid[tech_name] = np.zeros(c.timestep_number)
                        for i in range(len(balances[loc]['electricity'][tech_name])):
                            if balances[loc]['electricity'][tech_name][i] < 0:
                                pvsc[tech_name][i] = np.minimum(pv_still_available[i],-balances[loc]['electricity'][tech_name][i])
                                windsc[tech_name][i] = np.minimum(wind_still_available[i],-balances[loc]['electricity'][tech_name][i]-pvsc[tech_name][i])
                                from_grid[tech_name][i] = -balances[loc]['electricity'][tech_name][i] - pvsc[tech_name][i] - windsc[tech_name][i]
                                pv_still_available[i] = pv_still_available[i] - pvsc[tech_name][i]
                                wind_still_available[i] = wind_still_available[i] - windsc[tech_name][i]
                                pv_autoconsumption[i] += pvsc[tech_name][i]
                                wind_autoconsumption[i] += windsc[tech_name][i]
                        if print_ == True:
                            print('\npvsc_'+str(tech_name)+'      = '+str(round((sum(pvsc[tech_name]*(c.timestep/60))/1e3),2))+' MWh')
                            print('windsc_'+str(tech_name)+'    = '+str(round((sum(windsc[tech_name]*(c.timestep/60))/1e3),2))+' MWh') 
                            print('from_grid_'+str(tech_name)+' = '+str(round((sum(from_grid[tech_name]*(c.timestep/60))/1e3),2))+' MWh')
                             
            else:
                
                for tech_name in system: # (which is ordered py priority)
                    if tech_name == 'electricity demand':
                        tech_name = 'demand'
                    if tech_name in balances[loc]['electricity'] and all(num <= 0 for num in balances[loc]['electricity'][tech_name]) and tech_name != 'collective self consumption':
                        pvsc[tech_name] = np.zeros(c.timestep_number)  # initializing the array to store pv energy self consumption values for tech_name technology
                        windsc[tech_name] = np.zeros(c.timestep_number)  # initializing the array to store wind energy self consumption values for tech_name technology
                        from_grid[tech_name] = np.zeros(c.timestep_number)
                        for i in range(len(balances[loc]['electricity'][tech_name])):
                            if balances[loc]['electricity'][tech_name][i] < 0:
                                windsc[tech_name][i] = np.minimum(wind_still_available[i],-balances[loc]['electricity'][tech_name][i])
                                pvsc[tech_name][i] = np.minimum(pv_still_available[i],-balances[loc]['electricity'][tech_name][i]-windsc[tech_name][i])
                                from_grid[tech_name][i] = -balances[loc]['electricity'][tech_name][i] - windsc[tech_name][i] - pvsc[tech_name][i]
                                pv_still_available[i] = pv_still_available[i] - pvsc[tech_name][i]
                                wind_still_available[i] = wind_still_available[i] - windsc[tech_name][i]
                                pv_autoconsumption[i] += pvsc[tech_name][i]
                                wind_autoconsumption[i] += windsc[tech_name][i]
                        if print_ == True:
                            print('\nwindsc_'+str(tech_name)+'    = '+str(round((sum(windsc[tech_name]*(c.timestep/60))/1e3),2))+' MWh')
                            print('pvsc_'+str(tech_name)+'      = '+str(round((sum(pvsc[tech_name]*(c.timestep/60))/1e3),2))+' MWh') 
                            print('from_grid_'+str(tech_name)+' = '+str(round((sum(from_grid[tech_name]*(c.timestep/60))/1e3),2))+' MWh')

            wind_surplus = wind_still_available
            pv_surplus = pv_still_available
            if print_ == True:
                print('\nwind_autoconsumption     = '+str(round((sum(wind_autoconsumption*(c.timestep/60))/1e3),2))+' MWh')
                print('pv_autoconsumption       = '+str(round((sum(pv_autoconsumption*(c.timestep/60))/1e3),2))+' MWh')
                print('wind_surplus             = '+str(round((sum(wind_surplus*(c.timestep/60))/1e3),2))+' MWh')
                print('pv_surplus               = '+str(round((sum(pv_surplus*(c.timestep/60))/1e3),2))+' MWh')
            
  
        # Saving all technologies interacting with hydrogen
        tech_name1 = []
        tot_demand_for_hyd  = np.zeros(c.timestep_number)
        el_from_grid_hyd    = np.zeros(c.timestep_number)
        el_from_wind_hyd    = np.zeros(c.timestep_number)
        el_from_pv_hyd      = np.zeros(c.timestep_number)
        for tech_name in balances[loc]['hydrogen']:
            if tech_name != 'grid' and tech_name != 'demand':
                tech_name1.append(tech_name)
        
        # Calculating the electricity consumed to run the technologies interacting with hydrogen (only those having negative values, thus absorbing electricity)
        for tech_name in tech_name1:
            if tech_name in balances[loc]['electricity'] and all(num <= 0 for num in balances[loc]['electricity'][tech_name]):
                for i in range(len(balances[loc]['electricity'][tech_name])):
                    if balances[loc]['electricity'][tech_name][i] < 0:
                        tot_demand_for_hyd[i] += -balances[loc]['electricity'][tech_name][i]
                        el_from_grid_hyd[i] += from_grid[tech_name][i]
                        if tech_name in windsc:
                            el_from_wind_hyd[i] += windsc[tech_name][i]
                        if tech_name in pvsc:
                            el_from_pv_hyd[i] += pvsc[tech_name][i] 
        
        balances_pp = balances.copy()
        balances_pp[loc]['electricity']['hyd tot electricity']          = tot_demand_for_hyd
        balances_pp[loc]['electricity']['hyd grid electricity']         = el_from_grid_hyd
        balances_pp[loc]['electricity']['tech from grid']               = from_grid
        
        if 'wind' in balances[loc]['electricity']:
            balances_pp[loc]['electricity']['tech windsc']              = windsc
            balances_pp[loc]['electricity']['hyd wind electricity']     = el_from_wind_hyd
            balances_pp[loc]['electricity']['wind surplus']             = wind_surplus
            balances_pp[loc]['electricity']['wind autoconsumption']     = wind_autoconsumption
        if 'PV' in balances[loc]['electricity']:        
            balances_pp[loc]['electricity']['tech pvsc']                = pvsc
            balances_pp[loc]['electricity']['hyd pv electricity']       = el_from_pv_hyd
            balances_pp[loc]['electricity']['pv surplus']               = pv_surplus
            balances_pp[loc]['electricity']['pv autoconsumption']       = pv_autoconsumption
        
        if plot == True:
            k=1
            grid_cumulative_for_hyd     = np.cumsum(el_from_grid_hyd*(c.timestep/60))       # kW to kWh
            demand_cumulative_for_hyd   = np.cumsum(tot_demand_for_hyd*(c.timestep/60))     # kW to kWh
            if sum(el_from_wind_hyd) != 0:
                wind_sc_cumulative_for_hyd  = np.cumsum(el_from_wind_hyd*(c.timestep/60))   # kW to kWh
            if sum(el_from_pv_hyd) != 0:
                pv_sc_cumulative_for_hyd    = np.cumsum(el_from_pv_hyd*(c.timestep/60))     # kW to kWh
                
            # Cumulative plot
            fig, ax = plt.subplots(dpi=300)
            ax.plot(demand_cumulative_for_hyd, label = 'Electricity Dem for H2 chain', linewidth=k, zorder = 2)
            ax.plot(grid_cumulative_for_hyd, label = 'From grid', linewidth=k-0.2)
            if sum(el_from_wind_hyd) != 0:
                ax.plot(wind_sc_cumulative_for_hyd,  alpha = 1, zorder = 1, label = 'P$_{Wind SC}$', linewidth=k)
            if sum(el_from_pv_hyd) != 0:
                ax.plot(pv_sc_cumulative_for_hyd, label = 'P$_{PV SC}$', linewidth=k, zorder = 3)
            ax.grid(alpha = 0.3, zorder = 0)
            xticks = list(np.linspace(0, len(demand_cumulative_for_hyd) - 1, 13).astype(int))
            xticklabels = ['          Jan','         Feb','          Mar','         Apr','         May','          Jun','         Jul','          Aug','           Sep','          Oct','          Nov','           Dec','']
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticklabels)
            ax.set_ylabel('Energy [kWh]')
            ax.set_xlabel('Time [h]')
            ax.legend()
            plt.show()
            
            # Coloured areas
            fig, ax = plt.subplots(dpi=300)
            el_dem_H2_chain, = ax.plot(demand_cumulative_for_hyd, label='Electricity Dem for H2 chain', linewidth=k, zorder=2)
            if sum(el_from_wind_hyd) != 0 and sum(el_from_pv_hyd) == 0:
                ax.fill_between(range(len(demand_cumulative_for_hyd)), wind_sc_cumulative_for_hyd, color='green', alpha=0.3)
                ax.fill_between(range(len(demand_cumulative_for_hyd)), wind_sc_cumulative_for_hyd+grid_cumulative_for_hyd, wind_sc_cumulative_for_hyd, color='red', alpha=0.3)
                red_patch = mpatches.Patch(color='red', alpha=0.3, label='From grid')
                green_patch = mpatches.Patch(color='green', alpha=0.3, label='Wind$_{SC}$')
                ax.legend(handles=[el_dem_H2_chain, red_patch, green_patch])
            if sum(el_from_pv_hyd) != 0 and sum(el_from_wind_hyd) == 0:
                ax.fill_between(range(len(demand_cumulative_for_hyd)), pv_sc_cumulative_for_hyd, color='green', alpha=0.3)
                ax.fill_between(range(len(demand_cumulative_for_hyd)), pv_sc_cumulative_for_hyd+grid_cumulative_for_hyd, pv_sc_cumulative_for_hyd, color='red', alpha=0.3)
                red_patch = mpatches.Patch(color='red', alpha=0.3, label='From grid')
                green_patch = mpatches.Patch(color='green', alpha=0.3, label='PV$_{SC}$')
                ax.legend(handles=[el_dem_H2_chain, red_patch, green_patch])
            if sum(el_from_wind_hyd) != 0 and sum(el_from_pv_hyd) != 0:
                if list(system.keys()).index('PV') < list(system.keys()).index('wind'): # PV comes before wind
                    ax.fill_between(range(len(demand_cumulative_for_hyd)), pv_sc_cumulative_for_hyd, color='green', alpha=0.3)
                    ax.fill_between(range(len(demand_cumulative_for_hyd)), pv_sc_cumulative_for_hyd+wind_sc_cumulative_for_hyd, pv_sc_cumulative_for_hyd, color='orange', alpha=0.3)
                    ax.fill_between(range(len(demand_cumulative_for_hyd)), pv_sc_cumulative_for_hyd+wind_sc_cumulative_for_hyd+grid_cumulative_for_hyd, pv_sc_cumulative_for_hyd+wind_sc_cumulative_for_hyd, color='red', alpha=0.3)
                    red_patch = mpatches.Patch(color='red', alpha=0.3, label='From grid')
                    green_patch = mpatches.Patch(color='green', alpha=0.3, label='PV$_{SC}$')
                    orange_patch = mpatches.Patch(color='orange', alpha=0.3, label='Wind$_{SC}$')
                    ax.legend(handles=[el_dem_H2_chain, red_patch, green_patch, orange_patch])
                else:
                    ax.fill_between(range(len(demand_cumulative_for_hyd)), wind_sc_cumulative_for_hyd, color='green', alpha=0.3)
                    ax.fill_between(range(len(demand_cumulative_for_hyd)), wind_sc_cumulative_for_hyd+pv_sc_cumulative_for_hyd, wind_sc_cumulative_for_hyd, color='orange', alpha=0.3)
                    ax.fill_between(range(len(demand_cumulative_for_hyd)), wind_sc_cumulative_for_hyd+pv_sc_cumulative_for_hyd+grid_cumulative_for_hyd, wind_sc_cumulative_for_hyd+pv_sc_cumulative_for_hyd, color='red', alpha=0.3)
                    red_patch = mpatches.Patch(color='red', alpha=0.3, label='From grid')
                    green_patch = mpatches.Patch(color='green', alpha=0.3, label='Wind$_{SC}$')
                    orange_patch = mpatches.Patch(color='orange', alpha=0.3, label='PV$_{SC}$')
                    ax.legend(handles=[el_dem_H2_chain, red_patch, green_patch, orange_patch])
            
            ax.grid(alpha=0.3, zorder=0)
            xticks = list(np.linspace(0, len(demand_cumulative_for_hyd) - 1, 13).astype(int))
            xticklabels = ['          Jan','         Feb','          Mar','         Apr','         May','          Jun','         Jul','          Aug','           Sep','          Oct','          Nov','           Dec','']
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticklabels)
            ax.set_ylabel('Energy [kWh]')
            ax.set_xlabel('Time [h]')            
            plt.show()
                    
        # saving post process results
        with open('results/pkl/balances_pp_'+simulation_name+".pkl", 'wb') as f: pickle.dump(balances_pp, f)

    return balances_pp

def ghg_emissions(simulation_name,loc,energy_market, print_= False):
    """
    GHG emissions profile (green index) calculation for hydrogen production.
    Dependent on energy_balance_results function output.
    ----------
            
    simulation_name : str - name of energy_balances file .pkl where balances_pp results are stored
        
    loc : str - location_name

    emission_intensity : dict
    
    ref_year: str -> necessary to specify the reference year for grid intensity
                        
    output: float - ghg emission value [kgCO2/kgH2]
    """      

    with open('results/pkl/balances_pp_'+simulation_name+".pkl",'rb') as f: balances_pp = pickle.load(f) 

    el_from_grid_hyd_tot    = sum(balances_pp[loc]['electricity']['hyd grid electricity'])*c.timestep/60        # [kWh] total amount of electricity withdrown from grid to power hydrogen production chain
    emission_factor         = energy_market['electricity']['emission intensity']                                # [gCO2/kWh] electric grid intensity depending on the selected area for the analysis
    co2_tot                 = el_from_grid_hyd_tot*emission_factor/1000                                         # [kgCO2] total amount of carbon dioxide due to grid electricity utilization
    produced_hyd            = sum(balances_pp[loc]['hydrogen']['electrolyzer']*(c.timestep*60))                 # [kgH2] total amount of produced hydrogen via in situ electorlysis
    
    h2_ghg = round(co2_tot/produced_hyd,2)  # [kgCO2/kgH2] GHG intensity of the produced hydrogen
    if print_ == True:
        print(f"\nThe H2 GHG intensity calculated for the considered scenario results in {h2_ghg} kgCO2/kgH2")
    
    return h2_ghg

def hydrogen_production(simulation_name,loc,var=None):
    """
    Hydrogen production figures and plots
    
    simulationa_name : str 
    loc : str 
    var : str -> possibility of specifying **kargs
    
    """    
    with open('results/pkl/balances_'+simulation_name+'.pkl', 'rb') as f: balances = pickle.load(f)
    
    simulation_steps    = len(balances[loc]['hydrogen']['electrolyzer'])
    constantflow        = sum(balances[loc]['hydrogen']['electrolyzer'])/simulation_steps       # [kg/s] constant flow rate deliverable by the system if operated in supply-led mode
    demand_constant     = np.array([constantflow * (i + 1) for i in range(simulation_steps)])   # [kg/s] fictitious constant demand 
    # demand_variable     = -np.cumsum(balances[loc]['hydrogen']['demand'])
    production          = np.cumsum(balances[loc]['hydrogen']['electrolyzer'])                  # [kg/s] actual production form electorlysis stack 
    
    fig, ax = plt.subplots(dpi=1000)
    ax.plot(demand_constant*(c.timestep*60), label = 'constant demand')
    # ax.plot(demand_variable, label = 'variable demand')
    ax.plot(production*(c.timestep*60), label = 'production')
    ax.grid(alpha = 0.3, zorder = 0)
    ax.set_ylabel('Hydrogen [kg]')
    ax.set_xlabel('Time [h]')
    xticks = list(np.linspace(0, simulation_steps - 1, 13).astype(int))
    # xticklabels = [str(value) for value in xticks]
    xticklabels = ['         Jan','         Feb','          Mar','         Apr','         May','          Jun','        Jul','          Aug','           Sep','          Oct','          Nov','           Dec','']
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)
    ax.legend()
    plt.title('Cumulative H$_\mathregular{2}$ production and demand')
    plt.show()
   
    print("\nAnnual produced hydrogen is equal to "+str(round((sum(balances[loc]['hydrogen']['electrolyzer']*c.timestep*60)/1e3),2))+" t/y ensuring a deliverable hydrogen mass flow rate equal to "+str(round((constantflow),2))+" kg/s")
    
def plot_post_process(balances_pp,studycase,loc,first_day,last_day):
    
    """
    Total balances figures and graphs
    ----------
    balances_pp : dictionary cointaining energy balances elaborated during post process
            
    studycase : dictionary (all the inputs are optional)
        'location_1_name': inputs required to create a location object (see Location.py)
        'location_2_name': 
            ...
        'location_n_name':
    
    loc : str - location_name
    
    first_day : int - first day to start the balances bar plot
    
    last_day : int - last day to finish the balances bar plot

    """     
    
    if 'wind' not in balances_pp[loc]['electricity'] and 'PV' not in balances_pp[loc]['electricity']:    # if no renewables are included in the location        
        print("No renewable energy sources are present in the considered case study\nRES plots not available\n")
        return
    else:
        balances = balances_pp[loc]['electricity']
        x = np.arange(first_day*24*60/c.timestep,(last_day+1)*24*60/c.timestep)
        x = x*(c.timestep/60) # ensuring values are presented on a hourly basis
        width=0.9
        system = dict(sorted(studycase[loc].items(), key=lambda item: item[1]['priority'])) # ordered by priority
        hourly_steps = 60//c.timestep # number of simulation steps considered in 1 hour - depending on input parameters. If simulaton is on hourly basis, hourly_steps=1
        
    if 'wind' in balances and 'PV' not in balances:
        for tech_name_plot in balances['tech from grid']:                    
            fig = plt.figure(dpi=300)
            from mpl_toolkits.axisartist.axislines import SubplotZero
            ax = SubplotZero(fig, 1, 1, 1)
            fig.add_subplot(ax)
            ax.axis["xzero"].set_visible(True)
            ax.axis["xzero"].label.set_visible(False)
            ax.axis["xzero"].major_ticklabels.set_visible(False)
            for n in ["bottom","top", "right"]:ax.axis[n].set_visible(True)
            ax.grid(axis='y', alpha = 0.5, zorder = -4)
            ax.bar(x, balances['tech windsc'][tech_name_plot][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps], width,  label='Wind self consumption', color='yellowgreen')
            ax.bar(x, balances['tech from grid'][tech_name_plot][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps], width, bottom=balances['tech windsc'][tech_name_plot][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps], label='From grid', color='tomato')
            plt.title(tech_name_plot+' days '+str(first_day)+'-'+str(last_day))
            plt.plot(x,-balances[tech_name_plot][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps],'k',label='load')     
            plt.legend(ncol=2, bbox_to_anchor = (1.01,-0.11))
            plt.ylabel("Power [kW] ")
            plt.xlabel( "Time  [h] ")
            if (last_day-first_day) <= 10: # in order to better visualize daily behaviour if short timespans are selected
                plt.xticks(list(range(int(x[0]),int(np.ceil(x[-1]))+2,24)), [str(x) for x in list(range(first_day*24, ((last_day+1)*24)+1,24))], rotation=45) 
            plt.show()
            
                
    elif 'PV' in balances and 'wind' not in balances:
        for tech_name_plot in balances['tech from grid']:      
            fig = plt.figure(dpi=300)
            from mpl_toolkits.axisartist.axislines import SubplotZero
            ax = SubplotZero(fig, 1, 1, 1)
            fig.add_subplot(ax)
            ax.axis["xzero"].set_visible(True)
            ax.axis["xzero"].label.set_visible(False)
            ax.axis["xzero"].major_ticklabels.set_visible(False)
            for n in ["bottom","top", "right"]: ax.axis[n].set_visible(True)
            ax.grid(axis='y', alpha = 0.5, zorder = -4)
            ax.bar(x, balances['tech pvsc'][tech_name_plot][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps], width,  label='PV self consumption', color='yellowgreen')
            ax.bar(x, balances['tech from grid'][tech_name_plot][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps], width, bottom=balances['tech pvsc'][tech_name_plot][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps], label='From grid', color='tomato')
            plt.title(tech_name_plot+' days '+str(first_day)+'-'+str(last_day))
            plt.plot(x,-balances[tech_name_plot][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps],'k',label='load')     
            plt.legend(ncol=2, bbox_to_anchor = (1.01,-0.11))
            plt.ylabel("Power [kW] ")
            plt.xlabel( "Time  [h] ")
            if (last_day-first_day) <= 10: # in order to better visualize daily behaviour if short timespans are selected
                plt.xticks(list(range(int(x[0]),int(np.ceil(x[-1]))+2,24)), [str(x) for x in list(range(first_day*24, ((last_day+1)*24)+1,24))], rotation=45) 
            plt.show()


    elif 'PV' in balances and 'wind' in balances:
        
        if list(system.keys()).index('PV') < list(system.keys()).index('wind'): # PV comes before wind
            for tech_name_plot in balances['tech from grid']:
                fig = plt.figure(dpi=300)
                from mpl_toolkits.axisartist.axislines import SubplotZero
                ax = SubplotZero(fig, 1, 1, 1)
                fig.add_subplot(ax)
                ax.axis["xzero"].set_visible(True)
                ax.axis["xzero"].label.set_visible(False)
                ax.axis["xzero"].major_ticklabels.set_visible(False)
                for n in ["bottom","top", "right"]:
                    ax.axis[n].set_visible(True)
                ax.grid(axis='y', alpha = 0.5, zorder = -4)
                ax.bar(x, balances['tech pvsc'][tech_name_plot][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps], width,  label='PV self consumption', color='yellowgreen')
                ax.bar(x, balances['tech windsc'][tech_name_plot][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps], width, bottom=balances['tech pvsc'][tech_name_plot][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps], label='wind self consumption', color='orange')
                ax.bar(x, balances['tech from grid'][tech_name_plot][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps], width, bottom=balances['tech pvsc'][tech_name_plot][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps]+balances['tech windsc'][tech_name_plot][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps], label='From grid', color='tomato')
                plt.title(tech_name_plot+' days '+str(first_day)+'-'+str(last_day))
                plt.plot(x,-balances[tech_name_plot][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps],'k',label='load')     
                plt.legend(ncol=2, bbox_to_anchor = (1.01,-0.11))
                plt.ylabel("Power [kW]")
                plt.xlabel( "Time  [h] ")
                if (last_day-first_day) <= 10: # in order to better visualize daily behaviour if short timespans are selected
                    plt.xticks(list(range(int(x[0]),int(np.ceil(x[-1]))+2,24)), [str(x) for x in list(range(first_day*24, ((last_day+1)*24)+1,24))], rotation=45) 
                plt.show()
        
        else:
            for tech_name_plot in balances['tech from grid']:
                fig = plt.figure(dpi=300)
                from mpl_toolkits.axisartist.axislines import SubplotZero
                ax = SubplotZero(fig, 1, 1, 1)
                fig.add_subplot(ax)
                ax.axis["xzero"].set_visible(True)
                ax.axis["xzero"].label.set_visible(False)
                ax.axis["xzero"].major_ticklabels.set_visible(False)
                for n in ["bottom","top", "right"]:
                    ax.axis[n].set_visible(True)
                ax.grid(axis='y', alpha = 0.5, zorder = -4)
                ax.bar(x, balances['tech windsc'][tech_name_plot][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps], width,  label='wind self consumption', color='orange')
                ax.bar(x, balances['tech pvsc'][tech_name_plot][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps], width, bottom=balances['tech windsc'][tech_name_plot][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps], label='PV self consumption', color='yellowgreen')
                ax.bar(x, balances['tech from grid'][tech_name_plot][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps], width, bottom=balances['tech pvsc'][tech_name_plot][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps]+balances['tech windsc'][tech_name_plot][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps], label='From grid', color='tomato')
                plt.title(tech_name_plot+' days '+str(first_day)+'-'+str(last_day))
                plt.plot(x,-balances[tech_name_plot][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps],'k',label='load')     
                plt.legend(ncol=2, bbox_to_anchor = (1.01,-0.11))
                plt.ylabel("Power [kW]")
                plt.xlabel( "Time  [h] ")
                if (last_day-first_day) <= 10: # in order to better visualize daily behaviour if short timespans are selected
                    plt.xticks(list(range(int(x[0]),int(np.ceil(x[-1]))+2,24)), [str(x) for x in list(range(first_day*24, ((last_day+1)*24)+1,24))], rotation=45) 
                plt.show()
                
            
    # Hourly energy balances   
    if 'PV' in balances:
        pv = balances['PV'][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps]
    
    if 'wind' in balances:
        wind = balances['wind'][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps]
        
    if 'grid' in balances:                
        into_grid = np.zeros(24*(last_day-first_day+1)*hourly_steps)
        from_grid = np.zeros(24*(last_day-first_day+1)*hourly_steps)
        for i,e in enumerate(balances['grid'][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps]):
            if e > 0:
                from_grid[i] = e
            else:
                into_grid[i] = e 
    
    if 'battery' in balances:
        charge_battery = np.zeros(24*(last_day-first_day+1)*hourly_steps)
        discharge_battery = np.zeros(24*(last_day-first_day+1)*hourly_steps)
        for i,e in enumerate(balances['battery'][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps]):
            if e > 0:
                discharge_battery[i] = e
            else:
                charge_battery[i] = e 
    
    if 'electrolyzer' in balances:
        ele = balances['electrolyzer'][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps]
    
    if 'mechanical compressor' in balances:
        compressor = balances['mechanical compressor'][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps]
    
    if 'fuel cell' in balances:
        fc = balances['fuel cell'][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps]
        
    x = np.arange(first_day*24*60/c.timestep,(last_day+1)*24*60/c.timestep)        
    x = x*(c.timestep/60) # ensuring values are presented on a hourly basis

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
    
    if 'demand' in balances:
        load = -balances['demand'][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps]+balances['hyd tot electricity'][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps]
        el_demand = balances['demand'][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps]
        ax.bar(x, el_demand, width, label='To el-devices', color='grey')
    else:
        load = balances['hyd tot electricity'][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps] 
        el_demand = np.zeros(24*(last_day-first_day+1)*hourly_steps)
        
    if 'PV' in balances and 'wind' in balances:
        if list(system.keys()).index('PV') < list(system.keys()).index('wind'): # PV comes before wind
            prioritised_res_1 = pv
            label_res_1 = 'PV'
            prioritised_res_2 = wind
            label_res_2 = 'Wind'
        else:
            prioritised_res_1 = wind
            label_res_1 = 'Wind'
            prioritised_res_2 = pv
            label_res_2 = 'PV'
 
        ax.bar(x, np.minimum(load,prioritised_res_1), width,  label= label_res_1+' self consumption', color='yellowgreen')
        ax.bar(x, np.maximum(prioritised_res_1-load,np.zeros(len(prioritised_res_1))), width, bottom=load, label= label_res_1+' surplus', color='cornflowerblue')
        ax.bar(x, np.maximum(np.minimum(prioritised_res_1+prioritised_res_2,load)-prioritised_res_1,np.zeros(len(prioritised_res_1))), width, bottom= prioritised_res_1, label= label_res_2+' self consumption', color='tab:green')
        ax.bar(x, np.maximum(prioritised_res_1-load+prioritised_res_2,np.zeros(len(prioritised_res_1)))-np.maximum(prioritised_res_1-load,np.zeros(len(prioritised_res_1))), width, bottom=load+np.maximum(prioritised_res_1-load,np.zeros(len(prioritised_res_1))),label= label_res_2+' surplus', color='tab:blue')
        
        if not 'battery' in balances and not 'electrolyzer' in balances:
            ax.bar(x, from_grid, width ,bottom=prioritised_res_1+prioritised_res_2, label='from grid', color='tomato')
            ax.bar(x, into_grid, width, bottom=el_demand, label='into grid', color='orange')
            
        if 'battery' in balances:
            ax.bar(x, from_grid, width, bottom=prioritised_res_1+prioritised_res_2+discharge_battery, label='from grid', color='tomato')
            ax.bar(x, discharge_battery, width, bottom = prioritised_res_1+prioritised_res_2, label='discharge battery',  color='purple')
            ax.bar(x, into_grid, width, bottom=charge_battery+el_demand , label='into grid', color='orange')
            ax.bar(x, charge_battery, width,  bottom=el_demand, label='charge battery',  color='violet')
 
        if 'electrolyzer' in balances and 'fuel cell' in balances:
            ax.bar(x, from_grid, width ,bottom=prioritised_res_1+prioritised_res_2+fc, label='from grid', color='tomato')
            ax.bar(x, fc, width, bottom = prioritised_res_1+prioritised_res_2, label='from fuel cell',  color='peru')
            if 'mechanical compressor' in balances:
                ax.bar(x, ele, width, bottom = el_demand, label='to electrolyzer',  color='chocolate')
                ax.bar(x, compressor, width, bottom = ele+el_demand, label='to compressor',  color='maroon')
                ax.bar(x, into_grid, width, bottom= ele+compressor +el_demand, label='into grid', color='orange')
            else:
                ax.bar(x, into_grid, width,bottom=ele+el_demand , label='into grid', color='orange')
                ax.bar(x, ele, width, bottom = el_demand,  label='to electrolyzer',  color='chocolate')
                                
        if 'electrolyzer' in balances and not 'fuel cell' in balances:
            ax.bar(x, from_grid, width ,bottom=prioritised_res_1+prioritised_res_2, label='from grid', color='tomato')
            if 'mechanical compressor' in balances:
                ax.bar(x, ele, width, bottom = el_demand, label='to electrolyzer',  color='chocolate')
                ax.bar(x, compressor, width, bottom = ele+el_demand, label='to compressor',  color='maroon')
                ax.bar(x, into_grid, width, bottom= ele+compressor+el_demand , label='into grid', color='orange')
            else:
                ax.bar(x, into_grid, width,bottom=ele+el_demand , label='into grid', color='orange')
                ax.bar(x, ele, width, bottom = el_demand, label='to electrolyzer',  color='chocolate')
                
    elif 'PV' in balances and not 'wind' in balances:
        ax.bar(x, np.minimum(load,pv), width,  label='PV self consumption', color='yellowgreen')
        ax.bar(x, np.maximum(pv-load,np.zeros(len(pv))), width, bottom=load, label='PV surplus', color='cornflowerblue')
        
        if not 'battery' in balances and not 'electrolyzer' in balances:
            ax.bar(x, from_grid, width ,bottom=pv, label='from grid', color='tomato')
            ax.bar(x, into_grid, width, bottom=el_demand, label='into grid', color='orange')
            
        if 'battery' in balances:
            ax.bar(x, from_grid, width, bottom=pv+discharge_battery, label='from grid', color='tomato')
            ax.bar(x, discharge_battery, width, bottom = pv, label='discharge battery',  color='purple')
            ax.bar(x, into_grid, width, bottom=charge_battery+el_demand , label='into grid', color='orange')
            ax.bar(x, charge_battery, width, bottom=el_demand, label='charge battery',  color='violet')
       
        if 'electrolyzer' in balances and 'fuel cell' in balances:
            ax.bar(x, from_grid, width ,bottom=pv+fc, label='from grid', color='tomato')
            ax.bar(x, fc, width, bottom = pv, label='from fuel cell',  color='peru')
            if 'mechanical compressor' in balances:
                ax.bar(x, ele, width, bottom = el_demand, label='to electrolyzer',  color='chocolate')
                ax.bar(x, compressor, width, bottom = ele+el_demand, label='to compressor',  color='maroon')
                ax.bar(x, into_grid, width, bottom= ele+compressor+el_demand , label='into grid', color='orange')
            else:
                ax.bar(x, into_grid, width,bottom=ele+el_demand , label='into grid', color='orange')
                ax.bar(x, ele, width, bottom=el_demand, label='to electrolyzer',  color='chocolate')
                                
        if 'electrolyzer' in balances and not 'fuel cell' in balances:
            ax.bar(x, from_grid, width ,bottom=pv, label='from grid', color='tomato')
            if 'mechanical compressor' in balances:
                ax.bar(x, ele, width, bottom = el_demand,  label='to electrolyzer',  color='chocolate')
                ax.bar(x, compressor, width, bottom = ele+el_demand, label='to compressor',  color='maroon')
                ax.bar(x, into_grid, width, bottom= ele+compressor+el_demand , label='into grid', color='orange')
            else:
                ax.bar(x, into_grid, width,bottom=ele+el_demand, label='into grid', color='orange')
                ax.bar(x, ele, width, bottom=el_demand, label='to electrolyzer',  color='chocolate')
       
    elif 'wind' in balances and not 'PV' in balances:
        ax.bar(x, np.minimum(load,wind), width,  label='Wind self consumption', color='yellowgreen')
        ax.bar(x, np.maximum(wind-load,np.zeros(len(wind))), width, bottom=load, label='Wind surplus', color='cornflowerblue')
        
        if not 'battery' in balances and not 'electrolyzer' in balances:
            ax.bar(x, from_grid, width ,bottom=wind, label='from grid', color='tomato')
            ax.bar(x, into_grid, width,  bottom=el_demand,label='into grid', color='orange')
            
        if 'battery' in balances:
            ax.bar(x, from_grid, width, bottom=wind+discharge_battery, label='from grid', color='tomato')
            ax.bar(x, discharge_battery, width, bottom = wind, label='discharge battery',  color='purple')
            ax.bar(x, into_grid, width, bottom=charge_battery+el_demand , label='into grid', color='orange')
            ax.bar(x, charge_battery, width, bottom=el_demand, label='charge battery',  color='violet')
        
        if 'electrolyzer' in balances and 'fuel cell' in balances:
            ax.bar(x, from_grid, width ,bottom=wind+fc, label='from grid', color='tomato')
            ax.bar(x, fc, width, bottom = wind, label='from fuel cell',  color='peru')
            if 'mechanical compressor' in balances:
                ax.bar(x, ele, width, bottom = el_demand,  label='to electrolyzer',  color='chocolate')
                ax.bar(x, compressor, width, bottom = ele+el_demand, label='to compressor',  color='maroon')
                ax.bar(x, into_grid, width, bottom= ele+compressor+el_demand , label='into grid', color='orange')
            else:
                ax.bar(x, into_grid, width,bottom=ele+el_demand , label='into grid', color='orange')
                ax.bar(x, ele, width, bottom=el_demand, label='to electrolyzer',  color='chocolate')
                                
        if 'electrolyzer' in balances and not 'fuel cell' in balances:
            ax.bar(x, from_grid, width ,bottom=wind, label='from grid', color='tomato')
            if 'mechanical compressor' in balances:
                ax.bar(x, ele, width,  bottom = el_demand,  label='to electrolyzer',  color='chocolate')
                ax.bar(x, compressor, width, bottom = ele+el_demand, label='to compressor',  color='maroon')
                ax.bar(x, into_grid, width, bottom= ele+compressor+el_demand , label='into grid', color='orange')
            else:
                ax.bar(x, into_grid, width,bottom=ele+el_demand , label='into grid', color='orange')
                ax.bar(x, ele, width, bottom = el_demand, label='to electrolyzer',  color='chocolate')
       
       
    plt.title(loc+' days '+str(first_day)+'-'+str(last_day))
    plt.plot(x,load,'k',label='load')     
    plt.legend(ncol=2, bbox_to_anchor = (1.01,-0.11))
    plt.ylabel("Power [kW]")
    plt.xlabel("Time  [h]")
    if (last_day-first_day) <= 10: # in order to better visualize daily behaviour if short timespans are selected
        # plt.xticks(list(range(first_day*24*hourly_steps, ((last_day+1)*24*hourly_steps)+1,24*hourly_steps)), [str(x) for x in list(range(first_day*24, ((last_day+1)*24)+1,24))], rotation=45) 
        plt.xticks(list(range(int(x[0]),int(np.ceil(x[-1]))+2,24)), [str(x) for x in list(range(first_day*24, ((last_day+1)*24)+1,24))], rotation=45) 
    # ax.xaxis.set_tick_params(bottom=True,labelbottom=True)
    #plt.xticks([0,6,12,18,24],['0','6','12','18','24'],fontsize=10,color='g')
    #plt.xticks([0,6,12,18,24,30,36,42,48],['0','6','12','18','24','30','36','42','48'],fontsize=10,color='g')
    #plt.yticks([-2,-1,0,1,2],['-2','-1','0','1','2'])
    plt.show() 