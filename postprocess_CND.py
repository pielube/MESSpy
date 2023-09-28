"""
MESSpy - postprocessing
"""
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.io as pio
import plotly.graph_objects as go
from matplotlib.patches import Patch
from mpl_toolkits.axisartist.axislines import SubplotZero

   
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
    
    carriers = ['electricity','heating water','gas','hydrogen']
    units    = {'electricity': 'kWh', 'hydrogen': 'kg', 'gas':'kWh', 'heating water':'kWh'}
    
    if var: 
        
        carriers = [var]
        units    = {var : units[var]}

    for carrier in carriers:
        print('\nTotal '+carrier+' balances '+loc+':\n') 
        for b in balances[loc][carrier]:
            
            positiv=round(balances[loc][carrier][b][balances[loc][carrier][b]>0].sum(),1)
            negativ=round(balances[loc][carrier][b][balances[loc][carrier][b]<0].sum(),1)
            
            if positiv != 0:
                print(b+' '+str(round(positiv,1))+' '+units[carrier])
                
            if negativ != 0:
                print(b+' '+str(round(negativ,1))+' '+units[carrier])
                
    print('\n')
   
    
def NPV_scenarios(scenari,locations,title,devided=1):
    
    # scenari: list of simulation names
    # locations# list of locations
    
    plt.figure(dpi=1000)
    
    for simulation_name in scenari:
        with open('Results/economic_assessment_'+simulation_name+'.pkl', 'rb') as f:
            economic = pickle.load(f)
        
        y = np.zeros(31)
        for loc in locations:
        
            y += economic[loc]['NPV']/devided
        x = np.linspace(0,len(y)-1,len(y))
          #  plt.plot(x,y,label=loc+' '+simulation_name)
        plt.plot(x,y,label=simulation_name)
           
    plt.plot(x,np.zeros(31),color='k')
    plt.legend()
    plt.title(title)
    plt.grid()
    plt.ylabel('Net Present Value [€]')
    plt.xlabel('Time [years]')
    plt.ylabel('Valore attuale netto [€]')
    plt.xlabel('Tempo [anni]')
    #plt.xlim(0,len(y)-1)
    plt.xlim(0,30)
    #plt.ylim(-8000,18000)
    # plt.annotate(simulation_name,(15,-7800))
    plt.show()
    
def REC_electricity_balance(simulation_name, noprint=False, mounth=False):
    
    with open('Results/balances_'+simulation_name+'.pkl', 'rb') as f:
        balances = pickle.load(f)
        
    df = pd.DataFrame(0.00,columns=["Value [kWh]","Value / production [%]","Value / demand [%]"],
                      index=["SC","CSC","Into MV grid","From MV grid","Production","Demand"])
    
    csc = balances['REC']['electricity']['collective self consumption']
    into_grid = balances['REC']['electricity']['into grid']
    from_grid = balances['REC']['electricity']['from grid']
    
    dmc = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365] # duration of months: cumulate [days]             
    if mounth: # 1-12
        csc = csc [ dmc[mounth-1]*24 : dmc[mounth]*24]
        into_grid = into_grid [ dmc[mounth-1]*24 : dmc[mounth]*24]
        from_grid = from_grid [ dmc[mounth-1]*24 : dmc[mounth]*24]
        
    df.loc['CSC']['Value [kWh]'] = sum(csc)
    df.loc['Into MV grid']['Value [kWh]'] = -sum(into_grid)  -df.loc['CSC']['Value [kWh]']
    df.loc['From MV grid']['Value [kWh]'] = sum(from_grid) -df.loc['CSC']['Value [kWh]']
    
    for loc in balances:
        if loc != 'REC':   
            balance = balances[loc]['electricity']
            if 'demand' in balance:
                demand = balance['demand']
                if mounth: # 1-12
                    demand = demand [ dmc[mounth-1]*24 : dmc[mounth]*24]
                df.loc['Demand']['Value [kWh]'] += -sum(demand)
            if 'PV' in balance:
                pv = balance['PV']
                if mounth: # 1-12
                    pv = pv [ dmc[mounth-1]*24 : dmc[mounth]*24]
                df.loc['Production']['Value [kWh]'] += sum(pv)
                
    df.loc['SC']['Value [kWh]'] = df.loc['Demand']['Value [kWh]'] - df.loc['From MV grid']['Value [kWh]'] - df.loc['CSC']['Value [kWh]']
    
    for b in df.index:
        df.loc[b]['Value / production [%]'] = df.loc[b]['Value [kWh]'] / df.loc['Production']['Value [kWh]'] *100
        df.loc[b]['Value / demand [%]'] = df.loc[b]['Value [kWh]'] / df.loc['Demand']['Value [kWh]'] *100
      
    if not noprint:
        if mounth:
            print(f"mounth: {mounth}")
        print('\n REC electricity balances: '+simulation_name+'\n') 
        print(df.round(decimals=2))
    return(df.astype(int))
        
def hist_12_balances_pc(simulation_name,ymax):
    
    sc = []
    csc = []
    into_mv_grid = []
    from_mv_grid = []
    for m in np.arange(12):
        balances = REC_electricity_balance(simulation_name, noprint=True, mounth=m+1)
        sc.append(balances['Value [kWh]']['SC'])
        csc.append(balances['Value [kWh]']['CSC'])
        into_mv_grid.append(balances['Value [kWh]']['Into MV grid'])
        from_mv_grid.append(balances['Value [kWh]']['From MV grid'])
    sc = np.array(sc)
    csc = np.array(csc)
    x = np.arange(12)  # the label locations
    width = 0.8  # the width of the bars
    

    # Create a figure with two subplots arranged horizontally (1 row, 2 columns)
    fig, axes = plt.subplots(1, 2, dpi=1000, figsize=(12, 4))
    
    ax1 = axes[0]
    ax2 = axes[1]
    
    # Plot the data for the first subplot
    ax1.bar(x, sc, width, label='SC', color='yellowgreen')
    ax1.bar(x, csc, width, bottom=sc, label='CSC', color='gold')
    ax1.bar(x, into_mv_grid, width, bottom=sc + csc, label='Into grid', color='tab:blue')
    ax1.set_ylabel('Energia prodotta [kWh/mese]')
    ax1.set_xticks(np.arange(12))
    ax1.set_xticklabels(['gen', 'feb', 'mar', 'apr', 'mag', 'giu', 'lug', 'ago', 'set', 'ott', 'nov', 'dic'])
    colors = ["tab:blue", "yellowgreen"]
    labels = ["Immissione", "Autoconsumo"]
    leg = [Patch(facecolor=c) for c in colors]
    ax1.legend(leg, labels, ncol=1)  # Create the legend within the first subplot
    ax1.set_ylim(0, ymax)
    ax1.set_title('Produzione '+simulation_name)
    ax1.grid(axis='y')
    
    # Plot the data for the second subplot
    ax2.bar(x, sc, width, label='SC', color='yellowgreen')
    ax2.bar(x, csc, width, bottom=sc, label='CSC', color='gold')
    ax2.bar(x, from_mv_grid, width, bottom=sc + csc, label='From grid', color='tomato')
    ax2.set_ylabel('Energia consumata [kWh/mese]')
    ax2.set_xticks(np.arange(12))
    ax2.set_xticklabels(['gen', 'feb', 'mar', 'apr', 'mag', 'giu', 'lug', 'ago', 'set', 'ott', 'nov', 'dic'])
    colors = ["tomato", "yellowgreen"]
    labels = ["Prelievo", "Autosufficienza"]
    leg = [Patch(facecolor=c) for c in colors]
    ax2.legend(leg, labels, ncol=1)  # Create the legend within the second subplot
    ax2.set_ylim(0, ymax)
    ax2.set_title('Consumi '+simulation_name)
    ax2.grid(axis='y')
    
    # Adjust the layout to prevent overlapping of titles and labels
    plt.tight_layout()
    
    plt.show()
    
     
def csc_allocation_sum(simulation_name):
    with open('Results/balances_'+simulation_name+'.pkl', 'rb') as f:
        balances = pickle.load(f)
        
    for location_name in balances:
        csc = balances[location_name]['electricity']['collective self consumption']
        from_csc = csc.sum(where=csc>0)
        to_csc = csc.sum(where=csc<0)
    
        print(f"{location_name} {int(from_csc)}  {int(to_csc)}")
   

def hist12(simulation_name):
    
    demand = pd.DataFrame(np.tile(np.zeros(3),(12,1)),index=np.arange(12),columns=(1,2,3)) # dataframe 12x3
    
    time_slots = pd.read_csv(r'input_test\energy_price\time_slots.csv')['0']
    months = pd.read_csv(r'input_test\energy_price\months.csv')['0']
        
    with open('results/balances_'+simulation_name+'.pkl', 'rb') as f: balances = pickle.load(f)
    
    for h in time_slots.index:
        
        ts = time_slots[h]
        m = months[h]
        
        for loc in balances:
            
            if 'demand' in balances[loc]['electricity']:
                if loc != 'REC':
                    demand[ts][m] += - balances[loc]['electricity']['demand'][h]
                    
    fig, ax = plt.subplots(dpi=1000)
    width = 0.8
    x = [0,1,2,3,4,5,6,7,8,9,10,11]
    ax.bar(x, demand[1], width,  label='F1', edgecolor='k')
    ax.bar(x, demand[2], width, bottom=demand[1], label='F2', edgecolor='k')
    ax.bar(x, demand[3], width, bottom=demand[1]+demand[2], label='F3', edgecolor='k')
    
    plt.legend()
    
    ax.set_ylabel('Energia elettrica [kWh/mese]')
    
    #plt.xticks([0,1,2,3,4,5,6,7,8,9,10,11],['jan','feb','mar','apr','mar','jun','jul','ago','sep','opt','nov','dic'])
    plt.xticks([0,1,2,3,4,5,6,7,8,9,10,11],['gen','feb','mar','apr','mag','giu','lug','ago','set','ott','nov','dic'])

    plt.show()
    

def typedays(name_studycase,days,titles,ymax):
    
    with open('results/balances_'+name_studycase+'.pkl', 'rb') as f: balances = pickle.load(f)
    load = np.zeros(8760)
    production = np.zeros(8760)
    csc = np.zeros(8760)
    sc = np.zeros(8760)
    for loc in balances:
        if loc != 'REC':
            if 'demand' in balances[loc]['electricity']:
                load += - balances[loc]['electricity']['demand']   
            if 'collective self consumption' in balances[loc]['electricity']:
                csc += np.clip(balances[loc]['electricity']['collective self consumption'],0,None)
            if 'PV' in balances[loc]['electricity']:
                production += balances[loc]['electricity']['PV']
                if 'demand' in balances[loc]['electricity']: 
                    from_grid = np.clip(balances[loc]['electricity']['grid'],0,None)
                    sc +=  - balances[loc]['electricity']['demand'] - from_grid

    fig, axes = plt.subplots(2, 2, dpi=1000, figsize=(12, 8))
    x = np.arange(25)
    width = 0.8
    ax0 = axes[0][0]
    ax1 = axes[1][0]
    ax2 = axes[0][1]
    ax3 = axes[1][1]
    
    ax_dict = {0:ax0, 1:ax1, 2:ax2, 3:ax3}
    
    for d in range(4):
        
        load_d = load[days[d]*24:days[d]*24+25]
        production_d = production[days[d]*24:days[d]*24+25]
        csc_d = csc[days[d]*24:days[d]*24+25]
        sc_d = sc[days[d]*24:days[d]*24+25]
        
        ax_dict[d].bar(x, production_d, width,  label='Immessa', color='tab:blue')
        ax_dict[d].bar(x, sc_d, width,  label='Autoconsumo', color='yellowgreen')
        #ax_dict[d].bar(x, csc_d, width, bottom=scd,  label='Autoconsumo collettivo', color='gold')
        ax_dict[d].plot(x,load_d,'k',label='Consumo') 
        ax_dict[d].plot(x,production_d,'tab:blue',label='Produzione') 
        
        ax_dict[d].grid(axis='y', alpha = 0.5)
        ax_dict[d].set_xlim(0,23)
        ax_dict[d].set_ylim(0,ymax)
       # ax_dict[d].set_xticks([0,6,12,18,24],['0','6','12','18','24'])
        ax_dict[d].legend()
        ax_dict[d].set_xlabel('ora')
        ax_dict[d].set_ylabel('kWh')
        ax_dict[d].set_title(titles[d])

    plt.tight_layout()
    plt.show()
    
    
def bilanci_scenari(scenari):
    df = pd.DataFrame(index=['Consumo','Produzione','Prelievo','Immissione','Autoconsumo','Indice di autocosnumo','Indice di autosufficienza'],columns=scenari)
    for name_studycase in scenari:
        a = REC_electricity_balance(name_studycase,noprint=True)
        df.loc['Consumo',name_studycase] = a.loc['Demand','Value [kWh]']
        df.loc['Produzione',name_studycase] = a.loc['Production','Value [kWh]']
        df.loc['Immissione',name_studycase] = a.loc['Into MV grid','Value [kWh]']
        df.loc['Prelievo',name_studycase] = a.loc['From MV grid','Value [kWh]']
        df.loc['Autoconsumo',name_studycase] = a.loc['SC','Value [kWh]']
        df.loc['Indice di autocosnumo',name_studycase] = a.loc['SC','Value / production [%]']
        df.loc['Indice di autosufficienza',name_studycase] = a.loc['SC','Value / demand [%]']   
    return(df)

def flussi_di_cassa_scenari(scenari):
    scenari0 = ['Adesso'] + scenari
    df = pd.DataFrame(index=['Bolletta elettrica','Vendita energia','Totale','Risparmio vs adesso'],columns=scenari0)
    with open('results/economic_assessment_'+scenari[0]+'.pkl', 'rb') as f: economic = pickle.load(f)
    df.loc['Bolletta elettrica','Adesso'] = economic['GEC']['CF_refcase']['Purchase']['electricity'][0]
    df.loc['Vendita energia','Adesso'] = economic['GEC']['CF_refcase']['Sale']['electricity'][0]
    df.loc['Totale','Adesso'] = df.loc['Bolletta elettrica','Adesso'] + df.loc['Vendita energia','Adesso']
    df.loc['Risparmio vs adesso','Adesso'] = 0
    for name_economic in scenari:
        with open('results/economic_assessment_'+name_economic+'.pkl', 'rb') as f: economic = pickle.load(f)
        df.loc['Bolletta elettrica',name_economic] = economic['GEC']['CF_studycase']['Purchase']['electricity'][0]
        df.loc['Vendita energia',name_economic] = economic['GEC']['CF_studycase']['Sale']['electricity'][0]
        df.loc['Totale',name_economic] = df.loc['Bolletta elettrica',name_economic] + df.loc['Vendita energia',name_economic]
        df.loc['Risparmio vs adesso',name_economic] = df.loc['Totale',name_economic] - df.loc['Totale','Adesso']
    df.astype(int)
    return(df)

def flussi_di_cassa_scenari_con_gas(scenari):
    scenari0 = ['Adesso'] + scenari
    df = pd.DataFrame(index=['Bolletta elettrica','Bolletta GPL','Vendita energia','Totale','Risparmio vs adesso'],columns=scenari0)
    with open('results/economic_assessment_'+scenari[0]+'.pkl', 'rb') as f: economic = pickle.load(f)
    df.loc['Bolletta elettrica','Adesso'] = economic['GEC']['CF_refcase']['Purchase']['electricity'][0]
    df.loc['Bolletta GPL','Adesso'] = economic['GEC']['CF_refcase']['Purchase']['gas'][0]
    df.loc['Vendita energia','Adesso'] = economic['GEC']['CF_refcase']['Sale']['electricity'][0]
    df.loc['Totale','Adesso'] = df.loc['Bolletta elettrica','Adesso'] + df.loc['Vendita energia','Adesso'] + df.loc['Bolletta GPL','Adesso']
    df.loc['Risparmio vs adesso','Adesso'] = 0
    for name_economic in scenari:
        with open('results/economic_assessment_'+name_economic+'.pkl', 'rb') as f: economic = pickle.load(f)
        df.loc['Bolletta elettrica',name_economic] = economic['GEC']['CF_studycase']['Purchase']['electricity'][0]
        df.loc['Bolletta GPL',name_economic] = economic['GEC']['CF_studycase']['Purchase']['gas'][0]
        df.loc['Vendita energia',name_economic] = economic['GEC']['CF_studycase']['Sale']['electricity'][0]
        df.loc['Totale',name_economic] = df.loc['Bolletta elettrica',name_economic] + df.loc['Vendita energia',name_economic] + df.loc['Bolletta GPL',name_economic]
        df.loc['Risparmio vs adesso',name_economic] = df.loc['Totale',name_economic] - df.loc['Totale','Adesso']
    df.astype(int)
    return(df)
    
     