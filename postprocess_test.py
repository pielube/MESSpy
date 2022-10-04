"""
MESSpy - postprocessing
"""
import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def total_balances(simulation_name,loc):

    with open('Results/balances_'+simulation_name+'.pkl', 'rb') as f:
        balances = pickle.load(f)
    
    ###### load analysys
    
    ##### total electricity balances
    carriers = ['electricity']
    units = {'electricity': 'kWh', 'hydrogen': 'kg'}
    
    for carrier in carriers:
        print('\nTotal '+carrier+' balances:\n')
        print('\n'+loc)   
        for b in balances[loc][carrier]:
            
            positiv=balances[loc][carrier][b][balances[loc][carrier][b]>0].sum()
            negativ=balances[loc][carrier][b][balances[loc][carrier][b]<0].sum()
            
            if positiv != 0:
                print(b+' '+str(round(positiv,1))+' '+units[carrier])
                
            if negativ != 0:
                print(b+' '+str(round(negativ,1))+' '+units[carrier])
    print('\n')
   
    
def NPV_plot(study_case):
    ##### economic
    simulation_name = 'economic_assessment_'
    with open('Results/'+simulation_name+study_case+'.pkl', 'rb') as f:
        economic = pickle.load(f)
    
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
    
    with open('Results/balances_'+simulation_name+'.pkl', 'rb') as f:
        balances = pickle.load(f)
        
    df = pd.DataFrame(0.00,columns=["Value [kWh]","Value / production [%]","Value / demand [%]"],
                      index=["SC","CSC","Into MV grid","From MV grid","Production","Demand"])
    
    df.loc['CSC']['Value [kWh]'] = sum(balances['REC']['electricity']['collective self consumption'])
    df.loc['Into MV grid']['Value [kWh]'] = -sum(balances['REC']['electricity']['into grid'])  -df.loc['CSC']['Value [kWh]']
    df.loc['From MV grid']['Value [kWh]'] = sum(balances['REC']['electricity']['from grid']) -df.loc['CSC']['Value [kWh]']
    
    for loc in balances:
        if loc != 'REC':   
            balance = balances[loc]['electricity']
            df.loc['Demand']['Value [kWh]'] += -sum(balance['demand'])
            if 'PV' in balance:
                df.loc['Production']['Value [kWh]'] += sum(balance['PV'])
                
    df.loc['SC']['Value [kWh]'] = df.loc['Demand']['Value [kWh]'] - df.loc['From MV grid']['Value [kWh]'] - df.loc['CSC']['Value [kWh]']
    
    for b in df.index:
        df.loc[b]['Value / production [%]'] = df.loc[b]['Value [kWh]'] / df.loc['Production']['Value [kWh]'] *100
        df.loc[b]['Value / demand [%]'] = df.loc[b]['Value [kWh]'] / df.loc['Demand']['Value [kWh]'] *100
        
    print('\n REC electricity balances: '+simulation_name+'\n') 
    print(df.round(decimals=2))
    return(df.round(decimals=2))
    
def LOC_plot(simulation_name):
      
    with open('Results/LOC_'+simulation_name+'.pkl', 'rb') as f:
        LOC = pickle.load(f)
           
    unit = {'H tank': '[kg]', 'battery': '[kWh]'}
    
    for location_name in LOC:
        for tech in LOC[location_name]:
            
            y = LOC[location_name][tech]
            x = np.linspace(0,len(y)-1,len(y))        
            plt.plot(x,y,label=location_name)
            plt.grid()
            plt.ylabel('LOC '+unit[tech])
            plt.xlabel('Time [hours]')
            plt.title(location_name+' '+tech)
            plt.show()
            
def storage_control(simulation_name,e_cost=0.30,H_cost=0.05):
    with open('Results/LOC_'+simulation_name+'.pkl', 'rb') as f:
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
            
    
def hourly_balances(simulation_name,location_name,first_day,last_day,carrier='electricity',width=0.9,collective=0):
    
        with open('Results/balances_'+simulation_name+'.pkl', 'rb') as f:
            balances = pickle.load(f)
            
        balances = balances[location_name][carrier]
        
        if 'demand' in balances:
            load = -balances['demand'][first_day*24:last_day*24+24]
        
        if 'PV' in balances:
            pv = balances['PV'][first_day*24:last_day*24+24]
            
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
                    
        if 'electrolyzer' in balances and 'fuel cell' in balances:
            fc = balances['fuel cell'][first_day*24:last_day*24+24]
            ele = balances['electrolyzer'][first_day*24:last_day*24+24]
                     
        to_csc = np.zeros(24*(last_day-first_day+1))
        from_csc = np.zeros(24*(last_day-first_day+1))
        for i,e in enumerate(balances['collective self consumption'][first_day*24:last_day*24+24]):
            if e > 0:
                from_csc[i] = e
            else:
                to_csc[i] = e 
            
        x = np.linspace(first_day*24+1,(last_day+1)*24,(last_day-first_day+1)*24)
        x = np.arange((last_day-first_day+1)*24)
        
        fig = plt.figure(dpi=1000)
        from mpl_toolkits.axisartist.axislines import SubplotZero
        ax = SubplotZero(fig, 1, 1, 1)
        fig.add_subplot(ax)
        ax.axis["xzero"].set_visible(True)
        for n in ["bottom", "top", "right"]:
            ax.axis[n].set_visible(False)
        ax.grid(axis='y')
        
        if 'PV' in balances:
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
            
           # plt.ylim(-4,4)
            
        else:
            ax.bar(x, from_grid-from_csc, width ,bottom=from_csc, label='from grid', color='tomato')
            ax.bar(x, np.array(from_csc), width, label='collective self consumption',  color='gold')
           # plt.ylim(0,3.3)
        
        plt.title(location_name+' days '+str(first_day)+'-'+str(last_day))
        plt.plot(x,load,'k',label='load')     
        plt.legend(ncol=2, bbox_to_anchor=(1.2, 0))
        plt.ylabel("Hourly energy [kWh/h] ")
        #plt.xticks([0,6,12,18,24],['0','6','12','18','24'],fontsize=10,color='g')
        #plt.xticks([0,6,12,18,24,30,36,42,48],['0','6','12','18','24','30','36','42','48'],fontsize=10,color='g')
        #plt.yticks([-2,-1,0,1,2],['-2','-1','0','1','2'])
        plt.show()
        
            
def csc_allocation_sum(simulation_name):
    with open('Results/balances_'+simulation_name+'.pkl', 'rb') as f:
        balances = pickle.load(f)
        
    for location_name in balances:
        csc = balances[location_name]['electricity']['collective self consumption']
        from_csc = csc.sum(where=csc>0)
        to_csc = csc.sum(where=csc<0)
    
        print(f"{location_name} {int(from_csc)}  {int(to_csc)}")
   
    
   
        
