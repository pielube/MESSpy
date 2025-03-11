"""
MESSpy - postprocessing
"""
import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from core import constants as c
import matplotlib.patches as mpatches
from matplotlib import cm                         

        
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
    into_grid = -balances['REC']['electricity']['into electricity grid']*c.P2E/c.kWh2kJ
    from_grid = balances['REC']['electricity']['from electricity grid']*c.P2E/c.kWh2kJ
    
    
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
            if 'electricity demand' in balance:
                demand = balance['electricity demand']
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
                

def RES_plot(simulation_name,location_name):
      
    with open('results/pkl/balances_'+simulation_name+'.pkl', 'rb')  as f: balances  = pickle.load(f) 
    balances_RES = balances[location_name]['electricity']

    for tech in balances_RES:
        if tech in ['PV','wind']:
            plt.figure(dpi=600)                               
            y = balances_RES[tech]
            x = np.linspace(0,len(y),len(y))        
            plt.plot(x,y,label=location_name)
            plt.grid()
            plt.ylabel('RES production [kW]')
            if c.simulation_years == 1 and c.timestep == 60:
                xticks = list(np.linspace(0, len(x) - 1, 13).astype(int))
                xticklabels = ['          Jan','         Feb','          Mar','         Apr','         May','          Jun','         Jul','          Aug','           Sep','          Oct','          Nov','           Dec','']
                plt.xticks(xticks,xticklabels,rotation=45)
                plt.xlabel('Time [hours]')
            elif c.simulation_years != 1 and c.timestep == 60:
                plt.xlabel('Time [hours]')
            else:
                plt.xlabel('Timestep')
            plt.title(location_name+' '+tech)
            plt.xlim(0,x[-1])
            plt.show()

def demand_plot(simulation_name,location_name):
      
    with open('results/pkl/balances_'+simulation_name+'.pkl', 'rb')  as f: balances  = pickle.load(f) 

    UM = {
        'electricity': 'kW',
        'heating water': 'kW',
        'cooling water': 'kW',
        'process heat': 'kW',
        'process hot water': 'kW',
        'process cold water': 'kW',
        'process chilled water': 'kW',
        'hydrogen': 'kg/s',
        'LP hydrogen': 'kg/s',
        'HP hydrogen': 'kg/s',
        'oxygen': 'kg/s',
        'process steam': 'kg/s',
        'gas': 'Sm^3/s',
        'water': 'm^3/s'
    }

    for carrier in balances[location_name]:
        if carrier + ' demand' in balances[location_name][carrier]:
            plt.figure(dpi=600)                               
            y = - balances[location_name][carrier][carrier + ' demand']
            x = np.linspace(0,len(y),len(y))        
            plt.plot(x,y,label=location_name)
            plt.grid()
            plt.ylabel(f"{carrier} demand [{UM[carrier]}]")
            plt.xlabel('Time [hours]')
            if c.simulation_years == 1 and c.timestep == 60:
                xticks = list(np.linspace(0, len(x) - 1, 13).astype(int))
                xticklabels = ['          Jan','         Feb','          Mar','         Apr','         May','          Jun','         Jul','          Aug','           Sep','          Oct','          Nov','           Dec','']
                plt.xticks(xticks,xticklabels,rotation=45)
                plt.xlabel('Time [hours]')
            elif c.simulation_years != 1 and c.timestep == 60:
                plt.xlabel('Time [hours]')
            else:
                plt.xlabel('Timestep')
            plt.title(location_name+' '+carrier+' demand')
            plt.xlim(0,x[-1])
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
            if c.simulation_years == 1 and c.timestep == 60:
                xticks = list(np.linspace(0, len(x) - 1, 13).astype(int))
                xticklabels = ['          Jan','         Feb','          Mar','         Apr','         May','          Jun','         Jul','          Aug','           Sep','          Oct','          Nov','           Dec','']
                plt.xticks(xticks,xticklabels,rotation=45)
                plt.xlabel('Time [hours]')
            elif c.simulation_years != 1 and c.timestep == 60:
                plt.xlabel('Time [hours]')
            else:
                plt.xlabel('Timestep')                          
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


def hydrogen_chain_curves(studycase,simulation_name,loc,print_=False,plot=False):
    
    """
    Plot hydrogen supply chain cumulative electricity consumption curves and hydrogen production figures and plots
    ----------
    """  
    
    # Load the production and consumption data
    with open('results/pkl/production_'+simulation_name+'.pkl', 'rb') as f: production  = pickle.load(f)
    with open('results/pkl/consumption_'+simulation_name+'.pkl', 'rb') as f: consumption  = pickle.load(f)
    with open('results/pkl/balances_'+simulation_name+'.pkl', 'rb') as f: balances = pickle.load(f)
    
    # ---------- PART 1: Hydrogen Production Plot ----------
    
    # Hydrogen production details (previously hydrogen_production function)
    simulation_steps    = len(balances[loc]['hydrogen']['electrolyzer'])
    constantflow        = sum(balances[loc]['hydrogen']['electrolyzer']) / simulation_steps       # [kg/s] constant flow rate deliverable by the system
    demand_constant     = np.array([constantflow * (i + 1) for i in range(simulation_steps)])     # [kg/s] fictitious constant demand 
    production          = np.cumsum(balances[loc]['hydrogen']['electrolyzer'])                    # [kg/s] actual production from electrolyzer stack 
    
    if plot:  # Only plot if plot=True
        fig, ax = plt.subplots(dpi=1000)
        ax.plot(demand_constant * (c.timestep * 60), label='constant demand')
        ax.plot(production * (c.timestep * 60), label='production')
        ax.grid(alpha=0.3, zorder=0)
        ax.set_ylabel('Hydrogen [kg]')
        if c.simulation_years == 1 and c.timestep == 60:
            ax.set_xlabel('Time [hours]')
            xticks = list(np.linspace(0, simulation_steps - 1, 13).astype(int))
            xticklabels = ['         Jan', '         Feb', '          Mar', '         Apr', '         May', '          Jun',
                           '        Jul', '          Aug', '           Sep', '          Oct', '          Nov', '           Dec', '']
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticklabels)
        elif c.simulation_years != 1 and c.timestep == 60:
            ax.set_xlabel('Time [hours]')
        else:
            ax.set_xlabel('Timestep')
        ax.legend()
        plt.title('Cumulative H$_\mathregular{2}$ production and demand')
        plt.show()
    
    if print_:  # Only print if print_=True
        print("\nAnnual produced hydrogen is equal to " + 
              str(round((sum(balances[loc]['hydrogen']['electrolyzer'] * c.timestep * 60) / 1e3), 2)) +
              " t/y ensuring a deliverable hydrogen mass flow rate equal to " + str(round((constantflow), 2)) + " kg/s")

    # ---------- PART 2: Hydrogen Chain Curves (electricity consumption) ----------
    
    with open('results/pkl/production_'+simulation_name+'.pkl', 'rb')            as f: production  = pickle.load(f)
    with open('results/pkl/consumption_'+simulation_name+'.pkl', 'rb')           as f: consumption  = pickle.load(f)
    
    if 'wind' not in production[loc]['electricity'] and 'PV' not in production[loc]['electricity']:    # if no renewables are included in the location        
        print("No renewable energy sources are present in the considered case study\nAutoconsumption data calculation not available\n")
        return

    el_consumption = consumption[loc]['electricity']
 
    # Iterate over the el_consumption dictionary
    for tech_name, consumption_data in el_consumption.items():
        # Only process 'electrolyzer' and 'mechanical compressor'
        if tech_name in ['electrolyzer', 'mechanical compressor']:
                        
            plt.figure(figsize=(10, 6), dpi=600)
            # Loop through each supplying component (tech_name1)
            for tech_name1, energy_supplied in consumption_data.items():
                # Calculate cumulative sum for this component
                cumulative_energy = np.cumsum(energy_supplied*c.timestep/60)
                plt.plot(cumulative_energy, label=f'{tech_name1} to {tech_name}')
            
            plt.grid(alpha = 0.3, zorder = 0)
            plt.xlabel('Timestep',fontsize=14)
            plt.ylabel('Cumulative Energy [kWh]',fontsize=14)
            plt.title(f'Cumulative Energy Supplied to {tech_name}',fontsize=14)
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
            plt.legend(loc='upper left',fontsize=14)
            plt.show()
            
   
    # Function to extract priority from the studycase and sort tech_name1 components
    def get_sorted_components(studycase, tech_name):
        # Extract the system for the current location (assuming a generic loc)
        system = studycase.get("industrial_facility", {})
    
        # Create a sorted list of tech_name1 based on priority for the given tech_name
        sorted_tech_name1 = [
            key for key, value in sorted(system.items(), key=lambda item: item[1].get('priority', float('inf')))
            if key in el_consumption.get(tech_name, {})
        ]
        
        return sorted_tech_name1

    # stacked cumulative energy curves
    # Iterate over the el_consumption dictionary
    for tech_name, consumption_data in el_consumption.items():
        # Only process 'electrolyzer' and 'mechanical compressor'
        if tech_name in ['electrolyzer', 'mechanical compressor']:
            plt.figure(figsize=(10, 6), dpi=600)

            # Get sorted tech_name1 components based on priority from the studycase
            sorted_tech_name1 = get_sorted_components(studycase, tech_name)

            # Initialize an array to hold the sum of previous components for stacking
            previous_cumulative = np.zeros(len(next(iter(consumption_data.values()))))

            # Loop through each tech_name1 in priority order and plot
            for tech_name1 in sorted_tech_name1:
                # Calculate the cumulative energy for this component (tech_name1)
                cumulative_energy = np.cumsum(consumption_data[tech_name1]*c.timestep/60)

                # Plot the current cumulative energy stacked on the previous
                plt.fill_between(
                    range(len(cumulative_energy)),
                    previous_cumulative,
                    previous_cumulative + cumulative_energy,
                    label=f'{tech_name1} to {tech_name}'
                )

                # Update the previous cumulative sum for stacking
                previous_cumulative += cumulative_energy

            # Set plot labels and title
            plt.grid(alpha = 0.3, zorder = 0)
            plt.xlabel('Timestep',fontsize=14)
            plt.ylabel('Cumulative Energy [kWh]',fontsize=14)
            plt.title(f'Stacked Cumulative Energy Supplied to {tech_name}',fontsize=14)
            plt.legend(loc='upper left',fontsize=14)
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
            plt.show()
                   
    return 

def ghg_emissions(simulation_name,path,loc,energy_market, print_= False):
    """
    GHG emissions profile (green index) calculation for hydrogen production.
    Dependent on energy_balance_results function output.
    ----------
            
    simulation_name : str - name of energy_balances file .pkl where balances_pp results are stored
        
    loc : str - location_name                            
    
    energy_market : dictionary
                        
    output: float - ghg emission value [kgCO2/kgH2]
    """      

    with open('results/pkl/consumption_'+simulation_name+'.pkl', 'rb')           as f: consumption  = pickle.load(f) 
    with open('results/pkl/production_'+simulation_name+'.pkl', 'rb')           as f: production  = pickle.load(f)

    el_consumption = consumption[loc]['electricity']
    
    # Number of steps per hour based on the timestep
    steps_per_hour = 60 // c.timestep
    # Determine the total number of hours in the simulation
    total_hours = 8760 * c.simulation_years
    # Prepare an array to store hourly grid electricity consumption
    el_from_grid_hyd_tot = np.zeros(total_hours)
    
    # Load emission intensity, either a constant or a time series
    if isinstance(energy_market['electricity']['emission intensity'], str):  # Check if it's a file path (string)
        em_intensity = (pd.read_csv(path + '/grid_emission_intensity/' + energy_market['electricity']['emission intensity'])['0'].to_numpy()) / 1e3
        if c.simulation_years > 1:
            # Repeat the emission intensity series if the simulation lasts multiple years
            em_intensity = np.tile(em_intensity, int(c.simulation_years))
    else:
        # If the emission intensity is constant (numeric), use that value for the whole simulation period
        em_intensity = np.full(total_hours, energy_market['electricity']['emission intensity'] / 1e3)
    
    # Iterate over the electricity consumption dictionary
    for tech_name, consumption_data in el_consumption.items():
        # Only process 'electrolyzer' and 'mechanical compressor'
        if tech_name in ['electrolyzer', 'mechanical compressor']:
            for tech_name1, energy_supplied in consumption_data.items():
                if tech_name1 == 'electricity grid':
                    # Reshape and sum the consumption data to get hourly totals
                    if c.timestep != 60:
                        # Aggregate data into hourly intervals
                        energy_supplied_hourly = energy_supplied.reshape(-1, steps_per_hour).sum(axis=1)*c.timestep/60
                    else:
                        # If timestep is 60 minutes, no aggregation is needed
                        energy_supplied_hourly = energy_supplied
    
                    # Add the hourly grid consumption to the total
                    el_from_grid_hyd_tot[:len(energy_supplied_hourly)] += energy_supplied_hourly  # [kWh]
    
    # Calculate the CO2 emissions by multiplying electricity consumption with emission intensity
    co2_emissions = el_from_grid_hyd_tot * em_intensity
    
    # Total CO2 emissions over the entire simulation period
    co2_tot = np.sum(co2_emissions)  # [kgCO2] total amount of carbon dioxide due to grid electricity utilization                
    produced_hyd = sum(production[loc]['hydrogen']['electrolyzer']['Tot'])*(c.timestep*60)      # [kgH2] total amount of produced hydrogen via in situ electorlysis
    h2_ghg = round(co2_tot/produced_hyd,2)  # [kgCO2/kgH2] GHG intensity of the produced hydrogen
    if print_ == True:
        print(f"\nThe H2 GHG intensity calculated for the considered scenario results in {h2_ghg} kgCO2/kgH2")

    return h2_ghg
def plot_energy_balances(simulation_name,loc,first_day,last_day,carrier,width=0.9):
    
    with open('results/pkl/consumption_'+simulation_name+'.pkl', 'rb')           as f: consumption  = pickle.load(f) 
    with open('results/pkl/production_'+simulation_name+'.pkl', 'rb')           as f: production  = pickle.load(f)
    with open('results/pkl/balances_'+simulation_name+'.pkl', 'rb')             as f: balances = pickle.load(f)
     
    UM = {
        'electricity': 'kW',
        'heating water': 'kW',
        'cooling water': 'kW',
        'process heat': 'kW',
        'process hot water': 'kW',
        'process cold water': 'kW',
        'process chilled water': 'kW',
        'hydrogen': 'kg/s',
        'LP hydrogen': 'kg/s',
        'HP hydrogen': 'kg/s',
        'oxygen': 'kg/s',
        'process steam': 'kg/s',
        'gas': 'Sm^3/s',
        'water': 'm^3/s'
    }
    
    x = np.arange(first_day*24*60/c.timestep,(last_day+1)*24*60/c.timestep)
    hourly_steps = 60//c.timestep # number of simulation steps considered in 1 hour - depending on input parameters. If simulaton is on hourly basis, hourly_steps=1

    self_consumption = {}
    surplus = {}
    for tech_name in production[loc][carrier]:
        self_consumption[tech_name] = np.zeros(len(x))
        surplus[tech_name] = np.zeros(len(x))
        for tech in production[loc][carrier][tech_name]:
            if tech not in [f'{carrier} grid','Tot']:
                self_consumption[tech_name] += production[loc][carrier][tech_name][tech][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps]
            elif tech in [f'{carrier} grid']:
                surplus[tech_name] += production[loc][carrier][tech_name][tech][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps]
    
    consumption_tech = {}    
    for tech_name in consumption[loc][carrier]:
        consumption_tech[tech_name] = {}
        for tech in consumption[loc][carrier][tech_name]:
            if tech not in ['Tot']:
                consumption_tech[tech_name][tech] = np.zeros(len(x))
                consumption_tech[tech_name][tech] += -consumption[loc][carrier][tech_name][tech][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps]
            
    tot_consumption_tech = {}
    for tech_name in consumption[loc][carrier]:
        tot_consumption_tech[tech_name] = np.zeros(len(x))
        tot_consumption_tech[tech_name] +=  -consumption[loc][carrier][tech_name]['Tot'][first_day*24*hourly_steps:last_day*24*hourly_steps+24*hourly_steps]
        
    if 'collective self consumption' in balances[loc][carrier]:
        to_csc = np.zeros(int(24*(last_day-first_day+1)*60//c.timestep))
        from_csc = np.zeros(int(24*(last_day-first_day+1)*60//c.timestep))
        for i,e in enumerate(balances[loc][carrier]['collective self consumption'][first_day*24*60//c.timestep:last_day*24*60//c.timestep+24*60//c.timestep]):
            if e > 0:
                from_csc[i] = e
            else:
                to_csc[i] = e 
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
    
    custom_colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:gray',
                     'tab:olive','tab:cyan','lightgreen','tomato','deeppink','chocolate','gold','darkblue',
                     'bisque','darkgreen','maroon','magenta']
    ax.set_prop_cycle(color=custom_colors)

    bottom = np.zeros(len(x))
    for tech_name in self_consumption:
        if np.sum(self_consumption[tech_name]) != 0:
            if tech_name != f'{carrier} grid':
                label=f'{tech_name} sc'
            else:
                label = 'from grid'
            ax.bar(x,self_consumption[tech_name],bottom=bottom,width=width,label=label)
            bottom += self_consumption[tech_name]
    load = bottom.copy()
    for tech_name in surplus:
        if np.sum(surplus[tech_name]) != 0:  
            ax.bar(x,surplus[tech_name],bottom=bottom,width=width,label=f'{tech_name} surplus')
            bottom += surplus[tech_name]
            
    if 'collective self consumption' in balances[loc][carrier]:
        if np.sum(from_csc) != 0:
            bottom = sum(value for key,value in self_consumption.items() if key != 'electricity grid')
            ax.bar(x,from_csc,bottom=bottom,width=width,label='Demand from CSC')

    bottom = np.zeros(len(x))
    for tech_name in consumption_tech:
        for tech in consumption_tech[tech_name]:
            if np.sum(consumption_tech[tech_name][tech]) != 0:  
                if tech_name != f'{carrier} grid':
                    if tech_name != 'mechanical compressor':
                        label=f'{tech_name} from {tech}'
                    else:
                        label=f'compressor from {tech}'
                else:
                    label = f'{tech} into grid'
                ax.bar(x,consumption_tech[tech_name][tech],bottom=bottom,width=width,label=label)
                bottom += consumption_tech[tech_name][tech]
                
    if 'collective self consumption' in balances[loc][carrier]:
        if np.sum(to_csc) != 0:
            ax.bar(x,to_csc,bottom=-load,width=width,label='Surplus to CSC')
      
    plt.title(loc+' '+carrier+' balances days '+str(first_day)+'-'+str(last_day))
    plt.plot(x,load,'k',label='load')     
    plt.legend(ncol=2, bbox_to_anchor = (1.01,-0.11))
    plt.ylabel(f"[{UM[carrier]}]")
    if (last_day-first_day) <= 10 and c.timestep == 60: # in order to better visualize daily behaviour if short timespans are selected
        # plt.xticks(list(range(first_day*24*hourly_steps, ((last_day+1)*24*hourly_steps)+1,24*hourly_steps)), [str(x) for x in list(range(first_day*24, ((last_day+1)*24)+1,24))], rotation=45) 
        plt.xticks(list(range(int(x[0]),int(np.ceil(x[-1]))+2,24)), [str(x) for x in list(range(first_day*24, ((last_day+1)*24)+1,24))], rotation=45) 
        plt.xlabel("Time  [h]")
    else:
        plt.xlabel("Timestep")
    # ax.xaxis.set_tick_params(bottom=True,labelbottom=True)
    #plt.xticks([0,6,12,18,24],['0','6','12','18','24'],fontsize=10,color='g')
    #plt.xticks([0,6,12,18,24,30,36,42,48],['0','6','12','18','24','30','36','42','48'],fontsize=10,color='g')
    #plt.yticks([-2,-1,0,1,2],['-2','-1','0','1','2'])
    plt.show() 
    
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
    
    custom_colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:gray',
                     'tab:olive','tab:cyan','lightgreen','tomato','deeppink','chocolate','gold','darkblue',
                     'bisque','darkgreen','maroon','magenta']
    ax.set_prop_cycle(color=custom_colors)

    bottom = np.zeros(len(x))
    for tech_name in self_consumption:
        if np.sum(self_consumption[tech_name]) != 0:
            if tech_name != f'{carrier} grid':
                label=f'{tech_name} sc'
            else:
                label = 'from grid'
            ax.bar(x,self_consumption[tech_name],bottom=bottom,width=width,label=label)
            bottom += self_consumption[tech_name]
    load = bottom.copy()
    for tech_name in surplus:
        if np.sum(surplus[tech_name]) != 0:  
            ax.bar(x,surplus[tech_name],bottom=bottom,width=width,label=f'{tech_name} surplus')
            bottom += surplus[tech_name]
            
    if 'collective self consumption' in balances[loc][carrier]:
        if np.sum(from_csc) != 0:
            bottom = sum(value for key,value in self_consumption.items() if key != 'electricity grid')
            ax.bar(x,from_csc,bottom=bottom,width=width,label='Demand from CSC')
              
    bottom = np.zeros(len(x))
    for tech_name in tot_consumption_tech:
        if np.sum(tot_consumption_tech[tech_name]) != 0:  
            if tech_name != f'{carrier} grid':
                if tech_name != 'mechanical compressor':
                    label=f'{tech_name}'
                else:
                    label='compressor'
            else:
                label = 'into grid'
            ax.bar(x,tot_consumption_tech[tech_name],bottom=bottom,width=width,label=label)
            bottom += tot_consumption_tech[tech_name]
            
    if 'collective self consumption' in balances[loc][carrier]:
        if np.sum(to_csc) != 0:
            ax.bar(x,to_csc,bottom=-load,width=width,label='Surplus to CSC')
    plt.title(loc+' '+carrier+' balances days '+str(first_day)+'-'+str(last_day))
    plt.plot(x,load,'k',label='load')     
    plt.legend(ncol=2, bbox_to_anchor = (1.01,-0.11))
    plt.ylabel(f"[{UM[carrier]}]")
    if (last_day-first_day) <= 10 and c.timestep == 60: # in order to better visualize daily behaviour if short timespans are selected
        # plt.xticks(list(range(first_day*24*hourly_steps, ((last_day+1)*24*hourly_steps)+1,24*hourly_steps)), [str(x) for x in list(range(first_day*24, ((last_day+1)*24)+1,24))], rotation=45) 
        plt.xticks(list(range(int(x[0]),int(np.ceil(x[-1]))+2,24)), [str(x) for x in list(range(first_day*24, ((last_day+1)*24)+1,24))], rotation=45) 
        plt.xlabel("Time  [h]")
    else:
        plt.xlabel("Timestep")
    # ax.xaxis.set_tick_params(bottom=True,labelbottom=True)
    #plt.xticks([0,6,12,18,24],['0','6','12','18','24'],fontsize=10,color='g')
    #plt.xticks([0,6,12,18,24,30,36,42,48],['0','6','12','18','24','30','36','42','48'],fontsize=10,color='g')
    #plt.yticks([-2,-1,0,1,2],['-2','-1','0','1','2'])
    plt.show() 
    


def print_and_plot_annual_energy_balances(simulation_name, loc, print_= False):

    # Load consumption and production data from pickle files
    with open('results/pkl/consumption_'+simulation_name+'.pkl', 'rb') as f:
        consumption = pickle.load(f)
    with open('results/pkl/production_'+simulation_name+'.pkl', 'rb') as f:
        production = pickle.load(f)
    with open('results/pkl/balances_'+simulation_name+'.pkl', 'rb') as f:
        balances = pickle.load(f)
      
    UM = {
        'electricity': 'MWh',
        'heating water': 'MWh',
        'cooling water': 'MWh',
        'process heat': 'MWh',
        'process hot water': 'MWh',
        'process cold water': 'MWh',
        'process chilled water': 'MWh',
        'hydrogen': 't',
        'LP hydrogen': 't',
        'HP hydrogen': 't',
        'oxygen': 't',
        'process steam': 't',
        'gas': 'Sm^3',
        'water': 'm^3'
    }

    # Plot production data
    for carrier in production[loc]:
        if print_ == True:
            print(f"\n\nProduction: Annual balances {carrier} [{UM[carrier]}]")

        tech_names = []
        tech_totals = []

        for tech_name in production[loc][carrier]:
            for tech in production[loc][carrier][tech_name]:
                # Calculate the total production for each tech_name and carrier
                if carrier in ['hydrogen', 'LP hydrogen', 'HP hydrogen', 'oxygen', 'process steam']:
                    tot = ((sum(production[loc][carrier][tech_name][tech]) * c.timestep / 60) * 3600) / c.simulation_years
                else:
                    tot = (sum(production[loc][carrier][tech_name][tech]) * c.timestep / 60) / c.simulation_years

                if tot != 0:
                    if carrier not in ['gas', 'water']:
                        tot = round(tot / 1e3, 2)
                    else:
                        tot = round(tot, 2)

                    if tech != 'Tot':
                        if print_ == True:
                            print(f"{tech_name} to {tech} = {tot} {UM[carrier]}")
                        tech_names.append(f"{tech_name} to {tech}")
                    else:
                        if print_ == True:
                            print(f"{tech_name} production = {tot} {UM[carrier]}")
                        tech_names.append(f"{tech_name} production")
                    
                    tech_totals.append(tot)
         
        # Initialize CSC-related variables to avoid UnboundLocalError
        csc_value = []
        csc_index = []
        csc_names = []
        
        if 'collective self consumption' in balances[loc][carrier]:
            csc = balances[loc][carrier]['collective self consumption']
            to_csc = round(-np.sum(csc, where=(csc < 0)) / 1e3, 2)
            
            if to_csc != 0:
    
                if print_ == True:
                    print(f"surplus to CSC = {to_csc} {UM[carrier]}")
                
                remaining_csc= to_csc
                csc_index = []
                csc_value = []
                csc_names = []
                for tech_name in production[loc][carrier]:
                    while remaining_csc > 0:
                        index = tech_names.index(f"{tech_name} to electricity grid")
                        
                        tech_name_csc = np.min([tech_totals[index],remaining_csc])
                        
                        csc_index.append(index)
                        csc_value.append(tech_name_csc)
                        csc_names.append(f"{tech_name} to CSC")
                        
                        remaining_csc -= tech_name_csc
        # Plotting the bar chart with updated configurations
        if len(tech_totals) > 2 :  # Only plot if there are more than two bars
            plt.figure(figsize=(10, 8),dpi=600)
            
            # Generate a color map for the bars
            colors = cm.get_cmap('tab20', len(tech_totals)).colors
            
            # Plot bars with different colors
            bars = plt.bar(np.arange(len(tech_totals)), tech_totals, color=colors, edgecolor='k')
        
            # Modify the legend labels to remove 'production' or 'consumption'
            clean_tech_names = [name.replace(' production', 'tot').replace(' consumption', 'tot') for name in tech_names]
            
            if len(csc_value) > 0 and np.any(csc_value != 0):
                csc_colors = np.array([colors[i % len(colors)] for i in csc_index])
                csc_bars = plt.bar(csc_index, csc_value, color=csc_colors, edgecolor='k', hatch = '/')
                plt.legend(bars + csc_bars, tech_names + csc_names, loc='upper center', bbox_to_anchor=(0.5, -0.15),
                   fancybox=True, shadow=False, ncol=min(2, len(tech_names) + len(csc_names)), fontsize=14)
            else:
                plt.legend(bars, tech_names, loc='upper center', bbox_to_anchor=(0.5, -0.15),
                           fancybox=True, shadow=False, ncol=min(2, len(tech_names)), fontsize=14)

            # Remove xticks and set the fontsize for y-label and title
            plt.xticks([])  # Remove xticks
            plt.ylabel(f'Production [{UM[carrier]}]', fontsize=14)
            plt.title(f'Production - {carrier} - {loc}', fontsize=14)
            plt.yticks(fontsize=14)  # Set fontsize for yticks
        
            # Show the plot
            plt.show()

    # Plot consumption data
    for carrier in consumption[loc]:
        if print_ == True:
            print(f"\n\nConsumption: Annual balances {carrier} [{UM[carrier]}]")

        tech_names = []
        tech_totals = []

        for tech_name in consumption[loc][carrier]:
            for tech in consumption[loc][carrier][tech_name]:
                # Calculate the total consumption for each tech_name and carrier
                if carrier in ['hydrogen', 'LP hydrogen', 'HP hydrogen', 'oxygen', 'process steam']:
                    tot = ((sum(consumption[loc][carrier][tech_name][tech]) * c.timestep / 60) * 3600) / c.simulation_years
                else:
                    tot = (sum(consumption[loc][carrier][tech_name][tech]) * c.timestep / 60) / c.simulation_years

                if tot != 0:
                    if carrier not in ['gas', 'water']:
                        tot = round(tot / 1e3, 2)
                    else:
                        tot = round(tot, 2)

                    if tech != 'Tot':
                        if print_ == True:
                            print(f"{tech_name} from {tech} = {tot} {UM[carrier]}")
                        tech_names.append(f"{tech_name} from {tech}")
                    else:
                        if print_ == True:
                            print(f"{tech_name} consumption = {tot} {UM[carrier]}")
                        tech_names.append(f"{tech_name} consumption")
                    
                    tech_totals.append(tot)
          
        if 'collective self consumption' in balances[loc][carrier]:
            csc = balances[loc][carrier]['collective self consumption']
            from_csc = round(np.sum(csc, where=(csc > 0)) / 1e3, 2)
            
            if from_csc != 0:
                if print_ == True:
                    print(f"CSC to demand = {from_csc} {UM[carrier]}")
                               
                grid_index = tech_names.index("electricity demand from electricity grid")
                csc_name = ['electricity demand from csc']

        if len(tech_totals) > 2 or from_csc != 0:  # Only plot if there are more than two bars
            plt.figure(figsize=(10, 8),dpi=600)
            
            # Generate a color map for the bars
            colors = cm.get_cmap('tab20', len(tech_totals)).colors
            
            # Plot bars with different colors
            bars = plt.bar(np.arange(len(tech_totals)), tech_totals, color=colors, edgecolor='k')
        
            # Modify the legend labels to remove 'production' or 'consumption'
            clean_tech_names = [name.replace(' production', 'tot').replace(' consumption', 'tot') for name in tech_names]
        
            # Add the legend below the plot with up to 4 columns
            if from_csc != 0:
                csc_bars = plt.bar(grid_index, from_csc, color=colors[grid_index], edgecolor='k', hatch='/')
                plt.legend(bars + csc_bars, tech_names + csc_name, loc='upper center', bbox_to_anchor=(0.5, -0.15),
                           fancybox=True, shadow=False, ncol=min(2, len(tech_names) + len(csc_name)), fontsize=14)
            else:
                plt.legend(bars, tech_names, loc='upper center', bbox_to_anchor=(0.5, -0.15),
                           fancybox=True, shadow=False, ncol=min(2, len(tech_names)), fontsize=14)

            # Remove xticks and set the fontsize for y-label and title
            plt.xticks([])  # Remove xticks
            plt.ylabel(f'Consumption [{UM[carrier]}]', fontsize=14)
            plt.title(f'Consumption - {carrier} - {loc}', fontsize=14)
            plt.yticks(fontsize=14)  # Set fontsize for yticks
        
            # Show the plot
            plt.show()
                                                                                  