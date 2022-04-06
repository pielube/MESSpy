"""
Post processing
"""
import pickle
import numpy as np
import matplotlib.pyplot as plt
import plotly.io as pio
import plotly.graph_objects as go


def total_balances(simulation_name):

    with open('Results/balances_'+simulation_name+'.pkl', 'rb') as f:
        balances = pickle.load(f)
    
    ###### load analysys
    
    ##### total electricity balances
    carriers = ['electricity']
    units = {'electricity': 'kWh', 'hydrogen': 'kg'}
    
    for carrier in carriers:
        print('\nTotal '+carrier+' balances:\n')
        for loc in balances:
            print('\n'+loc)   
            for b in balances[loc][carrier]:
                
                positiv=balances[loc][carrier][b][balances[loc][carrier][b]>0].sum()
                negativ=balances[loc][carrier][b][balances[loc][carrier][b]<0].sum()
                
                if positiv != 0:
                    print(b+' '+str(round(positiv,1))+' '+units[carrier])
                    
                if negativ != 0:
                    print(b+' '+str(round(negativ,1))+' '+units[carrier])
    print('\n')
    
def NPV_plot():
    ##### economic
    simulation_name = 'economic_assessment'
    with open('Results/'+simulation_name+'.pkl', 'rb') as f:
        economic = pickle.load(f)
    
    plt.figure(dpi=1000)
    
    for loc in economic:    
        y = economic[loc]['NPV']
        x = np.linspace(0,len(y)-1,len(y))
        plt.plot(x,y,label=loc)
        
    plt.legend()
    plt.grid()
    plt.ylabel('Net Present Value [€]')
    plt.xlabel('Time [years]')
    plt.xlim(0,len(y)-1)
    #plt.ylim(-12000,12000)
    plt.show()
    
def SOC_plot(simulation_name):
      
    with open('Results/SOC_'+simulation_name+'.pkl', 'rb') as f:
        SOC = pickle.load(f)
           
    unit = {'H tank': '[m3]', 'battery': '[kWh]'}
    
    for loc in SOC:
        for tech in SOC[loc]:
            
            y = SOC[loc][tech]
            x = np.linspace(0,len(y)-1,len(y))        
            plt.plot(x,y,label=loc)
            plt.grid()
            plt.ylabel('SOC '+unit[tech])
            plt.xlabel('Time [hours]')
            plt.title(loc+' '+tech)
            plt.show()
        


#%%

def Flows(simulation_name,carrier='electricity'):

    pio.renderers.default='browser'
    with open('Results/balances_'+simulation_name+'.pkl', 'rb') as f:
        balances = pickle.load(f)
    
    ### data prepearing
    node_label = ["transformer","transformer","grid"]
    node_color = ['silver','silver','gray']
    source = []
    target = []
    value = []
    link_color = []
    link_label = []
        
    ### add from grid
    node_number = 3
    for loc in balances:
        if loc != 'REC':
            
            node_color.append('peru')
            node_label.append('load '+loc)
            
            source.append(0)
            target.append(node_number)
            value.append(np.sum(balances[loc][carrier]['grid'], where = balances[loc][carrier]['grid'] > 0))
            link_color.append('tomato')
            link_label.append('from grid')
            node_number += 1
    value.append(sum(value)-sum(balances['REC'][carrier]['collective self consumption']))
    source.append(2)
    target.append(0)
    link_color.append('tomato')
    link_label.append('from grid')
    
    ### add to grid and self consumption
    node_number2 = 3
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
                link_label.append('to grid')
                
                # self consumption
                source.append(node_number)
                target.append(node_number2)
                value.append(-np.sum(balances[loc][carrier]['demand']) - np.sum(balances[loc][carrier]['grid'], where = balances[loc][carrier]['grid'] > 0))
                link_color.append('yellowgreen')
                link_label.append('self consumption')
                
                node_number += 1
                node_number2 += 1
                
    ### LV to MV
    source.append(1)
    source.append(1)
    target.append(0)
    target.append(2)
    value.append(sum(balances['REC'][carrier]['collective self consumption'])) 
    value.append(- np.sum(balances['REC'][carrier]['into grid'])-sum(balances['REC'][carrier]['collective self consumption']))   
    link_color.append('gold')
    link_label.append('collective self consumption')
    link_color.append('orange')
    link_label.append('to grid')
                
    fig = go.Figure(data=[go.Sankey(
        valueformat = "0.f",
        valuesuffix = " kWh",
        arrangement = "snap",
        node = {
                'thickness': 30,
                #'line': dict(color = "black", width = 0.5),
                'label': node_label,
                'color': node_color,
                "x": [0.7,0.3,0.5],
                "y": [0.1,0.1,0.25], # stranamente è ribaltata rispetto ad annotations
                'pad': 100
                },
        link = {
                'source': source, # indices correspond to labels, eg A1, A2, A1, B1, ...
                'target': target,
                'value': value,
                'color': link_color,
                'label': link_label
                })])
    
    fig.update_layout(title_text="Energy flows", font_size=25,
                      annotations=[dict(
                          x=0.,
                          y=0.,
                          text='',
                          showarrow=False
                          )] )
    fig.show()


    
        