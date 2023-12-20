# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 14:55:54 2022

@author: mati
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os.path
import matplotlib.dates as mdates
from pandas.tseries.frequencies import to_offset

path = r'./NaturalGas'         # Folder contaning natural gas prices
os.chdir(path)

#%%
thermal_year = ['2018_2019','2019_2020','2020_2021','2021_2022','2022_2023']
year = []
plots = {}

for count,item in enumerate(thermal_year):       # data reading and handling
    
    if count < 4:
    
        price_file  = thermal_year[count] +'.xlsx'
        price_file1 = thermal_year[count+1] +'.xlsx'
        price_data  = pd.read_excel(price_file, sheet_name = 1, usecols = "A:B", skiprows=5)
        price_data1 = pd.read_excel(price_file1, sheet_name = 1, usecols = "A:B", skiprows=5)
        
        price_data['data']  = pd.to_datetime(price_data['data'].astype(str), format='%Y%m%d')
        price_data1['data'] = pd.to_datetime(price_data1['data'].astype(str), format='%Y%m%d')
        
        slices = str(price_data1['data'][0])[:4]
        year.append(slices)
        
        mask  = (price_data['data'] >= slices+'-01-01') & (price_data['data'] <= slices+'-09-30')
        mask1 = (price_data1['data'] >= slices+'-10-01') & (price_data1['data'] <= slices+'-12-31')

        df2   = price_data.loc[mask]
        df3   = price_data1.loc[mask1]
        
        couple =[df2,df3]
        
        df = pd.concat(couple)                   # creating annual series Jan-Dec 
        df = df.reset_index(drop=True)
        
        df['€/kWh'] = df['€/MWh']/1000           # defining costs in €/kWh 
        
        hours = len(df)*24                       # creating hourly series needed as input in MESS
        ng_cost = np.zeros(hours)
        
        for i in df.index:
            ng_cost[i*24:i*24+24] = df['€/kWh'][i]
            
        plots[slices]= ng_cost                   # dictionary where all hourly series are saved
        
        if len(ng_cost) > 8760:                  # check on leap-year
            ng_cost = ng_cost[0:8760]
            
        # Monthly mean for natural gas prices 
       
        a = df.groupby(pd.PeriodIndex(df['data'], freq="M"))['€/MWh'].mean()
        a = a.resample('M').mean()
        a.index = a.index.to_timestamp() + to_offset("14D")
        
        data = pd.DataFrame(ng_cost, columns=['0'])  
        data.to_csv('../NG_'+slices+'.csv',encoding='utf-8-sig')   # saving .csv files
        
# %%    Single Plots     

        fig,ax = plt.subplots(dpi=300, figsize = (8,5))
        
        ax.plot(df['data'],df['€/MWh'], label = 'Daily price', color = 'cornflowerblue', alpha = 0.8)
        ax.plot(df['data'],df['€/MWh'].rolling(5).mean(),label = 'Moving average 5 days')
        # ax.plot(a, label = 'Monthly Average')
        
        locator   = mdates.MonthLocator()
        formatter = mdates.ConciseDateFormatter(locator)
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(formatter)
        
        ax.set_xlim(df['data'][0],df['data'][len(df['data'])-1])
        
        mean_price = int(np.mean(df['€/MWh']))
        
        ax.axhline(mean_price, label = 'Mean price: %d €/MWh' %mean_price, color = 'forestgreen', alpha = 0.9)
        
        ax.set_xlabel('Time')
        ax.set_ylabel('MI-GAS [€/MWh]')
        ax.set_title('Natural Gas Price year '+slices)
        
        ax.grid(alpha = 0.3)
       
        ax.legend(loc = 'best')

# %%    Overall Plot - Comparison over years  
        
plt.figure(dpi=500)

for slices in plots:
    
    plt.plot(plots[slices],label=slices)     

plt.xlabel('Time [h]')
plt.ylabel('MI-GAS [€/kWh]')
plt.legend(loc='best',fontsize='small')
plt.grid(alpha = 0.3)
plt.title('Natural Gas prices - Italian market')
        
os.chdir('..')