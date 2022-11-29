# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 16:05:33 2022

@author: mati
"""

import pandas as pd
import os
import sys
path = os.getcwd()
sys.path.append(os.path.abspath(os.path.join(path, os.pardir)))   # temorarily adding constants module path 
import constants as c

#%%
LHV_NG=10.69 #[kWh/Sm3]  #Importare constants!!!  ed indirizzo relativo

def gas_price(p1,p2,p3,p4,filename):
    i=0
    ng_price=[]
    p1_kWh=p1/100/LHV_NG   #[€/kWh]
    p2_kWh=p2/100/LHV_NG   #[€/kWh]
    p3_kWh=p3/100/LHV_NG   #[€/kWh]
    p4_kWh=p4/100/LHV_NG   #[€/kWh]
    
    for i in range(0,int(8760/4)):
        ng_price.append(p1_kWh)
        i+=1
    for i in range(int(8760/4),int(8760/4*2)):
        ng_price.append(p2_kWh)
        i+=1
    for i in range(int(8760/4*2),int(8760/4*3)):
        ng_price.append(p3_kWh)
        i+=1
    for i in range(int(8760/4*3),8760):
        ng_price.append(p4_kWh)
        i+=1
        
    series_frame = pd.DataFrame(ng_price)
    series_frame.to_csv(f"{filename}.csv")

#%% ###########################################################################

if __name__ == "__main__":
    # energy hourly price
    p1 = 0.316
    p2 = 0.310
    p3 = 0.268
    p4 = 0.268
    
    gas_price(p1,p2,p3,p4,"naturalgas_price")