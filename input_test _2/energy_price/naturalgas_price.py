# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 16:05:33 2022

@author: mati
"""

import pandas as pd
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),os.path.pardir,os.path.pardir)))
import constants as c

#%%
LHV_NG=c.LHV_NGSVOL  #[kWh/Sm3]

def gas_price(p1,p2,p3,p4,filename):
    """
    Creates a general Natural Gas prices series 
    
    Parameters
    ----------
    p1 : float value [€/kWh]
         1st quarter value.
    p2 : float value [€/kWh]
         2nd quarter value.
    p3 : float value [€/kWh]
         3rd quarter value.
    p4 : float value [€/kWh]
         4th quarter value.
    filename: string value 
              name of the files where to store data

    Returns
    -------
    file.csv 

    """
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
        
    series_frame = pd.DataFrame(ng_price,columns=['0'])
    series_frame.to_csv(f"{filename}.csv",encoding='utf-8-sig')

#%% ###########################################################################

if __name__ == "__main__":
    # energy hourly price
    p1 = 0.316
    p2 = 0.310
    p3 = 0.268
    p4 = 0.268
    
    gas_price(p1,p2,p3,p4,"naturalgas_price")