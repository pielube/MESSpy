# -*- coding: utf-8 -*-
"""
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
"""
Valore PUN 
https://www.mercatoelettrico.org/It/Statistiche/ME/DatiSintesi.aspx
"""

### Il GSE fornisce i prezzi medi divisi per mese
### questo codice crea le serie oraria (8760 prezzi)
### Per adesso sono inseriti gli anni 2021, 2022 e 2023

def pun(year):
    # Input: year
    # Output: hourly serie €/MWh
   
    if year == 2021:
        #pun medio 2021
        pun_mean=[60.71, 56.57, 60.39, 69.02, 69.91, 84.80, 102.66, 112.40, 158.59, 217.63, 225.95, 281.24]
        
    if year == 2022:
        #pun medio 2022
        pun_mean = [224.50, 211.69, 308.07, 245.97, 230.06, 271.31, 441.65, 543.15, 429.92, 211.50, 224.51, 294.91]
       
    if year == 2023:
        #pun medio 2023
        pun_mean = [174.49, 161.07, 136.38, 134.97, 105.73, 105.34, 112.09, 111.89, 115.70, 134.26, 121.74, 115.46]
    
    
    energy_price = [] # idx = hour (0-8760), value = energy_price
    weight = np.genfromtxt('weight.csv')   #opening weight array generated in PUN_curve.py
    
    dm=[31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]  # duration of months [days]
    months_365 = [] # idx = day (0-365), value = month (0-11). Usefull below.
    for m in range(len(dm)):
        for d in range(dm[m]):    
            months_365.append(m)
    
    for d in range(365):   
    
        m = months_365[d]
        
        for h in range(24):
            
            p = pun_mean[m]
            
            energy_price.append(round(p/1000,3))
    
    energy_price = np.array(energy_price)
    energy_price = energy_price*weight
    
    series_frame = pd.DataFrame(energy_price)
    series_frame.to_csv(f"pun_{year}.csv")
    
    plt.plot(energy_price[3072:3240])
    
    return()
        

if __name__ == "__main__":
    pun(2021)
    pun(2022)
    pun(2023)
    
    