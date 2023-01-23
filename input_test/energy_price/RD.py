import numpy as np
import pandas as pd
"""
Prezzi energia venduta al Ritiro Dedicato
https://www.gse.it/servizi-per-te/fotovoltaico/ritiro-dedicato
https://www.gse.it/servizi-per-te/fotovoltaico/ritiro-dedicato/documenti
"""

### Il GSE fornisce i prezzi medi divisi per mese, fascia oraria e zona geografica
### questo codice crea le serie oraria (8760 prezzi)
### Per adesso sono inseriti gli anni 2019, 2020 e 2021 della zona 'Centro-Nord'

def RD(year):
    # Input: year
    # Output: hourly serie €/MWh
   
    if year == 2019:
        #prezzi 2019 Centro-Nord
        pf1=[75.45, 59.08, 55.15, 58.00, 52.72, 53.00, 56.69, 48.28, 54.41, 57.92, 56.59, 49.56]
        pf2=[68.21, 55.93, 46.72, 48.05, 51.12, 46.74, 48.49, 43.83, 45.10, 45.96, 50.58, 45.02]
        pf3=[57.98, 48.16, 48.17, 43.41, 43.94, 39.94, 41.21, 38.82, 42.07, 39.48, 39.76, 34.05]
        
    if year == 2020:
        #prezzi 2020 Centro-Nord
        pf1=[5.,70, 41.78, 32.73, 22.74, 22.69, 30.53, 42.76, 41.00, 53.29, 47.55, 53.18, 67.52]
        pf2=[48.02, 37.19, 33.70, 22.55, 21.40, 24.28, 35.37, 37.89, 41.93, 41.37, 47.41, 57.47]
        pf3=[39.42, 31.93, 23.25, 14.53, 14.55, 21.13, 25.81, 29.74, 35.47, 30.89, 39.70, 42.79]
       
    if year == 2021:
        #prezzi 2021 Centro-Nord
        pf1=[75.30, 60.78, 60.53, 70.41, 74.40, 88.61, 108.72, 113.63, 162.70, 230.79, 260.47, 328.68]
        pf2=[61.36, 57.63, 62.69, 67.80, 68.80, 78.89, 99.27, 103.12, 149.86, 206.17, 217.75, 291.38]
        pf3=[52.02, 40.85, 53.35, 55.64, 56.36, 70.35, 85.31, 92.72, 134.58, 181.99, 194.33, 245.34]
    
    
    ############################################ Time slot definition F1 F2 F3
    
    Weekday = np.ones(24)
    Saturday = np.zeros(24)
    Sunday = np.zeros(24)    
    Day_type = np.ones(365) # 1=Weekday 2=Saturday 3=Sunday
     
    ### F3
    Day_type[6:365:7] = 3 
    Sunday[0:24] = 3 
    Weekday[0:7] = 3
    Weekday[23] = 3
    Saturday[0:7] = 3
    Saturday[23] = 3    
    
    ### F2
    Day_type[5:365:7] = 2 
    Saturday[7:23] = 2
    Weekday[7] = 2
    Weekday[19:23] = 2     
    
    ### F1 where the value 1 remains
    
    energy_price = [] # idx = hour (0-8760), value = energy_price
    
    dm=[31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]  # duration of months [days]
    months_365 = [] # idx = day (0-365), value = month (0-11). Usefull below.
    for m in range(len(dm)):
        for d in range(dm[m]):    
            months_365.append(m)
    
    for d in range(365):   
    
        dt = Day_type[d] # 1=Weekday 2=Saturday 3=Sunday 
        m = months_365[d]
        
        for h in range(24):
            
            if dt == 1:
                f = Weekday[h]
            if dt == 2:
                f = Saturday[h]
            if dt == 3:
                f = Sunday[h]
                       
            if f == 1:
                p = pf1[m]
            if f == 2:
                p = pf2[m]
            if f == 3:
                p = pf3[m]
            
            energy_price.append(round(p/1000,3))
    
    
    series_frame = pd.DataFrame(energy_price)
    series_frame.to_csv(f"RD_{year}.csv")
    
    return()
        

if __name__ == "__main__":
    RD(2019)
    RD(2020)
    RD(2021)



