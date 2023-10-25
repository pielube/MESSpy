"""
time_slot

"""
import numpy as np
import pandas as pd

def hourly_serie(f1p,f2p,f3p,filename):
    """"
    f1p: energy price in time slot F1
    f2p: F2
    f3p: F3
    name: str name of the .csv file that is going to be written
    
    output:
        hourly serie of energy_price filename.csv    
    """
    
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
    
    for d in range(365):        
    
        dt = Day_type[d] # 1=Weekday 2=Saturday 3=Sunday 
        
        for h in range(24):
            
            if dt == 1:
                f = Weekday[h]
            if dt == 2:
                f = Saturday[h]
            if dt == 3:
                f = Sunday[h]
                       
            if f == 1:
                p = f1p
            if f == 2:
                p = f2p
            if f == 3:
                p = f3p
            
            energy_price.append(p)
    
    series_frame = pd.DataFrame(energy_price)
    series_frame.to_csv(f"{filename}.csv")
        

if __name__ == "__main__":
    # energy hourly price
    f1p = 0.316
    f2p = 0.310
    f3p = 0.268
    
    hourly_serie(f1p,f2p,f3p,"electricity_price")