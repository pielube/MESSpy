import pandas as pd
import numpy as np


lista = ['heating-dhw_test.csv','hydrogen_demand_constant.csv','hydrogen_demand_variable.csv','load_test_1.csv','load_test_2.csv','load_test_3.csv','load_test_4.csv']

for l in lista:
    
    try:
        a = pd.read_csv('E:\\mattia\\desktop\\GIT\\MESSpy\\input_test\\loads\\'+l)['kW']
    except:
        a = pd.read_csv('E:\\mattia\\desktop\\GIT\\MESSpy\\input_test\\loads\\'+l)['kg/s']
    
    b =  np.repeat(a, 4)
    
    b = pd.DataFrame(b)
    b.to_csv('E:\\mattia\\desktop\\GIT\\MESSpy\\input_test\\loads\\'+'4'+l)