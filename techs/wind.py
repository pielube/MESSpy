import pandas as pd
import numpy as np
import os
import pickle
import math

class wind:    
    
    def __init__(self,params,weather,ts):
        """
        Create a wind object based on the specified model
    
        params : dictionary
            'model': str type of model to be used for wind
            
            'area': float swept area [m2] e.g. 39.6 m^2 (Aircon 10/10 kW)
            'efficiency': float total efficiency = Betz*efficiency [-] default: 0.45 (ca. 0.593*0.76, i.e. Betz*efficiency)
            'Prated': float rated power [kW]
            'WSrated': float rated wind speed [m/s] e.g. 11.0 m/s (Aircon 10/10 kW)
            'WScutin': float cut in wind speed [m/s] e.g.  2.5 m/s (Aircon 10/10 kW)
            'WScutoff': float cut off wind speed [m/s] e.g. 32.0 m/s (Aircon 10/10 kW)
            
            'omega_min': [rpm] e.g. 50  rpm (Aircon 10/10 kW)
            'omega_max': [rpm] e.g. 130 rpm (Aircon 10/10 kW)
            'beta': [°] e.g. 0°
            'cp_max': [-] e.g. 0.44 ;values from 0.4 to 0.5
            'idx': [-] e.g. 5; values from 0 to 5
            
            'z_i': [m] wind turbine height, (?)
            'z_hub': [m] hub height, e.g. 30 m (Aircon 10/10 kW)
            'alpha': [-] Hellman or shear coefficient, values from 0 to 0.4
            'Vu': [°/m] Veer coefficient, values from 0 to 0.75 °/m
            'Nbands': [-] Number of horizontal bands
            
            
        general: dictionary
            see rec.py

        output : wind object able to:
            produce electricity .use(h)
        """
        
        if params['model'] == 'availability':
            print('tobeadded')
        elif params['model'] == 'betz':
            
            ws_turbine = weather['wind'].copy()
            
            ws_turbine[ws_turbine<params['WScutin']] = 0.
            ws_turbine[(ws_turbine>params['WSrated']) & (ws_turbine<params['WScutout'])] = params['WSrated']
            ws_turbine[ws_turbine>params['WScutout']] = 0.
            
            self.energy = 0.5*weather['rho']*params['area']*ws_turbine**3*params['efficiency']/1000.*ts # kWh
            
        elif params['model'] == 'detailed':
            
            ws_turbine = weather['wind'].copy()
            ws_turbine[ws_turbine<params['WScutin']] = 0.
            ws_turbine[(ws_turbine>params['WSrated']) & (ws_turbine<params['WScutout'])] = params['WSrated']
            ws_turbine[ws_turbine>params['WScutout']] = 0.
            
            omega_min = 0.
            omega_max = 0.
            cp_max = 0.44
            
            if 'omega_min' in params:
                omega_min = params['omega_min']
            if 'omega_max' in params:
                omega_max = params['omega_max']
            if 'cp_max' in params:
                cp_max = params['cp_max']
                
            powercoeff = self.cpfunc(ws_turbine,params['area'],params['beta'],params['idx'],omega_min=omega_min,omega_max=omega_max,cp_max=cp_max)
            wseq = self.eqspeed(ws_turbine,params['z_i'],params['z_hub'],params['alpha'],params['area'],params['Vu'],params['Nbands'])
            
            self.energy = 0.5*weather['rho']*params['area']*wseq**3*powercoeff/1000.*ts # kWh

        elif params['model'] == 'atlite':
            print('tobeadded')
        else:
            print('Error: wrong wind turbine model')
        

        
    def use(self,h):
        """
        Produce electricity
        
        h : int hour to be simulated
    
        output : float electricity produced that hour [kWh]    
        """
        
        return(self.energy[h])
    
    def cpfunc(self,ws,area,beta,idx,omega_min=0,omega_max=0,cp_max=0.44):
        
        c1 = [0.73, 0.5, 0.5176, 0.77, 0.5,  0.22]
        c2 = [151., 116.,116.,151.,116.,120.]
        c3 = [0.58,0.4,0.4,0.,0.,0.4]
        c4 = [0.0,0.0,0.0,0.0,0.4,0.]
        c5 = [0.002,0.,0.,0.,0.,0.]        
        x =  [2.14,0.0,0.0,0.0,0.0,0.0]
        c6 = [13.2,5.0,5.0,13.65,5.0,5.0]
        c7 = [18.4,21.,21.,18.4,21.,12.5]
        c8 = [0.0,0.0,0.006795,0.0,0.0,0.0]
        c9 = [-0.02,0.089,0.089,0.0,0.08,0.08]
        c10 = [0.03,0.035,0.035,0.0,0.035,0.035]
        
        diam = (4*area/math.pi)**(1/2)
        
        # cp_max and lambda_opt for beta = 0
        lambdas_0 = np.arange(0.1,12.,0.1)
        cps_0 = c1[idx]*((c2[idx]/(1.0/((1.0/lambdas_0)-c10[idx])))-c6[idx])*np.exp(-c7[idx]/(1.0/((1.0/lambdas_0)-c10[idx])))+c8[idx]*lambdas_0 # Eq. 2(I) + 2(II) (beta=0)
        cp_max_0 = np.max(cps_0)
        lambda_opt_0 = lambdas_0[np.argmax(cps_0)]
                        
        # if omega_min and omega_max are not given as inputs:     
        if omega_min == 0:
          omega_min = 188.8*diam**(-0.7081)*2*math.pi/60  # originally: omega_min = 1046.558*diam**(-1.0911)*2*math.pi/60        
        if omega_max == 0:
          omega_max = 793.7*diam**(-0.8504)*2*math.pi/60 # originally: omega_max = 705.406*diam**(-0.8349)*2*math.pi/60    

        # omega
        omega_min_arr = np.ones(len(ws))*omega_min
        omega_max_arr = np.ones(len(ws))*omega_max
        omega = np.minimum(omega_max_arr, np.maximum(omega_min_arr, lambda_opt_0/(diam/2)*ws)) # Eq. 5, all speeds in [rad/s]

        # cp and lambda for the actual beta
        lambdas = np.zeros(len(ws))
        cps = np.zeros(len(ws))
        lambdas[ws>0] = omega[ws>0]*(diam/2)/ws[ws>0] # Eq. 3
        cps[ws>0] = c1[idx]*(c2[idx]/(1.0/((lambdas[ws>0]+c9[idx]*beta)**(-1)-c10[idx]*((beta**3)+1.0)**(-1)))-c3[idx]*beta-c4[idx]*(1.0/((lambdas[ws>0]+c9[idx]*beta)**(-1)-c10[idx]*((beta**3)+1.0)**(-1)))*beta-c5[idx]*beta**x[idx]-c6[idx])*np.exp(-c7[idx]/(1.0/((lambdas[ws>0]+c9[idx]*beta)**(-1)-c10[idx]*((beta**3)+1.0)**(-1))))+c8[idx]*lambdas[ws>0] # Eq. 2 (I) + (II)
      
        # scaling cps wrt cp_max
        cps = cp_max/cp_max_0*cps
        
        # # cp cannot be < 0.
        # cps[cps<0.] = 0.

        return cps
    
    def hbandareas(self,area,Nbands):
        
        radius = (area/math.pi)**(1/2)
        
        area_bands = np.zeros(Nbands)
        area_tmp = 0.
        
        for i in range(int(Nbands/2)):
            
            dist_chord_center = radius-i*(2*radius/ Nbands)

            if Nbands % 2 == 0:
                theta = 2*math.acos(dist_chord_center/radius)        
                area_bands[i] = theta*180/math.pi/360*math.pi*radius**2-(radius*dist_chord_center*(math.sin(theta/2)))-area_tmp  # formula per il calcolo dell'area di un segmento circolare a due basi  
                area_bands[Nbands-1-i] = area_bands[i]
                area_tmp += area_bands[i]
            else:
                if dist_chord_center>0:
                    theta = 2*math.acos(dist_chord_center/radius)        
                    area_bands[i] = theta*180/math.pi/360*math.pi*radius**2-(radius*dist_chord_center*(math.sin(theta/2)))-area_tmp # formula per il calcolo dell'area di un segmento circolare a due basi  
                    area_bands[Nbands-1-i] = area_bands[i]
                    area_tmp += area_bands[i]
                else:
                    area_bands[i] = area-2*(sum(area_bands[0:i-1]))

        return area_bands
    
    def windshear(self,ws,z_hub,z_i,alpha):
        
        ws_shear_i = ws*(z_i/z_hub)**alpha # Eq. 7 wind speed corrected considering vertical wind profile for each h band
        
        return ws_shear_i
    
    def windveer(self,z_hub,z_i,Vu):
        
        deltaphi_i = Vu*(z_i-z_hub)  # Eq. 8
        
        return deltaphi_i
    
    def eqspeed(self,ws,z_i,z_hub,alpha,area,Vu,Nbands):
        
        diam = (4*area/math.pi)**(1/2)
        h = np.zeros(Nbands)
        ws_eq = np.zeros(len(ws))
        
        for i in range(Nbands):
            h[i] = z_hub-diam/2+diam/(Nbands*2)+i*diam/Nbands # height of the barycenter of the various h bands at which compute U_i e DeltaPhi_i
            
        area_bands = self.hbandareas(area,Nbands)
        
        for i in range(Nbands):
            ws_shear_i = self.windshear(ws,z_hub,z_i,alpha)
            deltaphi_i = self.windveer(z_hub,z_i,Vu)
            ws_eq = ws_eq+(area_bands[i]/area*(ws_shear_i*math.cos(deltaphi_i*2*math.pi/360))**3) # Eq. 9 (senza radice cubica)
           
        ws_eq = ws_eq**(1./3.)
        
        return ws_eq


###########################################################################################################################################################

if __name__ == "__main__":
    
    """
    Test
    """

    params = {
                'model': 'detailed',
                'area':  39.6 ,
                'efficiency': 0.45,
                'Prated': 10,
                'WSrated': 11.0,
                'WScutin': 2.5,
                'WScutout': 32.0,
                'omega_min': 50,
                'omega_max': 130,
                'beta': 0.,
                'cp_max': 0.44,
                'idx': 5,
                'z_i': 40.,
                'z_hub': 30.,
                'alpha': 0.25,
                'Vu': 0.5,
                'Nbands': 10
                }

    rho = np.array([1.2,1.2,1.2,1.2])
    windsp = np.array([0.5,3.0,12.0,35.0])

    weather = {
               'rho': rho,
               'wind': windsp
              }
    
    ts = 1.
    
    wind = wind(params,weather,ts)
    
    for h in range(4):
        print(wind.use(h))    
    