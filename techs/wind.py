import numpy as np
import math
import os
import sys 
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),os.path.pardir)))   # temorarily adding constants module path 
import constants as c

class wind:    
    
    def __init__(self,params,ts=1):
        """
        Create a wind object based on the specified model
    
        params : dictionary
            'model': str type of model to be used for wind
                    betz -> simple model based on Betz theory
                    detailed -> more detailed model based on
                    Saint-Drenan, Yves-Marie, et al. 
                    "A parametric model for wind turbine power curves incorporating environmental conditions." 
                    Renewable Energy 157 (2020): 754-768.
            
            'area': float swept area [m2] e.g. 39.6 m^2 (Aircon 10/10 kW)
            'efficiency': float total efficiency = Betz*efficiency [-] default: 0.45 (ca. 0.593*0.76, i.e. Betz*efficiency)
            'Prated': float rated power [kW] # NOTE: useless
            'WSrated': float rated wind speed [m/s] e.g. 11.0 m/s (Aircon 10/10 kW)
            'WScutin': float cut in wind speed [m/s] e.g.  2.5 m/s (Aircon 10/10 kW)
            'WScutoff': float cut off wind speed [m/s] e.g. 32.0 m/s (Aircon 10/10 kW)
            
            'omega_min': [rpm] OPTIONAL default = from eq.
            'omega_max': [rpm] OPTIONAL default = from eq.
            'beta': [°] e.g. 0°
            'cp_max': [-] OPTIONAL default = 0.44; values from 0.4 to 0.5
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
        
        self.params = params
        self.ts = ts
        
        if params['model'] not in ['betz','detailed']:
            print('Error: wrong wind turbine model')
        

        
    def use(self,h,ws_input,rho=1.225,ts=1):
        """
        Produce electricity
        
        h : int hour to be simulated
        windspeed: float wind speed [m/s]
        rho: density [kg/m3]
    
        output : float electricity produced that hour [kWh]
    
        """
        
        ws_turbine=ws_input
        
        if self.params['model'] == 'betz':
                        
            if ws_turbine<self.params['WScutin']:
                ws_turbine = 0.
            elif self.params['WSrated']<ws_turbine<self.params['WScutout']:
                ws_turbine = self.params['WSrated']
            elif ws_turbine>self.params['WScutout']:
               ws_turbine = 0.
            
            self.energy = 0.5*rho*self.params['area']*ws_turbine**3*self.params['efficiency']/1000.*ts # kWh
            
        elif self.params['model'] == 'detailed':
            
            if ws_turbine<self.params['WScutin']:
                ws_turbine = 0.
            elif self.params['WSrated']<ws_turbine<self.params['WScutout']:
                ws_turbine = self.params['WSrated']
            elif ws_turbine>self.params['WScutout']:
               ws_turbine = 0.
            
            cp_max = 0.44
            
            if 'cp_max' in self.params:
                cp_max = self.params['cp_max']
                
            powercoeff = self.cpfunc(ws_turbine,self.params['area'],self.params['beta'],self.params['idx'],cp_max=cp_max)
            wseq = self.eqspeed(ws_turbine,self.params['z_i'],self.params['z_hub'],self.params['alpha'],self.params['area'],self.params['Vu'],self.params['Nbands'])
            
            self.energy = 0.5*rho*self.params['area']*wseq**3*powercoeff/1000.*ts # kWh
            
        return(self.energy)


    def cpfunc(self,ws,area,beta,idx,cp_max=0.44):
        
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
        if 'omega_min' in self.params:
            omega_min = self.params['omega_min']
        else:
            omega_min = 1046.558*diam**(-1.0911)*2*math.pi/60
            # omega_min = 188.8*diam**(-0.7081)*2*math.pi/60 # alternative equation proposed by Niccolò Baldi  
        
        if 'omega_max' in self.params:
            omega_max = self.params['omega_max']
        else:
            omega_max = 705.406*diam**(-0.8349)*2*math.pi/60
            # omega_max = 793.7*diam**(-0.8504)*2*math.pi/60 # alternative equation proposed by Niccolò Baldi

        omega = np.minimum(omega_max, np.maximum(omega_min, lambda_opt_0/(diam/2)*ws)) # Eq. 5, all speeds in [rad/s]

        # cp and lambda for the actual beta
        lambdaparam = 0.
        cp = 0
        if not(ws == 0.):
            lambdaparam = omega*(diam/2)/ws # Eq. 3
            cp = c1[idx]*(c2[idx]/(1.0/((lambdaparam+c9[idx]*beta)**(-1)-c10[idx]*((beta**3)+1.0)**(-1)))-c3[idx]*beta-c4[idx]*(1.0/((lambdaparam+c9[idx]*beta)**(-1)-c10[idx]*((beta**3)+1.0)**(-1)))*beta-c5[idx]*beta**x[idx]-c6[idx])*np.exp(-c7[idx]/(1.0/((lambdaparam+c9[idx]*beta)**(-1)-c10[idx]*((beta**3)+1.0)**(-1))))+c8[idx]*lambdaparam # Eq. 2 (I) + (II)
        
        # scaling cps wrt cp_max
        cp = cp_max/cp_max_0*cp

        return cp
    
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
        ws_eq = 0.
        
        for i in range(Nbands):
            h[i] = z_hub-diam/2+diam/(Nbands*2)+i*diam/Nbands # height of the barycenter of the various h bands at which compute U_i e DeltaPhi_i
            
        area_bands = self.hbandareas(area,Nbands)
        
        for i in range(Nbands):
            ws_shear_i = self.windshear(ws,z_hub,z_i,alpha)
            deltaphi_i = self.windveer(z_hub,z_i,Vu)
            ws_eq = ws_eq+(area_bands[i]/area*(ws_shear_i*math.cos(deltaphi_i*2*math.pi/360))**3) # Eq. 9 (senza radice cubica)
           
        ws_eq = ws_eq**(1./3.)
        
        return ws_eq

    def tech_cost(self,tech_cost):
        """
        Parameters
        ----------
        tech_cost : dict
            'cost per unit': float [€/kW]
            'OeM': float, percentage on initial investment [%]
            'refud': dict
                'rate': float, percentage of initial investment which will be rimbursed [%]
                'years': int, years for reimbursment
            'replacement': dict
                'rate': float, replacement cost as a percentage of the initial investment [%]
                'years': int, after how many years it will be replaced

        Returns
        -------
        self.cost: dict
            'total cost': float [€]
            'OeM': float, percentage on initial investment [%]
            'refud': dict
                'rate': float, percentage of initial investment which will be rimbursed [%]
                'years': int, years for reimbursment
            'replacement': dict
                'rate': float, replacement cost as a percentage of the initial investment [%]
                'years': int, after how many years it will be replaced
        """
        tech_cost = {key: value for key, value in tech_cost.items()}

        size = self.param['Prated'] # KW
        
        if tech_cost['cost per unit'] == 'default price correlation':
            C0 = 1270 # €/kW
            scale_factor = 0.8 # 0:1
            C = size * C0 **  scale_factor
        else:
            C = size * tech_cost['cost per unit']

        tech_cost['total cost'] = tech_cost.pop('cost per unit')
        tech_cost['total cost'] = C
        tech_cost['OeM'] = tech_cost['OeM'] *C /100 # €


        self.cost = tech_cost

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
                'beta': 0.,
                'idx': 5,
                'z_i': 40.,
                'z_hub': 30.,
                'alpha': 0.25,
                'Vu': 0.5,
                'Nbands': 10
                }

    # windsp = np.array([0.5,3.0,12.0,35.0])
    windsp = np.array([11.0])    
    
    wind = wind(params)
    
    for h in range(len(windsp)):
        print(wind.use(h,windsp[h]))    
    