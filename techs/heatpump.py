import numpy as np
import constants as c
from scipy.interpolate import interp1d

class heatpump:
    
        def __init__(self,parameters,simulation_hours):
            
            """
            Create heat-pump object
            
            Parameters
            ----------
            parameters : dict
            
                "type": 1 = air-water (other types not yet implemented...)
                
                "nom Pth": float [kW] nominal condition: t_amb=5° t_out=35° 6000 rpm
                
                "t rad heat": float [C°] temperature radiant system in heating mode
                "t rad cool": float [C°] temperatura radiant system in cooling mode
                
                "inertial TES volume": thermal energy storage float [lt]
                "inertial TES dispersion": float [W/m2K]
                               
                "PV surplus": bool # allow to use PV surplus to charge inertial_TES
                "REC surplus": bool # allow to use REC PV surplus to charhe intertial_TES

                
            Returns
            -------
            Air-water HP object able to:
                ...
                
    
            """
        
            self.type = parameters['type'] 
            
            self.nom_Pth = parameters['nom Pth']
            
            self.t_rad_h = parameters['t rad heat']
            #self.t_rad_c = parameters['t rad cool']
                    
            self.PV_surplus = parameters['PV surplus']
            self.REC_surplus = parameters['REC surplus']
            if self.REC_surplus:
                self.PV_surplus = True
            
            self.mode = 1 # 1 = "heat" initial mode, during MESS simulation it is changed to 2 = "cool" when cooling is required
            
            #### inertial TES#################################################
            self.i_TES_volume = parameters['inertial TES volume']
            self.i_TES_dispersion = parameters['inertial TES dispersion'] 
            self.i_TES_mass = self.i_TES_volume # lt -> kg
            self.i_TES_surface = 6 * (self.i_TES_volume/c.H2OADENSITY)**(2/3) # cube surface [m2]
            self.cp = c.CP_WATER  # J/kgK     
            self.cp_kWh = self.cp/3600000 # kWh/kgK        
            self.i_TES_t = self.t_rad_h # initial temperature C°
            
            ### stories #######################################################
            self.i_TES_story = np.zeros(simulation_hours) # T_inertial_TES
            self.satisfaction_story = np.zeros(simulation_hours) # 0 no demand, 1 demand satisfied by iTES, 2 demand satisfied, 3 demand satisfied and iTES_T raised, 4 damand satisfied and iTES_T reaches maximum, -1 unsatisfied demand, -2 unsatisfied demand and iTES under minimum -3 t_amb too cold
            
            self.cop_story = np.zeros(simulation_hours) # coefficient of performance
            self.surplus_story = np.zeros(simulation_hours) # 0 no surplus is used by HP. 1 HP use PV surplus to charge iTES.
            
            #### HP MODEL GU' #################################################
            # danfoss coolselector software available at https://www.danfoss.com/it-it/service-and-support/downloads/dcs/coolselector-2/
            # name='VZH028CH' fluid='R410A
            
            # 6000 rpm. Polynomial model coefficients (electrical and cooling power coefficients)
            self.C_Pele6000 = [0.559950351	,-0.100388322,	0.112208669	,-0.002576625	,0.003712633	,-0.001175774	,-1.77E-05	,3.70E-05	,-2.81E-05	,1.13E-05]
            self.C_Pq6000 = [18.13823604	,0.646313462,	-0.12747235,	0.008771572,	-0.003863774	,-3.45E-05	,4.25E-05	,-5.80E-05	,-1.48E-05	,-5.46E-06]
            #self.C_p3000 = [1.726798023	,-74.34272731	,74.46650928	,-2.079362567	,2.915533504,	-0.934003915,	-0.018207035,	0.036409513,	-0.024701877,	0.008844176]
            #self.C_q3000 = [10442.69118	,374.0860539,	-71.42510762,	4.970490623	,-2.285688038,	-0.161289545	,0.023442599	,-0.029741993	,-0.005992321,	-0.000337462]

            self.tc_max=interp1d([-32,-25,-8,21,27],[35,47,65,65,60]) # t_amb and t_water_out working range
        
            self.pinch_air=10  # pinch air
            self.overH=5       # overheating air
            self.dT_eva= self.pinch_air+ self.overH # dT evaporator
            self.pinch_water=3 #pinch water
            
            ### Regulation [Fahlén P. Capacity control of heat pumps. REHVA J 2012:28–31.]
            x=np.array([0.15 , 0.2, 0.4,  0.5, 0.6, 0.8, 1])
            y=np.array([2.87, 3.2 , 4 ,  4.1, 4 ,  3.65, 3.1])
            p = np.polyfit(x, y/3.1, 4)
            self.f_regulation_Pele = np.poly1d(p)            
            X=np.linspace(0.1,1,1000)
            Y=self.f_regulation_Pele (X)
            self.f_regulation_Pth  = interp1d(X * Y, Y)            
            self.Pth_min_regulation= 0.15*self.f_regulation_Pele(0.15)
            
            # size factor
            self.c_t=0.85  # correcting factor to consider less efficiency at 6000 rpm
            Pth_7_35 = 13.197816888999997*self.c_t #  Pth at nominal condition of the HP used as reference model
            self.size_factor=  self.nom_Pth/ Pth_7_35 
            
        def output(self,C,Te,Tc):
            # polynomial model
            Y = C[0] + C[1]*Te + C[2]*Tc + C[3]*Te**2 + C[4]*Te*Tc + C[5]*Tc**2 + C[6]*Te**3 + C[7]*Tc*Te**2 + C[8]*Te*Tc**2 + C[9]*Tc**3
            return Y
        
        def nominal_performance(self,t_amb,t_w):    
            t_w_eff = t_w # t water out
            Te = t_amb - self.dT_eva # evaporator temperature
            Tc = t_w + self.pinch_water # condenser temperature
            Tc_max = self.tc_max(Te) # max condenser temperature
            if Tc > Tc_max:
                Tc = Tc_max
                t_w_eff = Tc_max - self.pinch_water
            Pele= self.output(self.C_Pele6000,Te,Tc)
            Pth= self.output(self.C_Pq6000,Te,Tc) + Pele
            Pth=Pth*self.c_t  # correcting factor to consider less efficiency at 6000 rpm
            cop=Pth/Pele
            Pele= Pele*self.size_factor
            Pth= Pth*self.size_factor
            return cop,Pth,Pele, t_w_eff
        
        def HP_follows_thermal(self,t_amb,t_w,e_th):
            # heatpump follows thermal demand
            
            cop,Pth,Pele,t_w_eff = self.nominal_performance(t_amb,t_w) # nominal working
            rf = e_th / Pth # regulation factor
            if rf<1:
                if rf<self.Pth_min_regulation:
                    rf = self.Pth_min_regulation
                    Pth= Pth*rf
                else:
                    Pth = e_th
                cop = cop * self.f_regulation_Pth(rf)
                Pele = Pth/cop
            #else nominal            
            
            return cop,Pth,Pele,t_w_eff  
        
        def HP_follows_electricity(self,t_amb,t_w,e_ele):
            # heatpump follows electricyty available insted of thermal demand
            
            cop,Pth,Pele, t_w_eff = self.nominal_performance(t_amb,t_w)
            rf = e_ele / Pele # regulation factor
            if rf<1:
                if rf<0.15:
                    rf = 0.15
                    Pele = Pele*rf
                else:
                    Pele = e_ele
                cop = cop * self.f_regulation_Pth(rf)
                Pth = Pele*cop
            #else nominal       
            
            return cop,Pth,Pele,t_w_eff
             
        def use(self,t_amb,e_th,e_ele,h):
            """
            heat (e_th<0) or cool (e_th>0) required by radiation system to heatpump system

            Parameters
            ----------
            t_amb : float [C°]
            e_th : float [kWh]
                <0 heat
                >0 cool
            e_ele: float [kWh]
                if e_ele > 0 and PV_surplus or REC_surplus == True the HP can use the suprlus of electicity to heats the inertial_TES

            Returns
            -------
            electricity used (-), heat supply (+) or absorbed (-) by the HP, heat supply (+) or absorbed (-) by the inertial_TES

            """
            
            # check mode
            if e_th < 0:
                self.mode = 1 # heat
                
            if e_th > 0:
                self.mode = 2 # cool  
                
            if e_th == 0 and self.mode != 0:
                if np.count_nonzero(self.satisfaction_story[h-48:h]== 0) == 48: 
                    self.mode = 0 # off after 24 hours of inactivity
                
            # initialise
            e_th_i_TES = 0
            e_ele_hp = 0
            e_th_hp = 0
        
            # inertial_TES dispersion (one hour)
            self.i_TES_story[h] = self.i_TES_t
            self.i_TES_t += self.i_TES_dispersion * self.i_TES_surface * (20-self.i_TES_t) * 3600 / (self.i_TES_mass*self.cp)
                
            e_overcharging_tot = 0 # parameter used in PV_surplus working
            e_charging = 0 # parameter that have to be initialise = 0 to be used in PV_surplus working
            ### normal working: heatpump follows thermal demand
            
            if e_th < 0:
                
                # inertial_TES heats radiation system
                if self.i_TES_t > self.t_rad_h:
                    e_th_i_TES = min(-e_th, self.i_TES_mass*self.cp_kWh*(self.i_TES_t-self.t_rad_h))
                    self.i_TES_t += - e_th_i_TES/(self.i_TES_mass*self.cp_kWh)
                    e_th += e_th_i_TES  
                    self.satisfaction_story[h] = 1
                    
                    
                if e_th < 0: # energy inside i_TES is not enough
       
                    ### HP switch-on     
                    
                    ### heat to recharge the i_TES
                    e_charging = self.i_TES_mass*self.cp_kWh*(self.t_rad_h-self.i_TES_t)      
                    cop,Pth,Pele,t_w_eff = self.HP_follows_thermal(t_amb, self.t_rad_h, e_charging)
                    
                    if t_w_eff < self.t_rad_h: # the air temperature is too low to generate water at the required temperature
                        self.satisfaction_story[h] = -3
                        return(0, 0, 0)
                                     
                    if Pth < e_charging: # i_TES can't be charged less than one hour
                        self.satisfaction_story[h] = -2 
                        self.i_TES_t += Pth/(self.i_TES_mass*self.cp_kWh)
                        e_th_i_TES += - Pth
                        e_th_hp = Pth
                        e_ele_hp = Pele
                        
                    
                    else: # i_TES is charge, HP can heats the radiation system
                    
                        self.i_TES_t = self.t_rad_h
                        e_th_i_TES += -e_charging
                        e_th += - e_charging
                    
                        cop,Pth,Pele,t_w_eff = self.HP_follows_thermal(t_amb, self.i_TES_t, -e_th)
                        
                        if Pth < -e_th:
                            self.satisfaction_story[h] = -1
                            e_th_hp = Pth
                            e_ele_hp = Pele
                            
                        if Pth == -e_th:
                            self.satisfaction_story[h] = 2
                            e_th_hp = Pth
                            e_ele_hp = Pele
                            
                        if Pth > -e_th:
                            
                            # HP heats i_TES   
                            self.satisfaction_story[h] = 3
                            j = 0 # used to calculate number of switch on number
                            e_th_r = e_th # demand ramained to satisfied, used for the while cycle
                            while e_th_r < -0.00000001:
                                
                                cop,Pth,Pele,t_w_eff = self.HP_follows_thermal(t_amb, self.i_TES_t, -e_th)
                                
                                if t_w_eff == self.i_TES_t: # maximum temperature not still reached 
                                    e_th_hp += Pth/60
                                    e_ele_hp += Pele/60
                                    e_overcharging = (Pth+e_th)/60
                                    e_th_r += -e_th/60
                                
                                    if e_overcharging > 0: # charge iTES if there is energy to do it
                                        e_th_i_TES += -e_overcharging
                                        e_overcharging_tot += e_overcharging # used in PV_surplus working
                                        self.i_TES_t += e_overcharging/(self.i_TES_mass*self.cp_kWh)
                                        
                                else: # maximum temperaeture reached: used energy inside iTES to satisfied demand
                                    self.satisfaction_story[h] = 4+j 
                                    j += 1
                                    e_th_available = min(-e_th_r, self.i_TES_mass*self.cp_kWh*(self.i_TES_t-self.t_rad_h))
                                    e_th_i_TES += e_th_available
                                    self.i_TES_t += - e_th_available/(self.i_TES_mass*self.cp_kWh) 
                                    e_th_r += e_th_available
                                    if e_th_r < 0:
                                        self.satisfaction_story[h] = 4+j
                                        e_th = e_th_r
                        
            ### PV_surplus working: heatpump follows available electricity insted of thermal demand
            if self.PV_surplus and e_ele > e_ele_hp and self.satisfaction_story[h] in [0,1,2,3] and self.mode != 0:
                
                self.surplus_story[h] = 1
                
                # reinitialise balances                
                e_ele_hp = 0
                e_th_hp = 0
                e_th_i_TES += e_overcharging_tot
                self.i_TES_t += - e_overcharging_tot/(self.i_TES_mass*self.cp_kWh)
                
                # HP heats i_TES more than the noral working
                
                if e_th < 0: # demand
                    j = 0 # used to calculate number of switch on number
                    e_th_r = e_th # demand ramained to satisfied, used for the while cycle
                    while e_th_r < -0.00000001:
                        
                        cop,Pth,Pele,t_w_eff = self.HP_follows_electricity(t_amb, self.i_TES_t, e_ele)
                        
                        if t_w_eff == self.i_TES_t: # maximum temperature not still reached 
                            
                            e_th_hp += Pth/60
                            e_ele_hp += Pele/60
                            e_overcharging = (Pth+e_th)/60
                            e_th_r += -e_th/60
                        
                            if e_overcharging > 0: # charge iTES if there is energy to do it
                                e_th_i_TES += -e_overcharging
                                self.i_TES_t += e_overcharging/(self.i_TES_mass*self.cp_kWh)
                                               
                        else: # maximum temperaeture reached: used energy inside iTES to satisfied demand
                            self.satisfaction_story[h] = 4+j 
                            j += 1
                            e_th_available = min(-e_th_r, self.i_TES_mass*self.cp_kWh*(self.i_TES_t-self.t_rad_h))
                            e_th_i_TES += e_th_available
                            self.i_TES_t += - e_th_available/(self.i_TES_mass*self.cp_kWh) 
                            e_th_r += e_th_available
                            if e_th_r < 0:
                                self.satisfaction_story[h] = 4+j
                                e_th = e_th_r
                                
                elif e_th == 0:
                    for m in np.arange(60):
                        cop,Pth,Pele,t_w_eff = self.HP_follows_electricity(t_amb, self.i_TES_t, e_ele)
                        
                        if t_w_eff == self.i_TES_t: # maximum temperature not still reached 
                            
                            e_th_hp += Pth/60
                            e_ele_hp += Pele/60
                            e_overcharging = (Pth+e_th)/60
                        
                            if e_overcharging > 0: # charge iTES if there is energy to do it
                                e_th_i_TES += -e_overcharging
                                self.i_TES_t += e_overcharging/(self.i_TES_mass*self.cp_kWh)
                                               
                        else: # maximum temperaeture reached: used energy inside iTES to satisfied demand
                            break
                             
            if e_th_hp > 0:
                self.cop_story[h] = e_th_hp/e_ele_hp
                            
            return(-e_ele_hp, e_th_hp, e_th_i_TES)
                                        
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

            size = self.nom_Pth
            
            if tech_cost['cost per unit'] == 'default price correlation':
                C0 = 1400 # €/kW
                scale_factor = 0.8 # 0:1
                C = size * C0 **  scale_factor
            else:
                C = size * tech_cost['cost per unit']
    
            tech_cost['total cost'] = tech_cost.pop('cost per unit')
            tech_cost['total cost'] = C
            tech_cost['OeM'] = tech_cost['OeM'] *C /100 # €

    
            self.cost = tech_cost    