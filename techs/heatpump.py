import numpy as np
import matplotlib.pyplot as plt
from hplib import hplib as hpl # https://github.com/RE-Lab-Projects/hplib

class heatpump:
    
        def __init__(self,parameters):
            """
            Create heat-pump object
            
            Parameters
            ----------
            parameters : dict
            
                "type": 1 = air-water (other types not yet implemented...)
                
                "usage": 1 = heat and cool, 2 = heat and dhw, 3 = dhw 
                    n.b "cool" and "dhw" can't be supplied by the same hp
                    if usage = 2 an 'heatpump_boiler' (usage=3) can be used to supply dhw 
                        
                "nom Pth": float [kW]
                "nom Tamb": float [C°]
                "nom Tout" float [C°]
                
                "t max hp": float [C°]
                "t max res": float [C°]
                "t min hp": float [C°]
                
                "t min rad heat": float [C°]
                "t max rad cool": float [C°]
                
                "tank volume": float [lt]
                "tank dispersion": float [W/m2K]
                
                "regulation": bool 
                "variable set point": false or 24h list
                
            Returns
            -------
            Air-water HP object able to:
                supply (heating) or extract (cooling) thermal energy abrosrbing electricity .use()
    
            """
            
            # create a new hp object based on hplib database
            
            par = hpl.get_parameters(model='Generic',  # generic model created on the basis of a database of real models
                                     group_id=parameters['type'],       # air-water
                                     t_in=parameters["nom Tamb"],
                                     t_out=parameters["nom Tout"], 
                                     p_th=parameters["nom Pth"]*1000) # kW -> W
            
            self.hp = hpl.HeatPump(par)       
            self.usage = parameters['usage']
            self.mode = 1 # 1 = "heat" initial mode, during MESS simulation it is changed to 2 = "cool" when cooling is required
            self.switch = "stand by" # on / off / standy by
            
            self.t_max_hp = parameters['t max hp']
            self.t_max_res = parameters['t max res']
            self.t_min_hp = parameters['t min hp']
            self.t_min_rad_h = parameters['t min rad heat']
            self.t_max_rad_c = parameters['t max rad cool']
            
            self.regulation = parameters['regulation']
            self.set_point = parameters['set point']*np.ones(8760)
            self.start_point = self.t_min_rad_h*np.ones(8760)
            self.PV_surplus = parameters['PV surplus'] # False or float °C
            
            ######### partial load regulation https://www.rehva.eu/rehva-journal/chapter/capacity-control-of-heat-pumps-full-version
            if self.regulation:
                # modern
                #x=np.array([0.1 , 0.2,  0.5, 0.735,  1])
                #y=np.array([4.85 , 5.3,  4.9,  4  ,  3.1])
                
                # state of art
                x=np.array([0.1 , 0.2, 0.4,  0.5, 0.6, 0.8, 1]) # sperimental profiles
                y=np.array([2.45, 3.2 , 4 ,  4.1, 4 ,  3.65, 3.1]) # sperimental profiles
                
                p = np.polyfit(x, y/3.1, 4)
                self.regulate = np.poly1d(p)
        
            ######################################################## inertial tank
            self.tank_volume = parameters['tank volume'] # lt
            self.tank_dispersion = parameters['tank dispersion']
            self.tank_mass = self.tank_volume # lt -> kg
            self.tank_surface = 6 * (self.tank_volume/1000)**(2/3) # cube surface [m2]
            self.cp = 4187 # J/kgK                    
            self.tank_t = self.t_min_rad_h # initial temperature C°
                      
            self.hp_story = [] # P_th
            self.tank_story = [] # T_tank
            self.satisfied_control = []
        
        def partial_load_graph(self):
            """
            This function is not used in MESS simulations
    
            Returns
            -------
            Graph of COP variation at partial loads
    
            """
            X=np.linspace(0.1,1)
            Y=self.regulation(X)
            plt.figure(dpi=200)
            plt.plot(X,Y)
            plt.ylim(0,2)
            plt.grid()
            plt.xlabel('$ \dfrac{load}{load_{nominal}}$')
            plt.ylabel('$ \dfrac{COP}{COP_{nominal}}$')
            plt.title("Partial load regulation")
            plt.show()
            
        
        def use(self,t_amb,e_th,h):
            """
            heat (e_th<0) or cool (e_th>0) required by radiation system to heatpump system

            Parameters
            ----------
            t_amb : float [C°]
            e_th : float [kWh]
                <0 heat
                >0 cool

            Returns
            -------
            electricity used (-), heat supply (+) or absorbed (-) by the tank, heat supply (+) or absorbed (-) by the hp, 

            """
            
            # check mode
            if e_th < 0:
                self.mode = 1 # heat
                if self.switch == "off": 
                    self.switch = "stand by" # ready to be used after tank temperature control
                
            if e_th > 0:
                self.mode = 2 # cool  
                if self.switch == "off":
                    self.switch = "stand by" # ready to be used after tank temperature control
                
            # after 24 hours of inactivity, the machine shuts down until heat or cool is again required from the radiant system
            if e_th == 0 and self.hp_story[len(self.hp_story)-24*60:len(self.hp_story)] == list(np.zeros(24*60)):
                self.switch = "off" 
 
            e_th = e_th*3600000 # kWh/h -> J/h
            p_th = e_th/60 # J/h = J/min
                           
            # initialise hourly balances    
            e_ele_tot = 0
            e_th_hp_tot = 0
            e_th_tank_tot = 0
            
            # satisfaction control
            if e_th != 0: # if energy is required
                if self.t_max_rad_c+5 < self.tank_t < self.t_min_rad_h-5:
                    self.satisfied_control.append(1)
                else: 
                    self.satisfied_control.append(0)
            else:
                self.satisfied_control.append(0)
            
            # time step 1 minute   
            for ts in np.arange(60):
                
                # tank dispersion Q = UAdT
                self.tank_t += self.tank_dispersion * self.tank_surface * (20-self.tank_t) * 60 / (self.tank_mass*self.cp)
                                             
                # hp heats/cools tank
                if self.switch == "on":         
                    res = self.hp.simulate(t_in_primary=t_amb, t_in_secondary=self.tank_t, t_amb=t_amb, mode=self.mode)
                    cop = res['COP']
                    eer = res['EER']
                    e_th_hp = res['P_th']*60 # J/min heat (+) cool (-)
                    e_el_hp = res['P_el']*60 # J/min (+)
                    
                    # regulation (cannot adjust immediately after starting and when the tank temperature is at the limit)
                    # 0.5 is the point of max cop
                    if self.regulation and self.tank_t > self.t_min_rad_h+5:
                        cop = cop*self.regulate(0.5) 
                        e_el_hp = e_el_hp*0.5
                        e_th_hp = e_el_hp * cop
                        
                    self.tank_t += e_th_hp / (self.tank_mass*self.cp) # J / kg * (J / kg K)
                    
                else:
                    e_th_hp = 0
                    e_el_hp = 0
                    
                # tank heats/cools radiation system
                if (e_th < 0 and self.mode == 1 and self.tank_t > self.t_min_rad_h-5) or (e_th > 0 and self.mode == 2 and self.tank_t < self.t_max_rad_c+5): 
                    self.tank_t += p_th / (self.tank_mass*self.cp)
                    e_th += - p_th
                    
                # save balances (W)
                self.hp_story.append(e_th_hp/60)
                e_ele_tot += e_el_hp
                e_th_hp_tot += e_th_hp
                e_th_tank_tot += - p_th
                
                # tank temperature control
                if self.switch != "off":
                    if self.mode == 1: # heating
                        if self.tank_t <= self.start_point[h]:
                            self.switch = "on"
                        if self.tank_t >= self.set_point[h]:
                            self.switch = "stand by"
                            
                    if self.mode == 2: # cooling
                        if self.tank_t >= self.t_max_rad_c:
                            self.switch = "on"
                        if self.tank_t <= self.t_min_hp:
                            self.switch = "stand by"  
            
            # J -> kWh
            self.tank_story.append(self.tank_t)
            e_ele_tot = e_ele_tot/3600000
            e_th_hp_tot = e_th_hp_tot/3600000
            e_th_tank_tot = e_th_tank_tot/3600000
            
            return(-e_ele_tot, e_th_hp_tot, e_th_tank_tot)
                                        
            
        def use_surplus(self,t_amb,e_th,h):
            #print('surplus switch')
            # clean hystory
            self.tank_story = self.tank_story[:-1]
            self.tank_t = self.tank_story[-1]
            self.hp_story = self.hp_story[:-60]
            self.satisfied_control = self.satisfied_control[:-1]
            
            # change set point 
            self.set_point[h] = self.PV_surplus
            #self.switch = "on"
            
            # resimulate
            a = self.use(t_amb,e_th,h)
            return(a)
                
                
                
            
            
#%% ##########################################################################################################################################################

if __name__ == "__main__":
    
    """
    Functional test
    """
    

    
