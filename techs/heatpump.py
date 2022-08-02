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
                "nom Pth": float [kW]
                "nom Tamb": float [C°]
                "nom Tout" float [C°]
                
                "t max": float [C°]
                "t min": float [C°]
                "t rad heat": float [C°]
                "t rad cool": float [C°]
                
                "tank volume": float [lt]
                
            Returns
            -------
            Air-water HP object able to:
                supply (heating) or extract (cooling) thermal energy abrosrbing electricity .use()
                absorb electricity surplus heating or cooling its tank .use2()
    
            """
            
            # create a new hp object based on hplib database
            
            par = hpl.get_parameters(model='Generic',  # generic model created on the basis of a database of real models
                                     group_id=parameters['type'],       # air-water
                                     t_in=parameters["nom Tamb"],
                                     t_out=parameters["nom Tout"], 
                                     p_th=parameters["nom Pth"]*1000) # kW -> W
            
            self.hp = hpl.HeatPump(par)       
            self.mode = "heat" # initial mode, during MESS simulation it is changed to "cool" when cooling is required
           
            self.t_max = parameters['t max'] 
            self.t_min = parameters['t min']
            self.t_rad_heat = parameters['t rad heat'] # C°
            self.t_rad_cool = parameters['t rad cool'] # C°
            # self.hp.delta_t = 5 # delta_t = 5 is the ref value: it could be changed
            
            ######### partial load regulation https://www.rehva.eu/rehva-journal/chapter/capacity-control-of-heat-pumps-full-version
        
            # modern
            #x=np.array([0.1 , 0.2,  0.5, 0.735,  1])
            #y=np.array([4.85 , 5.3,  4.9,  4  ,  3.1])
            
            # state of art
            x=np.array([0.1 , 0.2, 0.4,  0.5, 0.6, 0.8, 1]) # sperimental profiles
            y=np.array([2.45, 3.2 , 4 ,  4.1, 4 ,  3.65, 3.1]) # sperimental profiles
            
            p = np.polyfit(x, y/3.1, 4)
            self.regulation = np.poly1d(p)
            self.max_regulation = 0.1
        
            ######################################################## inertial tank
            self.tank_volume = parameters['tank volume'] # lt
            self.tank_mass = self.tank_volume # kg
            self.cp = 4200 # J/kgK                    
            self.tank_t = parameters['t rad heat'] # initial temperature C°
            
            self.story = [] # on, off, on/off
        
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
        
        def use(self,t_amb,e_th):
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
            
            ##########################################################################################
            if e_th < 0: # heating (-)
                
                self.mode = "heat"
                
                ### first use heat inside the tank
                t = self.tank_t + (e_th*3600*1000) / (self.tank_mass*self.cp) # J / kg 
                
                ### energy in the tank is sufficient to met demand, hp can stay off
                if t > self.t_rad_heat: 
                    self.tank_t = t
                    e_tank = -e_th # (+) tank to radiation
                    e_el = 0
                    e_hp = 0
                    self.story.append('off')
                    
                ### tank riched the minimum temperature (irradiation system temperature), so heat pump must be turned on
                else: 
                    e_tank = self.tank_mass*self.cp*(self.tank_t-self.t_rad_heat)/3600000 # heat that can be provided by the tank (+)
                    self.tank_t = self.t_rad_heat # minimum temperature riched
                    
                    ### switch on hp
                    e_hp = e_th + e_tank # (-)
                    
                    res = self.hp.simulate(t_in_primary=t_amb, t_in_secondary=self.tank_t, t_amb=t_amb, mode=1)
                    cop = res['COP']
                    nom_e_th = res['P_th']/1000 # (+)
                    nom_e_el = res['P_el']/1000 # (+)
                    
                    # unfulfilled demand (hp nominal)
                    if -e_hp > nom_e_th:
                        e_hp = nom_e_th # (+)
                        e_el = - nom_e_el # (-)
                        print("unfulfilled demand")
                        self.story.append('on')
                        
                    # partial load: heat to the radiation system
                    else: 
                        if (-e_hp/nom_e_th) > self.max_regulation:
                            cop = cop*self.regulation(-e_hp/nom_e_th) 
                            self.story.append('on')
                        else:
                            cop = cop*self.regulation(self.max_regulation)
                            self.story.append('on/off')
                        e_el = e_hp / cop # (-)
                        e_hp = - e_hp # (+) hp to radiaiton
                        
            ######################################################################################
            elif e_th > 0: # cooling (+)
                self.mode = "cool"
                
                ### first use the cold water inside the tank (heat the tank)
                t = self.tank_t + (e_th*3600*1000) / (self.tank_mass*self.cp) # J / kg 
                
                ### energy in the tank is sufficient to met demand, hp can stay off
                if t < self.t_rad_cool: 
                    self.tank_t = t
                    e_tank = -e_th # (-) # tank to radiation
                    e_el = 0
                    e_hp = 0
                    self.story.append('off')
                    
                ### tank riched the maximum temperature (irradiation system temperature), so heat pump must be turned on
                else: 
                    e_tank = self.tank_mass*self.cp*(self.tank_t-self.t_rad_cool)/3600000 # heat that can be absorbed by the tank (-)
                    self.tank_t = self.t_rad_cool # maximum temperature riched
                    
                    ### switch on hp
                    e_hp = e_th + e_tank # (+)
                    
                    res = self.hp.simulate(t_in_primary=t_amb, t_in_secondary=self.tank_t, t_amb=t_amb, mode=2)
                    eer = res['EER']
                    nom_e_th = - res['P_th']/1000 # (-) -> (+)
                    nom_e_el = res['P_el']/1000 # (+)
                    
                    # unfulfilled demand (hp nominal)
                    if e_hp > nom_e_th:
                        e_hp = nom_e_th # (+)
                        e_el = - nom_e_el # (-)
                        print("unfulfilled demand")
                        self.story.append('on')
                        
                    # partial load: heat to the radiation system
                    else: 
                        if (e_hp/nom_e_th) > self.max_regulation:
                            eer = eer*self.regulation(e_hp/nom_e_th) 
                            self.story.append('on')
                        else: 
                            eer = eer*self.regulation(self.max_regulation)
                            self.story.append('on/off')
                        e_el = - e_hp / eer # (-)
                        e_hp = - e_hp # (-) hp to radiaiton
                                    
                    
            elif e_th == 0: # switched off
                self.story.append('off')
                e_el = 0
                e_tank = 0
                e_hp = 0
            
            return(e_el, e_tank, e_hp)
        
        
        def use2(self,t_amb,e_el):
            """
            there is electricity surplus so hp can be turn on to heat or cool the tank

            Parameters
            ----------
            t_amb : float [C°]
            e_el : float [kWh] (+)

            Returns
            -------
            electricity used (-), heat supply (+) or absorbed (-) by the tank, heat supply (+) or absorbed (-) by the hp, 

            """
            
            if e_el > 0:
            
                if self.mode == "heat":
                    
                    res = self.hp.simulate(t_in_primary=t_amb, t_in_secondary=self.tank_t, t_amb=t_amb, mode=1)
                    nom_cop = res['COP']
                    nom_e_el = res['P_el']/1000 # (+)
                    nom_e_th = res['P_th']/1000 # (+)
                    
                    
                    
                    ## quanto regolo???
                    # onoff < max_regulation*nom
                    # on == e_el
                    # on == tank_full
                    
                    
                    ### can hp absorb all that electricity?
                    if e_el > nom_e_el: # no
                        e_el = nom_e_el
                       
                    ### partial load (include nominal case)
                    if e_el < nom_e_el*self.max_regulation:
                        cop = nom_cop * self.regulation(self.max_regulation)
                        self.story.append("on/off")
                    else:       
                         cop = nom_cop * self.regulation (e_el/nom_e_el)    
                         self.story.append("on")
                       
                    e_th = e_el * cop
                    
                    ### can tank absorb all that energy?
                    t = self.tank_t + (e_th*3600*1000) / (self.tank_mass*self.cp) # J / kg 
                    
                    if t < self.t_max: # yes
                        self.tank_t = t
                        
                    else: # no
                        e_th = self.tank_mass*self.cp*(self.t_max-self.tank_t)/3600000 # heat that can be absorbed by the tank (+)
                        self.tank_t = self.t_max
                        
                        if e_th < nom_e_th*self.max_regulation:
                            cop = nom_cop*self.regulation(self.max_regulation)
                            self.story[-1] = "on/off"
                        else:
                            cop = nom_cop * self.regulation (e_th/nom_e_th)
                        e_el = e_th / cop
                        ## controlla self.story poi fai cool, poi location
                    return(-e_el,-e_th,e_th)
                
                ###################################################################################
                if self.mode == "cool":
                    
                    res = self.hp.simulate(t_in_primary=t_amb, t_in_secondary=self.tank_t, t_amb=t_amb, mode=2)
                    nom_eer = res['EER']
                    nom_e_el = res['P_el']/1000 # (+)
                    nom_e_th = - res['P_th']/1000 # (+)
                    
                    ### can hp absorb all that electricity?
                    if e_el > nom_e_el: # no
                        e_el = nom_e_el
                        
                    ### partial load (include nominal case)
                    if e_el < nom_e_el*self.max_regulation:
                        eer = nom_eer * self.regulation(self.max_regulation)
                        self.story.append("on/off")
                    else:       
                         eer = nom_eer * self.regulation (e_el/nom_e_el)    
                         self.story.append("on")
                        
                    e_th = e_el * eer
                    
                    ### can tank absorb all that energy?
                    t = self.tank_t - (e_th*3600*1000) / (self.tank_mass*self.cp) # J / kg 
                    
                    if t > self.t_min: # yes
                        self.tank_t = t
                        
                    else: # no
                        e_th = self.tank_mass*self.cp*(self.tank_t-self.t_min)/3600000 # heat that can be absorbed by the tank (+)
                        self.tank_t = self.t_min
                        
                        if e_th/nom_e_th < self.max_regulation:
                              eer = nom_eer * self.regulation (self.max_regulation)
                              self.story[-1] = "on/off"
                        else:
                            eer = nom_eer * self.regulation (e_th/nom_e_th)
                        e_el = e_th / eer
                    
                    return(-e_el,e_th,-e_th)
                
            else:
                self.story.append('off')
                return(0,0,0)

              
                
###########################################################################################################################################################

if __name__ == "__main__":
    
    """
    Functional test
    """
    
    
    parameters = {"type":1,
                  "nom Pth":10,
                  "nom Tamb":10,
                  "nom Tout":60,
                  "t max": 70,
                  "t min": 5,
                  "t rad heat": 45,
                  "t rad cool": 15,                
                  "tank volume": 300
                  }
    
    hp = heatpump(parameters)
    
    hp.partial_load_graph()
    
    hp.tank_t = 70
    
    hp.tank_t = 50
    hp.mode = "cool"
    
    hp.tank_t = 15
    hp.use2(35,1)
    hp.tank_t
    
    

    

    

    
    


