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
        self.mode = "heat" # initial mode, during MESS simulation it is changed to "cool" when cooling is required
        
        # self.hp.delta_t = 5 # delta_t = 5 is the ref value: it could be changed
        
        ### partial load regulation https://www.rehva.eu/rehva-journal/chapter/capacity-control-of-heat-pumps-full-version
        
# =============================================================================
#         'modern'
#         x=np.array([0.1 , 0.2,  0.5, 0.735,  1])
#         y=np.array([4.85 , 5.3,  4.9,  4  ,  3.1])
# =============================================================================

        "state of art"
        x=np.array([0.1 , 0.2, 0.4,  0.5, 0.6, 0.8, 1]) # sperimental profiles
        y=np.array([2.45, 3.2 , 4 ,  4.1, 4 ,  3.65, 3.1]) # sperimental profiles
        p = np.polyfit(x, y/3.1, 4)
        self.regulation = np.poly1d(p)
        self.max_regulation = 0.1 ### max partial load 
        
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
        
    def use(self,t_water_in,t_amb,mod2,e):
        """
        

        Parameters
        ----------
        t_water_in : TYPE
            DESCRIPTION.
        t_amb : TYPE
            DESCRIPTION.
        mode : TYPE
            DESCRIPTION.
        mod2 : TYPE
            DESCRIPTION.
        e : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """

        if self.mode == "heat":
            
            res = self.hp.simulate(t_in_primary=t_amb, t_in_secondary=t_water_in, t_amb=t_amb, mode=1)
            
            cop = res['COP']
            nom_P_th = res['P_th']/1000 
            nom_P_el = res['P_el']/1000
            
            if mod2 == "request": # heat requested, houw much electricity is needed? (heat -> ele)
                P_th = -e
                
                if P_th < nom_P_th * self.max_regulation: # swith off
                    P_th = 0
                    P_el = 0
                
                elif P_th > nom_P_th: # nominal
                    P_th = nom_P_th
                    P_el = nom_P_el
                    
                else: # partial load
                    cop = cop*self.regulation(P_th/nom_P_th) 
                    P_el = P_th / cop
                    
                
            if mod2 == "surplus": # surplus of electricity can be used to produce heat (ele -> heat)
                P_el = e
                
                if P_el < nom_P_el*self.max_regulation: # switch off
                    P_th = 0
                    P_el = 0
                
                elif P_el > nom_P_el: # nominal
                    P_th = nom_P_th
                    P_el = nom_P_el
                    
                else: # partial load
                    cop = cop*self.regulation(P_el/nom_P_el)
                    P_th = P_el * cop  
                          
            return(-P_el,P_th,cop)
            
            
        if self.mode == "cool":
            
            res = self.hp.simulate(t_in_primary=t_amb, t_in_secondary=t_water_in, t_amb=t_amb, mode=2)
            
            eer = res['EER']
            nom_P_th = res['P_th']/1000
            nom_P_el = res['P_el']/1000
            
            if mod2 == "request": # cool requested, houw much electricity is needed? (cool -> ele)             
                P_th = e
                
                if -P_th < nom_P_th*self.max_regulation: # switch off
                    P_th = 0
                    P_el = 0
                    
                elif -P_th > nom_P_th: # nominal
                    P_th = nom_P_th
                    P_el = nom_P_el
                    
                else: # partial load
                    eer = eer*self.regulation(-P_th/nom_P_th)
                    P_el = P_th / eer
                    

            if mod2 == "surplus":  # surplus of electricity can be used to produce cool (ele -> cool)
                P_el = e
                
                if P_el < nom_P_el*self.max_regulation: # switch off
                    P_th = 0
                    P_el = 0
                    
                elif P_el > nom_P_el: # nominal
                    P_th = nom_P_th
                    P_el = nom_P_el
                    
                else: # partial load
                    eer = eer*self.regulation(P_el/nom_P_el)
                    P_th = P_el * eer     
                
            return(-P_el,-P_th,eer)
            
            


###########################################################################################################################################################

if __name__ == "__main__":
    
    """
    Functional test
    """
    
    
    parameters = {"type":1,
                  "nom Pth":10,
                  "nom Tamb":5,
                  "nom Tout":60}
    
    hp = heatpump(parameters)
    
    print(hp.use(35,10,"request",-6))
    print(hp.use(36,10,"request",-6))
    
    print(hp.use(35,10,"surplus",2))
    print(hp.use(36,10,"surplus",2))
    
    hp.mode = "cool"

    print(hp.use(12,30,"request",6))
    print(hp.use(10,30,"request",6))
    
    print(hp.use(12,30,"surplus",2))
    print(hp.use(10,30,"surplus",1))
            
    hp.partial_load_graph()
    
            