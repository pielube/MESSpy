import numpy as np
from scipy.interpolate import interp1d

class inverter:

    def __init__(self,parameters,simulation_hours):
        """
        

        Parameters
        ----------
        max efficiency : float 0-1
        number : int number of inverters
        peakP : float peak power of a single inverter

        Returns
        -------
        inverter object

        """
        
        self.eta_max = parameters['max efficiency']
        self.n = parameters['number']
        self.peakP = parameters['peakP']
        self.eta_story = np.zeros(simulation_hours)
        
        # efficiency curve
        x=np.array([0,2.5,5,10,20,30,50,100])/100
        y_ref=np.array([0,70,82,89,92,94,96.5,92])/100
        y = y_ref * self.eta_max/max(y_ref) 
        self.eta = interp1d(x, y)
        
        
    def use(self,h,e):
        """

        h: int hour to be simulated
        e: electricity provided (e>0) [kWh]
      
        output : e_lost of electricity that hour (e<0) [kWh]
        """
        
        if e <= 0:
            e_lost = 0

        else:
            e = e/self.n
            
            if e > self.peakP:
                self.eta_story[h] = self.eta(1)
                e_lost = e-self.peakP + self.peakP*(1-self.eta_story[h])
                
            else:  
                self.eta_story[h] = self.eta(e/self.peakP)
                e_lost = e*(1-self.eta_story[h])
                
            e_lost = e_lost*self.n
                          
        return(-e_lost)
    
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

        size = self.peakP
        
        if tech_cost['cost per unit'] == 'default price correlation':
            C0 = 500 # €/kW
            scale_factor = 0.6 # 0:1
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
    Functional test
    """
    parameters = {"max efficiency":0.95, 
                  "peakP":5,
                  "number":1}
    inv = inverter(parameters,8760)
    inv.use(0,4)
    inv.eta_story
