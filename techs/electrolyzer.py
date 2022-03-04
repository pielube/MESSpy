import numpy as np

class electrolyzer:
    
    def __init__(self,parameters,simulation_hours):
        """
        Create an electrolyzer object
    
        parameters : dictionary
            'Npower': float nominal power [kW]
            'stack model': str 'Enapter 2.1' or 'McLyzer 800' are aviable
                      
        output : electrolyzer object able to:
            abrosrb electricity and produce hydrogen .use(e)
        """
        
        self.Npower = parameters['Npower'] # float nominal power [kW]
        H_N_density = 0.08988237638480538 # PropsSI('D', 'T', 273.15, 'P', 101325, 'H2') # hydrogen density under normal condition

        if parameters['stack model'] == 'Enapter 2.1': # https://www.enapter.com/it/newsroom/enapter-reveals-new-electrolyser-el-2-1
            stack_operative_power_consumption = 2.4  # [kW]
            stack_production_rate_L = 500  # [NL/hr]
            stack_production_rate = stack_production_rate_L / 1000 * H_N_density # [kg/hr]
            
        if parameters['stack model'] == 'McLyzer 800': # https://mcphy.com/it/apparecchiature-e-servizi/elettrolizzatori/large/
            stack_operative_power_consumption = 4000 #[kW]
            stack_production_rate_m3 = 800 # [Nm3/hr]
            stack_production_rate = stack_production_rate_m3 * H_N_density # [kg/hr]
    
        self.production_rate = stack_production_rate * self.Npower / stack_operative_power_consumption
        
        self.energy_balance = {'electricity': {'in': np.zeros(simulation_hours)}, 'hydrogen': {'out': np.zeros(simulation_hours)}}
        
    def use(self,h,e,storable_hydrogen):
        """
        Produce hydrogen
        
        e : float > 0 electricity provided in one hour [kWh]
        storable_hydrogen: float storable hydrogen H tank max_capacity - SOC[h-1] [kg] or maximum absorbable production if there is no tank
    
        output : 
            float hydrogen produced that hour [kg]    
            float electricity absorbed that hour [kWh]
        """
        
        e_absorbed = min(self.Npower,e) 
        hyd = self.production_rate * e_absorbed / self.Npower
        
        if hyd > storable_hydrogen: # if the is not enough space in the H tank to store the hydrogen (H tank is nearly full)
            hyd = 0
            e_absorbed = 0
             # turn off the electrolyzer
            # this behavior could be solved with more advanced models, necessary inverse production functions.
            
        self.energy_balance['electricity']['in'][h] = e_absorbed
        self.energy_balance['hydrogen']['out'][h] = hyd
        return(hyd,-e_absorbed) # return hydrogen supplied and electricity absorbed
        
    
##########################################################################################

if __name__ == "__main__":
    
    """
    Functional test
    """
    
    inp_test = {'Npower': 5,
                'stack model': 'Enapter 2.1'}
    
    el = electrolyzer(inp_test)

    flow = [1,2,3,3,5,7,7,7]
    for h in range(len(flow)):
        print(el.use(flow[h]))
