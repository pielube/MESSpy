from scipy.interpolate import interp1d

class fuel_cell:
    
    def __init__(self,parameters):
        """
        Create a Fuel ell object
    
        parameters : dictionary
            'Npower': float nominal power [kW]
                      
        output : Fuel cell object able to:
            abrosrb hydrogen and produce electricity .use(e)
        """
        
        # this model is based on FCS-C5000 characteristic curves https://www.horizonfuelcell.com/hseries
        
        self.Npower = parameters['Npower'] * 1000 # W
        self.nc = 120 * self.Npower / 5000 # # number of cells (120 cells are required to have a nominal power of 5000 W)
        
        # characteristic curves
        I=[0, 10, 20, 30, 40, 50, 60, 70, 80]
        v=[0.96, 0.83, 0.796, 0.762, 0.728, 0.694, 0.66, 0.6, 0.5]
        V=[]
        P=[]
        for i in range(len(I)):
            volt=self.nc*v[i]
            V.append(volt)
            pot=I[i]*volt
            P.append(pot)
        
        self.IP=interp1d(I,P,kind='cubic',bounds_error=False,fill_value='extrapolate')  
        self.Pmax=max(P)
       
        # PI inverse curve that will be used to find the operating point
        I=[]
        P=[]
        for a in range(800):
            i=a/10
            pot=self.IP(i)
            if pot>self.Pmax:
                break  
            I.append(i)
            P.append(pot)
        self.PI=interp1d(P,I,kind='cubic',bounds_error=False,fill_value='extrapolate')
        
        
    def use(self,e,available_hyd):
        """
        The Fuel celle can absorb hydrogen and produce electricity: H2 --> 2H+ + 2e
    
        e: float < 0 energy required  [kWh]
        available_hyd: float available hydrogen H tank SOC[h-1] [kg]
      
        output : hydrogen absorbed and electricity supplied that hour [kg]
        """
        
        p_required = -e*1000  # kWh --> W
        
        p = min(p_required,self.Npower) # how much electricity can Fuel Cell absorb?
        
        # find the operating point on the characteristic curves
        I=self.PI(p)
       
        # calculate the hydrogen consumed by each single cell
        qe=I/96485 #[mol_e/s] Faraday
        qH=qe*0.5*3600 #[mol_H2/h] the moles of H2 are 1/2 the moles of electrons. 3600s in one hour.
        QH=qH*2.058 #[g_H2/h] H2 molar mass = 2.058
        
        hyd = QH*self.nc/1000 # total stack hydrogen [kg/hr]
        
        if hyd > available_hyd: # if not enough hydrogen is available to meet demand (H tank is nearly empty)
            hyd = 0
            p = 0
            # turn off the fuel cell
            # this behavior could be solved with more advanced models, necessary inverse production functions.
        
        return (-hyd,p/1000) # return hydrogen absorbed [kg] and electricity required [kWh]
        
    
    
##########################################################################################

if __name__ == "__main__":
    
    """
    Functional test
    """
    
    inp_test = {'Npower': 5}
    
    fc = fuel_cell(inp_test)

    flow = [-1,-2,0,-6,-5,-8,-1]
    for h in range(len(flow)):
        print(fc.use(flow[h],100))



