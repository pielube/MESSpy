# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 15:45:44 2023

@author: AndreaAde
"""

import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
from numpy import log as ln
from numpy import e
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression
import pandas as pd
from sklearn.preprocessing import PolynomialFeatures
from scipy.optimize import curve_fit
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.getcwd(),os.path.pardir)))   # temporarily adding constants module path
from core import constants as c

class mhhc_compressor:

    def __init__(self,parameters, timestep_number,timestep=False):

        """
        Create a Hydride Hydrogen Compressor object

        parameters : dictionary
            'Beta': float compression ratio [-]
            
            timestep_number : int number of timesteps considered in the simulation

        output : Hydride Hydrogen Compressor able to:
            absosrb hydrogen at a certain level of pressure and temperature
            and desorbs it at an higher level of pressure and temperature .use(e)
        """
        
        self.ABS_Temp               = 20+273.15     # [K] Isoterma di assorbimento, T_low temperatura a cui entra l'idrogeno   REF:doi:10.1016/j.jallcom.2009.05.069
        self.DES_Temp               = 100+273.15    # [K] Isoterma di desorbimento T_high temperatura a cui esce l'idrogeno    REF:https://doi.org/10.1016/j.ijhydene.2020.09.249
        self.PressioneLowABS        = 0.1090        # [bar] pressione minima ABS @ 20°C        ref:doi:10.1016/j.jallcom.2009.05.069
        self.PressioneHigh_ABS      = 56.3946       # [bar] pressione massima ABS @ 20°C      ref:doi:10.1016/j.jallcom.2009.05.069
        self.PressioneLow_DES       = 76.5490       # [bar] pressione minima DES @ 100°C      ref: https://doi.org/10.1016/j.ijhydene.2020.09.249
        self.PressioneHigh_DES      = 221.9176      # [bar] pressione massima DES @ 100°C    ref: https://doi.org/10.1016/j.ijhydene.2020.09.249
        self.Conc_Low_ABS           = 0.035         # [wt%] percentuale MINIMA di idrogeno assorbito su massa di TiZrV    ref:doi:10.1016/j.jallcom.2009.05.069
        self.Conc_High_ABS          = 1.798         # [wt%] percentuale MASSIMA di idrogeno assorbito su massa di TiZrV   ref:doi:10.1016/j.jallcom.2009.05.069
        self.Conc_Low_DES           = 0.171         # [wt%] percentaule MINIMA di idrogeno desorbito su massa di TiZrV    ref: https://doi.org/10.1016/j.ijhydene.2020.09.249
        self.Conc_High_DES          = 1.585         # [wt%] percentuale MASSIMA di idrogeno desorbito su massa di TiZrV   ref: https://doi.org/10.1016/j.ijhydene.2020.09.249
        self.CycleTime              = 1000          # [s] tempo necessario per assorbire e desorbire                          ref: https://doi.org/10.1016/j.ijhydene.2020.09.249
        self.MetalTankMass          = 30            # [kg] massa tank cilindrico in cui viene alloggiato il letto di MH  --
        self.H2AbsAlloyMass         = 20            # [kg] massa lega che assorbe l'idrogeno, cioè massa del letto di MH
        self.CvMH                   = 0.50          # [kJ/(kg*K)] calore specifico a v costante del MH  doi:10.1016/j.ijhydene.2009.12.027
        self.HeatCapacityTank       = 0.50          # [kJ/(kg*K)] "               " tank  doi:10.1016/j.ijhydene.2012.04.088
        self.PolytropicCoeff        = 1.3           # [-] ref:doi:10.1016/j.ijhydene.2012.04.088
        #COEFFICIENTI FASE ALFA ABS
        self.A_alfaABS              = 1.5e-03       # [-]
        self.Gamma_alfaABS          = 0.96          # [-]
        self.V_alfaABS              = -31           # [m^3/mol]
        self.Enthalpy_alfa_ABS      = -9490         # [J/mol]
        self.Conc_media             = 0.97          # [wt%]
        self.SlopeFactor            = 0.36          # [-]
        #COEFFICIENTI FASE BETA ABS
        self.A_beta_ABS             = 0.07          # [-]
        self.Gamma_betaABS          = 0.966         # [-]
        self.V_betaABS              = 20.0          # [m^3/mol]
        self.Enthalpy_beta_ABS      = -4537         # [J/mol]
        #COEFFICIENTI DI FORMAZIONE ABS
        self.DeltaH_formazione_ABS  = 20120         # [J/mol]
        self.DeltaS_formazione_ABS  = 97.4          # [J/mol*K]
        #COEFFICIENTI FASE ALFA DES
        self.A_alfa_DES             = 6.5e-03       # [-]
        self.Gamma_alfaDES          = 1.55          # [-]
        self.V_alfaDES              = 28            # [m^3/mol]
        self.Enthalpy_alfa_DES      = -4000         # [J/mol]
        self.SlopeFactor_DES        = 0.34
        #COEFFICIENTI FASE BETA DES
        self.A_beta_DES             = 0.78          # [-]
        self.Gamma_betaDES          = 0.213         # [-]
        self.V_betaDES              = 18.7          # [m^3/mol]
        self.Enthalpy_beta_DES      = -5042         # [J/mol]
        #COEFFICIENTI DI FORMAZIONE DES
        self.DeltaH_formazione_DES  = 24700         # [J/mol]
        self.DeltaS_formazione_DES  = 107.5         # [J/mol]
        ##################################
        self.Conc_media = 0.97
        coeff = 0.02                                # correction coefficient added to make the representation of the curve smoother, near the points where the phase changes the model does not perform well and needs this experimental coefficient.
        self.H2MolMass  = 2.01588e-3                # [kg/mol] hydrogen molar mass
        
        if __name__ == "__main__":                  # if code is being executed from chp_gt.py script
            self.timestep   = timestep            # [min]       simulation timestep if launched from main
        else:
            self.timestep   = c.timestep              # [min]       simulation timestep if launched from fuelcell.py
        
        self.n_compressor = parameters['compressor number']   #[-] number of compressors working at the same time
        self.Q         = 3                                    #[kWh] Heat requested at design point--->equivalent to kW at the equivalent timestep
        self.n_compressors_used = np.zeros(timestep_number)
        self.ETA_Polytropic = np.zeros(timestep_number)
        
        'abs and des curves are divided into three parts to represent the absorption, transition and desorption phases'

        self.pressione_abs_data = np.array((0.109038103,0.203484305,0.492970063,1.076464064,2.324669654,4.329380067,7.661159279,14.00865151,16.29187226,20.36902385,22.50415024,23.96121874,25.89608039,27.42838985,28.85380274,29.19111169,29.72353623,30.0783013,31.7531471,33.24368276,35.03720615,36.82811925,39.09345154,40.776238,43.53396474,43.96834197,46.62616884,50.32246938,51.54261613,56.39462298), float)
        self.pressione_des_data = np.array((76.54907473,90.62949413,101.4176079,108.120863,114.132208,117.9774571,120.7459488,123.4458739,126.1428977,127.3642353,127.6032787,127.6132562,129.5009164,130.3620339,132.2435345,135.1168243,138.3496002,142.6933504,144.5821304,145.1823102,150.0186411,156.5416732,161.0380261,168.9237622,175.6874213,180.0136742,192.4804189,203.6347936,209.5552888,221.9176829),float)
        self.conc_abs_data      = np.array((0.035355127,0.042985636,0.057000856,0.070237452,0.110102966,0.148099784,0.199582733,0.271393583,0.294948142,0.349633766,0.387427954,0.412614073,0.441402498,0.504382805,0.606958453,0.727545629,0.81210975,0.869686599,0.977665706,1.094620557,1.211605427,1.330361431,1.440141691,1.533711575,1.596691883,1.618275696,1.688460615,1.724423631,1.746007445,1.798210855), float)
        self.conc_des_data      = np.array((0.170884537,0.23782002,0.304755504,0.371690988,0.438626472,0.505561956,0.572125576,0.639432924,0.661744752,0.68405658,0.70636841,0.72868024,0.75099206,0.79561572,0.840239375,0.90717486,0.974110343,1.08566948,1.130293139,1.15260497,1.22107655,1.33109959,1.398035074,1.464970558,1.509594214,1.53344214,1.567441116,1.576036866,1.579109063,1.585253456), float)
        
        
        'Absorption'

        self.Npoints = 30

        self.conc_abs = np.zeros(self.Npoints)

        for i in range(0,self.Npoints):

            if 0<=i<=9:                     # ALPHA PHASE
                self.conc_abs[i]=self.A_alfaABS*((self.pressione_abs_data[i])**(self.Gamma_alfaABS/2))*e**((-self.Gamma_alfaABS*self.V_alfaABS*self.pressione_abs_data[i])/(c.R_UNIVERSAL*self.ABS_Temp))*e**((-self.Gamma_alfaABS*self.Enthalpy_alfa_ABS)/(c.R_UNIVERSAL*self.ABS_Temp))
            elif 10<=i<=25:                 # BETA + ALPHA PHASE
                self.conc_abs[i]=self.Conc_media+(1/(self.SlopeFactor))*(ln(self.pressione_abs_data[i])+(self.DeltaH_formazione_ABS/(c.R_UNIVERSAL*self.ABS_Temp))-(self.DeltaS_formazione_ABS/(c.R_UNIVERSAL)))
            elif 26<=i<=30:                 # BETA PHASE
                self.conc_abs[i]=self.A_beta_ABS*((self.pressione_abs_data[i])**(self.Gamma_betaABS/2))*e**((-self.Gamma_betaABS*self.V_betaABS*self.pressione_abs_data[i])/(c.R_UNIVERSAL*self.ABS_Temp))*e**((-self.Gamma_betaABS*self.Enthalpy_beta_ABS)/(c.R_UNIVERSAL*self.ABS_Temp))

        self.conc_abs=np.sort(self.conc_abs)

        'Curve smooth visualization'

        self.conc_abs_smooth= np.zeros(self.Npoints)  # aggiunta del coefficiente correttivo coeff

        for i in range(0,self.Npoints):
            if 0<=i<=7:
                self.conc_abs_smooth[i]=self.conc_abs[i]
            elif 8<=i<=9:
                self.conc_abs_smooth[8]=self.conc_abs[8]+coeff
                self.conc_abs_smooth[9]=self.conc_abs[9]
            elif i>=10:
                self.conc_abs_smooth[i]=self.conc_abs[i]

        'Desorption'

        self.conc_des=np.zeros(self.Npoints)

        for i in range(0, self.Npoints):

            if 0<=i<=6:
                self.conc_des[i]=self.A_alfa_DES*((self.pressione_des_data[i])**(self.Gamma_alfaDES/2))*e**((-self.Gamma_alfaDES*self.V_alfaDES*self.pressione_des_data[i])/(c.R_UNIVERSAL*self.DES_Temp))*e**((-self.Gamma_alfaDES*self.Enthalpy_alfa_DES)/(c.R_UNIVERSAL*self.DES_Temp))
            elif 7<=i<=22:
                self.conc_des[i]=self.Conc_media+(1/(self.SlopeFactor))*(ln(self.pressione_des_data[i])+(self.DeltaH_formazione_DES/(c.R_UNIVERSAL*self.DES_Temp))-(self.DeltaS_formazione_DES/(c.R_UNIVERSAL)))
            elif 23<=i<=30:
                self.conc_des[i]=self.A_beta_DES*((self.pressione_des_data[i])**(self.Gamma_betaDES/2))*e**((-self.Gamma_betaDES*self.V_betaDES*self.pressione_des_data[i])/(c.R_UNIVERSAL*self.DES_Temp))*e**((-self.Gamma_betaDES*self.Enthalpy_beta_DES)/(c.R_UNIVERSAL*self.DES_Temp))


        self.conc_des=np.sort(self.conc_des)

        'Definition of interpolation methods for absorption-desorption curves'

        self.interp_abs=interp1d(self.pressione_abs_data, self.conc_abs_smooth, kind='cubic', bounds_error=False, fill_value='extrapolate')
        self.interp_des=interp1d(self.pressione_des_data, self.conc_des, kind='cubic', bounds_error=False, fill_value='extrapolate')
        
        Beta=[2.706,2.500,2.371,2.247,2.128,2.013,1.903,1.747,1.599,1.369]       # [-]
        DeltaConc=[1.399,1.415,1.426,1.437,1.448,1.459,1.470,1.487,1.503,1.531]  # [wt%]
        x=np.linspace(2.706,1.369,100)                                           # [-]

        self.interp_beta=interp1d(Beta, DeltaConc, kind='cubic', bounds_error=None, fill_value='extrapolate')

        # plt.figure(dpi=1000)
        # plt.plot(self.conc_des, self.pressione_des_data, linewidth=3)
        # plt.plot(self.interp_des(self.pressione_des_data), self.pressione_des_data, linewidth=3, linestyle='--')
        # plt.yscale('log')

        # plt.figure(dpi=1000)
        # plt.plot(self.conc_abs, self.pressione_abs_data, linewidth=3)
        # plt.plot(self.interp_abs(self.pressione_abs_data), self.pressione_abs_data, linewidth=3, linestyle='--')
        # plt.yscale('log')
        # plt.show()

    def abs_validationplot(self):

        print('---------------absorption process validation-----------------')
        plt.figure(dpi=1000, figsize = (6,5))
        plt.plot(self.conc_abs_data, self.pressione_abs_data, label='Experimental', linewidth=1.9,color = '#eb4034', linestyle = '--', zorder =2)
        plt.scatter(self.conc_abs_data, self.pressione_abs_data, edgecolors='k', color ='#eb4034', zorder = 3, s=20)
        plt.plot(self.conc_abs_smooth, self.pressione_abs_data, label='Model', linewidth=1.9, color ='cornflowerblue', zorder = 1)
        plt.legend(loc='lower center', fontsize=10)
        plt.grid(which='both', zorder = 0, alpha = 0.3, )
        plt.xlim(0,2)
        plt.xlabel('Concentration [wt%]')
        plt.ylim(0.1,100)
        plt.ylabel('P$_{abs}$ [bar]')
        plt.yscale('log')
        plt.axvline(self.conc_abs_data[11], ymax =0.8, linestyle = '--', c='k')
        # plt.axvline(self.conc_abs_data[11], ymax =self.pressione_abs_data[11]/max(self.pressione_abs_data), linestyle = '--', c='k')
        # plt.axvline(self.conc_abs_data[25], ymax =self.pressione_abs_data[26]/max(self.pressione_abs_data), linestyle = '--', c='k')
        plt.axvline(1.6, ymax = 0.86, linestyle = '--', c='k')
        plt.text(0.25, 2, r'$\mathit{\alpha}$', fontsize=18)
        plt.text(0.85, 2, r'$\mathit{\alpha + \beta}$', fontsize=18)
        plt.text(1.65, 2, r'$\mathit{\beta}$', fontsize=18)
        # plt.savefig('AbsorptionValidation.png')
        plt.show()

        MSE = np.zeros(self.Npoints)                    # mean square error
        SS_res = np.zeros(self.Npoints)                 # residual value
        Media_Conc_real = np.mean(self.conc_abs_data)
        SS_tot = np.zeros(self.Npoints)

        for i in range(0,self.Npoints):

            MSE[i]=((self.conc_abs_smooth[i]-self.conc_abs_data[i])**2)/self.Npoints
            SS_res[i]=((self.conc_abs_smooth[i]-self.conc_abs_data[i])**2)
            SS_tot[i]=(self.conc_abs_data[i]-Media_Conc_real)**2


        DEVST = np.std(MSE)
        MSE=sum(MSE)
        RMSE =MSE**0.5
        SS_res=sum(SS_res)
        SS_tot=sum(SS_tot)

        R2=1-SS_res/SS_tot


        print('\nMean Square Error_ABS:', MSE)
        print('Standard Deviation MSE_ABS:', DEVST)
        print('Root Mean Square Error_ABS:', RMSE)
        print('Coefficient of Determination_ABS:', R2)
        print('\n---------------end of absorption process validation----------------')

    def des_validationplot(self):

        print('\n--------------desorption process validation-------------')

        plt.figure(dpi=1000, figsize = (6,5))
        plt.plot(self.conc_des_data, self.pressione_des_data, label='Experimental', linewidth=1.9,color = '#eb4034', linestyle = '--', zorder =2)
        plt.scatter(self.conc_des_data, self.pressione_des_data, edgecolors='k', color ='#eb4034', zorder = 3, s=20)
        plt.plot(self.conc_des, self.pressione_des_data, label='Model', linewidth=1.9, color ='cornflowerblue', zorder = 1)
        plt.legend(loc='lower center', fontsize=10)
        plt.grid(which='both', zorder = 0, alpha = 0.3, )
        plt.xlim(0,1.75)
        plt.xticks(np.arange(0,2, 0.25))
        plt.xlabel('Concentration [wt%]')
        # plt.ylim(70, 250)
        plt.axvline(0.455, ymax = 0.42, linestyle = '--', c='k')
        plt.axvline(1.46, ymax = 0.7, linestyle = '--', c='k')
        plt.text(0.25, 80, r'$\mathit{\alpha}$', fontsize=18)
        plt.text(0.85, 100, r'$\mathit{\alpha + \beta}$', fontsize=18)
        plt.text(1.65, 100, r'$\mathit{\beta}$', fontsize=18)
        plt.ylabel('P$_{des}$ [bar]')
        plt.yscale('log')
        # plt.title('DESORPTION CURVE')
        # plt.savefig('DesorptionValidation.png')
        plt.show()

        MSE = np.zeros(self.Npoints)                     # mean square error
        SS_res = np.zeros(self.Npoints)                  # residual values
        Media_Conc_real = np.mean(self.conc_des_data)
        SS_tot = np.zeros(self.Npoints)

        for i in range(0,self.Npoints):

            MSE[i]=((self.conc_des[i]-self.conc_des_data[i])**2)/self.Npoints
            SS_res[i]=((self.conc_des[i]-self.conc_des_data[i])**2)
            SS_tot[i]=(self.conc_des_data[i]-Media_Conc_real)**2


        DEVST = np.std(MSE)
        MSE=sum(MSE)
        RMSE =MSE**0.5
        SS_res=sum(SS_res)
        SS_tot=sum(SS_tot)

        R2=1-SS_res/SS_tot


        print('\nMean Square Error_DES:', MSE)
        print('Standard Deviation MSE_DES:', DEVST)
        print('Root Mean Square Error_DES:', RMSE)
        print('Coefficient of Determination_DES:', R2)
        print('\n----------------end desorption process validation-----------------')
        
    def plot_absdesplot(self):

        'Absorption'
        plt.figure(dpi=1000)
        plt.plot(self.conc_abs_smooth,self.pressione_abs_data,label='ABS-data', color='b',marker='.', linestyle='None', mec='r', markersize=7, markerfacecolor='white', zorder=0)
        plt.plot(self.interp_abs(self.pressione_abs_data),self.pressione_abs_data,label='cubic', linestyle='--')
        plt.yscale('log')
        plt.grid()
        plt.legend(fontsize=8)
        plt.xlabel('Concentrazione [wt%]')
        plt.ylabel('Pressione [bar]')
        plt.title('Absorption Curve' )

        'Desorption'
        plt.figure(dpi=1000)
        plt.plot(self.conc_des,self.pressione_des_data,label='DES-data', color='b',marker='.', linestyle='None', mec='r', markersize=7, markerfacecolor='white', zorder=0)
        plt.plot(self.interp_des(self.pressione_des_data),self.pressione_des_data,label='cubic', linestyle='--')
        plt.yscale('log')
        plt.grid()
        plt.legend(fontsize=8)
        plt.xlabel('Concentrazione [wt%]')
        plt.ylabel('Pressione [bar]')
        plt.title('Desorption Curve' )
        plt.show()
        
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

        size = self.Q
            
        if tech_cost['cost per unit'] == 'default price correlation':
            mh_price = 200   #[euros/kg]
            steel_inox_price = 3.52 #[euros/kg]
            C = (mh_price*self.H2AbsAlloyMass + steel_inox_price*self.MetalTankMass)*self.n_compressor
        else:
            C = size * tech_cost['cost per unit']

        tech_cost['total cost'] = tech_cost.pop('cost per unit')
        tech_cost['total cost'] = C
        tech_cost['OeM'] = tech_cost['OeM'] *C /100 # €


        self.cost = tech_cost    
       

    def plot_performancemhhc(self):
        

        self.conc_abs_beta=np.split(self.conc_abs,6)   # split the concentration vector during absorption to take only the beta phase points, which are 5

        self.conc_abs_betaS=self.conc_abs_beta[5]

        self.Npunti=25

        self.conc_abs_betaplus=np.zeros(self.Npunti)   # add points to the beta phase, then do the same with the alpha phase of desorption
        for i in range(0,self.Npunti):
            if i == 0:
                self.conc_abs_betaplus[i] = min(self.conc_abs_betaS)
            else:
                self.conc_abs_betaplus[i] = min(self.conc_abs_betaS)+((max(self.conc_abs_betaS)-min(self.conc_abs_betaS))/self.Npunti)*i

        self.conc_des_alfa=np.zeros(7)
        for i in range(0,7):
            self.conc_des_alfa[i]=self.conc_des[i]

        self.conc_des_alfaplus=np.zeros(self.Npunti)    #add points to the alpha phase of desorption
        for i in range(0,self.Npunti):
            if i ==0:
                self.conc_des_alfaplus[i]=min(self.conc_des_alfa)
            else:
                self.conc_des_alfaplus[i]=min(self.conc_des_alfa)+((max(self.conc_des_alfa)-min(self.conc_des_alfa))/self.Npunti)*i

        self.conc_des_alfaplus = sorted(self.conc_des_alfaplus,reverse=True) #invert to create the vector delta conc and then beta so that the components vary between the maximum and minimum values
        self.DeltaConc=np.zeros(self.Npunti)  #vector corresponding to the possible combinations of concentration changes, ranging from maximum (which will then correspond to the minimum beta) to minimum

        for i in range(0,self.Npunti):
            self.DeltaConc[i]=self.conc_abs_betaplus[i]-self.conc_des_alfaplus[i]

        self.pressione_abs_beta=np.split(self.pressione_abs_data,6)  # the same procedure adopted for the concentration is applied for the pressure, so as to derive beta, i.e. compression ratio
        self.pressione_abs_betaS=self.pressione_abs_beta[5]

        self.pressione_abs_betaplus=np.zeros(self.Npunti)
        for i in range(0,self.Npunti):
            if i==0:
                self.pressione_abs_betaplus[i]=min(self.pressione_abs_betaS)
            else:
                self.pressione_abs_betaplus[i]=min(self.pressione_abs_betaS)+((max(self.pressione_abs_betaS)-min(self.pressione_abs_betaS))/self.Npunti)*i

        self.pressione_des_alfa=np.zeros(7)
        for i in range(0,7):
            self.pressione_des_alfa[i]=self.pressione_des_data[i]

        self.pressione_des_alfaplus=np.zeros(self.Npunti)
        for i in range(0,self.Npunti):
            if i==0:
                self.pressione_des_alfaplus[i]=min(self.pressione_des_alfa)
            else:
                self.pressione_des_alfaplus[i]=min(self.pressione_des_alfa)+((max(self.pressione_des_alfa)-min(self.pressione_des_alfa))/self.Npunti)*i
        self.pressione_des_alfaplus=sorted(self.pressione_des_alfaplus, reverse=True)


        self.pressione_absMAX=max(self.pressione_abs_betaplus)   # [bar] max value of pressure during absorption
        

        self.pressione_desMAX=max(self.pressione_des_alfaplus)   # [bar] max value of pressure during desorption
       

        self.Beta=np.zeros(self.Npunti)                          # [-] compressione ratio

        for i in range(0,self.Npunti):
            self.Beta[i]=self.pressione_des_alfaplus[i]/self.pressione_abs_betaplus[i]

        #specific volume of hydrogen at 20°C (absorption temperature) between minimum and maximum absorption pressure. This is necessary to calculate the work supplied to the hydrogen
        self.H2SpecificVolume=np.array((0.294419, 0.291031, 0.287722, 0.284439, 0.281330, 0.278242, 0.275223, 0.272271, 0.269383, 0.266557, 0.263792, 0.261085, 0.258435, 0.255840, 0.253297, 0.250806, 0.248366, 0.245973, 0.243628, 0.241329, 0.239073, 0.236861, 0.234691, 0.232561, 0.230471), float) #    [m^3/kg] NIST:https://webbook.nist.gov/cgi/fluid.cgi?T=303K&PLow=43.53&PHigh=55.88&PInc=0.52&Applet=on&Digits=10&ID=C1333740&Action=Load&Type=IsoTherm&TUnit=K&PUnit=bar&DUnit=kg%2Fm3&HUnit=kJ%2Fmol&WUnit=m%2Fs&VisUnit=Pa*s&STUnit=N%2Fm&RefState=DEF
        #Work
        self.Work_Polytropic = ((((self.PolytropicCoeff/(self.PolytropicCoeff-1))*self.pressione_abs_betaplus*self.H2SpecificVolume*((self.Beta)**((self.PolytropicCoeff-1)/self.PolytropicCoeff)-1))/10)*(2*self.H2AbsAlloyMass*self.DeltaConc/100)/3600)*1000  #[kW] ref:doi:10.1016/j.ijhydene.2012.04.088
        #Heat Requested
        self.Heat_Required = (((((self.H2AbsAlloyMass*self.DeltaConc/100)/self.H2MolMass)*self.DeltaH_formazione_DES/1000 + (self.DES_Temp-self.ABS_Temp)*self.CvMH*(self.MetalTankMass+self.H2AbsAlloyMass))/500)/3600)*1000 #[kW] ref:doi:10.1016/j.ijhydene.2012.04.088 #sarebbe diviso 1000 per esprimere in MJ ma visto che il calore richiesto è necessario sia per riscaldare che per raffreddare, 2/1000 = 500

        #interpolation to obtain the alpha- and beta-phase concentration to be used in the next function

        self.interp_abs_beta=interp1d(self.pressione_abs_betaplus, self.conc_abs_betaplus, kind='cubic', bounds_error=None, fill_value='extrapolate')
        self.interp_des_alfa=interp1d(self.pressione_des_alfaplus, self.conc_des_alfaplus, kind='cubic', bounds_error=None, fill_value='extrapolate')

        Beta=[2.706,2.500,2.371,2.247,2.128,2.013,1.903,1.747,1.599,1.369]       # [-]
        DeltaConc=[1.399,1.415,1.426,1.437,1.448,1.459,1.470,1.487,1.503,1.531]  # [wt%]
        x=np.linspace(2.706,1.369,100)                                           # [-]
        Beta1 = [1.3647, 1.4972, 1.6478, 1.7984, 1.949, 2.0996, 2.2482, 2.4008, 2.5474, 2.698] # valori estratti a mano da grafico tesi Niccolò 
        Efficienza = [0.0222, 0.0291, 0.0362, 0.0428, 0.049, 0.0548, 0.0602, 0.0654, 0.070, 0.0748]

        self.interp_beta=interp1d(Beta, DeltaConc, kind='cubic', bounds_error=None, fill_value='extrapolate')

        # plt.figure(dpi=1000, figsize = (6,5))
        fig,ax = plt.subplots(dpi=1000, figsize = (6,5))
        ax2 = ax.twinx()
        ax3 = ax.twinx() 
        ax3.spines["right"].set_position(("axes", 1.15)) # posizione dell'asse
        ax3.spines["right"].set_visible(True) # rendi visibili le spine
        # ax3.spines["right"].set_color('#FF8C00') # colore delle spine
        # plt.title('MHHC Performance')
        ax.plot(self.Beta, self.DeltaConc, label=r'$\Delta_{conc}$', c='#F55072')
        ax2.plot(self.Beta, self.Work_Polytropic, label=r'Power', c='#3BBDD4')
        ax3.plot(Beta1, Efficienza, label=r'$\eta$', c='#4D6AEB')
        ax3.set_ylim(0,0.1)
        # plt.plot(DeltaConc, Beta, linestyle='None', marker='.',mec='r',markersize=10)
        # plt.plot(self.interp_beta(x), x, label=r'$\beta_{data}$', marker='.', markersize=2)
        ax.set_ylabel(r'$\Delta$Concentration [wt%]')
        ax2.set_ylabel(r'Power [kW]')
        ax3.set_ylabel(r'Efficiency [-]')
        ax.set_xlabel(r'$\beta$ [-]')
        ax.grid(alpha= 0.3)
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        h3, l3 = ax3.get_legend_handles_labels()
        ax.legend(h1+h2+h3, l1+l2+l3, loc='upper center', ncols= 3)
        # plt.savefig('PerformancePlot1.png')
        plt.show()


    def use(self,step,hyd,storable_hydrogen):

        p_in=45

        p_out=120
                  
        beta=p_out/p_in  #[-]   compression ratio

        Q_requested=((((((self.H2AbsAlloyMass*(self.interp_beta(beta))/100)/self.H2MolMass)*self.DeltaH_formazione_DES/1000 + (self.DES_Temp-self.ABS_Temp)*self.CvMH*(self.MetalTankMass+self.H2AbsAlloyMass))/500)/3600)*1000)#*self.n_compressor #[kW] heat requested
        # print(Q_requested)
        H2percycle_h = ((self.H2AbsAlloyMass*(self.interp_beta(beta))/100)/(c.H2SDENSITY*self.CycleTime/3600))#*self.n_compressor #- (Q_requested/(c.LHVH2*1000))*3600/c.H2SDENSITY  #[Sm^3/h]   flow rate that can be desorbed in the time interval considered (one hour) NET VALUE
        
        Work_Polytropic = (((((self.PolytropicCoeff/(self.PolytropicCoeff-1))*p_in*0.27593*((beta)**((self.PolytropicCoeff-1)/self.PolytropicCoeff)-1))/10)*(2*self.H2AbsAlloyMass*(self.interp_beta(beta))/100)/3600)*1000)#*self.n_compressor  #[kW] work supplied to hydrogen in the form of pressure
        # print(Work_Polytropic)
        ETA_Polytropic = Work_Polytropic/Q_requested #[-] compressor efficiency
        # print(ETA_Polytropic)
        self.ETA_Polytropic[step] = ETA_Polytropic

        self.H2_kg = H2percycle_h*c.H2SDENSITY
       
        if hyd == 0:
            n_compressor_used = 0
        else:
            n_compressor_used = hyd/self.H2_kg
        if n_compressor_used > self.n_compressor:
            print('Warning: The number of Methal Hydride Hydrogen Compressors is not sufficient \n')
            self.n_compressors_used[step] = self.n_compressor
            hyd_compressed = self.H2_kg*self.n_compressor
            Q_requested = Q_requested*self.n_compressor
            Sm3_requested = Q_requested / (c.LHV_H2*3600)
        else:
            self.n_compressors_used[step] = int(n_compressor_used)+1
            hyd_compressed = hyd
            Q_requested = Q_requested*n_compressor_used
            Sm3_requested = Q_requested / (c.LHV_H2*3600)

        return (hyd_compressed,-Sm3_requested)

##########################################################################################

if __name__ == "__main__":

    """
    Functional test
    """

    inp_test = {'compressor number': 90,
                'compressor model':'Ti-V-Zr',
                }

    sim_hours=36                             # [h] simulated period of time - usually it's 8760 hours

    mhhc=mhhc_compressor(inp_test, sim_hours)
    mhhc.plot_absdesplot()                     # mhhc abs-des curve
    # mhhc.plot_performancemhhc()                # mhhc performance curve
    mhhc.abs_validationplot()                  # mhhc absorption validation curve
    mhhc.des_validationplot()                  # mhhc desorption validation curve
    storable_hydrogen = 10000

    hyd_to_be_compressed = np.linspace(0.5,2.5,sim_hours)                   # Hydrogen to be compressed by the MHHC
    hyd_compressed = np.zeros(sim_hours)
    Sm3_requested = np.zeros(sim_hours)
    Q_requested = np.zeros(sim_hours)
    # n_compressors_used = []
    eff = np.zeros(sim_hours)
    for i in range (len(hyd_to_be_compressed)):
        hyd_compressed[i],Sm3_requested[i] = mhhc.use(i,hyd_to_be_compressed[i], storable_hydrogen)
        Q_requested[i] = Sm3_requested[i] * c.LHV_H2*3600
        n_compressors_used = mhhc.n_compressors_used
    
    '-------Plots-----'
    fig,ax = plt.subplots(dpi=1000, figsize = (6,5))
    ax2 = ax.twinx()
    # ax2.plot(hyd_compressed, eff, label=r'$\eta$', c='#3BBDD4' )
    ax.plot(hyd_compressed, -Q_requested, label=r'Heat', c='#F55072', linewidth= 1.5)
    ax2.bar(hyd_compressed,n_compressors_used,width=0.03,zorder=3,edgecolor='k',label='Modules', alpha =0.8) 
    ax.set_ylim(0,None)
    ax2.set_ylim(0,6)
    ax.set_ylabel(r'P$_{th,req}$ [kW]')
    ax2.set_ylabel(r'N$_{modules}$ [-]')
    ax.set_xlabel('ṁ$_{H2}$ [kg/h]')
    ax.grid(alpha= 0.3)
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1+h2, l1+l2, loc='upper center', ncol=2)
    # plt.savefig('PerformancePlot2.png')
    # plt.show()
    

    # # create figure and axis objects with subplots()
    # fig,ax = plt.subplots(dpi=600)
    # ax.plot(hyd_to_be_compressed,n_compressors_used,color="red",marker="o")
    # ax.set_xlabel("Hydrogen to be compressed [kg]", fontsize = 14)
    # ax.set_ylabel("N. compressors used",color="red",fontsize=14)
    # # twin object for two different y-axis on the sample plot
    # props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    # plt.text(hyd_to_be_compressed[-1]/20,n_compressors_used[-1]-n_compressors_used[-1]/10,'Hyd compressed by one MHHC=%.3f kg' % (mhhc.H2_kg),fontsize=10,va='bottom',backgroundcolor='none', bbox=props)
    # ax2 = ax.twinx()
    # ax2.plot(hyd_to_be_compressed,hyd_compressed,color="blue",marker="o")
    # ax2.set_ylabel("Hydrogen compressed [kg]",color="blue",fontsize=14)
    # ax.grid()
    # plt.show()
    
    # fig,ax = plt.subplots(dpi=600)
    # ax.plot(hyd_compressed,-Q_requested,color="red",marker="o")
    # ax.set_xlabel("Hydrogen compressed [kg]", fontsize = 14)
    # ax.set_ylabel("Necessary heat [kW]",color="red",fontsize=14)
    # ax.grid()
    # plt.show()