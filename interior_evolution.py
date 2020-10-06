#################################################################
import numpy as np
from scipy.optimize import fsolve, fixed_point, brentq, root
from matplotlib import pyplot as plt

from get_input import *
from get_parameters import *
import support_functions as suppf
import plotting_functions as plotf
#################################################################

class interior_evolution:

####################################################################################
    def __init__(self, body='Earth', inpfile='input.json'):        
        """Set initial conditions and model parameters"""
####################################################################################

        # Set planet parameters and initial conditions
        get_input(self, body=body, inpfile=inpfile)

        # Set model parameters
        get_parameters(self)

######################################################################
    def calculate_evolution(self, outfile=None):
        """Compute evolution of mantle and core temperatures"""
######################################################################
            
        # Initialize arrays
        suppf.initialize_arrays(self)
        
        # Set initial conditions
        suppf.set_initial_conditions(self)

        # Set/compute some parameters
        D  = self.Rp - self.Rc                                   # Mantle thickness                                
        M  = 4./3.*np.pi*((self.Rp**3.-self.Rc**3)*self.rhom + (self.Rc**3)*self.rhoc) # Planet mass
        Mm = 4./3.*np.pi*((self.Rp**3.-self.Rc**3)*self.rhom)    # Mantle mass
        Mc = 4./3.*np.pi*(self.Rc**3*self.rhoc)                  # Core mass
        Ap = 4.*np.pi*self.Rp**2                                 # Surface area             
        Ac = 4.*np.pi*self.Rc**2                                 # Core area             
        kappa = self.km/self.rhom/self.cm                        # Thermal diffusivity  

        if ((self.Q0 > 0.) and (self.lam > 0)):
            pass
        else:             
            # Scale heat production back in time
            X0_U238 = self.X_U*suppf.initialize_heatproduction(self.maxtime,self.tau_U238)
            X0_U235 = self.X_U*suppf.initialize_heatproduction(self.maxtime,self.tau_U235)
            X0_Th232 = self.X_Th*suppf.initialize_heatproduction(self.maxtime,self.tau_Th232)
            X0_K40 = self.X_K*suppf.initialize_heatproduction(self.maxtime,self.tau_K40)

        # Time loop
        for i in np.arange(0,self.n_steps): 
            
            # Mantle Viscosity
            if(i == 0):
               Pm = self.rhom*self.g*self.delta_s[i]    # Pressure at the base of the lithosphere
            else:   
               Pm = self.rhom*self.g*self.delta_s[i-1]  # Pressure at the base of the lithosphere from previous time step
            self.etam[i] = suppf.calculate_viscosity(self, self.Tm[i], Pm)             
            #self.etam[i] = suppf.calculate_viscosity(self, self.Tm[i], 0)             

            # Rayleigh number
            self.Ra[i] = self.rhom*self.g*self.alpha*(self.Tm[i] - self.Ts)*D**3./(kappa*self.etam[i])    

            # Compute surface heat flux based on scaling laws for 
            # mobile lid convection
            if self.tectonics == 'ML':
                self.delta_s[i] = D*(self.Racrit/self.Ra[i])**self.beta        # Thickness of the upper TBL
                self.qs[i] = self.km*(self.Tm[i] - self.Ts) / self.delta_s[i]  # Surface heat flux

            # stagnant lid convection
            elif self.tectonics == 'SL':
                gamma = self.E*(self.Tm[i] - self.Ts)/(self.Rg*self.Tm[i]**2)     # Frank-Kamentzkii parameter                                  
                if(i == 0):                                                       
                    ds_i = self.delta_s0
                else: 
                    ds_i = self.delta_s[i-1]                
                dssol = root( suppf.calculate_ds_root, ds_i, args = (self, self.Tm[i], gamma), method=self.rootmethod )                 
                self.delta_s[i] = dssol.x                                         # Thickness of the upper TBL
                self.qs[i] = self.km*(self.Tm[i] - self.Ts) / self.delta_s[i]     # Surface heat flux
                
            # Consider core cooling 
            if self.core_cooling == 'yes':
                # Compute bottom heat flux with iterations
                if(i == 0):
                    dc_i = self.delta_c0
                else: 
                    dc_i = self.delta_c[i-1]
                # Determine bottom TBL thickness iteratively    
                dcsol = root( suppf.calculate_dc_root, dc_i, args = (self, self.delta_s[i], self.Tm[i], self.Tc[i]), method=self.rootmethod )                 
                dc = dcsol.x

                # Temperature at the top of the lower TBL going adiabatically down from Tm
                zb = self.Rp - (self.Rc + dc)           # depth of the lower TBL
            
                # Pressure at the top of the bottom TBL (Pb) and at the CMB (Pc)
                Pb = self.rhom*self.g*zb
                Pc = self.rhom*self.g*(self.Rp - self.Rc)

                # Temperature at the top of the lower TBL going adiabatically down from Tm
                self.Tprofile[i,:] = suppf.calculate_adiabat(self, self.n_layers, self.Tm[i], Pm, Pb)  
                self.Tb[i] = self.Tprofile[i][-1]

                # Internal and critical Rayleigh number
                deltaTm = (self.Tm[i] - self.Ts)
                deltaTc = (self.Tc[i] - self.Tb[i])
                Ra_int = self.rhom*self.g*self.alpha*( deltaTm + deltaTc )*D**3./(kappa*self.etam[i])    
                Racrit_int = 0.28*Ra_int**0.21

                # Average viscosity near the bottom TBL (etab) and at the CMB (etac)
                Tbmean = (self.Tb[i] + self.Tc[i])/2
                Pbmean = (Pb + Pc)/2
                self.etab[i] = suppf.calculate_viscosity(self, Tbmean, Pbmean)            
                self.etac[i] = suppf.calculate_viscosity(self, self.Tc[i], Pc)

                # Update bottom TBL thickness
                self.delta_c[i] = (kappa*self.etab[i]*Racrit_int/(self.rhom*self.alpha*self.g*np.abs(self.Tc[i] - self.Tb[i])))**self.beta  

                # CMB heat flux
                self.qc[i] = self.km*(self.Tc[i] - self.Tb[i])/self.delta_c[i]    
                
            # Neglect core cooling 
            else: 
                zb = self.Rp - self.Rc
                Pb = self.rhom*self.g*zb  
                #self.Tb[i] = suppf.calculate_adiabat(self, self.Tm[i], zb)                
                self.Tprofile[i,:] = suppf.calculate_adiabat(self, self.n_layers, self.Tm[i], Pm, Pb)  
                self.Tb[i] = self.Tprofile[i][-1]
                self.delta_c[i] = 0.0                
                self.qc[i] = 0.0                
                self.etab[i] = suppf.calculate_viscosity(self, self.Tb[i], Pb)            
                self.etac[i] = suppf.calculate_viscosity(self, self.Tb[i], Pb)
                
            # Radioactive decay
            if ((self.Q0 > 0.) and (self.lam > 0)):
                # based on initial heat production and a single decay constant
                self.Q_tot[i]   = suppf.calculate_radiodecay_simple(self.Q0, self.lam, self.t[i])
            else: 
                # based on concentration and heat production of U, Th and K isotopes
                self.Q_U238[i]  = suppf.calculate_radiodecay(self.f_U238, X0_U238, self.H_U238, self.tau_U238, self.t[i])
                self.Q_U235[i]  = suppf.calculate_radiodecay(self.f_U235, X0_U235, self.H_U235, self.tau_U235, self.t[i])
                self.Q_Th232[i] = suppf.calculate_radiodecay(self.f_Th232, X0_Th232, self.H_Th232, self.tau_Th232, self.t[i])
                self.Q_K40[i]   = suppf.calculate_radiodecay(self.f_K40, X0_K40, self.H_K40, self.tau_K40, self.t[i])
                self.Q_tot[i]   = self.Q_U238[i] + self.Q_U235[i] + self.Q_Th232[i] + self.Q_K40[i]  

            # Urey ratio
            self.Ur[i] = self.Q_tot[i]*Mm / (self.qs[i]*Ap)

            # Advance in time mantle and CMB temperature        
            Mcm = 4./3.*np.pi*( (self.Rp - self.delta_s[i])**3 - self.Rc**3 )*self.rhom # Mass and surface area of the convecting mantle
            Acm = 4.*np.pi*(self.Rp - self.delta_s[i])**2
            Mcm = Mm
            Acm = Ap
            self.Tm[i+1] = self.Tm[i] + self.dt*(self.Q_tot[i]/self.cm + self.Qtidal/self.cm - (Acm*self.qs[i])/(Mcm*self.cm) + (Ac*self.qc[i])/(Mc*self.cm))
            self.Tc[i+1] = self.Tc[i] - ((self.dt*Ac*self.qc[i])/(Mc*self.cc))                                                        
            
        # Write timeseries on file
        if not(outfile):
            print('no output file written')
            pass
        else:
            print('output written in ' + outfile)
            suppf.write_output_file(self, outfile)
        