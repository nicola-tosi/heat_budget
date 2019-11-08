######################################################
import numpy as np
from scipy.optimize import fsolve, fixed_point
from matplotlib import pyplot as plt

from get_input import *
from get_parameters import *
import support_functions as suppf
import plotting_functions as plotf
######################################################

class interior_evolution:

#####################################################
    def __init__(self, body='Earth', inpfile='input.json'):        
        """Set initial conditions and model parameters"""
#####################################################

        # Set planet parameters and initial conditions
        get_input(self, body=body, inpfile=inpfile)

        # Set model parameters
        get_parameters(self)

######################################################################
    def calculate_evolution(self, tecmode):
        """Compute evolution of mantle and core temperatures"""
######################################################################
        
        # Number of timesteps
        n_steps = np.int(self.maxtime/self.dt)       
        
        # Initialize arrays
        self.t       = np.linspace(0, self.maxtime, n_steps+1)
        self.Tm      = np.zeros((n_steps+1))
        self.Tc      = np.zeros((n_steps+1))
        self.Tb      = np.zeros((n_steps+1))
        self.etam    = np.zeros((n_steps+1))
        self.etac    = np.zeros((n_steps+1))
        self.etab    = np.zeros((n_steps+1))
        self.delta_s = np.zeros((n_steps+1))
        self.delta_c = np.zeros((n_steps+1))
        self.qs      = np.zeros((n_steps+1))
        self.qc      = np.zeros((n_steps+1))
        self.Q_U238  = np.zeros((n_steps+1))
        self.Q_U235  = np.zeros((n_steps+1))
        self.Q_Th232 = np.zeros((n_steps+1))
        self.Q_K40   = np.zeros((n_steps+1))
        self.Q_tot   = np.zeros((n_steps+1))
        self.Ur      = np.zeros((n_steps+1))
        
        # Set initial conditions
        self.Tm[0]   = self.Tm0           # Upper mantle temperature
        self.Tc[0]   = self.Tc0           # Core temperature
        self.delta_c[0] = self.delta_c0   # Initial thickness of the bottom TBL
        self.delta_s[0] = self.delta_s0   # Initial thickness of the upper TBL

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
        for i in np.arange(0,n_steps): 

            # Mantle Viscosity
            if(i == 0):
               Pm = self.rhom*self.g*self.delta_s[i]    # Pressure at the base of the lithosphere
            else:   
               Pm = self.rhom*self.g*self.delta_s[i-1]  # Pressure at the base of the lithosphere from previous time step
            self.etam[i] = suppf.calculate_viscosity(self, self.Tm[i], Pm)             

            # Rayleigh number
            Ra = self.rhom*self.g*self.alpha*(self.Tm[i] - self.Ts)*D**3./(kappa*self.etam[i])    

            # Compute surface heat flux based on scaling laws for 
            # mobile lid convection
            if tecmode == 'ML':
                self.delta_s[i] = D*(self.Racrit/Ra)**self.beta                # Thickness of the upper TBL
                self.qs[i] = self.km*(self.Tm[i] - self.Ts) / self.delta_s[i]  # Surface heat flux

            # stagnant lid convection
            elif tecmode == 'SL':
                gamma = self.E*(self.Tm[i] - self.Ts)/(self.Rg*self.Tm[i]**2)                                 # Frank-Kamentzkii parameter      
                self.qs[i] = self.aa * self.km*(self.Tm[i] - self.Ts) / D * gamma**(-4./3.) * Ra**self.beta  # Surface heat flux
                self.delta_s[i] = self.km*(self.Tm[i] - self.Ts)/self.qs[i]                                  # Thickness of the upper TBL
            
 
            # Compute bottom heat flux with iterations
            if(i == 0):
                dc_i = self.delta_c0
            else: 
                dc_i = self.delta_c[i-1]

            # Determine bottom TBL thickness iteratively    
            dc = fixed_point(self.calculate_dc, dc_i, args = (self.delta_s[i], self.Tm[i], self.Tc[i]) ) 

            # Temperature at the top of the lower TBL going adiabatically down from Tm
            zb = self.Rp - (self.Rc + dc)           # depth of the lower TBL
            self.Tb[i] = suppf.calculate_adiabat(self, self.Tm[i], zb)

            # Internal and critical Rayleigh number
            deltaTm = (self.Tm[i] - self.Ts)
            deltaTc = (self.Tc[i] - self.Tb[i])
            Ra_int = self.rhom*self.g*self.alpha*( deltaTm + deltaTc )*D**3./(kappa*self.etam[i])    
            Racrit_int = 0.28*Ra_int**0.21

            # Pressure at the top of the bottom TBL (Pb) and at the CMB (Pc)
            Pb = self.rhom*self.g*zb
            Pc = self.rhom*self.g*(self.Rp - self.Rc)

            # Viscosity at the top of the bottom TBL (etab) and at the CMB (etac)
            self.etab[i] = suppf.calculate_viscosity(self, self.Tb[i], Pb)
            self.etac[i] = suppf.calculate_viscosity(self, self.Tc[i], Pc)

            # Update bottom TBL thickness
            self.delta_c[i] = (kappa*self.etab[i]*Racrit_int/(self.rhom*self.alpha*self.g*np.abs(self.Tc[i] - self.Tb[i])))**self.beta  # Thickenss of the lower TBL
            #self.delta_c[i] = min(self.delta_c[i], self.delta_s[i])

            # CMB heat flux
            self.qc[i] = self.km*(self.Tc[i] - self.Tb[i])/self.delta_c[i]    

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
            self.Tm[i+1] = self.Tm[i] + self.dt*(self.Q_tot[i]/self.cm + self.Qtidal/self.cm - (Ap*self.qs[i])/(Mm*self.cm) + (Ac*self.qc[i])/(Mm*self.cm))
            self.Tc[i+1] = self.Tc[i] - ((self.dt*Ac*self.qc[i])/(Mc*self.cc))                                                        
            
        # Write timeseries on file
        suppf.write_output_file(self)
            
            
#####################################################
    def calculate_dc(self, x, ds, Tm, Tc): 
        """"""
#####################################################

        D = self.Rp - self.Rc  
        kappa = self.km/self.rhom/self.cm                     
        Pm = self.rhom*self.g*ds
        etam = suppf.calculate_viscosity(self, Tm, Pm)   
        zb = self.Rp - (self.Rc + x)
        Tb = suppf.calculate_adiabat(self, Tm, zb)
        deltaTc = Tc - Tb
        deltaTm = Tm - self.Ts
        Ra_int = self.rhom*self.g*self.alpha*(deltaTm + deltaTc)*D**3./(kappa*etam)    
        Racrit_int = 0.28*Ra_int**0.21
        PPb = self.rhom*self.g*zb
        etab = self.etaref*np.exp( (self.E + PPb*self.V)/(self.Rg*Tb) - (self.E + self.Pref*self.V)/(self.Rg*self.Tref) )
        dc = ( kappa*etab*Racrit_int / (self.rhom*self.alpha*self.g*np.abs(deltaTc)) )**self.beta 
        #f = dc - x
        f = dc
        
        return f

