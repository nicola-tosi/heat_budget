import numpy as np
from matplotlib import pyplot as plt

class interior_evolution:

#####################################################
    def __init__(self, body):        
        """Set parameters and initial conditions"""
#####################################################
        
        if body == 'Mars':
            self.name = 'Mars'    
            self.Rp = 3400e3
            self.Rc = 1700e3      
            self.g  = 3.7         
            self.Ts = 220.0       
            self.rhom = 3500.0          
            self.rhoc = 7200.0          
            self.H = 23e-12      
            self.Tm0 = 1700.            
            self.Tc0 = 2500.            
            # Dreibus & Waenke composition model
            self.X_U  = 16e-9 
            self.X_Th = 56e-9
            self.X_K  = 305e-6 
    	
        elif body == 'Mercury':
            self.name = 'Mercury'
            self.Rp = 2440e3
            self.Rc = 2020e3
            self.g = 3.7         
            self.Ts = 440.0       
            self.rhom = 3500.0          
            self.rhoc = 7200.0          
            self.H = 23e-12      
            self.Tm0 = 1700.            
            self.Tc0 = 1900.            
            # To be changed!
            self.X_U  = 16e-9    
            self.X_Th = 56e-9
            self.X_K  = 305e-6 
    	
        elif body == 'Moon':
            self.name = 'Moon'
            self.Rp = 1740e3
            self.Rc = 390e3
            self.g = 1.6         
            self.Ts = 250.0       
            self.rhom = 3500.0          
            self.rhoc = 7200.0          
            self.H = 23e-12      
            self.Tm0 = 1700.            
            self.Tc0 = 1900.            
            # To be changed!
            self.X_U  = 16e-9 
            self.X_Th = 56e-9
            self.X_K  = 305e-6 
            	
        elif body == 'Earth':
            self.name = 'Earth'
            self.Rp = 6370e3
            self.Rc = 3485e3
            self.g = 9.8         
            self.Ts = 280.0       
            self.rhom = 3500.0          
            self.rhoc = 7200.0          
            self.H = 23e-12      
            self.Tm0 = 1700.            
            self.Tc0 = 3500.            
            # To be changed!
            self.X_U  = 16e-9 
            self.X_Th = 56e-9
            self.X_K  = 305e-6 
              	
        elif body == 'Venus':
            self.name = 'Venus'
            self.Rp = 6050e3
            self.Rc = 3110e3
            self.g = 9.0         
            self.Ts = 735.0       
            self.rhom = 3500.0          
            self.rhoc = 7200.0          
            self.H = 23e-12      
            self.Tm0 = 1700.            
            self.Tc0 = 3500.            
            # To be changed!
            self.X_U  = 16e-9 
            self.X_Th = 56e-9
            self.X_K  = 305e-6 

        #################################
        # Radiogenic elements 
        #################################
        self.Racrit = 1100. 
        self.cc = 840.0    
        self.cm = 1142.0          
        self.alpha = 3e-5       
        self.lamb = 1.42e-17
        self.E = 3e5    
        self.R = 8.3144 
        self.etaref = 1e21         
        self.Tref = 1600.0
        self.km = 4. 	                 

        #################################
        # Radiogenic elements 
        #################################
        # Half lives
        yr = 365.0*24.0*60.0*60.0   # 1 year in seconds
        self.tau_U238  = 4.47e9*yr
        self.tau_U235  = 0.704e9*yr
        self.tau_Th232 = 14.5e9*yr
        self.tau_K40   = 1.2483e9*yr
        # Isotopes abundaces
        self.f_U238  = 0.9928
        self.f_U235  = 7.2e-3
        self.f_Th232 = 1.0
        self.f_K40   = 1.19e-4
        # Decay energies (J)
        self.H_U238  = 9.46e-5
        self.H_U235  = 5.69e-4
        self.H_K40   = 2.92e-5
        self.H_Th232 = 2.64e-5
            
######################################################################
    def calculate_evolution(self, tecmode):
        """Compute evolution of mantle and core temperatures"""
######################################################################
        
        # 
        yr = 365.0*24.0*60.0*60.0   # 1 year in seconds
        time = 4.5e9*yr             # age of the Solar System 4.5 Gyr
        Myr = 1e6*yr                # 1 Myr in seconds
        dt = 1e7*yr                 # fixed timestep of 10 Myr
        n_steps = np.int(time/dt)   # number of timesteps    
        
        # Initialize arrays
        self.t       = np.linspace(0,time,n_steps)
        self.Tm      = np.linspace(0,time,n_steps)
        self.Tc      = np.linspace(0,time,n_steps)
        self.etam    = np.linspace(0,time,n_steps)
        self.etac    = np.linspace(0,time,n_steps)
        self.delta_s = np.linspace(0,time,n_steps)
        self.delta_c = np.linspace(0,time,n_steps)
        self.qs      = np.linspace(0,time,n_steps)
        self.qc      = np.linspace(0,time,n_steps)
        self.Q_U238  = np.linspace(0,time,n_steps)
        self.Q_U235  = np.linspace(0,time,n_steps)
        self.Q_Th232 = np.linspace(0,time,n_steps)
        self.Q_K40   = np.linspace(0,time,n_steps)
        self.Q_tot   = np.linspace(0,time,n_steps)
        
        # Set initial conditions
        self.Tm[0]   = self.Tm0   # Upper mantle temperature
        self.Tc[0]   = self.Tc0   # Core temperature
        self.etam[0] = self.etaref*np.exp(self.E/self.R*(self.Tref-self.Tm0)/(self.Tref*self.Tm0)) # Upper mantle viscosity
        self.etac[0] = self.etaref*np.exp(self.E/self.R*(self.Tref-self.Tc0)/(self.Tref*self.Tc0)) # Upper mantle viscosity
        
        # Set/compute some parameters
        beta = 1./3.    # Nu-Ra scaling exponent
        c  = 0.5        # Prefactor for Nu-Ra stagnant-lid scaling
        D  = self.Rp - self.Rc   # Mantle thickness                                
        M  = 4./3.*np.pi*((self.Rp**3.-self.Rc**3)*self.rhom + (self.Rc**3)*self.rhoc) # planet mass
        Mm = 4./3.*np.pi*((self.Rp**3.-self.Rc**3)*self.rhom) # mantle mass
        Mc = 4./3.*np.pi*(self.Rc**3*self.rhoc)               # core mass
        Ap = 4.*np.pi*self.Rp**2                              # surface area             
        Ac = 4.*np.pi*self.Rc**2                              # core area             
        kappa = self.km/self.rhom/self.cm                     # thermal diffusivity  

        # Scale initial concentrations of radiogenic isotpes back in time
        ee = lambda tbp, tau : np.exp(np.log(2.)*tbp/tau)
        X0_U238 = self.X_U*ee(time,self.tau_U238)
        X0_U235 = self.X_U*ee(time,self.tau_U235)
        X0_Th232 = self.X_Th*ee(time,self.tau_Th232)
        X0_K40 = self.X_K*ee(time,self.tau_K40)

        # Radiodecay law
        qrad = lambda f, X0, H, tau, t: f*X0*H*np.exp(-np.log(2.)*t/tau)

        for i in np.arange(1,n_steps): 

            # Radioactive decay
            self.Q_U238[i]  = qrad(self.f_U238, X0_U238, self.H_U238, self.tau_U238, self.t[i])
            self.Q_U235[i]  = qrad(self.f_U235, X0_U235, self.H_U235, self.tau_U235, self.t[i])
            self.Q_Th232[i] = qrad(self.f_Th232, X0_Th232, self.H_Th232, self.tau_Th232, self.t[i])
            self.Q_K40[i]   = qrad(self.f_K40, X0_K40, self.H_K40, self.tau_K40, self.t[i])
            self.Q_tot[i]   = self.Q_U238[i] + self.Q_U235[i] + self.Q_Th232[i] + self.Q_K40[i]  

            # Rayleigh number
            Ra = self.rhom*self.g*self.alpha*(self.Tm[i-1] - self.Ts)*D**3./(kappa*self.etam[i-1])          

            # Compute surface heat flux based on scaling laws for 
            # mobile lid convection
            if tecmode == 'ML':
                self.delta_s[i] = D*(self.Racrit/Ra)**beta  # Thickness of the upper TBL
                self.qs[i] = self.km*(self.Tm[i-1] - self.Ts) / self.delta_s[i]  # Surface heat flux
            # stagnant lid convection
            elif tecmode == 'SL':
                gamma = self.E*(self.Tm[i-1] - self.Ts)/(self.R*self.Tm[i-1]**2)               # Frank-Kamentzkii parameter      
                self.qs[i] = c * self.km*(self.Tm[i-1] - self.Ts) / D *gamma**(-4/3)*Ra**beta  # Surface heat flux
                self.delta_s[i] = self.km*(self.Tm[i-1] - self.Ts)/self.qs[i]                  # Thickness of the upper TBL

            
            # Compute bottom heat flux
            self.delta_c[i] = (kappa*self.etac[i-1]*self.Racrit/(self.rhom*self.alpha*self.g*np.abs(self.Tc[i-1] - self.Tm[i-1])))**beta  # Thickenss of the lower TBL
            self.qc[i] = self.km*(self.Tc[i-1] - self.Tm[i-1])/self.delta_c[i]                                       # CMB heat flux

            # Advance in time mantle and CMB temperature
            self.Tm[i] = self.Tm[i-1] + dt*(self.Q_tot[i]/self.cm - (Ap*self.qs[i])/(Mm*self.cm) + (Ac*self.qc[i])/(Mm*self.cm)) # Mantle temperature
            self.Tc[i] = self.Tc[i-1] - ((dt*Ac*self.qc[i])/(Mc*self.cc))                                    # Core temperature

            # Compute the corresponding viscosities
            self.etam[i] = self.etaref*np.exp(self.E/self.R*(self.Tref-self.Tm[i])/(self.Tref*self.Tm[i]))   # Mantle viscosity
            self.etac[i] = self.etaref*np.exp(self.E/self.R*(self.Tref-self.Tc[i])/(self.Tref*self.Tc[i]))   # CMB viscosity



#####################################################
    def plot_evolution(self): 
        """Make some plots"""
#####################################################

        yr = 365.0*24.0*60.0*60.0   # 1 year in seconds
            
        fig = plt.figure(figsize=(12,12))    
        small_size = 10
        medium_size = 15
        bigger_size = 20
        
        plt.rc('font', size=small_size)          # controls default text sizes
        plt.rc('axes', titlesize=small_size)     # fontsize of the axes title
        plt.rc('axes', labelsize=medium_size)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=medium_size)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=medium_size)    # fontsize of the tick labels
        plt.rc('legend', fontsize=medium_size)    # legend fontsize
        plt.rc('figure', titlesize=bigger_size)  # fontsize of the figure title


        ############################################
        # Heat production
        ############################################
        ax = fig.add_subplot(3,2,1)
        ax.plot(self.t[1:]/yr/1e6, self.Q_U238[1:]*1e12, label='$^{238}U$')
        ax.plot(self.t[1:]/yr/1e6, self.Q_U235[1:]*1e12, label='$^{235}U$')
        ax.plot(self.t[1:]/yr/1e6, self.Q_Th232[1:]*1e12, label='$^{232}Th$')
        ax.plot(self.t[1:]/yr/1e6, self.Q_K40[1:]*1e12, label='$^{40}K$')
        ax.plot(self.t[1:]/yr/1e6, self.Q_tot[1:]*1e12, label='Total')
        ax.set_xlabel('Time [Myr]')
        ax.set_ylabel('Heat production [pW/kg]')
        ax.grid()
        ax.legend()

        ############################################
        # Mantle and CMB temperatures
        ############################################
        ax = fig.add_subplot(3,2,3)
        ax.plot(self.t/yr/1e6, self.Tm, label='$T_m$')
        ax.plot(self.t/yr/1e6, self.Tc, label='$T_c$')
        ax.set_xlabel('Time [Myr]')
        ax.set_ylabel('Temperature [K]')
        ax.grid()
        ax.legend()

        ############################################
        # Mantle and CMB viscosities
        ############################################
        ax = fig.add_subplot(3,2,4)
        ax.plot(self.t/yr/1e6, self.etam, label='$\eta_m$')
        ax.plot(self.t/yr/1e6, self.etac, label='$\eta_c$')
        ax.set_xlabel('Time [Myr]')
        ax.set_ylabel('Viscosity [Pa s]')
        ax.set_yscale('log')
        ax.grid()
        ax.legend()

        ############################################
        # Surface and CMB heat flux
        ############################################
        ax = fig.add_subplot(3,2,5)
        ax.plot(self.t[1:]/yr/1e6, self.qs[1:]*1e3, label='$q_s$')
        ax.plot(self.t[1:]/yr/1e6, self.qc[1:]*1e3, label='$q_c$')
        ax.set_xlabel('Time [Myr]')
        ax.set_ylabel('Heat flux [mW/m$^2$]')
        ax.grid()
        ax.legend()

        ############################################
        # Top and bottom boundary layer thicknesses
        ############################################
        ax = fig.add_subplot(3,2,6)
        ax.plot(self.t[1:]/yr/1e6, self.delta_s[1:]/1e3, label='$\delta_s$')
        ax.plot(self.t[1:]/yr/1e6, self.delta_c[1:]/1e3, label='$\delta_c$')
        ax.set_xlabel('Time [Myr]')
        ax.set_ylabel('Boundary layer thickness [km]')
        ax.grid()
        ax.legend()

 
        plt.tight_layout()

