#######################################################################
def get_parameters(self):
    """Set some basic model constants and parameters"""
#######################################################################

    #################################
    # Time stuff
    #################################
    yr = 365.0*24.0*60.0*60.0   # 1 year in seconds
    self.maxtime = 4.55e9*yr    # Simulation time
    self.dt = 5e7*yr            # Time stepping
    
    #################################
    # Various constants
    #################################
    self.Racrit = 1100.       # Critical Rayleigh number
    self.beta = 1./3.         # Nu-Ra scaling exponent
    self.aa  = 0.5            # Prefactor for Nu-Ra stagnant-lid scaling
    self.delta_s0 = 5e3       # Initial thickness of top boundary layer     
    self.delta_c0 = 5e3       # Initial thickness of bottom boundary layer     
    self.cc = 840.0           # Core heat capacity (J/(kg K))
    self.cm = 1200.0          # Mantle heat capacity (J/(kg K))    
    self.alpha = 3e-5         # Mantle thermal expansivity (1/K)   
    self.km = 3.0 	          # Mantle thermal conductivity (W/(mK))         
    self.Tref = 1600.0        # Reference temperature (K)
    self.Pref = 3e9           # Reference pressure (Pa)
    self.E = 3e5              # Activation energy (J/mol)
    self.V = 4e-6             # Activation volume (m^3/mol)
    self.Rg = 8.3144          # Gas constant (J/(mol K))
    
    #################################
    # Radioactive elements 
    #################################
    # Half lives (s)
    self.tau_U238  = 4.47e9*yr    
    self.tau_U235  = 0.704e9*yr   
    self.tau_Th232 = 14.5e9*yr    
    self.tau_K40   = 1.2483e9*yr  
    # Isotope abundaces
    self.f_U238  = 0.9928
    self.f_U235  = 7.2e-3
    self.f_Th232 = 1.0
    self.f_K40   = 1.19e-4
    # Present-day heat productions (W/kg)
    self.H_U238  = 9.46e-5
    self.H_U235  = 5.69e-4
    self.H_K40   = 2.92e-5
    self.H_Th232 = 2.54e-5
