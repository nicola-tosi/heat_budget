
import numpy as np
#####################################################


#####################################################
def initialize_arrays(self):
    """Initialize arrays for thermal evolution"""
#####################################################

    self.t       = np.linspace(0, self.maxtime, self.n_steps+1)
    self.Tprofile= np.zeros((self.n_steps+1,self.n_layers))
    self.Tm      = np.zeros((self.n_steps+1))
    self.Tc      = np.zeros((self.n_steps+1))
    self.Tb      = np.zeros((self.n_steps+1))
    self.etam    = np.zeros((self.n_steps+1))
    self.etac    = np.zeros((self.n_steps+1))
    self.etab    = np.zeros((self.n_steps+1))
    self.delta_s = np.zeros((self.n_steps+1))
    self.delta_c = np.zeros((self.n_steps+1))
    self.qs      = np.zeros((self.n_steps+1))
    self.qc      = np.zeros((self.n_steps+1))
    self.Q_U238  = np.zeros((self.n_steps+1))
    self.Q_U235  = np.zeros((self.n_steps+1))
    self.Q_Th232 = np.zeros((self.n_steps+1))
    self.Q_K40   = np.zeros((self.n_steps+1))
    self.Q_tot   = np.zeros((self.n_steps+1))
    self.Ra      = np.zeros((self.n_steps+1))
    self.Ur      = np.zeros((self.n_steps+1))
    
    return

#####################################################
def set_initial_conditions(self):
    """Set initial
       upper mantle temperature Tm0
       core temperature Tc0
       thickness of the upper thermal boundary layer delta_s0
       thickness of the lower thermal boundary layer delta_c0
    """
#####################################################

    self.Tm[0]   = self.Tm0           # Initial upper mantle temperature
    self.Tc[0]   = self.Tc0           # Initial core temperature
    self.delta_c[0] = self.delta_c0   # Initial thickness of the bottom TBL
    self.delta_s[0] = self.delta_s0   # Initial thickness of the upper TBL
    
    return

#####################################################
def calculate_viscosity(s, T, P):
    """Calculate T- and P-dependent viscosity based on Arrhenius law for diffusion creep"""
#####################################################

    eta = s.etaref*np.exp( (s.E + P*s.V)/(s.Rg*T) - (s.E + s.Pref*s.V)/(s.Rg*s.Tref) )

    return eta 

#####################################################################################
def calculate_viscosity_karato(s, T, P):
    """Calculate T- and P-dependent viscosity according to Karato & Wu (1993)"""
#####################################################################################

    B = 6.1e-19
    m = 2.5
    d = 2e-3
    E = 3e5
    V = 2.5e-6
    eta = 1/(2*B) * d**m * np.exp( (E+P*V)/(s.Rg*T) )

    return eta 

###################################################################################################
def calculate_thermal_expansivity(T, P):
    """Calculate P- and T-dependent thermal expansivity using the lower mantle parametrization of 
    Tosi et al. (PEPI, 2013).
    T: temperature in K
    P: pressure in GPa
    alpha: thermal expansivity in 1/K
    """
###################################################################################################

    a0 = 2.68e-5  # 1/K
    a1 = 2.77e-9  # 1/K**2
    a2 = -1.21    # K
    a3 = 8.63e-3  # 1/GPa
    
    alpha = (a0 + a1*T + a2/T**2) * np.exp(-a3*P)

    return alpha

###################################################################################################
def calculate_adiabat(s, n_layers, Tm, Pm, Pb):
    """"""
###################################################################################################
     
    Tprof = np.zeros((n_layers))
    Pprof = np.linspace(Pm, Pb, n_layers)
    dP = (np.diff(Pprof))[0]
    dP = Pprof[1]-Pprof[0]
    
    Tprof[0] = Tm
    Pprof[0] = Pm

    if (s.var_alpha == 'yes'):        
        for i in np.arange(0, n_layers-1):         
            alpha = calculate_thermal_expansivity(Tprof[i], Pprof[i]/1e9)
            Tprof[i+1] = Tprof[i] + alpha*Tprof[i]/s.rhom/s.cm * dP
    else:
        Tprof = Tm * np.exp(s.alpha * Pprof / s.rhom /s.cm)
            
    return Tprof

#####################################################
def calculate_dc(x, s, ds, Tm, Tc): 
    """"""
#####################################################

    D = s.Rp - s.Rc  
    kappa = s.km/s.rhom/s.cm                     
    Pm = s.rhom*s.g*ds
    etam = calculate_viscosity(s, Tm, Pm)   
    zb = s.Rp - (s.Rc + x)
    Tb = calculate_adiabat_old(s, Tm, zb)
    deltaTc = Tc - Tb
    deltaTm = Tm - s.Ts
    Ra_int = s.rhom*s.g*s.alpha*(deltaTm + deltaTc)*D**3./(kappa*etam)    
    Racrit_int = 0.28*Ra_int**0.21
    # Racrit_int = s.Racrit           
    PPb = s.rhom*s.g*zb
    Tbmean = (Tb + Tc)/2
    etab = s.etaref*np.exp( (s.E + PPb*s.V)/(s.Rg*Tbmean) - (s.E + s.Pref*s.V)/(s.Rg*s.Tref) )        
    dc = ( kappa*etab*Racrit_int / (s.rhom*s.alpha*s.g*np.abs(deltaTc)) )**s.beta 
        
    return dc

#####################################################
def calculate_dc_root(x, s, ds, Tm, Tc): 
    """"""
#####################################################

    D = s.Rp - s.Rc  
    kappa = s.km/s.rhom/s.cm                     
    Pm = s.rhom*s.g*ds
    etam = calculate_viscosity(s, Tm, Pm)   
    zb = s.Rp - (s.Rc + x)
    Pb = s.rhom*s.g*zb
    Pc = s.rhom*s.g*(s.Rp - s.Rc)
    Tbp = calculate_adiabat(s, s.n_layers, Tm, Pm, Pb) 
    Tb = Tbp[-1]
    deltaTc = Tc - Tb
    deltaTm = Tm - s.Ts
    Ra_int = s.rhom*s.g*s.alpha*(deltaTm + deltaTc)*D**3./(kappa*etam)    
    Racrit_int = 0.28*Ra_int**0.21
    Tbmean = (Tb + Tc)/2
    Pbmean = (Pb + Pc)/2
    etab = s.etaref*np.exp( (s.E + Pbmean*s.V)/(s.Rg*Tbmean) - (s.E + s.Pref*s.V)/(s.Rg*s.Tref) )        
    dc = ( kappa*etab*Racrit_int / (s.rhom*s.alpha*s.g*np.abs(deltaTc)) )**s.beta 
    f = dc - x
        
    return f

#####################################################
def calculate_ds(x, s, Tm, gamma):
    """"""
#####################################################

    D = s.Rp - s.Rc
    kappa = s.km/s.rhom/s.cm        
    Pm = s.rhom*s.g*x
    etam = calculate_viscosity(s, Tm, Pm)
    Ra = s.rhom*s.g*s.alpha*(Tm - s.Ts)*D**3/(kappa*etam)
    ds = 1. / ( s.aa/D * gamma**(-4./3.) * Ra**s.beta )

    return ds

#####################################################
def calculate_ds_root(x, s, Tm, gamma): 
    """"""
#####################################################

    D = s.Rp - s.Rc  
    kappa = s.km/s.rhom/s.cm                       
    Pm = s.rhom*s.g*x
#     etam = calculate_viscosity(s, Tm, Pm)
    etam = calculate_viscosity(s, Tm, 0)    
    Ra = s.rhom*s.g*s.alpha*(Tm - s.Ts)*D**3/(kappa*etam)    
    ds = 1. / ( s.aa/D * gamma**(-4./3.) * Ra**s.beta )
    f = ds - x
    
    return f

    
#####################################################
def calculate_ds_fsolve(x, s, Tm, gamma): 
    """"""
#####################################################

    D = s.Rp - s.Rc  
    kappa = s.km/s.rhom/s.cm                       
    Pm = s.rhom*s.g*x
    etam = calculate_viscosity(s, Tm, Pm)             
    Ra = s.rhom*s.g*s.alpha*(Tm - s.Ts)*D**3/(kappa*etam)    
    ds = 1. / ( s.aa/D * gamma**(-4./3.) * Ra**s.beta )
    f = ds - x
    
    return f

#####################################################
def initialize_heatproduction(tbp, tau):
    """Scale present-day heat production back by tbp"""
#####################################################

    f = np.exp(np.log(2.)*tbp/tau)

    return f

##############################################################
def calculate_radiodecay(f, X0, H, tau, t):
    """Calculate radioactive decay for a specific isotope"""
##############################################################

    g = f*X0*H*np.exp(-np.log(2.)*t/tau)

    return g

########################################################################
def calculate_radiodecay_simple(Q0, lam, t):
    """Calculate radioactive decay based on a single decay constant"""
########################################################################

    g = Q0*np.exp(-lam*t)

    return g

######################################################################################
def calculate_dry_solidus(P):
    """Calculate dry solidus of Katz et al. (Gcubed, 2003).
    Solidus in K, input pressure in GPa"""
######################################################################################

    a1 = 1358.85
    a2 = 132.9
    a3 = -5.1
    Tsol = a1 + a2*P + a3*P*P

    return Tsol

######################################################################################
def calculate_dry_liquidus(P):
    """Calculate dry liquidus of Katz et al. (Gcubed, 2003).
    Liquidus in K, input pressure in GPa"""
######################################################################################

    b1 = 2053.15
    b2 = 45.0
    b3 = -2.0
    Tliq = b1 + b2*P + b3*P*P

    return Tliq

######################################################################################
def melting_idx(s):
    """Determine the indeces of the timeseries where the mantle temperature is above
       and below the solidus"""
######################################################################################    
    rho = s.rhom
    g = s.g
    ds = s.delta_s
    # Pressure (in Pa) at the base of the lid
    P = rho*g*ds
    # Indices where Tm >= Tsol and Tm < Tsol
    idx_above_solidus = np.where( s.Tm[:-1] >= calculate_dry_solidus(P[:-1]/1e9) )
    idx_below_solidus = np.where( s.Tm[:-1] < calculate_dry_solidus(P[:-1]/1e9) )    
    
    return idx_above_solidus, idx_below_solidus

######################################################################################
def melting_range(s):
    """"""
######################################################################################    
    rho = s.rhom
    g = s.g
    ds = s.delta_s
    # Pressure (in Pa) at the base of the lid
    P = rho*g*ds
    
    # Arrays containing temperature above and below the solidus
    # For Stagnant lid, consider solidus at the base of the lithosphere
    if (s.tectonics == 'SL'):
        T_below_solidus = np.ma.masked_where(s.Tm[:-1] >= calculate_dry_solidus(P[:-1]/1e9), s.Tm[:-1])
        T_above_solidus = np.ma.masked_where(s.Tm[:-1] < calculate_dry_solidus(P[:-1]/1e9), s.Tm[:-1])
    # For Mobile lid, consider surface solidus
    elif (s.tectonics == 'ML'):
        T_below_solidus = np.ma.masked_where(s.Tm[:-1] >= calculate_dry_solidus(0), s.Tm[:-1])
        T_above_solidus = np.ma.masked_where(s.Tm[:-1] < calculate_dry_solidus(0), s.Tm[:-1])

    return T_below_solidus, T_above_solidus

######################################################################################
def calculate_initial_CMBtemperature(Pcmb):
    """ """
######################################################################################
    X_Fe = 0.2
    Tcmb = 5400*(Pcmb/140)**0.48 / (1 - np.log(1. - X_Fe))

    return Tcmb

######################################################################################
def write_output_file(s, outfile):
    """ """
######################################################################################

    yrs = 365.0*24.0*60.0*60.0   # 1 year in seconds
#     outfile = s.name + ".out"
    outdata = np.array([s.t/yrs/1e6, s.Q_tot, s.Ur, s.Tm, s.Tb, s.Tc, s.etam, s.etab, s.etac, s.qs, s.qc, s.delta_s, s.delta_c])
    outdata = outdata.T
    with open(outfile, 'w+') as outfile_id:
        np.savetxt(outfile_id, outdata, fmt=['%.4e', '%.4e', '%.4e', '%.4e', '%.4e', '%.4e', '%.4e', '%.4e', '%.4e', '%.4e', '%.4e', '%.4e', '%.4e'])
    
    return

######################################################################################
def mass_radius_relations_withFe(Mr, X_Fe):
    """Given the mass planetary mass to Earth mass ratio Mr = Mp/M_E, and iron
    mass fraction, calculate planetary radius (Rp), 
    core radius (Rc), surface gravity (g), mantle density (rhom) and core density (rhoc) 
    using mass-radius relations from L. Noack (unpublished)
    """
######################################################################################

    M_E = 5.972e24  # Earth mass
    G = 6.67e-11    # Gravitational constant

    Rp = (7e3 - 1.8e3*X_Fe)*Mr**0.3
    rhoc = 12300.*Mr**0.2
    Rc = 1e-3*( (X_Fe*Mr*M_E)/ (4./3*np.pi*rhoc) )**(1./3) 
    rhom = (1 - X_Fe)*Mr*M_E/ (4./3*np.pi*((Rp*1e3)**3 - (Rc*1e3)**3))
    g = G*M_E*Mr / (Rp*1e3)**2
    
    print('Rp = ', Rp, 'km')
    print('Rc = ', Rc, 'km')
    print('g  = ', g, 'm/s^2')
    print('rhom = ', rhom, 'kg/m^3')
    print('rhoc = ', rhoc, 'kg/m^3')

    return Rp, Rc, g, rhom, rhoc

######################################################################################
def mass_radius_relations(Mr):
    """Given the mass planetary mass to Earth mass ratio Mr = Mp/M_E,
    calculate planetary radius (Rp), core radius (Rc), surface gravity (g),
    mantle density (rhom) and core density (rhoc) using mass-radius relations
    from Valencia et al. (Icarus, 2006) 
    """
######################################################################################

    Rp_E = 6371e3   # Earth raidus
    Rc_E = 3480e3   # Earth's core radius
    M_E = 5.972e24  # Earth mass
    g_E = 9.81      # Earth surface gravity
    X_Fe = 0.326    # Earth core mass fraction

    Rp = Rp_E*Mr**0.27
    Rc = Rc_E*Mr**0.247
    g = g_E*Mr**0.46
    rhom = Mr*M_E*(1-X_Fe) / (4*np.pi/3*(Rp**3 - Rc**3))
    rhoc = Mr*M_E*X_Fe / (4*np.pi/3*Rc**3)
    
    print('Rp = ', Rp/1e3, 'km')
    print('Rc = ', Rc/1e3, 'km')
    print('g  = ', g, 'm/s^2')
    print('rhom = ', rhom, 'kg/m^3')
    print('rhoc = ', rhoc, 'kg/m^3')

    return Rp, Rc, g, rhom, rhoc