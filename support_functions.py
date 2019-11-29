
import numpy as np

#####################################################
def calculate_viscosity(s, T, P):
    """Calculate T- and P-dependent viscosity based on Arrhenius law for diffusion creep"""
#####################################################

    eta = s.etaref*np.exp( (s.E + P*s.V)/(s.Rg*T) - (s.E + s.Pref*s.V)/(s.Rg*s.Tref) )

    return eta 

#####################################################
def calculate_adiabat(s, Tm, z):
    """Calculate adiabatic temperature at depth z starting from temperature Tm"""
#####################################################

    Tad = Tm * np.exp(s.alpha * s.g * z / s.cm)

    return Tad

#####################################################
def calculate_dc(x, s, ds, Tm, Tc): 
    """"""
#####################################################

    D = s.Rp - s.Rc  
    kappa = s.km/s.rhom/s.cm                     
    Pm = s.rhom*s.g*ds
    etam = calculate_viscosity(s, Tm, Pm)   
    zb = s.Rp - (s.Rc + x)
    Tb = calculate_adiabat(s, Tm, zb)
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
    idx_above_solidus = np.where( s.Tm[:-1] >= calculate_dry_solidus(P[:-1]/1e9))
    idx_below_solidus = np.where( s.Tm[:-1] < calculate_dry_solidus(P[:-1]/1e9))    
    
    return idx_above_solidus, idx_below_solidus

######################################################################################
def calculate_initial_CMBtemperature(Pcmb):
    """ """
######################################################################################
    X_Fe = 0.2
    Tcmb = 5400*(Pcmb/140)**0.48 / (1 - np.log(1. - X_Fe))

    return Tcmb

######################################################################################
def write_output_file(s):
    """ """
######################################################################################

    yrs = 365.0*24.0*60.0*60.0   # 1 year in seconds
    outfile = s.name + ".out"
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
    using mass-radius relations from Noack (unpublished)
    """
######################################################################################

    M_E = 5.972e24  # Earth mass
    G = 6.67e-11    # Gravitational constant

    Rp = (7000 - 18*X_Fe)*Mr**0.3
    rhoc = 12300.*Mr**0.2
    Rc = (X_Fe*Mr*M_E)**(1./3) / (4./3*np.pi*rhoc)
    rhom = (1 - X_Fe)*Mr*M_E/ (4./3*np.pi*((Rp*1e3)**3 - (Rc*1e3)**3))
    g = G*M_E*Mr / (Rp*1e3)**2

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

    return Rp, Rc, g, rhom, rhoc
