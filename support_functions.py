
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
    fc = 0.326      # Earth core mass fraction

    Rp = Rp_E*Mr**0.27
    Rc = Rc_E*Mr**0.247
    g = g_E*Mr**0.46
    rhom = Mr*M_E*(1-fc) / (4*np.pi/3*(Rp**3 - Rc**3))
    rhoc = Mr*M_E*fc / (4*np.pi/3*Rc**3)

    return Rp, Rc, g, rhom, rhoc
