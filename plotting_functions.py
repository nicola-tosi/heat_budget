#####################################################
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import ticker, cm, colors
import scipy as scp

import support_functions as suppf
#####################################################

#####################################################
def plot_evolution(s): 
    """Plot time evolution of main quantities"""
#####################################################

    yr = 365.0*24.0*60.0*60.0   # 1 year in seconds
        
    fig = plt.figure(figsize=(12,12))    
    small_size = 10
    medium_size = 15
    bigger_size = 20

    lw = 2.5 

    mcolor = 'tab:blue'
    bcolor = 'tab:orange'
    ccolor = 'tab:red'
    
    plt.rc('font', size=small_size)          # controls default text sizes
    plt.rc('axes', titlesize=small_size)     # fontsize of the axes title
    plt.rc('axes', labelsize=medium_size)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=medium_size)   # fontsize of the tick labels
    plt.rc('ytick', labelsize=medium_size)   # fontsize of the tick labels
    plt.rc('legend', fontsize=medium_size)   # legend fontsize
    plt.rc('figure', titlesize=bigger_size)  # fontsize of the figure title

    nx_panels = 3
    ny_panels = 2
    n = 0

    ############################################
    # Heat production
    ############################################
    n = n + 1
    ax = fig.add_subplot(nx_panels,ny_panels,n)
    ax.plot(s.t[:-1]/yr/1e6,  s.Q_U238[:-1]*1e12, label='$^{238}$U', lw=lw)
    ax.plot(s.t[:-1]/yr/1e6,  s.Q_U235[:-1]*1e12, label='$^{235}$U', lw=lw)
    ax.plot(s.t[:-1]/yr/1e6, s.Q_Th232[:-1]*1e12, label='$^{232}$Th', lw=lw)
    ax.plot(s.t[:-1]/yr/1e6,   s.Q_K40[:-1]*1e12, label='$^{40}$K', lw=lw)
    ax.plot(s.t[:-1]/yr/1e6,   s.Q_tot[:-1]*1e12, label='Total', lw=lw)
    ax.set_xlabel('Time [Myr]')
    ax.set_ylabel('Heat production [pW/kg]')
    ax.grid()
    ax.legend(loc=1)

    ############################################
    # Urey ratio
    ############################################
    n = n + 1
    ax = fig.add_subplot(nx_panels,ny_panels,n)
    ax.plot(s.t[1:-1]/yr/1e6,  s.Ur[1:-1], label='', lw=lw)
    ax.set_xlabel('Time [Myr]')
    ax.set_ylabel('Urey ratio')
    ax.grid()

    ##################################################################################
    # Mantle temperature, CMB temperature and temperature at the top of the lower TBL
    ##################################################################################
    n = n + 1
    ax = fig.add_subplot(nx_panels,ny_panels,n)
    ax.plot(s.t[:-1]/yr/1e6, s.Tm[:-1], label='$T_m$', color=mcolor, lw=lw)
    ax.plot(s.t[:-1]/yr/1e6, s.Tb[:-1], label='$T_b$', color=bcolor, lw=lw)
    ax.plot(s.t[:-1]/yr/1e6, s.Tc[:-1], label='$T_c$', color=ccolor, lw=lw)
    ax.set_xlabel('Time [Myr]')
    ax.set_ylabel('Temperature [K]')
    ax.grid()
    ax.legend(loc=1)

    ############################################
    # Mantle and CMB viscosities
    ############################################
    n = n + 1
    ax = fig.add_subplot(nx_panels,ny_panels,n)
    ax.plot(s.t[1:-1]/yr/1e6, s.etam[1:-1], label='$\eta_m$', color=mcolor, lw=lw)
    ax.plot(s.t[1:-1]/yr/1e6, s.etab[1:-1], label='$\eta_b$', color=bcolor, lw=lw)
    ax.plot(s.t[1:-1]/yr/1e6, s.etac[1:-1], label='$\eta_c$', color=ccolor, lw=lw)
    ax.set_xlabel('Time [Myr]')
    ax.set_ylabel('Viscosity [Pa s]')
    ax.set_yscale('log')
    ax.grid()
    ax.legend(loc=1)

    ############################################
    # Surface and CMB heat flux
    ############################################
    n = n + 1
    ax = fig.add_subplot(nx_panels,ny_panels,n)
    ax.plot(s.t[1:-1]/yr/1e6, s.qs[1:-1]*1e3, label='$q_s$', color=mcolor, lw=lw)
    ax.plot(s.t[1:-1]/yr/1e6, s.qc[1:-1]*1e3, label='$q_c$', color=ccolor, lw=lw)
    
    ax.set_xlabel('Time [Myr]')
    ax.set_ylabel('Heat flux [mW/m$^2$]')
    ax.grid()
    ax.legend(loc=1)

    ############################################
    # Top and bottom boundary layer thicknesses
    ############################################
    n = n + 1
    ax = fig.add_subplot(nx_panels,ny_panels,n)
    ax = fig.add_subplot(3,2,6)
    ax.plot(s.t[1:-1]/yr/1e6, s.delta_s[1:-1]/1e3, label='$\delta_s$', color=mcolor, lw=lw)
    ax.plot(s.t[1:-1]/yr/1e6, s.delta_c[1:-1]/1e3, label='$\delta_c$', color=ccolor, lw=lw)
    ax.set_xlabel('Time [Myr]')
    ax.set_ylabel('Boundary layer thickness [km]')
    ax.grid()
    ax.legend(loc=1)
 
    plt.tight_layout()


####################################################################
def plot_profiles(s, time, plot_solidus=0, plot_liquidus=0):
    """Plot temperature and viscosity profiles at a given time"""
####################################################################
    
    yr = 365.0*24.0*60.0*60.0
    etamax = 1e26
    Pmax = 10e9
    # Index of the time array closest to input time
    i = np.abs(s.t/yr/1e6 - time).argmin()
    
    fig = plt.figure(figsize=(10,5))    
    small_size = 10
    medium_size = 15
    bigger_size = 20

    Tcolor = 'tab:blue'
    Tsolcolor = 'tab:red'
    Tliqcolor = 'tab:red'
    etacolor = 'tab:blue'

    lw = 2.5 
        
    plt.rc('font', size=small_size)          # controls default text sizes
    plt.rc('axes', titlesize=small_size)     # fontsize of the axes title
    plt.rc('axes', labelsize=medium_size)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=medium_size)   # fontsize of the tick labels
    plt.rc('ytick', labelsize=medium_size)   # fontsize of the tick labels
    plt.rc('legend', fontsize=medium_size)   # legend fontsize
    plt.rc('figure', titlesize=bigger_size)  # fontsize of the figure title
 
    nx_panels = 1
    ny_panels = 2
    n = 0
    
    nptb = 10
    npta = 50
    npts = 10
    npt = nptb + npta + npts

    ########################
    # Temperature
    ########################
    n = n + 1
    ax = fig.add_subplot(nx_panels,ny_panels,n)

    # Adiabatic mantle
    ra = np.linspace(s.Rc + s.delta_c[i], s.Rp - s.delta_s[i],  npta)
    Ta = s.Tm[i]*np.exp(s.alpha * s.g * (s.Rp - ra) / s.cm )    
    ax.plot(Ta, ra/1e3, '-', color=Tcolor, lw=lw)

    # Lower boundary layer
    rb = np.linspace(s.Rc, s.Rc + s.delta_c[i], nptb)
    Trb = np.linspace(s.Tc[i], Ta[0], nptb)
    ax.plot(Trb, rb/1e3, '-', color=Tcolor, lw=lw)    

    # Upper boundary layer
    rs = np.linspace(s.Rp - s.delta_s[i], s.Rp, npts)
    Trs = np.linspace(Ta[npta-1], s.Ts, npts)
    ax.plot(Trs, rs/1e3, '-', color=Tcolor, lw=lw)    
    
    # Solidus
    if plot_solidus == 'yes':
        Psol = np.linspace(Pmax,0,npt)
        rsol = s.Rp - Psol/(s.rhom*s.g)
        Tsol = suppf.calculate_dry_solidus(Psol/1e9)
        ax.plot(Tsol, rsol/1e3, '--', color=Tsolcolor, lw=lw)    
        
    # Liquidus
    if plot_liquidus == 'yes':
        Pliq = np.linspace(Pmax,0,npt)
        rliq = s.Rp - Pliq/(s.rhom*s.g)
        Tliq = suppf.calculate_dry_liquidus(Pliq/1e9)
        ax.plot(Tliq, rsol/1e3, '--', color=Tliqcolor, lw=lw)    
             
    ax.grid()    
    ax.set_ylabel('Radius [km]')
    ax.set_xlabel('Temperautre [K]')

    ax.text(0.05, 0.1, 'Time =' + str(time) + ' Myr', horizontalalignment='left', verticalalignment='center',
            transform=ax.transAxes, fontsize = medium_size)
    
    ########################
    #  Viscosity 
    ########################
    n = n + 1
    ax = fig.add_subplot(nx_panels,ny_panels,n)

    # Lower boundary layer
    Prb = s.rhom*s.g*(s.Rp - rb)
    etab = suppf.calculate_viscosity(s, Trb, Prb)
    ax.plot(etab, rb/1e3, '-', color=etacolor, lw=lw)

    # Adiabatic mantle  
    Pa = s.rhom*s.g*(s.Rp  - ra)             
    etaa = suppf.calculate_viscosity(s, Ta, Pa)
    ax.plot(etaa, ra/1e3, '-', color=etacolor, lw=lw)
    
    # Upper boundary layer
    Prs = s.rhom*s.g*(s.Rp - rs)   
    etas = suppf.calculate_viscosity(s, Trs, Prs)
    ax.plot(etas, rs/1e3, '-', color=etacolor, lw=lw)
     
    ax.grid()    
    ax.set_xscale('log')
    ax.set_xlim(right=etamax)
    ax.set_ylabel('Radius [km]')
    ax.set_xlabel('Viscosity [Pa s]')
        
    plt.tight_layout()    


####################################################################
def plot_profiles_evolution(s):
    """Plot evolution of the temperature and viscosity profiles"""
####################################################################

    yr = 365.0*24.0*60.0*60.0   # 1 year in seconds
    etamax = 1e25
    colormap_1 = plt.cm.plasma
    colormap_2 = plt.cm.plasma_r
            
    fig = plt.figure(figsize=(15,4))  
    n = 0
    nx_panels = 1
    ny_panels = 2
    lw = 1
    
    small_size = 10
    medium_size = 15
    bigger_size = 20
    plt.rc('font', size=small_size)          # controls default text sizes
    plt.rc('axes', titlesize=small_size)     # fontsize of the axes title
    plt.rc('axes', labelsize=medium_size)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=medium_size)   # fontsize of the tick labels
    plt.rc('ytick', labelsize=medium_size)   # fontsize of the tick labels
    plt.rc('legend', fontsize=medium_size)   # legend fontsize
    plt.rc('figure', titlesize=bigger_size)  # fontsize of the figure title
    
    nt = np.size(s.t)-1
    nptb = 10
    npta = 150
    npts = 10
    npt = nptb + npta + npts
    Tprof = np.zeros((npt,nt))
    etaprof = np.zeros((npt,nt))
    meltzone_top = np.zeros((nt))
    meltzone_bot = np.zeros((nt))
    
    r = np.linspace(s.Rc, s.Rp, npt)
    Pr = np.linspace(s.rhom*s.g*(s.Rp - s.Rc)/1e9, 0., npt)
    Tsolr = suppf.calculate_dry_solidus(Pr)
    
    # Prepare arrays for contour plots
    for i in np.arange(0, nt):

        # Adiabatic mantle
        ra = np.linspace(s.Rc + s.delta_c[i], s.Rp - s.delta_s[i], npta)
        Ta = s.Tm[i]*np.exp(s.alpha * s.g * (s.Rp - s.delta_s[i] - ra) / s.cm )  
        Pa = s.rhom*s.g*(s.Rp - ra)                     
        etaa = suppf.calculate_viscosity(s, Ta, Pa)    
                    
        # Lower boundary layer
        rb = np.linspace(s.Rc, s.Rc + s.delta_c[i], nptb)
        Trb = np.linspace(s.Tc[i], Ta[0], nptb)
        Prb = s.rhom*s.g*(s.Rp - (s.Rc + s.delta_c[i]))
        etab = suppf.calculate_viscosity(s, Trb, Prb)
        
        # Upper boundary layer
        rs = np.linspace(s.Rp - s.delta_s[i], s.Rp, npts)
        Trs = np.linspace(Ta[npts-1], s.Ts, npts)
        Prs = s.rhom*s.g*(s.Rp - rs)
        etas = suppf.calculate_viscosity(s, Trs, Prs)    
                
        Ttemp = np.concatenate((Trb, Ta, Trs))
        etatemp = np.concatenate((etab, etaa, etas))
        rtemp = np.concatenate((rb, ra, rs))
        
        Tprof[:,i] = np.interp(r, rtemp, Ttemp) 
        etaprof[:,i] = np.interp(r, rtemp, etatemp) 
        
        # Partial melt zone
        idx_meltzone = np.argwhere(np.diff(np.sign(Tprof[:,i] - Tsolr))).flatten()
        if (np.size(idx_meltzone) >= 2) :
            meltzone_top[i] = r[idx_meltzone[-1]]
            meltzone_bot[i] = r[idx_meltzone[-2]]
        else:    
            meltzone_top[i] = None
            meltzone_bot[i] = None
    
    ############################################
    # Temperature and melt zone
    ############################################
    n = n + 1
    ax = fig.add_subplot(nx_panels,ny_panels,n)
    
    levsT_cont = np.arange(s.Ts, np.max(Tprof), 10)
    levsT_disc = np.arange(s.Ts, np.max(Tprof), 300)
    
    cf = ax.contourf(s.t[:-1]/yr/1e6, r/1e3, Tprof, levsT_cont, cmap=colormap_1)
    cb = plt.colorbar(cf, extend='both', ticks = levsT_disc)
    
    ax.plot(s.t[:-1]/yr/1e6, meltzone_top/1e3, '--', color='black', lw=lw)
    ax.plot(s.t[:-1]/yr/1e6, meltzone_bot/1e3, '--', color='black', lw=lw)
    
    cb.ax.set_ylabel('Temperature [K]')
    ax.set_xlabel('Time [Myr]')
    ax.set_ylabel('Radius [km]')
    ax.grid()
    
    ############################################
    # Viscosity
    ############################################
    n = n + 1
    ax = fig.add_subplot(nx_panels,ny_panels,n)
    
    # Continous colorscale for contour plot
    levse_exp_cont = np.arange(np.floor(np.log10(etaprof.min())-1), np.ceil(np.log10(etamax)+1), 0.02)
    levse_cont = np.power(10, levse_exp_cont)
    # Discrete intervals for colorbar ticks
    levse_exp_disc = np.arange(np.floor(np.log10(etaprof.min())-1), np.ceil(np.log10(etamax)+1), 2)
    levse_disc = np.power(10, levse_exp_disc)
    
    cf = ax.contourf(s.t[:-1]/yr/1e6, r/1e3, etaprof, levse_cont, norm=colors.LogNorm(vmin = etaprof.min(), vmax = etamax), extend='max', cmap=colormap_2)
    cb = plt.colorbar(cf, extend='max', ticks=levse_disc)
    
    cb.ax.set_ylabel('Viscosity [Pa s]')
    ax.set_xlabel('Time [Myr]')
    ax.set_ylabel('Radius [km]')
    ax.grid()
    
    plt.tight_layout()
    
