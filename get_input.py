import json

#####################################################
def get_input(self, body='Earth', inpfile='input.json'):
        """Set input parameters: 

        name: planet name
        core_cooling: consider core cooling ("yes") or neglect it ("no")
        Rp: planet radius (m) 
        Rc: core radius (m)
        g: surface gravity (m/s^2)  
        rhom: mantle density (kg/m^3)
        rhoc: core density (kg/m^3)
        Ts: surface temperature (K) 
        Tm0: initial mantle temperature (K)
        Tc0: initial core temperature (K)
        etaref: reference viscosity at Tref and Pref set in get_parameters.py (Pa s)
        X_U: bulk abundance of U 
        X_Th: bulk abundance of Th
        X_K: bulk abundance of K
            or 
        Q0: initial specific heat production (W/kg)
        lam: decay constant (s)

        """
#####################################################

        if body == 'Mercury':

            self.name = 'Mercury'
            self.tectonics = 'SL'
            self.core_cooling = 'yes'
            self.var_alpha = 'yes'
            # Core radius, core density and mantle density from the representative model of Margot et al. (2019). In: Mercury - The view after MESSENGER. Ch. 4, 85-113, CUP.
            self.Rp   = 2440e3
            self.Rc   = 2024e3   
            self.g    = 3.7         
            self.rhom = 3295.0          
            self.rhoc = 7034.0          
            self.Ts   = 440.0       
            self.Tm0  = 1700.            
            self.Tc0  = 1900.  
            self.etaref = 1e19
            # Heat sources from Hauck et al. (2019). In: Mercury - The view after MESSENGER. Ch. 19, 516–543, CUP.
            self.X_U  = 25e-9    
            self.X_Th = 44e-9
            self.X_K  = 368e-6 
            self.Q0   = 0
            self.lam  = 0
            self.Qtidal = 0            
            
        elif body == 'Venus':
            
            self.name = 'Venus'
            self.tectonics = 'SL'
            self.core_cooling = 'yes'
            self.var_alpha = 'no'
            self.Rp   = 6050e3
            self.Rc   = 3186e3
            self.g    = 8.87         
            self.rhom = 4400.0          
            self.rhoc = 10100.0          
            self.Ts   = 730.0       
            self.Tm0  = 1700.            
            self.Tc0  = 4200.
            self.etaref = 1e19           
            # Heat sources from Kaula (Icarus, 1999)
            self.X_U  = 21e-9
            self.X_Th = 86e-9 
            self.X_K  = 153e-6 
            self.Q0   = 0
            self.lam  = 0
            self.Qtidal = 0
    
        elif body == 'Earth':

            self.name = 'Earth'
            self.tectonics = 'PT'
            self.core_cooling = 'yes'
            self.var_alpha = 'yes'
            self.Rp   = 6370e3
            self.Rc   = 3480e3
            self.g    = 9.8         
            self.rhom = 4460.0          
            self.rhoc = 10640.0          
            self.Ts   = 288.0       
            self.Tm0  = 1700.            
            self.Tc0  = 4000.
            self.etaref = 1e21
            # Heat sources from McDonough & Sun (1995), Chem. Geol., 223-253.
            self.X_U  = 20e-9
            self.X_Th = 80e-9 
            self.X_K  = 240e-6 
            self.Q0   = 0
            self.lam  = 0
            self.Qtidal = 0
            
        elif body == 'Moon':

            self.name = 'Moon'
            self.tectonics = 'SL'
            self.core_cooling = 'yes'
            self.var_alpha = 'yes'
            self.Rp   = 1740e3
            self.Rc   = 330e3
            self.g    = 1.6         
            self.rhom = 3400.0          
            self.rhoc = 7400.0          
            self.Ts   = 270.0       
            self.Tm0  = 1700.            
            self.Tc0  = 1998. 
            self.etaref = 1e19
            # Heat sources from Taylor (1982). Planetary Science: A Lunar Perspective, LPI.
            self.X_U  = 33e-9
            self.X_Th = 125e-9
            self.X_K  = 82.5e-6 
            self.Q0   = 0
            self.lam  = 0
            self.Qtidal = 0                        
            
        elif body == 'Mars':

            self.name = 'Mars'    
            self.tectonics = 'SL'
            self.core_cooling = 'yes'
            self.var_alpha = 'yes'
            self.Rp   = 3390e3
            self.Rc   = 1850e3   # Core radius from Plesa et al. (2018). Geophys. Res. Lett., 45(22), 12198-12209.       
            self.g    = 3.7         
            self.rhom = 3500.0          
            self.rhoc = 7200.0          
            self.Ts   = 220.0       
            self.Tm0  = 1700.            
            self.Tc0  = 2160.     
            self.etaref = 1e20
            # Heat sources from Waenke & Dreibus (1994). Philos. Trans. R. Soc. London, A349, 2134–2137.
            self.X_U  = 16e-9 
            self.X_Th = 56e-9
            self.X_K  = 305e-6 
            self.Q0   = 0
            self.lam  = 0
            self.Qtidal = 0
                               
        else:
            
            #Read input from input file
            with open(inpfile, 'r') as inp:
                data = inp.read()
            inpvar = json.loads(data)

            self.name         = inpvar['name']
            self.tectonics    = inpvar['tectonics']
            self.core_cooling = inpvar['core_cooling']
            self.var_alpha    = inpvar['var_alpha']
            self.Rp           = inpvar['Rp']
            self.Rc           = inpvar['Rc']
            self.g            = inpvar['g']
            self.rhom         = inpvar['rhom']
            self.rhoc         = inpvar['rhoc']
            self.Ts           = inpvar['Ts'] 
            self.Tm0          = inpvar['Tm0']
            self.Tc0          = inpvar['Tc0'] 
            self.etaref       = inpvar['etaref'] 
            self.Q0           = inpvar['Q0']
            self.lam          = inpvar['lam']            
            self.X_U          = inpvar['X_U']
            self.X_Th         = inpvar['X_Th']
            self.X_K          = inpvar['X_K']
            self.Qtidal       = inpvar['Qtidal']            
        
        return