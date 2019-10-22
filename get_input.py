import json

#####################################################
def get_input(self, body='Earth', inpfile='input.json'):
        """Set input parameters: 

        Rp: planet radius (m) 
        Rc: core radius (m)
        g: surface gravity (m/s^2)  
        Ts: surface temperature (K) 
        rhom: mantle density (kg/m^3)
        rhoc: core density (kg/m^3)
        Tm0: initial mantle temperature (K)
        Tc0: initial core temperature (K)

        X_U: bulk abundance of U 
        X_Th: bulk abundance of Th
        X_K: bulk abundance of K
            or 
        Q0: initial specific heat production
        lam: decay constant

        """
#####################################################

        if body == 'Mars':

            self.name = 'Mars'    
            self.Rp   = 3400e3
            self.Rc   = 1800e3      
            self.g    = 3.7         
            self.rhom = 3500.0          
            self.rhoc = 7200.0          
            self.Ts   = 220.0       
            self.Tm0  = 1800.            
            self.Tc0  = 2300.            
            self.delta_s0 = 50e3            
            self.delta_c0 = 50e3            
            # Heat sources from Waenke & Dreibus (1994). Philos. Trans. R. Soc. London, A349, 2134–2137.
            self.X_U  = 16e-9 
            self.X_Th = 56e-9
            self.X_K  = 305e-6 
            self.Q0   = 0
            self.lam  = 0

        elif body == 'Mercury':

            self.name = 'Mercury'
            self.Rp   = 2440e3
            self.Rc   = 2020e3
            self.g    = 3.7         
            self.rhom = 3500.0          
            self.rhoc = 7200.0          
            self.Ts   = 440.0       
            self.Tm0  = 1700.            
            self.Tc0  = 2000.            
            self.delta_s0 = 5e3            
            self.delta_c0 = 5e3            
            # Heat sources from Hauck et al. (2019). In: Mercury - The view after MESSENGER. Ch. 19, 516–543, CUP.
            self.X_U  = 25e-9    
            self.X_Th = 44e-9
            self.X_K  = 368e-6 
            self.Q0   = 0
            self.lam  = 0
    	
        elif body == 'Moon':

            self.name = 'Moon'
            self.Rp   = 1740e3
            self.Rc   = 390e3
            self.g    = 1.6         
            self.rhom = 3500.0          
            self.rhoc = 7200.0          
            self.Ts   = 250.0       
            self.Tm0  = 1700.            
            self.Tc0  = 1900.            
            self.delta_s0 = 5e3            
            self.delta_c0 = 5e3            
            # Heat sources from Taylor (1982). Planetary Science: A Lunar Perspective, LPI.
            self.X_U  = 33e-9
            self.X_Th = 125e-9
            self.X_K  = 82.5e-6 
            self.Q0   = 0
            self.lam  = 0

        elif body == 'Earth':

            self.name = 'Earth'
            self.Rp   = 6370e3
            self.Rc   = 3485e3
            self.g    = 9.8         
            self.rhom = 4400.0          
            self.rhoc = 7200.0          
            self.Ts   = 280.0       
            self.Tm0  = 1700.            
            self.Tc0  = 3500.            
            self.delta_s0 = 5e3            
            self.delta_c0 = 5e3            
            # Heat sources from McDonough & Sun (1995), Chem. Geol., 223-253.
            self.X_U  = 20e-9
            self.X_Th = 80e-9 
            self.X_K  = 240e-6 
            self.Q0   = 0
            self.lam  = 0

        elif body == 'Venus':

            self.name = 'Venus'
            self.Rp   = 6050e3
            self.Rc   = 3110e3
            self.g    = 9.0         
            self.rhom = 4400.0          
            self.rhoc = 7200.0          
            self.Ts   = 735.0       
            self.Tm0  = 1700.            
            self.Tc0  = 3500.            
            self.delta_s0 = 5e3            
            self.delta_c0 = 5e3            
            # Heat sources from McDonough & Sun (1995), Chem. Geol., 223-253.
            self.X_U  = 20e-9
            self.X_Th = 80e-9 
            self.X_K  = 240e-6 
            self.Q0   = 0
            self.lam  = 0
         
        else:
            
            #Read input from file
            with open(inpfile, 'r') as inp:
                data = inp.read()

            inpvar = json.loads(data)

            self.name = inpvar['name']
            self.Rp   = inpvar['Rp']
            self.Rc   = inpvar['Rc']
            self.g    = inpvar['g']
            self.rhom = inpvar['rhom']
            self.rhoc = inpvar['rhoc']
            self.Ts   = inpvar['Ts'] 
            self.Tm0  = inpvar['Tm0']
            self.Tc0  = inpvar['Tc0'] 
            self.delta_s0 = inpvar['delta_s0']            
            self.delta_c0 = inpvar['delta_c0']            
            self.Q0   = inpvar['Q0']
            self.lam  = inpvar['lam']
            self.X_U  = inpvar['X_U']
            self.X_Th = inpvar['X_Th']
            self.X_K  = inpvar['X_K'] 
