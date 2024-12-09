import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution
import matplotlib.pyplot as plt

def time_to_threshold(time, PAH, threshold):
    below = np.where(PAH < threshold)[0]
    return time[below[0]] if len(below) > 0 else np.inf

class BaseDegradationModel:
    def __init__(self, params):
        """
        params is a dictionary containing parameters:
        {
            'kdeg': float, # degradation rate constant for PAH by E
            'decay_IclR': float,
            'decay_HNS': float,
            'decay_E': float,
            'alpha_pscpA': float,
            'alpha_pompC': float,
            'K_pompC': float,
            'K_I': float,
            'K_H': float,
            'n': float,
            'm': float,
            'threshold': float, # PAH safe threshold
            ...
        }
        """
        self.params = params
    
    def degradation_rate(self, E, PAH):
        # Simple linear degradation rate
        kdeg = self.params['kdeg']
        return kdeg * E * PAH
    
    def simulate(self, y0, t_span, t_eval):
        # Must be implemented by subclass: defines the ODE system
        pass

class BistableModel(BaseDegradationModel):
    def __init__(self, params):
        super().__init__(params)

    def ode_system(self, t, y):
        # y = [IclR, HNS, E, PAH]
        IclR, HNS, E, PAH = y
        
        # Extract parameters
        p = self.params
        alpha_pompC = p['alpha_pompC']
        alpha_pscpA = p['alpha_pscpA']
        K_pompC = p['K_pompC']
        K_I = p['K_I']
        K_H = p['K_H']
        n = p['n']
        m = p['m']
        
        decay_IclR = p['decay_IclR']
        decay_HNS = p['decay_HNS']
        decay_E = p['decay_E']

        # promoter activities
        # pscpA activity (Dormant promoter) = alpha_pscpA / (1 + (HNS/K_H)^m)
        pscpA_activity = alpha_pscpA / (1.0 + (HNS/K_H)**m)
        
        # pompC activity (Active promoter) = (alpha_pompC * PAH^n / (K_pompC^n + PAH^n)) * 1/(1+(IclR/K_I)^m)
        pompC_activity = alpha_pompC * (PAH**n/(K_pompC**n + PAH**n)) * (1.0/(1.0+(IclR/K_I)**m))
        
        # ODEs:
        dIclR_dt = pscpA_activity - decay_IclR*IclR
        dHNS_dt = pompC_activity - decay_HNS*HNS
        
        # In the bistable model, E is produced by the active promoter
        dE_dt = pompC_activity - decay_E*E
        
        # PAH degradation
        dPAH_dt = -self.degradation_rate(E, PAH)
        
        return [dIclR_dt, dHNS_dt, dE_dt, dPAH_dt]

    def simulate(self, y0, t_span, t_eval):
        sol = solve_ivp(self.ode_system, t_span, y0, t_eval=t_eval, max_step=0.1)
        return sol


class InducibleModel(BaseDegradationModel):
    def __init__(self, params, k_E=0.1):
        super().__init__(params)
        self.k_E = k_E  # proportionality constant for enzyme production from PAH

    def ode_system(self, t, y):
        # y = [E, PAH]
        E, PAH = y
        
        # In the inducible model, E production is directly proportional to PAH
        # No bistable switch, no IclR/HNS needed
        decay_E = self.params['decay_E']
        dE_dt = self.k_E * PAH - decay_E*E
        
        # PAH degradation
        dPAH_dt = -self.degradation_rate(E, PAH)
        
        return [dE_dt, dPAH_dt]

    def simulate(self, y0, t_span, t_eval):
        sol = solve_ivp(self.ode_system, t_span, y0, t_eval=t_eval, max_step=0.1)
        return sol