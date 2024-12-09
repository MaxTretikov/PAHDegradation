import numpy as np
from scipy.optimize import differential_evolution

from pahdegradation.models import BistableModel, InducibleModel, time_to_threshold

class DegradationOptimizer:
    def __init__(self, params, t_span, t_eval):
        self.params = params
        self.t_span = t_span
        self.t_eval = t_eval
    
    def objective_func(self, x, y0, model_type='bistable'):
        p = self.params.copy()
        p['kdeg'] = x[0]
        p['decay_E'] = x[1]
        if model_type == 'bistable':
            p['decay_IclR'] = x[2]
            p['decay_HNS'] = x[3]
            p['K_pompC'] = x[4]
            p['K_I'] = x[5]
            p['K_H'] = x[6]
            # p['n'] = x[7]
            # p['m'] = x[8]

        # Choose which model to optimize
        if model_type == 'bistable':
            model = BistableModel(p)
            sol = model.simulate(y0, self.t_span, self.t_eval)
            PAH_sol = sol.y[-1, :]
        else:
            model = InducibleModel(p, k_E=0.5)
            sol = model.simulate(y0, self.t_span, self.t_eval)
            PAH_sol = sol.y[-1, :]

        tt = time_to_threshold(sol.t, PAH_sol, p['threshold'])
        return tt

    def optimize_parameters(self, y0, bounds, model_type):
        result = differential_evolution(
            lambda x: self.objective_func(x, y0, model_type), bounds, maxiter=20, disp=True
        )
        return result.x

    def bistable(self, y0, bounds):
        return self.optimize_parameters(y0, bounds, 'bistable')

    def inducible(self, y0, bounds):
        return self.optimize_parameters(y0, bounds, 'inducible')