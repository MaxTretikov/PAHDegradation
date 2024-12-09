import numpy as np

from pahdegradation.models import BistableModel, InducibleModel, time_to_threshold
from pahdegradation.optimize import DegradationOptimizer
from pahdegradation.plot import plot_results

def run(optimization=False, steps=50):
    # Set parameters (default values)
    params = {
        'kdeg': 0.01,     # rate constant for PAH degradation
        'decay_IclR': 0.1,
        'decay_HNS': 0.1,
        'decay_E': 0.05,
        'alpha_pscpA': 1.0,
        'alpha_pompC': 5.0,
        'K_pompC': 1.0,
        'K_I': 1.0,
        'K_H': 1.0,
        'n': 2.0,
        'm': 2.0,
        'threshold': 0.1
    }

    # Initial conditions
    IclR0 = 2.0
    HNS0 = 0.1
    E0 = 0.0
    PAH0 = 5.0

    y0_bistable = [IclR0, HNS0, E0, PAH0]
    y0_inducible = [E0, PAH0]

    t_span = (0, steps)
    t_eval = np.linspace(0, steps, 1001)
    threshold = params['threshold']

    if optimization:
        opt = DegradationOptimizer(params, t_span, t_eval)

        # Optimization bounds
        bounds_bistable = [(0.001, 0.1), (0.01, 0.1), (0.01, 0.1), (0.01, 0.1), (0.01, 0.1), (0.1, 10.0), (0.1, 10.0), (0.1, 10.0)]
        bounds_inducible = [(0.001, 0.1), (0.01, 0.1)]

        # Optimize parameters
        best_params_bistable = opt.bistable(y0_bistable, bounds_bistable)
        print("Optimized parameters for bistable model:", best_params_bistable)

        best_params_inducible = opt.inducible(y0_inducible, bounds_inducible)
        print("Optimized parameters for inducible model:", best_params_inducible)

        # Update parameters with optimized values
        params.update({
            'kdeg': best_params_bistable[0],
            'decay_E': best_params_bistable[1],
            'decay_IclR': best_params_bistable[2],
            'decay_HNS': best_params_bistable[3],
            'K_pompC': best_params_bistable[4],
            'K_I': best_params_bistable[5],
            'K_H': best_params_bistable[6]
        })

    # Run bistable model
    bistable_model = BistableModel(params)
    sol_bistable = bistable_model.simulate(y0_bistable, t_span, t_eval)

    # Run inducible model
    inducible_model = InducibleModel(params, k_E=0.5)
    sol_inducible = inducible_model.simulate(y0_inducible, t_span, t_eval)

    # Extract solutions
    PAH_bistable = sol_bistable.y[-1, :]
    PAH_inducible = sol_inducible.y[-1, :]

    # Find times to reach threshold
    t_thresh_bistable = time_to_threshold(sol_bistable.t, PAH_bistable, threshold)
    t_thresh_inducible = time_to_threshold(sol_inducible.t, PAH_inducible, threshold)

    print(f"Time for bistable model to reach threshold: {t_thresh_bistable}")
    print(f"Time for inducible model to reach threshold: {t_thresh_inducible}")

    # Plot results
    plot_results(sol_bistable, sol_inducible, PAH_bistable, PAH_inducible, threshold, optimization)