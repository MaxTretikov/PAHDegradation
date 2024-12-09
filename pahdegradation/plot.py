import matplotlib.pyplot as plt

def plot_results(sol_bistable, sol_inducible, PAH_bistable, PAH_inducible, threshold, optimization=False):
    fig, ax = plt.subplots()
    ax.plot(sol_bistable.t, PAH_bistable, label='Bistable Model PAH' + (' (Optimized)' if optimization else ''))
    ax.plot(sol_inducible.t, PAH_inducible, label='Inducible Model PAH' + (' (Optimized)' if optimization else ''))

    ax.axhline(threshold, color='r', linestyle='--', label='Threshold')
    ax.set_xlabel("Time")
    ax.set_ylabel("PAH Concentration")
    ax.legend()
    plt.title("PAH Degradation Comparison" + (' (Optimized)' if optimization else ''))
    plt.show()