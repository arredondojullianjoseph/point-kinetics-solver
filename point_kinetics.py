"""
point_kinetics.py
This script simulates how a nuclear reactor's power changes when we suddenly
introduce reactivity (a step insertion) using the six-group point kinetics model.
We're using scipy's Radau method to solve the ODEs because these equations are
stiff (The prompt neutron term decays on a ~1e-4 s timescale while the precursors evolve over seconds).
The script then plots normalized power vs. time.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# From Lamarsh & Baratta, "Introduction to Nuclear Engineering," 3rd ed., Table 3.5
lambda_decay = np.array([0.0124, 0.0305, 0.111, 0.301, 1.14, 3.01])  # How fast the precursors decay (in s^-1)
beta = np.array(
    [0.000215, 0.001424, 0.001274, 0.002568, 0.000748, 0.000273]
)  # The fraction of neutrons that are delayed in each group
beta_total = beta.sum()  # The total delayed neutron fraction (~0.0065 for U-235)
gen_time = 1e-4  # The prompt neutron generation time (Lambda, in seconds).


def reactivity(t):
    """
    Step insertion: reactor sits at delayed critical until t = 1 s,
    then a 200 pcm step is inserted and maintained.
    """
    if t < 1.0:
        return 0.0
    else:
        return 0.002 # ~30% of beta (650 pcm for U-235), so we stay safely below prompt critical 


def kinetics_odes(t, y):
    """
    y[0] is the reactor power (neutron population, n).
    y[1:] are the concentrations of our delayed neutron precursors (C_i).
    """
    n = y[0]
    C = y[1:]
    rho = reactivity(t)

    dydt = np.zeros_like(y)

    # Prompt term (rho - beta)/Lambda * n, plus the delayed source from precursors decaying back into neutrons.
    dydt[0] = (rho - beta_total) / gen_time * n + np.sum(lambda_decay * C)

    # Each precursor group is grown by fission at a rate of beta_i/Lambda * n, and depletes by its own decay constant.
    for i in range(6):
        dydt[i + 1] = (beta[i] / gen_time) * n - lambda_decay[i] * C[i]

    return dydt


def main():
    # Assumes the reactor has been running at a steady power level of 1 for a while.
    # This means everything is in balance (the derivatives are zero),
    # so the precursors are in equilibrium (dC_i/dt = 0 => C_i = beta_i * n0 / (Lambda * lambda_i)).
    n0 = 1.0
    C0 = beta * n0 / (gen_time * lambda_decay)
    y0 = np.concatenate(([n0], C0))

    # Generates what happens over 10 seconds, grabbing 1000 data points
    t_span = (0.0, 10.0)
    t_eval = np.linspace(*t_span, 1000)

    # Solving the ODEs
    solution = solve_ivp(
        fun=kinetics_odes,
        t_span=t_span,
        y0=y0,
        t_eval=t_eval,
        method="Radau",  # Specifies the implicit solver: Radau
        max_step=1e-3,
        atol=1e-10,
        rtol=1e-8,
    )

    # Pulls out the results for the graph
    time = solution.t
    power = solution.y[0]

    # Plots
    plt.figure(figsize=(10, 6))
    plt.plot(time, power, "b-", linewidth=2, label="Relative reactor power (n)")

    # Draws a line where the 200 pcm reactivity step is inserted
    plt.axvline(x=1.0, color="r", linestyle="--", alpha=0.5, label="Reactivity step insertion")

    # Labels graph
    plt.title("Point Reactor Kinetics Transient (Step Reactivity Insertion)", fontsize=14)
    plt.xlabel("Time (seconds)", fontsize=12)
    plt.ylabel("Normalized power", fontsize=12)
    plt.grid(True, which="both", linestyle="--", alpha=0.7)

    # Uses log scale due to reactor power being able to swing by orders of magnitude
    plt.yscale("log")
    plt.legend(fontsize=12)
    plt.tight_layout()

    # Saves plot
    plt.savefig("step_response.png", dpi=150)
    plt.show()


if __name__ == "__main__":
    main()
