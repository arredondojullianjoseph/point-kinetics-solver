"""
Prompt_jump_verification.py

Verifies the point-kinetics ODE solver using the Prompt Jump Approximation. 
It runs a short simulation to capture the immediate power spike right after a step reactivity insertion, then compares results to the Prompt Jump Approximation. 
"""
import numpy as np
from scipy.integrate import solve_ivp
from point_kinetics import beta_total, reactivity, kinetics_odes, steady_state_y0

def prompt_jump_ratio(rho):
    """
    Solves the Prompt Jump Approximation for n(0+)/n0,
    which gives us the power ratio right after a step insertion.
    """
    return beta_total / (beta_total - rho)
    
def main():
    """
    Sets up the initial steady-state conditions, runs the numerical integration,
    and compares the resulting power jump against the theoretical value.
    """
    t_sample = 1.1  # Sampled shortly after the step, once the jump has settled
    
    y0 = steady_state_y0()
    
    # Generates what happens over 2 seconds, grabbing 4000 data points
    t_span = (0.0, 2.0)
    t_eval = np.linspace(*t_span, 4000)
    
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
    
    # Pulls out the results
    time = solution.t
    power = solution.y[0]
    
    # Finds the grid point closest to t_sample and grabs the reactivity there
    idx = np.searchsorted(time, t_sample)
    rho = reactivity(time[idx])
    
    # Calculate the numerical jump from our simulation and the theoretical jump
    n_numerical = power[idx]
    n_analytical = prompt_jump_ratio(rho)
    
    # Calculate how far off our simulation is from the analytical result
    error = 100 * abs(n_numerical - n_analytical) / n_analytical
    
    # Print the results
    print("Step insertion: {:.1f} pcm".format(rho * 1e5))
    print("Analytical n(0+)/n0:               {:.4f}".format(n_analytical))
    print("Numerical  n(t={:.2f}s)/n0:          {:.4f}".format(t_sample, n_numerical))
    print("Percent difference from analytical value:    {:.3f} %".format(error))
    
if __name__ == "__main__":
    main()
