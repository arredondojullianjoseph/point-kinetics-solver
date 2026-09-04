"""
Ramp_period_verification.py

Checks the accuracy of a ramp reactivity simulation by looking at its final 
growth rate. Once the ramp ends and reactivity stays flat at rho_final, the 
system should settle into the same stable period as step insertion. 
"""
import numpy as np
from scipy.integrate import solve_ivp

from point_kinetics import reactivity_ramp, kinetics_odes, steady_state_y0
from Inhour_verification import inhour_omega

def run_ramp(t_start, t_end, rho_final, t_span, n_points, max_step):
    """
    Runs the six-group solver with a ramp insertion and returns the
    time and power arrays.
    """
    y0 = steady_state_y0()
    
    def ramp_fn(t):
        return reactivity_ramp(t, t_start, t_end, rho_final)
        
    t_eval = np.linspace(*t_span, n_points)
    solution = solve_ivp(
        fun=kinetics_odes,
        t_span=t_span,
        y0=y0,
        t_eval=t_eval,
        method="Radau",  # Specifies the implicit solver: Radau
        max_step=max_step,
        atol=1e-10,
        rtol=1e-8,
        args=(ramp_fn,),  # passed through to kinetics_odes as reactivity_fn
    )
    return solution.t, solution.y[0]
    
def main():
    t_start, t_end, rho_final = 1.0, 3.0, 0.002
    t_span = (0.0, 60.0)
    t_fit_start = 40.0  # Late enough that only the slowest-decaying precursor group is still contributing
    
    time, power = run_ramp(t_start, t_end, rho_final, t_span, 6000, max_step=1e-3)
    
    # Keeps only late-term points
    mask = time >= t_fit_start
    
    # Linear fit of ln(power) vs time
    slope, _ = np.polyfit(time[mask], np.log(power[mask]), 1)
    
  
    period_numerical = 1.0 / slope
    period_analytical = 1.0 / inhour_omega(rho_final)
    
    # Calculate how far off our simulation is from the analytical result
    error = 100 * abs(period_numerical - period_analytical) / period_analytical
    
    # Print the results
    print("Ramp: {:.1f} pcm from t={:.1f}s to t={:.1f}s".format(rho_final * 1e5, t_start, t_end))
    print("Asymptotic period check (t>={:.0f}s, after ramp settles to rho_final):".format(t_fit_start))
    print("Analytical period (inhour root at rho_final): {:.3f} s".format(period_analytical))
    print("Numerical period (fit):                       {:.3f} s".format(period_numerical))
    print("Percent difference from analytical value:     {:.3f} %".format(error))
    
if __name__ == "__main__":
    main()
