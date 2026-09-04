"""
Ramp_convergence_verification.py

Verifies the ramp reactivity insertion by checking numerical convergence:
halves the max step size and confirms the solution barely changes, showing
the integration is converged for the ramp case rather than just settled on
some step-size-dependent answer. Reuses run_ramp from Ramp_period_verification.py.
"""

import numpy as np

from ramp_period_verification import run_ramp

def main():
    t_start, t_end, rho_final = 1.0, 3.0, 0.002
    t_span = (0.0, 60.0)
    
    # Runs two simulations, one with normal step size and another with half, then compares.
    time_a, power_a = run_ramp(t_start, t_end, rho_final, t_span, 6000, max_step=1e-3)
    time_b, power_b = run_ramp(t_start, t_end, rho_final, t_span, 6000, max_step=5e-4)
    max_rel_diff = np.max(np.abs(power_a - power_b) / power_b)
    
    # Print the results
    print("Ramp: {:.1f} pcm from t={:.1f}s to t={:.1f}s".format(rho_final * 1e5, t_start, t_end))
    print("Convergence check (max_step 1e-3 vs 5e-4):")
    print("Max relative difference in power:            {:.2e}".format(max_rel_diff))
    
if __name__ == "__main__":
    main()
