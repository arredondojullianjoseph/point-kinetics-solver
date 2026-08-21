"""
Inhour_verification.py

This script verifies the six-group point-kinetics ODE solver by comparing it
with the analytical inhour equation. The transient is allowed to evolve until
the faster precursor groups have decayed away, and the asymptotic growth rate
is extracted from a fit to the power curve. That rate is then compared with
the corresponding root of the inhour equation.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from point_kinetics import lambda_decay, beta, beta_total, gen_time, reactivity, kinetics_odes


def inhour_omega(rho):
    """
    Solves the inhour equation for w (omega),
    which gives us the inverse of the asymptotic reactor period.
    """
    def f(omega):
        # The inhour equation rearranged to equal 0:
        # Lambda*omega + sum(beta_i * omega / (omega + lambda_i)) - rho = 0
        return (gen_time * omega
                + np.sum(beta * omega / (omega + lambda_decay))
                - rho)

    # Finds the root of the equation
    return brentq(f, 1e-8, 100.0)


def main():
    """
    Sets up the initial steady-state conditions, runs the numerical integration,
    and compares the resulting reactor period against the theoretical value.
    """
    rho = 0.002  # The reactivity step (200 pcm)
    t_fit_start = 40.0  # Late enough that only the slowest-decaying precursor group is still contributing

    # Assumes the reactor has been running at a steady power level of 1 for a while.
    # This means everything is in balance (the derivatives are zero),
    # so the precursors are in equilibrium (dC_i/dt = 0 => C_i = beta_i * n0 / (Lambda * lambda_i)).
    n0 = 1.0
    C0 = beta * n0 / (gen_time * lambda_decay)
    y0 = np.concatenate(([n0], C0))

    # Generates what happens over 60 seconds, grabbing 4000 data points
    t_span = (0.0, 60.0)
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

    # Keeps only late-term points
    mask = time >= t_fit_start

    # Linear fit of ln(power) vs time gives us a slope corresponding to omega
    slope, _ = np.polyfit(time[mask], np.log(power[mask]), 1)

    # Calculate the numerical period from our fit and the theoretical period
    period_numerical = 1.0 / slope
    period_analytical = 1.0 / inhour_omega(rho)

    # Calculate how far off our simulation is from the analytical result
    error = 100 * abs(period_numerical - period_analytical) / period_analytical

    # Print the results
    print("Step insertion: {:.1f} pcm".format(rho * 1e5))
    print("Analytical period (1/omega):                 {:.3f} s".format(period_analytical))
    print("Numerical period (fit, t>={:.0f}s):              {:.3f} s".format(t_fit_start, period_numerical))
    print("Percent difference from analytical value:    {:.3f} %".format(error))


if __name__ == "__main__":
    main()
