"""
test_point_kinetics.py

Automated tests for the point kinetics project. Turns the same checks from
the verification scripts into pass/fail assertions.
"""
import numpy as np
from scipy.integrate import solve_ivp

from point_kinetics import reactivity, kinetics_odes, steady_state_y0
from Inhour_verification import inhour_omega
from Prompt_Jump_Approximation import prompt_jump_ratio
from Ramp_period_verification import run_ramp

def test_steady_state_y0_is_actually_steady():
    # At the reactor's starting steady state (rho = 0), the derivatives
    # should all be ~0 - nothing should be changing yet.
    y0 = steady_state_y0()
    dydt = kinetics_odes(0.0, y0, reactivity_fn=lambda t: 0.0)
    assert np.allclose(dydt, 0.0, atol=1e-8)

def test_inhour_period_matches_analytical():
    # Same setup as Inhour_verification.py, just as an assertion instead of a print
    rho = 0.002
    t_fit_start = 40.0
    y0 = steady_state_y0()
    t_span = (0.0, 60.0)
    t_eval = np.linspace(*t_span, 4000)
    solution = solve_ivp(
        fun=kinetics_odes,
        t_span=t_span,
        y0=y0,
        t_eval=t_eval,
        method="Radau",
        max_step=1e-3,
        atol=1e-10,
        rtol=1e-8,
    )
    time = solution.t
    power = solution.y[0]
    mask = time >= t_fit_start
    slope, _ = np.polyfit(time[mask], np.log(power[mask]), 1)
    period_numerical = 1.0 / slope
    period_analytical = 1.0 / inhour_omega(rho)
    error = 100 * abs(period_numerical - period_analytical) / period_analytical
    assert error < 1.0  

def test_prompt_jump_matches_analytical():
    # Same setup as Prompt_Jump_Approximation.py, just as an assertion instead of a print
    t_sample = 1.1
    y0 = steady_state_y0()
    t_span = (0.0, 2.0)
    t_eval = np.linspace(*t_span, 4000)
    solution = solve_ivp(
        fun=kinetics_odes,
        t_span=t_span,
        y0=y0,
        t_eval=t_eval,
        method="Radau",
        max_step=1e-3,
        atol=1e-10,
        rtol=1e-8,
    )
    time = solution.t
    power = solution.y[0]
    idx = np.searchsorted(time, t_sample)
    rho = reactivity(time[idx])
    n_numerical = power[idx]
    n_analytical = prompt_jump_ratio(rho)
    error = 100 * abs(n_numerical - n_analytical) / n_analytical
    assert error < 1.0  

def test_ramp_period_matches_analytical():
    # Same setup as Ramp_period_verification.py, just as an assertion instead of a print
    t_start, t_end, rho_final = 1.0, 3.0, 0.002
    t_span = (0.0, 60.0)
    t_fit_start = 40.0
    time, power = run_ramp(t_start, t_end, rho_final, t_span, 6000, max_step=1e-3)
    mask = time >= t_fit_start
    slope, _ = np.polyfit(time[mask], np.log(power[mask]), 1)
    period_numerical = 1.0 / slope
    period_analytical = 1.0 / inhour_omega(rho_final)
    error = 100 * abs(period_numerical - period_analytical) / period_analytical
    assert error < 1.0  

def test_ramp_is_converged():
    # Same setup as Ramp_convergence_verification.py, just as an assertion instead of a print
    t_start, t_end, rho_final = 1.0, 3.0, 0.002
    t_span = (0.0, 60.0)
    time_a, power_a = run_ramp(t_start, t_end, rho_final, t_span, 6000, max_step=1e-3)
    time_b, power_b = run_ramp(t_start, t_end, rho_final, t_span, 6000, max_step=5e-4)
    max_rel_diff = np.max(np.abs(power_a - power_b) / power_b)
    assert max_rel_diff < 1e-6 
