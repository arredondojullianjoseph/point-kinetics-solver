"""
point_kinetics.py
Solves the six-group point reactor kinetics equations for a step
reactivity insertion.

The script evaluates the transient using scipy's stiff ODE integrator (Radau) and plots normalized power vs time using matplotlib.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

Lambda = np.array([0.0124, 0.0305, 0.111, 0.301, 1.14, 3.01]) 
Beta = np.array([0.000215, 0.001424, 0.001274, 0.002568, 0.000748, 0.000273]) 
Beta_total = sum(Beta) 
Lambda_gen = 1e-4 

def reactivity(t):
    if t < 1:
        return 0
    else:
        return 0.002

def Kinetic_ODEs(t, y):
    n = y[0]
    C = y[1:]
    rho = reactivity(t)
    dydt = np.zeros_like(y)
    dydt[0] = ((rho - Beta_total) / Lambda_gen) * n + np.sum(Lambda * C)
    for i in range(6):
        dydt[i + 1] = (Beta[i] / Lambda_gen) * n - Lambda[i] * C[i]
    return dydt

def main():
    n_0 = 1
    C_0 = (Beta / (Lambda_gen * Lambda)) * n_0
    y0 = np.concatenate(([n_0], C_0))
    t_span = (0.0, 10.0)
    t_eval = np.linspace(t_span[0], t_span[1], 1000)
    solution = solve_ivp(
    fun=Kinetic_ODEs,
    t_span=t_span,
    y0=y0,
    t_eval=t_eval,
    method='Radau',
    max_step=1e-3, 
    atol=1e-10,
    rtol=1e-8) 
    time = solution.t
    power = solution.y[0]
    plt.figure(figsize=(10, 6))
    plt.plot(time, power, 'b-', linewidth=2, label='Relative Reactor Power (n)')
    plt.axvline(x=1.0, color='r', linestyle='--', alpha=0.5, label='Reactivity Step Insertion')
    plt.title('Point Reactor Kinetics Transient (Step Reactivity Insertion)', fontsize=14)
    plt.xlabel('Time (seconds)', fontsize=12)
    plt.ylabel('Normalized Power', fontsize=12)
    plt.grid(True, which='both', linestyle='--', alpha=0.7)
    plt.yscale('log')
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig('step_response.png', dpi=150)
    plt.show()

if __name__ == '__main__':
    main()
