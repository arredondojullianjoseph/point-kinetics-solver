"""
point_kinetics.py
Solves the six-group point reactor kinetics equations for a step
reactivity insertion using scipy's stiff ODE integrator (Radau).
"""

import numpy as np
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
    solution = solve_ivp(fun=Kinetic_ODEs, t_span=t_span, y0=y0, t_eval=t_eval, method='Radau', max_step=1e-3, atol=1e-10, rtol=1e-8)
    power = solution.y[0]
    print("Final relative power: {:.4f}".format(power[-1]))

if __name__ == '__main__':
    main()
