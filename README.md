# Point Reactor Kinetics Solver

Status: Work in Progress 🚧

Welcome! I'm currently building this point reactor kinetics solver from the
ground up. The core numerical solver has been added, and I am now working
on expanding its visualization, verification, and documentation.

## Mathematical Model

The script (`point_kinetics.py`) solves the standard point kinetics
equations for one prompt neutron group and six delayed neutron precursor
groups:

$$\frac{dn}{dt} = \frac{\rho(t) - \beta_{total}}{\Lambda} n(t) + \sum_{i=1}^{6} \lambda_i C_i(t)$$

$$\frac{dC_i}{dt} = \frac{\beta_i}{\Lambda} n(t) - \lambda_i C_i(t)$$

**Note:** in the code, the variable `Lambda` refers to the array of
precursor decay constants (λ_i above), and `Lambda_gen` refers to the
prompt neutron generation time (Λ above). This is the opposite of how the
symbols are typically laid out in textbook notation. The swap was done in order to avoid a python syntax error since 'lambda' is a reserved keyword. Keep that in mind
when comparing the equations to the code.

## Current Implementation

- **Core Solver:** Solves the kinetics equations in Python, using `numpy`
  and `scipy`.
- **Stiff ODE Integration:** Uses `scipy.integrate.solve_ivp` with the
  `Radau` method to handle the stiffness (prompt neutron lifetime is
  ≈10⁻⁴ s and the delayed neutron precursors are much slower).
- **Reactivity Insertion:** Currently models a step reactivity insertion
  (ρ = 0.002 at t = 1 s).
- **Steady-State Initialization:** Automatically sets the initial
  precursor concentrations so the system starts from a critical steady
  state (n₀ = 1).

## What's Next?

- **Plotting & Visualization:** Adding graphing tools using `matplotlib`
  to visualize the transient reactor power and delayed neutron precursor
  concentrations.
- **Analytical Verification:** Benchmarking against established models
  (Inhour Equation & Prompt Jump Approximation).
- **Documentation:** Adding inline comments, expanding setup instructions,
  usage examples, and a deeper breakdown of the underlying physics.

## Usage

To run the current solver, ensure you have `numpy` and `scipy` installed
and then run:
```bash
python point_kinetics.py
```

## Expected Output

Running the solver with the default parameters will output the final normalized neutron population at $t = 10$ seconds:
```text
Final relative power: 3.0305
```



