# Point Reactor Kinetics Solver

Status: Work in Progress 🚧

Welcome! I'm currently building this point reactor kinetics solver from the ground up. The core numerical solver, transient visualization, and inline documentation are now in place. Next up: expanding analytical verification!

## Mathematical Model

The script (`point_kinetics.py`) solves the standard point kinetics equations for one prompt neutron group and six delayed neutron precursor groups:

$$\frac{dn}{dt} = \frac{\rho(t) - \beta_{total}}{\Lambda} n(t) + \sum_{i=1}^{6} \lambda_i C_i(t)$$

$$\frac{dC_i}{dt} = \frac{\beta_i}{\Lambda} n(t) - \lambda_i C_i(t)$$

**Note:** In the Python code, the variable `lambda_decay` refers to the array of precursor decay constants ($\lambda_i$ above), and `gen_time` refers to the prompt neutron generation time ($\Lambda$ above). We use `lambda_decay` instead of standard notation because `lambda` is a reserved keyword in Python! Keep that in mind when comparing the textbook equations to the codebase.

## Current Implementation

- **Core Solver:** Solves the kinetics equations in Python, using `numpy` and `scipy`.
- **Stiff ODE Integration:** Uses `scipy.integrate.solve_ivp` with the `Radau` method to handle the stiffness (the prompt neutron generation time is ~10⁻⁴ s, while the delayed neutron precursors evolve over seconds).
- **Reactivity Insertion:** Models a step reactivity insertion ($\rho = 0.002$ at $t = 1$ s). This is roughly 30% of $\beta$ for U-235, simulating a controllable transient safely below the prompt critical threshold.
- **Steady-State Initialization:** Automatically sets the initial precursor concentrations so the system starts from a critical steady state ($n_0 = 1$).
- **Plotting & Visualization:** Uses `matplotlib` to graph the normalized reactor power versus time on a logarithmic scale, visually indicating the step insertion.
- **Detailed Documentation:** The code features conversational, beginner-friendly inline comments that explain the physics and mathematical reasoning behind the ODEs.

## What's Next?

- **Analytical Verification:** Benchmarking against established models (Inhour Equation & Prompt Jump Approximation).
- **Expanded Usage:** Adding setup instructions, varying step sizes, and a deeper breakdown of the underlying physics.

## Usage

To run the current solver, ensure you have `numpy`, `scipy`, and `matplotlib` installed and then run:
```bash
python point_kinetics.py
```

## Expected Output

Running the solver with the default parameters will generate a transient response plot and automatically save it to your directory as `step_response.png`.
![Step Reactivity Insertion Response](step_response.png)
