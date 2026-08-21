# Point Reactor Kinetics Solver
Status: Work in Progress 🚧

Welcome! I'm currently building this point reactor kinetics solver from the ground up. The core numerical solver, transient visualization, inline documentation, step and ramp insertions with two verifications for step insertion (Inhour and Prompt Jump) and one for Ramp insertion (Period Verification) are now in place. Next up: Ramp insertion convergance verification!

## Mathematical Model
The script (`point_kinetics.py`) solves the standard point kinetics equations for one prompt neutron group and six delayed neutron precursor groups:

$$\frac{dn}{dt} = \frac{\rho(t) - \beta_{total}}{\Lambda} n(t) + \sum_{i=1}^{6} \lambda_i C_i(t)$$

$$\frac{dC_i}{dt} = \frac{\beta_i}{\Lambda} n(t) - \lambda_i C_i(t)$$

**Note:** In the script, the variable `lambda_decay` refers to the array of precursor decay constants ($\lambda_i$ above), and `gen_time` refers to the prompt neutron generation time ($\Lambda$ above). We use `lambda_decay` instead of the single-letter textbook notation because `lambda` is a reserved keyword in Python. Keep that in mind when comparing the equations to the script.

## Current Implementation
- **Core Solver:** Solves the kinetics equations in Python, using `numpy` and `scipy`.
- **Stiff ODE Integration:** Uses `scipy.integrate.solve_ivp` with the `Radau` method to handle the stiffness (the prompt neutron generation time is ~10⁻⁴ s, while the delayed neutron precursors evolve over seconds).
- **Reactivity Insertion:** Models both step insertion ($\rho = 0.002$ at $t = 1$ s) and ramp insertion (linearly from $\rho = 0$ to $\rho = 0.002$ between $t = 1$ s and $t = 3$ s). 200 pcm is roughly 30% of $\beta$ for U-235, simulating a controllable transient safely below the prompt critical threshold.
- **Steady-State Initialization:** Automatically sets the initial precursor concentrations so the system starts from a critical steady state ($n_0 = 1$).
- **Plotting & Visualization:** Uses `matplotlib` to graph the normalized reactor power versus time on a logarithmic scale, for both the step case and the ramp case.
- **Analytical Verification (Inhour):** `Inhour_verification.py` solves the six-group model out to 60 s, finds the numerical reactor period, and calculates the percent error vs the root of the inhour equation (found via `scipy.optimize.brentq`). Currently agrees to within 0.21%.
- **Analytical Verification (Prompt Jump):** `Prompt_jump_verification.py` runs a short simulation to capture the almost instant power spike right after the step insertion and compares it against the analytical Prompt Jump Approximation. Currently agrees within 0.66%.
- **Ramp Verification (Period):** `Ramp_period_verification.py` Checks the accuracy of a ramp reactivity simulation by looking at its final 
growth rate. Once the ramp ends and reactivity stays flat at rho_final, the system should settle into the same stable period as step insertion. Currently the period agrees with the insertion within .237%. 
- **Detailed Documentation:** The script features conversational, beginner-friendly inline comments that explain the physics and mathematical reasoning behind the ODEs.

## Verification Results
### Step
After a step insertion, power doesn't ramp up smoothly
it jumps almost instantly (precursors can't respond on a 10⁻⁴ s timescale), then climbs
as an asymptotic exponential afterward on a period set by the precursor decay. You'll see this clearly on
the log-scale plot below.

#### Inhour Verification
$$\rho = \Lambda \omega + \sum_{i=1}^{6} \frac{\beta_i \omega}{\omega + \lambda_i}$$
 
$1/\omega$ is the asymptotic reactor period that the reactivity settles into once the fast transients have died out. The script solves it by rearranging the equation to equal zero and finding the root with `scipy.optimize.brentq`, then it compares the root to an exponential fit of the tail of the simulated power curve.

Running `Inhour_verification.py` with the default 200 pcm step gives:

| Quantity | Value |
|---|---|
| Reactivity step | 200.0 pcm |
| Analytical period (inhour root) | 17.404 s |
| Numerical period (fit, t ≥ 40 s) | 17.367 s |
| Percent difference | 0.213 % |

The six-group solver's asymptotic period agrees with the inhour equation to within 0.21%, confirming the ODE solver is behaving correctly in that timescale.

That remaining 0.21% isn't necessarily solver error. The inhour equation describes the *pure* asymptotic mode, but the actual transient is six exponential terms (one per precursor group), and the fastest still hasn't completely vanished by t = 40 s. Pushing the fit window later would lower the error even further.

#### Prompt Jump Approximation 
$$\frac{n(0^+)}{n_0} = \frac{\beta_{total}}{\beta_{total} - \rho}$$
 
This gives the near-instant power ratio right after a step insertion. The script computes this directly, then compares it against the simulated power ratio shortly (0.1 s) after the step is inserted.

Running `Prompt_jump_verification.py` with the default 200 pcm step gives:

| Quantity | Value |
|---|---|
| Reactivity step | 200.0 pcm |
| Analytical n(0+)/n0 | 1.4442 |
| Numerical n(t=1.1s)/n0 | 1.4538 |
| Percent difference | 0.661 % |

The solver's immediate post-insertion power ratio agrees with the Prompt Jump Approximation within .66%, confirming the fast-timescale behavior of the solver is correct in addition to its long-term asymptotic behavior.

### Ramp
For a ramp insertion, reactivity rises linearly instead of jumping all at once:
 
$$\rho(t) = \rho_{final} \cdot \frac{t - t_{start}}{t_{end} - t_{start}}, \quad t_{start} \le t \le t_{end}$$
 
Power accelerates smoothly through the ramp, then it settles into the same asymptotic exponential as step insertion once reactivity is constant. Since the reactivity itself is changing throughout the ramp window, there's no sharp jump the way there is for a step, which you'll see in the log-scale plot for ramp reactivity below. 

#### Inhour Verification
There's no simple closed-form solution for a ramp insertion in the six-group model, so once the ramp ends and reactivity remains constant, the system should settle into the same stable period as the step insertion above. 
 
Running `Ramp_period_verification.py` with the default ramp (0 to 200 pcm between t=1s and t=3s) gives:
 
| Quantity | Value |
|---|---|
| Analytical period (inhour root at rho_final) | 17.404 s |
| Numerical period (fit, t ≥ 40 s) | 17.362 s |
| Percent difference | 0.237 % |
 
The ramp case agrees with the inhour equation to within 0.24%, consistent with the step-case result above. This confirms the solver handles the ramp insertion correctly on the long timescale.
 

## What's next
- **Ramp Convergence Verification:** Will check the simulation's numerical convergence by running two simulations, one with normal step size and another with half, then comparing. If cutting the step size in half doesn't change the answer, that's good evidence the solver's converged on the correct answer.
- 
## Usage
To run the current solver, ensure you have `numpy`, `scipy`, and `matplotlib` installed and then run:
```bash
python point_kinetics.py
```
To run the inhour equation verification:
```bash
python Inhour_verification.py
```
To run the Prompt Jump Approximation verification:
```bash
python Prompt_Jump_Approximation.py
```
To run the Ramp period verification:
```bash
python Ramp_period_verification.py
```
## Expected Output
Running the solver with the default parameters will generate two transient response plots and automatically save them to your directory as `step_response.png` and `ramp_response.png`.
 
![Step Reactivity Insertion Response](step_response.png)
![Ramp Reactivity Insertion Response](ramp_response.png)

## References
- J. J. Duderstadt, L. J. Hamilton, *Nuclear Reactor Analysis*, John Wiley & Sons, 1976.
  - (Tabel 2-3) - six-group delayed neutron data used in `point_kinetics.py`.
  - Chapter 6 - Six-group point kinetics equations in the normalized reactivity/generation-time form used in `point_kinetics.py`, corresponding inhour equation in this same normalized form used in `Inhour_verification.py`, and Prompt Jump Approximation used in `Prompt_Jump_Approximation.py`.
