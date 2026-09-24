# Point Reactor Kinetics Solver

Six-group point reactor kinetics solver in Python. Solves the coupled prompt-neutron, delayed-precursor, and lumped-fuel-temperature ODEs with a stiff Radau integrator. Supports step and ramp reactivity insertions. Fuel temperature is integrated with Newton cooling and, when the linear Doppler coefficient is enabled, feeds back into reactivity so power peaks and settles instead of growing forever. Verified against the inhour equation and prompt jump approximation for the step insertion, and against period and convergence checks for the ramp insertion. With default insertions and steady-state conditions, all verifications agree to within 1%. Automated tests cover the same checks.

Scope: Intended as a self-study verification of the standard six-group equations, extended step by step toward Doppler reactivity feedback.

## Mathematical model

The script (`point_kinetics.py`) solves the standard point kinetics equations for one prompt neutron group and six delayed neutron precursor groups, plus a lumped fuel energy balance.

$$\frac{dn}{dt} = \frac{\rho(t) - \beta}{\Lambda} n(t) + \sum_{i} \lambda_i C_i(t)$$

$$\frac{dC_i}{dt} = \frac{\beta_i}{\Lambda} n(t) - \lambda_i C_i(t)$$

$$\frac{dT}{dt} = \frac{(T_0 - T_c)\, n(t) - (T - T_c)}{\tau}$$

$\rho(t)$ is the external insertion (step or ramp). $T$ feeds back into reactivity through an optional linear Doppler term, added inside `kinetics_odes`:

$$\rho_{\text{total}} = \rho(t) + \alpha_D (T - T_0)$$

`alpha_D` defaults to `0.0` everywhere it isn't explicitly passed, so $\rho_{\text{total}} = \rho(t)$ and every existing script/test is unaffected. The heat-balance form is constructed so $(n, T) = (1, T_0)$ is an equilibrium: heat generation at $n = 1$ exactly matches Newton cooling from $T_0$ down to coolant temperature $T_c$.

Hardcoded thermal constants (illustrative UO2-ish values, not a specific core design): $T_0 = 900$ K (`T0`), $T_c = 580$ K (`T_coolant`), $\tau = 5$ s (`tau_fuel`), $\alpha_D = -2 \times 10^{-5}$ /K, i.e. about -2 pcm/K (`ALPHA_D`). $\alpha_D$ is negative because hotter fuel adds negative reactivity (Doppler broadening of resonance absorption). At this value, a 200 pcm step is compensated once fuel heats by $\Delta T = -\rho / \alpha_D \approx 100$ K.

Note: in the script, `lambda_decay` refers to the array of precursor decay constants ($\lambda_i$ above), and `gen_time` refers to the prompt neutron generation time ($\Lambda$ above). This is the reverse of textbook notation, since `lambda` is a reserved keyword in Python. Keep that in mind when comparing the equations to the script.

## Implementation

- **Core solver:** solves the kinetics and fuel-temperature equations using numpy and scipy.
- **Stiff ODE integration:** uses `scipy.integrate.solve_ivp` with the Radau method (the prompt neutron generation time is ~10⁻⁴ s, while the delayed precursors evolve over seconds).
- **Reactivity insertion:** models step insertion ($\rho = 0.002$ at $t = 1$ s) and ramp insertion ($\rho$ from $0$ to $0.002$ linearly between $t = 1$ s and $t = 3$ s). 200 pcm is roughly 30% of $\beta$ for U-235. This keeps the transient controllable and below the prompt critical threshold. With Doppler off (`alpha_D = 0`, the default), these insertions produce unbounded exponential growth after the prompt jump, same as a constant-$\rho$ model.
- **Lumped fuel temperature:** always integrated as `y[7]` with Newton cooling; precursors occupy `y[1:7]`.
- **Linear Doppler feedback (opt-in):** `kinetics_odes(t, y, reactivity_fn, alpha_D)` adds $\alpha_D (T - T_0)$ to $\rho(t)$ when `alpha_D` is passed as nonzero. Every existing caller omits it, so `alpha_D` defaults to `0.0` and behaves exactly as before Doppler was added.
- **Steady-state initialization:** `steady_state_y0()` sets precursor concentrations so the system starts from a critical steady state ($n_0 = 1$) and fuel temperature $T_0$ so $dT/dt = 0$ at that power.
- **Plotting:** `main()`/`main_ramp()` graph normalized reactor power versus time on a log scale (Doppler off), for both the step and ramp cases. `main_step_doppler()`/`main_ramp_doppler()` run the same insertions with Doppler on and plot power (left axis, log scale) together with fuel temperature (right axis) on one figure.
- **Step verification (inhour):** `inhour_verification.py` solves the model out to 60 s, finds the numerical reactor period, and checks percent error against the root of the inhour equation (via `scipy.optimize.brentq`). Agrees to within 0.21%.
- **Step verification (prompt jump):** `prompt_jump_verification.py` captures the near-instant power spike after the step insertion and compares it against the analytical prompt jump approximation. Agrees to within 0.66%.
- **Ramp verification (period):** `ramp_period_verification.py` checks the ramp simulation's final growth rate against the step-insertion period, once reactivity holds flat at `rho_final`. Agrees to within 0.24%.
- **Ramp verification (convergence):** `ramp_convergence_verification.py` runs the simulation at normal step size and at half step size and compares. Agrees to within 8.58e-09.
- **Inline documentation:** comments explain the physics and mathematical reasoning behind the ODEs.
- **Automated tests:** `test_point_kinetics.py` turns the verification checks into pass/fail assertions using pytest.

## Verification results

### Step

After a step insertion, power jumps almost instantly (precursors can't respond on a 10⁻⁴ s timescale), then climbs as an asymptotic exponential set by the precursor decay. The log-scale plot below shows this.

### Inhour verification

$$\rho = \Lambda \omega + \sum_{i} \frac{\beta_i \omega}{\omega + \lambda_i}$$

$1/\omega$ is the asymptotic reactor period the reactivity settles into once the fast transients have died out. The script finds the root with `scipy.optimize.brentq`, then compares it to an exponential fit of the simulated power curve's tail.

Running `inhour_verification.py` with the default 200 pcm step:

| Quantity | Value |
| --- | --- |
| Reactivity step | 200.0 pcm |
| Analytical period (inhour root) | 17.404 s |
| Numerical period (fit, t ≥ 40 s) | 17.367 s |
| Percent difference | 0.213% |

The asymptotic period agrees with the inhour equation to within 0.21%. That remaining error isn't necessarily solver error: the inhour equation describes the pure asymptotic mode, but the actual transient is six exponential terms, and the fastest hasn't fully vanished by t = 40 s. A later fit window would lower the error further.

### Prompt jump approximation

$$\frac{n(0^+)}{n_0} = \frac{\beta}{\beta - \rho}$$

This gives the near-instant power ratio right after a step insertion. The script computes it directly and compares it against the simulated power ratio 0.1 s after the step.

Running `prompt_jump_verification.py` with the default 200 pcm step:

| Quantity | Value |
| --- | --- |
| Reactivity step | 200.0 pcm |
| Analytical $n(0^+)/n_0$ | 1.4442 |
| Numerical $n(t=1.1\text{s})/n_0$ | 1.4538 |
| Percent difference | 0.661% |

The solver's immediate post-insertion power ratio agrees with the prompt jump approximation to within 0.66%, confirming correct fast-timescale behavior.

### Ramp

Reactivity rises linearly instead of jumping:

$$\rho(t) = \rho_{\text{final}} \cdot \frac{t - t_{\text{start}}}{t_{\text{end}} - t_{\text{start}}}, \quad t_{\text{start}} \le t \le t_{\text{end}}$$

Power accelerates smoothly through the ramp, then settles into the same asymptotic exponential as step insertion once reactivity is constant. There's no sharp jump, since reactivity itself is changing through the ramp window. See the log-scale ramp plot below.

### Period verification

There's no simple closed-form solution for a ramp insertion in the six-group model. Once the ramp ends and reactivity holds constant, the system should settle into the same stable period as step insertion.

Running `ramp_period_verification.py` with the default ramp (0 to 200 pcm, t = 1 s to 3 s):

| Quantity | Value |
| --- | --- |
| Analytical period (inhour root at rho_final) | 17.404 s |
| Numerical period (fit, t ≥ 40 s) | 17.362 s |
| Percent difference | 0.237% |

Agrees with the inhour equation to within 0.24%, consistent with the step case.

### Convergence verification

Runs the simulation at normal step size and half step size, then compares. If halving the step size doesn't change the answer, that's evidence the solver has converged.

Running `ramp_convergence_verification.py` with the default ramp:

| Quantity | Value |
| --- | --- |
| Max relative difference (max_step 1e-3 vs 5e-4) | 8.58e-09 |

A difference this small confirms the ramp solution is converged.

## Limitations

- Point kinetics only; local power tilts and rod-position effects aren't represented.
- Six-group parameters are hardcoded and not configurable for fuels other than U-235.
- Fuel is a single lumped node with a fixed coolant temperature; there is no clad or coolant dynamics.
- Doppler feedback is off by default (`alpha_D = 0`); the default step/ramp plots and every verification script still show a positive insertion growing without bound. Doppler is only active in `main_step_doppler()`/`main_ramp_doppler()`, which pass `alpha_D` explicitly.
- Only a linear Doppler law ($\rho_D = \alpha_D (T - T_0)$) is implemented so far; no sqrt(T) or other functional form yet, and no analytical verification of the new equilibrium yet (planned for a later step).
- Thermal constants ($T_0$, $T_c$, $\tau$, $\alpha_D$) are hardcoded illustrative values, not a specific core design.
- Limited to step and ramp insertion types.

## Usage

```bash
python point_kinetics.py
python inhour_verification.py
python prompt_jump_verification.py
python ramp_period_verification.py
python ramp_convergence_verification.py
```

## Running tests

    pytest test_point_kinetics.py -v

## Expected output

Running `python point_kinetics.py` generates four transient response plots:

- `step_response.png`, `ramp_response.png`: normalized power only, Doppler off (unbounded growth after the prompt jump).
- `step_doppler_response.png`, `ramp_doppler_response.png`: normalized power (left axis, log scale) and fuel temperature (right axis), Doppler on. Power peaks and levels off around $n \approx 1.3$ as $T$ settles near $T_0 + 100\,K$.

## References

Duderstadt JJ, Hamilton LJ. 1976. Nuclear Reactor Analysis. New York: Wiley.
- Table 2-3: six-group delayed neutron data used in `point_kinetics.py`.
- Chapter 6: six-group point kinetics equations in the normalized reactivity/generation-time form used in `point_kinetics.py`, the corresponding inhour equation used in `inhour_verification.py`, and the prompt jump approximation used in `prompt_jump_verification.py`.
