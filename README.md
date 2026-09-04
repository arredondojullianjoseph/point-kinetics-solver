# Point Reactor Kinetics Solver

Six-group point reactor kinetics solver in Python. Solves the coupled prompt-neutron and delayed-precursor ODEs with a stiff Radau integrator. Supports step and ramp reactivity insertions. Verified against the inhour equation and prompt jump approximation for the step insertion, and against period and convergence checks for the ramp insertion. With default insertions and steady-state conditions, all verifications agree to within 1%. Automated tests cover the same checks.

Scope: Intended as a self-study verification of the standard six-group equations.

## Mathematical model

The script (`point_kinetics.py`) solves the standard point kinetics equations for one prompt neutron group and six delayed neutron precursor groups:

dn/dt = [(ρ(t) − β_total) / Λ] n(t) + Σᵢ λᵢ Cᵢ(t)

dCᵢ/dt = (βᵢ / Λ) n(t) − λᵢ Cᵢ(t)

Note: in the script, `lambda_decay` refers to the array of precursor decay constants (λᵢ above), and `gen_time` refers to the prompt neutron generation time (Λ above). This is the reverse of textbook notation, since `lambda` is a reserved keyword in Python. Keep that in mind when comparing the equations to the script.

## Implementation

- **Core solver:** solves the kinetics equations using numpy and scipy.
- **Stiff ODE integration:** uses `scipy.integrate.solve_ivp` with the Radau method (the prompt neutron generation time is ~10⁻⁴ s, while the delayed precursors evolve over seconds).
- **Reactivity insertion:** models step insertion (ρ = 0.002 at t = 1 s) and ramp insertion (ρ = 0 to 0.002 linearly between t = 1 s and t = 3 s). 200 pcm is roughly 30% of β for U-235. This keeps the transient controllable and below the prompt critical threshold.
- **Steady-state initialization:** sets initial precursor concentrations so the system starts from a critical steady state (n₀ = 1).
- **Plotting:** uses matplotlib to graph normalized reactor power versus time on a log scale, for both the step and ramp cases.
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

ρ = Λω + Σᵢ [βᵢω / (ω + λᵢ)]

1/ω is the asymptotic reactor period the reactivity settles into once the fast transients have died out. The script finds the root with `scipy.optimize.brentq`, then compares it to an exponential fit of the simulated power curve's tail.

Running `inhour_verification.py` with the default 200 pcm step:

| Quantity | Value |
|---|---|
| Reactivity step | 200.0 pcm |
| Analytical period (inhour root) | 17.404 s |
| Numerical period (fit, t ≥ 40 s) | 17.367 s |
| Percent difference | 0.213% |

The asymptotic period agrees with the inhour equation to within 0.21%. That remaining error isn't necessarily solver error: the inhour equation describes the pure asymptotic mode, but the actual transient is six exponential terms, and the fastest hasn't fully vanished by t = 40 s. A later fit window would lower the error further.

### Prompt jump approximation

n(0+)/n₀ = β_total / (β_total − ρ)

This gives the near-instant power ratio right after a step insertion. The script computes it directly and compares it against the simulated power ratio 0.1 s after the step.

Running `prompt_jump_verification.py` with the default 200 pcm step:

| Quantity | Value |
|---|---|
| Reactivity step | 200.0 pcm |
| Analytical n(0+)/n₀ | 1.4442 |
| Numerical n(t=1.1s)/n₀ | 1.4538 |
| Percent difference | 0.661% |

The solver's immediate post-insertion power ratio agrees with the prompt jump approximation to within 0.66%, confirming correct fast-timescale behavior.

### Ramp

Reactivity rises linearly instead of jumping:

ρ(t) = ρ_final · (t − t_start) / (t_end − t_start), for t_start ≤ t ≤ t_end

Power accelerates smoothly through the ramp, then settles into the same asymptotic exponential as step insertion once reactivity is constant. There's no sharp jump, since reactivity itself is changing through the ramp window. See the log-scale ramp plot below.

### Period verification

There's no simple closed-form solution for a ramp insertion in the six-group model. Once the ramp ends and reactivity holds constant, the system should settle into the same stable period as step insertion.

Running `ramp_period_verification.py` with the default ramp (0 to 200 pcm, t = 1 s to 3 s):

| Quantity | Value |
|---|---|
| Analytical period (inhour root at rho_final) | 17.404 s |
| Numerical period (fit, t ≥ 40 s) | 17.362 s |
| Percent difference | 0.237% |

Agrees with the inhour equation to within 0.24%, consistent with the step case.

### Convergence verification

Runs the simulation at normal step size and half step size, then compares. If halving the step size doesn't change the answer, that's evidence the solver has converged.

Running `ramp_convergence_verification.py` with the default ramp:

| Quantity | Value |
|---|---|
| Max relative difference (max_step 1e-3 vs 5e-4) | 8.58e-09 |

A difference this small confirms the ramp solution is converged.

## Limitations

- Point kinetics only; local power tilts and rod-position effects aren't represented.
- Six-group parameters are hardcoded and not configurable for fuels other than U-235.
- No reactivity feedback.
- Limited to step and ramp insertion types.

## Usage

To run the solver:

    python point_kinetics.py

To run the inhour equation verification:

    python inhour_verification.py

To run the prompt jump approximation verification:

    python prompt_jump_verification.py

To run the ramp period verification:

    python ramp_period_verification.py

To run the ramp convergence verification:

    python ramp_convergence_verification.py

## Running tests

    pytest test_point_kinetics.py -v

## Expected output

Running the solver with default parameters generates two transient response plots, saved as `step_response.png` and `ramp_response.png`.

## References

Duderstadt JJ, Hamilton LJ. 1976. Nuclear Reactor Analysis. New York: Wiley.
- Table 2-3: six-group delayed neutron data used in `point_kinetics.py`.
- Chapter 6: six-group point kinetics equations in the normalized reactivity/generation-time form used in `point_kinetics.py`, the corresponding inhour equation used in `inhour_verification.py`, and the prompt jump approximation used in `prompt_jump_verification.py`.
