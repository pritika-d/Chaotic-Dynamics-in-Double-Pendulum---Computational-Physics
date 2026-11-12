# Chaotic Dynamics of the Double Pendulum: A Computational Study

An interactive Python simulation of a planar double pendulum. It integrates the full nonlinear equations of motion with the classical fourth-order Runge-Kutta (RK4) method and uses the results to study how the motion changes from near-periodic at small angles to sensitive and irregular at large angles.

## Features

- **Interactive simulator** (Tkinter + Matplotlib): enter `m1`, `m2`, `l1`, `l2`, `θ1`, `θ2`, `ω1` and `ω2`, with Pause/Play and Restart controls. The animation shows a trail for the second bob.
- **Analysis plots**: angles vs. time, phase-space plots (`ω1` vs. `θ1`, `ω2` vs. `θ2`), relative energy error, and sensitivity to initial conditions. Each runs a fresh 30 s simulation with step `h = 5×10⁻³ s`.
- **Chaos diagnostics**: state-space separation `D(t)` between two nearby trajectories, and an estimate of the largest Lyapunov exponent by the perturb-and-renormalise method.
- **Input validation**: numbers only; masses and lengths must be positive. Angles are entered in degrees and converted to radians internally.

## Physics and Method

The system is two point masses on massless rigid rods, moving in a vertical plane under gravity (`g = 9.81 m/s²`), with a frictionless pivot and no air resistance.

The equations of motion come from the Euler-Lagrange equations. With `Δ = θ1 − θ2`:

```
θ1'' = [ -g(2m1 + m2) sin θ1 - m2 g sin(θ1 - 2θ2)
         - 2 sin Δ · m2 (ω2² l2 + ω1² l1 cos Δ) ] / [ l1 (2m1 + m2 - m2 cos 2Δ) ]

θ2'' = [ 2 sin Δ ( ω1² l1 (m1 + m2) + g (m1 + m2) cos θ1 + ω2² l2 m2 cos Δ ) ]
         / [ l2 (2m1 + m2 - m2 cos 2Δ) ]
```

The denominator is always positive (`≥ 2m1`), so the equations are well defined for all states. They are integrated as four first-order equations in `y = (θ1, ω1, θ2, ω2)` using RK4. The full derivation is in the project report.

Sensitivity is measured with a dimensionless state-space distance, with angle differences wrapped to (−π, π] and velocities scaled by `ω* = √(g/l1)`:

```
D(t) = sqrt( Δθ1² + Δθ2² + (Δω1/ω*)² + (Δω2/ω*)² )
```

## Test Cases

All cases use `m1 = m2 = 1 kg`, `l1 = l2 = 1 m` and `ω1(0) = ω2(0) = 0`.

| Case | θ1(0) | θ2(0) | λmax (s⁻¹) | Max relative energy error | Observation |
|------|-------|-------|------------|---------------------------|-------------|
| A: small angles | 10° | 15° | 1.9×10⁻² | 6.7×10⁻¹² | Near-periodic; a small separation does not grow |
| B: intermediate | 40° | 60° | 4.0×10⁻² | 9.8×10⁻¹⁰ | Weak growth; a 30 s run cannot decide |
| C: large angles | 70° | 100° | 4.9×10⁻¹ | 1.2×10⁻⁶ | Clear exponential separation |

## Key Results

- For Case C, a perturbation of `10⁻³°` (≈ 1.745×10⁻⁵ rad) in `θ1(0)` leaves the motion nearly unchanged for about 10 s, then grows to order one by about 20 s. After that `D(t)` only fluctuates, because the distance is bounded.
- The largest Lyapunov exponent for Case C is about **0.49 s⁻¹**, so a small error grows by a factor of *e* roughly every 2 s. It is positive in both halves of the run.
- In Case A, `D(t)` stays close to its starting value and never exceeds about 3.4×10⁻⁵.
- RK4 energy error decreases as the time step is reduced (tested for `h` from 10⁻² s down to 1.25×10⁻³ s), but it is largest for the most energetic motion.

## Limitations

- The Lyapunov exponent is a finite-time estimate from a single trajectory, so small positive values need cautious interpretation.
- RK4 does not conserve energy exactly, and drift may build up over longer runs.
- The time step is fixed, with no adaptive solver.
- The animation uses a larger step than the analysis runs, so it is for visualisation, not precision measurement.
- The model is idealised: no damping or friction.

## Getting Started

### Requirements

- Python 3
- NumPy
- Matplotlib
- Tkinter (included with most Python installations)

```bash
pip install numpy matplotlib
```

### Run

```bash
python double_pendulum_project.py
```

Enter your parameters in the control panel and use the buttons to start the animation or generate the analysis plots.


## Author

*Pritika D*
