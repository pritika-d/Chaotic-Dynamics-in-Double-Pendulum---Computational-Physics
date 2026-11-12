# Chaotic Dynamics of the Double Pendulum: A Computational Study

An interactive Python simulation of a planar double pendulum. It integrates the full nonlinear equations of motion with the classical fourth-order Runge-Kutta (RK4) method and uses the results to study how the motion changes from near-periodic at small angles to sensitive and irregular at large angles.

## Features

- **Interactive simulator** (Tkinter + Matplotlib): enter `m1`, `m2`, `l1`, `l2`, `θ1`, `θ2`, `ω1` and `ω2`, with Pause/Play and Restart controls. The animation shows a trail for the second bob.
- **Analysis plots**: angles vs. time, phase-space plots ($\omega_1$ vs. $\theta_1$, $\omega_2$ vs. $\theta_2$), relative energy error, and sensitivity to initial conditions. Each runs a fresh 30 s simulation with step $h = 5\times10^{-3}$ s.
- **Chaos diagnostics**: state-space separation $D(t)$ between two nearby trajectories, and an estimate of the largest Lyapunov exponent by the perturb-and-renormalise method.
- **Input validation**: numbers only; masses and lengths must be positive. Angles are entered in degrees and converted to radians internally.

## Physics and Method

The system is two point masses on massless rigid rods, moving in a vertical plane under gravity ($g = 9.81\ \text{m/s}^2$), with a frictionless pivot and no air resistance.

The equations of motion come from the Euler-Lagrange equations. With $\Delta = \theta_1 - \theta_2$:

$$\ddot\theta_1=\frac{-g(2m_1+m_2)\sin\theta_1-m_2 g\sin(\theta_1-2\theta_2)-2\sin\Delta\,m_2\left(\omega_2^2 l_2+\omega_1^2 l_1\cos\Delta\right)}{l_1\left(2m_1+m_2-m_2\cos 2\Delta\right)}$$

$$\ddot\theta_2=\frac{2\sin\Delta\left[\omega_1^2 l_1(m_1+m_2)+g(m_1+m_2)\cos\theta_1+\omega_2^2 l_2 m_2\cos\Delta\right]}{l_2\left(2m_1+m_2-m_2\cos 2\Delta\right)}$$

The denominator is always positive ($\geq 2m_1$), so the equations are well defined for all states. They are integrated as four first-order equations in $y = (\theta_1, \omega_1, \theta_2, \omega_2)$ using RK4. The full derivation is in the project report.

Sensitivity is measured with a dimensionless state-space distance, with angle differences wrapped to $(-\pi, \pi]$ and velocities scaled by $\omega^{\ast} = \sqrt{g/l_1}$:

$$D(t)=\sqrt{(\Delta\theta_1)^2+(\Delta\theta_2)^2+\left(\frac{\Delta\omega_1}{\omega^{\ast}}\right)^2+\left(\frac{\Delta\omega_2}{\omega^{\ast}}\right)^2}$$

## Test Cases

All cases use $m_1 = m_2 = 1$ kg, $l_1 = l_2 = 1$ m and $\omega_1(0) = \omega_2(0) = 0$.

| Case | θ1(0) | θ2(0) | λmax (s⁻¹) | Max relative energy error | Observation |
|------|-------|-------|------------|---------------------------|-------------|
| A: small angles | 10° | 15° | 1.9×10⁻² | 6.7×10⁻¹² | Near-periodic; a small separation does not grow |
| B: intermediate | 40° | 60° | 4.0×10⁻² | 9.8×10⁻¹⁰ | Weak growth; a 30 s run cannot decide |
| C: large angles | 70° | 100° | 4.9×10⁻¹ | 1.2×10⁻⁶ | Clear exponential separation |

## Results

**Angles vs. time.** Case A oscillates smoothly and almost periodically, while Case C is irregular and covers a much larger range of angles.

| Case A (10°, 15°) | Case C (70°, 100°) |
|:---:|:---:|
| ![Case A](images/02_angles_vs_time_case_A.png) | ![Case C](images/01_angles_vs_time_case_C.png) |

**Phase space.** Case A gives thin, closed bands, while Case C fills a much larger region without closing.

| Case A: ω₁ vs. θ₁ | Case C: ω₁ vs. θ₁ |
|:---:|:---:|
| ![Phase space A](images/05_phase_space_w1_theta1_case_A.png) | ![Phase space C](images/03_phase_space_w1_theta1_case_C.png) |

**Sensitivity to initial conditions.** Changing $\theta_1(0)$ by only $10^{-3}$° ($D(0) \approx 1.745\times10^{-5}$) has no visible effect on Case A. In Case C the separation stays tiny for about 10 s, then grows to order one by about 20 s.

| Case A | Case C |
|:---:|:---:|
| ![Sensitivity A](images/07_sensitivity_case_A.png) | ![Sensitivity C](images/08_sensitivity_case_C.png) |

**Largest Lyapunov exponent.** For Case C, $\lambda_{\max} \approx 0.49\ \text{s}^{-1}$, so a small error grows by a factor of $e$ about every 2 s. For Cases A and B the values are too small to distinguish from finite-time effects in a 30 s run.

**Accuracy.** The relative energy error decreases as the time step is reduced, and stays below about $10^{-6}$ even for Case C.

## Limitations

- RK4 does not conserve energy exactly, and the step size is fixed.
- The Lyapunov exponent is a finite-time estimate from a single trajectory.
- The animation uses a larger step than the analysis runs, so it is for visualisation only.
- The model is idealised: no friction or air resistance.

## Getting Started

**Requirements:** Python 3, NumPy, Matplotlib and Tkinter (included with most Python installations).

```bash
pip install numpy matplotlib
python double_pendulum_project.py
```

## Author

*Your name here*
