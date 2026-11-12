import tkinter as tk
from tkinter import ttk, messagebox

import numpy as np

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt

G = 9.81                      
FRAME_INTERVAL_MS = 33 # ~30 fps

PARAMS = [
    ("m1",     "m1 (kg)",        1.0),
    ("m2",     "m2 (kg)",        1.0),
    ("l1",     "l1 (m)",         1.0),
    ("l2",     "l2 (m)",         1.0),
    ("theta1", "θ1 (deg)",       0.0),
    ("theta2", "θ2 (deg)",       0.0),
    ("omega1", "ω1 (rad/s)",     0.0),
    ("omega2", "ω2 (rad/s)",     0.0),
]

#PHYSICS
class DoublePendulum:

    def __init__(self, m1, m2, l1, l2, theta1, theta2, omega1, omega2):
        self.m1, self.m2 = m1, m2
        self.l1, self.l2 = l1, l2
        # State vector: [theta1, omega1, theta2, omega2]
        self.state = np.array([theta1, omega1, theta2, omega2], dtype=float)

    def derivatives(self, t, y):
        #Lagrangian form
        th1, w1, th2, w2 = y
        m1, m2, l1, l2 = self.m1, self.m2, self.l1, self.l2

        delta = th1 - th2
        sin_d, cos_d = np.sin(delta), np.cos(delta)
        den = 2 * m1 + m2 - m2 * np.cos(2 * delta)

        a1 = (
            -G * (2 * m1 + m2) * np.sin(th1)
            - m2 * G * np.sin(th1 - 2 * th2)
            - 2 * sin_d * m2 * (w2 ** 2 * l2 + w1 ** 2 * l1 * cos_d)
        ) / (l1 * den)

        a2 = (
            2 * sin_d * (
                w1 ** 2 * l1 * (m1 + m2)
                + G * (m1 + m2) * np.cos(th1)
                + w2 ** 2 * l2 * m2 * cos_d
            )
        ) / (l2 * den)

        return [w1, a1, w2, a2]

    def step(self, dt):
        #Fourth-order Runge-Kutta (RK4)
        t = 0.0
        y = self.state
        h = dt

        k1 = np.asarray(self.derivatives(t, y), dtype=float)
        k2 = np.asarray(self.derivatives(t + h / 2, y + h * k1 / 2), dtype=float)
        k3 = np.asarray(self.derivatives(t + h / 2, y + h * k2 / 2), dtype=float)
        k4 = np.asarray(self.derivatives(t + h, y + h * k3), dtype=float)

        self.state = y + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

    def positions(self):
        #Cartesian positions of pivot, bob 1 and bob 2
        th1, _, th2, _ = self.state
        x1 = self.l1 * np.sin(th1)
        y1 = -self.l1 * np.cos(th1)
        x2 = x1 + self.l2 * np.sin(th2)
        y2 = y1 - self.l2 * np.cos(th2)
        return (0.0, 0.0), (x1, y1), (x2, y2)


# ----------------------------------------------------------------------------
# UI + animation loop
class PendulumApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Double Pendulum")

        self.running = True          
        self.pendulum = None
        self.trail_x, self.trail_y = [], []
        self.entries = {}             

        self._build_ui()
        self._build_plot()
        self.restart()               
        self._animate()              

    # ---------------------------- UI setup ----------------------------------
    def _build_ui(self):
        """Create the control panel with entry fields and buttons."""
        panel = ttk.Frame(self.root, padding=10)
        panel.pack(side=tk.LEFT, fill=tk.Y)

        ttk.Label(panel, text="Parameters", font=("TkDefaultFont", 12, "bold")
                  ).grid(row=0, column=0, columnspan=2, pady=(0, 8))

        # One labelled entry field per adjustable parameter
        for i, (key, label, default) in enumerate(PARAMS, start=1):
            ttk.Label(panel, text=label).grid(row=i, column=0, sticky="w", pady=2)
            var = tk.StringVar(value=str(default)) 
            ttk.Entry(panel, textvariable=var, width=10).grid(
                row=i, column=1, padx=(8, 0), pady=2)
            self.entries[key] = var

        # Playback buttons
        self.play_btn = ttk.Button(panel, text="Pause", command=self.toggle_play)
        self.play_btn.grid(row=len(PARAMS) + 1, column=0, columnspan=2,
                           sticky="ew", pady=(12, 3))
        ttk.Button(panel, text="Restart", command=self.restart).grid(
            row=len(PARAMS) + 2, column=0, columnspan=2, sticky="ew", pady=3)

        # Analysis buttons: each plot starts a fresh simulation from the values
        button_row = len(PARAMS) + 3
        ttk.Button(panel, text="Plot θ1 and θ2 vs Time",
                   command=self.plot_angles_vs_time).grid(
            row=button_row, column=0, columnspan=2, sticky="ew", pady=(10, 3))
        ttk.Button(panel, text="Phase Space: θ1 vs ω1",
                   command=self.plot_phase_space_1).grid(
            row=button_row + 1, column=0, columnspan=2, sticky="ew", pady=3)
        ttk.Button(panel, text="Phase Space: θ2 vs ω2",
                   command=self.plot_phase_space_2).grid(
            row=button_row + 2, column=0, columnspan=2, sticky="ew", pady=3)
        ttk.Button(panel, text="Plot Relative Energy Error",
                   command=self.plot_relative_energy_error).grid(
            row=button_row + 3, column=0, columnspan=2, sticky="ew", pady=3)
        ttk.Button(panel, text="Sensitivity to Initial Conditions",
                   command=self.plot_sensitivity).grid(
            row=button_row + 4, column=0, columnspan=2, sticky="ew", pady=3)
        
    def _build_plot(self):
        """Create the matplotlib figure embedded in the Tk window."""
        self.fig = Figure(figsize=(6, 6), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_aspect("equal")
        self.ax.grid(True, alpha=0.3)
        self.ax.set_xlabel("x (m)")
        self.ax.set_ylabel("y (m)")

    
        self.rods, = self.ax.plot([], [], "-", color="black", lw=2)
        self.bobs, = self.ax.plot([], [], "o", color="tab:red", ms=10)
        self.trail, = self.ax.plot([], [], "-", color="tab:blue", lw=1, alpha=0.6)
        self.time_text = self.ax.text(0.02, 0.96, "", transform=self.ax.transAxes,
                                      va="top")

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

    # ---------------------------- Controls ----------------------------------
    def _read_inputs(self):
        try:
            values = {k: float(v.get()) for k, v in self.entries.items()}
        except ValueError:
            messagebox.showerror("Invalid input", "All fields must be numbers.")
            return None
        if min(values["m1"], values["m2"], values["l1"], values["l2"]) <= 0:
            messagebox.showerror("Invalid input",
                                 "Masses and lengths must be greater than zero.")
            return None
        # The user enters angles in degrees but the physics equations use radians.
        values["theta1"] = np.deg2rad(values["theta1"])
        values["theta2"] = np.deg2rad(values["theta2"])
        return values

    def toggle_play(self):
        #Pause/Play
        self.running = not self.running
        self.play_btn.config(text="Pause" if self.running else "Play")

    def restart(self):
        #Reset
        values = self._read_inputs()
        if values is None:
            return
        self.pendulum = DoublePendulum(**values)
        self.current_values = values.copy()
        self.sim_time = 0.0
        self.trail_x.clear()
        self.trail_y.clear()

        # Rescale axes to fit the full reach of the pendulum
        reach = (values["l1"] + values["l2"]) * 1.1
        self.ax.set_xlim(-reach, reach)
        self.ax.set_ylim(-reach, reach)
        self._draw()


    # ---------------------------- Analysis plots ----------------------------
    def _simulate_from_inputs(self, duration=30.0, dt=0.005):
        
        values = self._read_inputs()
        if values is None:
            return None

        model = DoublePendulum(**values)
        steps = int(round(duration / dt))
        times = np.linspace(0.0, steps * dt, steps + 1)
        states = np.empty((steps + 1, 4), dtype=float)
        states[0] = model.state.copy()
        for i in range(1, steps + 1):
            model.step(dt)
            states[i] = model.state.copy()
        return times, states, values

    def plot_angles_vs_time(self):
        #Plot both angular displacements against time for current inputs
        result = self._simulate_from_inputs()
        if result is None:
            return
        times, states, values = result
        fig, ax = plt.subplots(figsize=(9, 5.5))
        ax.plot(times, np.rad2deg(states[:, 0]), color="tab:blue", label=r"$\theta_1(t)$")
        ax.plot(times, np.rad2deg(states[:, 2]), color="tab:orange", label=r"$\theta_2(t)$")
        ax.set_title("Angular Displacement vs Time (RK4)")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Angular displacement (degrees)")
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.show()

    def plot_phase_space_1(self):
        #Plot omega1 against theta1 for current initial conditions
        result = self._simulate_from_inputs()
        if result is None:
            return
        times, states, values = result
        fig, ax = plt.subplots(figsize=(7, 5.5))
        ax.plot(np.rad2deg(states[:, 0]), states[:, 1], color="tab:blue", lw=1)
        ax.set_title(r"Phase Space: $\omega_1$ vs $\theta_1$ (RK4)")
        ax.set_xlabel(r"Angular position $\theta_1$ (degrees)")
        ax.set_ylabel(r"Angular velocity $\omega_1$ (rad/s)")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.show()

    def plot_phase_space_2(self):
        #Plot omega2 against theta2 for current initial conditions
        result = self._simulate_from_inputs()
        if result is None:
            return
        times, states, values = result
        fig, ax = plt.subplots(figsize=(7, 5.5))
        ax.plot(np.rad2deg(states[:, 2]), states[:, 3], color="tab:green", lw=1)
        ax.set_title(r"Phase Space: $\omega_2$ vs $\theta_2$ (RK4)")
        ax.set_xlabel(r"Angular position $\theta_2$ (degrees)")
        ax.set_ylabel(r"Angular velocity $\omega_2$ (rad/s)")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.show()

    def plot_relative_energy_error(self):
        #Plot relative total-energy error for the current initial conditions
        result = self._simulate_from_inputs()
        if result is None:
            return
        times, states, values = result

        m1, m2 = values["m1"], values["m2"]
        l1, l2 = values["l1"], values["l2"]
        theta1, omega1 = states[:, 0], states[:, 1]
        theta2, omega2 = states[:, 2], states[:, 3]

        kinetic_energy = (
            0.5 * (m1 + m2) * l1**2 * omega1**2
            + 0.5 * m2 * l2**2 * omega2**2
            + m2 * l1 * l2 * omega1 * omega2 * np.cos(theta1 - theta2)
        )
        potential_energy = (
            -(m1 + m2) * G * l1 * np.cos(theta1)
            - m2 * G * l2 * np.cos(theta2)
        )
        total_energy = kinetic_energy + potential_energy
        initial_energy = total_energy[0]

        if np.isclose(initial_energy, 0.0):
            messagebox.showerror(
                "Relative energy error undefined",
                "The initial total energy is zero, so relative energy error "
                "cannot be calculated using |E(t) - E(0)| / |E(0)|."
            )
            return

        relative_error = np.abs(total_energy - initial_energy) / abs(initial_energy)

        fig, ax = plt.subplots(figsize=(9, 5.5))
        ax.plot(times, relative_error, color="tab:purple",
                label="Relative energy error")
        ax.set_title("Relative Energy Error vs Time (RK4)")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel(r"Relative energy error $\epsilon_E(t)$")
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.show()

    # ---------------------- Sensitivity / chaos analysis -------------------
    @staticmethod
    def _wrapped_angle_difference(a, b):
        #Shortest signed angular difference a-b, in radians
        return np.arctan2(np.sin(a - b), np.cos(a - b))

    @classmethod
    def _state_separation(cls, state_a, state_b, omega_scale):
        #Dimensionless distance in angle/angular-velocity state space
        dth1 = cls._wrapped_angle_difference(state_a[0], state_b[0])
        dth2 = cls._wrapped_angle_difference(state_a[2], state_b[2])
        dw1 = (state_a[1] - state_b[1]) / omega_scale
        dw2 = (state_a[3] - state_b[3]) / omega_scale
        return np.sqrt(dth1**2 + dth2**2 + dw1**2 + dw2**2)

    def _simulate_pair(self, values, duration=30.0, dt=0.005,
                       perturbation_deg=0.001):
        #Simulate reference and perturbed systems with identical RK4 steps
        reference = DoublePendulum(**values)
        perturbed_values = values.copy()
        perturbed_values["theta1"] += np.deg2rad(perturbation_deg)
        perturbed = DoublePendulum(**perturbed_values)

        steps = int(round(duration / dt))
        times = np.arange(steps + 1, dtype=float) * dt
        ref_states = np.empty((steps + 1, 4), dtype=float)
        pert_states = np.empty((steps + 1, 4), dtype=float)
        ref_states[0] = reference.state.copy()
        pert_states[0] = perturbed.state.copy()
        for i in range(1, steps + 1):
            reference.step(dt)
            perturbed.step(dt)
            ref_states[i] = reference.state.copy()
            pert_states[i] = perturbed.state.copy()
        return times, ref_states, pert_states, values

    def _get_sensitivity_data(self):
        values = self._read_inputs()
        if values is None:
            return None
        result = self._simulate_pair(values)
        times, ref_states, pert_states, values = result
        omega_scale = np.sqrt(G / values["l1"])
        dth1 = self._wrapped_angle_difference(ref_states[:, 0], pert_states[:, 0])
        dth2 = self._wrapped_angle_difference(ref_states[:, 2], pert_states[:, 2])
        angular_separation = np.sqrt(dth1**2 + dth2**2)
        separation = np.array([
            self._state_separation(a, b, omega_scale)
            for a, b in zip(ref_states, pert_states)
        ])
        return times, ref_states, pert_states, angular_separation, separation

    def plot_sensitivity(self):
        #Plot full state-space distance D(t), mark D(0), and print diagnostics
        result = self._get_sensitivity_data()
        if result is None:
            return
        times, ref_states, pert_states, _, separation = result
        initial_separation = separation[0]

        fig, ax = plt.subplots(figsize=(9, 5.5))
        ax.plot(times, separation, color="tab:blue", label=r"State-space distance $D(t)$")
        ax.scatter([times[0]], [initial_separation], color="tab:red", zorder=5,
                   label=rf"Initial separation $D(0)={initial_separation:.3e}$")
        ax.axhline(initial_separation, color="tab:red", linestyle="--", linewidth=1,
                   alpha=0.8, label=r"Reference level $D(0)$")
        ax.annotate(f"D(0) = {initial_separation:.3e}",
                    xy=(times[0], initial_separation),
                    xytext=(10, 18), textcoords="offset points",
                    fontsize=9, color="tab:red")
        ax.set_title("Sensitivity to Initial Conditions")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel(r"Dimensionless state-space distance $D(t)$")
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.show()

        # Print relative energy error for the reference trajectory
        values = self.current_values if hasattr(self, "current_values") else self._read_inputs()
        if values is None:
            return
        relative_error = self._relative_energy_error(ref_states, values)
        print("\nSensitivity to Initial Conditions")
        print("Initial perturbation in theta1: 0.001 degrees")
        print(f"Initial state-space distance D(0): {initial_separation:.6e}")
        print(f"Final state-space distance D(t_end): {separation[-1]:.6e}")
        print(f"Maximum state-space distance: {np.max(separation):.6e}")
        print("\nRelative Energy Error (RK4)")
        print(f"Final relative energy error:   {relative_error[-1]:.6e}")
        print(f"Maximum relative energy error: {np.max(relative_error):.6e}")
        self.estimate_lyapunov_exponent(values)

    @staticmethod
    def _relative_energy_error(states, values):
        #Return relative total-energy error for a stored trajectory
        m1, m2 = values["m1"], values["m2"]
        l1, l2 = values["l1"], values["l2"]
        th1, w1 = states[:, 0], states[:, 1]
        th2, w2 = states[:, 2], states[:, 3]
        kinetic = (0.5 * (m1 + m2) * l1**2 * w1**2
                   + 0.5 * m2 * l2**2 * w2**2
                   + m2 * l1 * l2 * w1 * w2 * np.cos(th1 - th2))
        potential = (-(m1 + m2) * G * l1 * np.cos(th1)
                     - m2 * G * l2 * np.cos(th2))
        energy = kinetic + potential
        initial_energy = energy[0]
        if np.isclose(initial_energy, 0.0):
            return np.full_like(energy, np.nan, dtype=float)
        return np.abs(energy - initial_energy) / abs(initial_energy)

    def estimate_lyapunov_exponent(self, values=None):
        #Estimate and print the largest Lyapunov exponent
        if values is None:
            values = self._read_inputs()
            if values is None:
                return None
        dt = 0.0025
        duration = 30.0
        renormalization_interval = 0.1
        steps_per_interval = max(1, int(round(renormalization_interval / dt)))
        actual_interval = steps_per_interval * dt
        number_intervals = int(duration / actual_interval)
        omega_scale = np.sqrt(G / values["l1"])
        delta0 = 1e-7

        reference = DoublePendulum(**values)
        perturbed = DoublePendulum(**values)
        perturbed.state[0] += delta0
        accumulated_log_growth = 0.0
        valid_intervals = 0

        for _ in range(number_intervals):
            for _ in range(steps_per_interval):
                reference.step(dt)
                perturbed.step(dt)
            delta = self._state_separation(reference.state, perturbed.state, omega_scale)
            if not np.isfinite(delta) or delta <= 1e-15:
                print("Largest Lyapunov exponent estimate failed: invalid trajectory separation.")
                return None
            accumulated_log_growth += np.log(delta / delta0)
            valid_intervals += 1

            dth1 = self._wrapped_angle_difference(perturbed.state[0], reference.state[0])
            dth2 = self._wrapped_angle_difference(perturbed.state[2], reference.state[2])
            dw1_scaled = (perturbed.state[1] - reference.state[1]) / omega_scale
            dw2_scaled = (perturbed.state[3] - reference.state[3]) / omega_scale
            norm = np.sqrt(dth1**2 + dth2**2 + dw1_scaled**2 + dw2_scaled**2)
            if norm <= 1e-15 or not np.isfinite(norm):
                print("Largest Lyapunov exponent estimate failed during renormalization.")
                return None
            factor = delta0 / norm
            perturbed.state[0] = reference.state[0] + factor * dth1
            perturbed.state[2] = reference.state[2] + factor * dth2
            perturbed.state[1] = reference.state[1] + factor * dw1_scaled * omega_scale
            perturbed.state[3] = reference.state[3] + factor * dw2_scaled * omega_scale

        if valid_intervals == 0:
            print("Largest Lyapunov exponent estimate failed: no valid intervals.")
            return None
        exponent = accumulated_log_growth / (valid_intervals * actual_interval)
        print("\nLargest Lyapunov Exponent Estimate")
        print(f"Integration step: {dt:g} s; renormalization interval: {actual_interval:g} s")
        print(f"Duration: {valid_intervals * actual_interval:g} s")
        print(f"Estimated lambda_max: {exponent:.6e} s^-1")
        return exponent

    # ---------------------------- Animation ---------------------------------
    def _draw(self):
        #Update to match the current pendulum state
        p0, p1, p2 = self.pendulum.positions()
        self.rods.set_data([p0[0], p1[0], p2[0]], [p0[1], p1[1], p2[1]])
        self.bobs.set_data([p1[0], p2[0]], [p1[1], p2[1]])
        self.trail.set_data(self.trail_x, self.trail_y)
        self.time_text.set_text(f"t = {self.sim_time:6.2f} s")
        self.canvas.draw_idle()

    def _animate(self):
        if self.running and self.pendulum is not None:
            dt = FRAME_INTERVAL_MS / 1000.0
            self.pendulum.step(dt)        
            self.sim_time += dt

            #Record the second bob's path for the trail
            _, _, p2 = self.pendulum.positions()
            self.trail_x.append(p2[0])
            self.trail_y.append(p2[1])
            self._draw()

        #Schedule the next frame (~30 fps)
        self.root.after(FRAME_INTERVAL_MS, self._animate)


# ----------------------------------------------------------------------------
# Entry point
def main():
    root = tk.Tk()
    PendulumApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()