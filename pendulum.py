import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ellipk
import matplotlib.pyplot as plt

# ---- Physical constants ----
g = 9.81      # gravitational acceleration (m/s^2)
L = 1.0       # length of the pendulum (m)
w0 = np.sqrt(g / L)          # natural angular frequency (small-angle)
T0 = 2 * np.pi / w0          # small-angle period (independent of amplitude)


def equation_of_motion(t, y):
    """
    The pendulum's equation of motion written as two first-order equations.
    y[0] = theta        (angle)
    y[1] = omega        (angular velocity = d theta/dt)
    Returns [d theta/dt, d omega/dt].
    """
    theta, omega = y
    dtheta_dt = omega
    domega_dt = -(g / L) * np.sin(theta)
    return [dtheta_dt, domega_dt]


def zero_crossing(t, y):
    """Event: the pendulum passes through the bottom (theta = 0)."""
    return y[0]
zero_crossing.terminal = True      # stop integrating at the first crossing
zero_crossing.direction = -1       # only when swinging from + to - (downward)


def measure_period(theta0):
    """
    Release the pendulum FROM REST at angle theta0 and measure its period.

    Trick: by symmetry, the time to swing from theta0 down to theta = 0 is
    exactly a QUARTER of the full period. So we integrate until the first
    zero-crossing and multiply that time by 4.
    """
    sol = solve_ivp(
        equation_of_motion,
        t_span=(0, 10 * T0),          # integrate long enough to catch the crossing
        y0=[theta0, 0.0],             # start at theta0, at rest (omega = 0)
        events=zero_crossing,
        rtol=1e-10, atol=1e-12,       # tight tolerances = accurate period
        dense_output=True,
    )
    quarter_period = sol.t_events[0][0]
    return 4 * quarter_period


def exact_period(theta0):
    """
    Exact period from theory, using the complete elliptic integral K.
      T = (4 / w0) * K( sin^2(theta0 / 2) )
    (SciPy's ellipk takes the parameter m = k^2.)
    """
    m = np.sin(theta0 / 2) ** 2
    return (4 / w0) * ellipk(m)


# ------------------------------------------------------------------
# 1) Period vs amplitude
# ------------------------------------------------------------------
amplitudes_deg = np.arange(5, 171, 5)          # 5 deg to 170 deg
amplitudes_rad = np.radians(amplitudes_deg)

T_numerical = np.array([measure_period(a) for a in amplitudes_rad])
T_exact = np.array([exact_period(a) for a in amplitudes_rad])

# Print a small table so you can SEE the anharmonicity
print(f"Small-angle period T0 = {T0:.4f} s\n")
print(f"{'amp (deg)':>9} | {'T_numerical':>12} | {'T_exact':>9} | {'T/T0':>6}")
print("-" * 48)
for deg, Tn, Te in zip(amplitudes_deg, T_numerical, T_exact):
    if deg % 20 == 0 or deg == 5:
        print(f"{deg:>9} | {Tn:>12.4f} | {Te:>9.4f} | {Tn/T0:>6.3f}")

# ------------------------------------------------------------------
# 2) Plot: period ratio vs amplitude
# ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5))
ax.axhline(1.0, color="gray", ls="--", lw=1.2, label="Small-angle prediction (constant)")
ax.plot(amplitudes_deg, T_exact / T0, "-", color="#1f3a5f", lw=2, label="Exact theory (elliptic K)")
ax.plot(amplitudes_deg, T_numerical / T0, "o", color="#e07a3c", ms=5, label="Numerical simulation")
ax.set_xlabel("Amplitude  $\\theta_0$  (degrees)")
ax.set_ylabel("Period ratio  $T / T_0$")
ax.set_title("Large-angle pendulum: the period grows with amplitude (anharmonicity)")
ax.legend()
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig("period_vs_amplitude.png", dpi=150)

# ------------------------------------------------------------------
# 3) Plot: the waveform distorts at large amplitude
# ------------------------------------------------------------------
fig2, ax2 = plt.subplots(figsize=(8, 5))
for theta0_deg in [10, 90, 160]:
    theta0 = np.radians(theta0_deg)
    T = measure_period(theta0)
    sol = solve_ivp(equation_of_motion, (0, 2 * T), [theta0, 0.0],
                    rtol=1e-10, atol=1e-12, dense_output=True)
    t = np.linspace(0, 2 * T, 1000)
    theta_t = sol.sol(t)[0]
    ax2.plot(t / T0, np.degrees(theta_t), lw=2, label=f"$\\theta_0 = {theta0_deg}°$")
ax2.set_xlabel("Time  $t / T_0$")
ax2.set_ylabel("Angle  $\\theta$  (degrees)")
ax2.set_title("At large amplitude the swing is slower and no longer a clean sine")
ax2.legend()
ax2.grid(alpha=0.3)
fig2.tight_layout()
fig2.savefig("waveform.png", dpi=150)

print("\nSaved plots: period_vs_amplitude.png, waveform.png")
