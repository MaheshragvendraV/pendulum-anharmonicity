# Large-Angle Pendulum: Anharmonicity

A numerical study of how a simple pendulum stops behaving like a simple
harmonic oscillator as its swing amplitude grows — its period increases and
its motion is no longer a clean sine wave. This is **anharmonicity**.

## What this does

- Integrates the exact pendulum equation of motion, `θ'' = -(g/L) sin θ`,
  with no small-angle approximation (using `scipy.integrate.solve_ivp`).
- Measures the oscillation period for amplitudes from 5° to 170°.
- Validates the simulation against the **exact theory** — the period in terms
  of the complete elliptic integral of the first kind,
  `T = (4/ω₀) · K(sin²(θ₀/2))`.

## Key result

The numerical simulation matches the exact theory to ~4 decimal places, and
shows the period growing with amplitude:

| Amplitude | T / T₀ |
|-----------|--------|
| 20°       | 1.008  |
| 60°       | 1.073  |
| 100°      | 1.232  |
| 160°      | 2.008  |

At 160° the pendulum takes **twice as long** per swing as the small-angle
formula predicts.

## Figures

- `period_vs_amplitude.png` — period ratio vs amplitude (numerical vs exact vs small-angle).
- `waveform.png` — the swing waveform distorting away from a sine at large amplitude.

## Run it

```bash
pip install numpy scipy matplotlib
python pendulum.py
```

## Possible extensions

- Add an **energy-conservation check** to verify the integrator's accuracy.
- Fit the small-amplitude data to the series `T/T₀ ≈ 1 + θ₀²/16 + ...`.
- Extend to **two coupled pendulums** and study energy exchange / normal modes.
- Add damping and a driving force to explore the route to **chaos**.

## Author

Maheshragavendra V — undergraduate physics, PSG College of Technology.
Built as a computational companion to a large-angle / coupled-pendulum project.
