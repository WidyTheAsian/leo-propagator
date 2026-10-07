"""Cartesian propagation with optional J2; the CLI runs milestone 1 in SI units."""

import argparse
import json
from pathlib import Path
import platform

import matplotlib
import numpy as np
import scipy
from scipy.integrate import solve_ivp


# Adopted model constants, rather than an Earth shape or gravity field model.
MU = 3.986004418e14  # Earth's gravitational parameter [m^3 / s^2]
EARTH_RADIUS = 6_371_000.0  # spherical mean radius [m]
# NASA Earth fact sheet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
EQUATORIAL_RADIUS = 6_378_137.0  # J2 gravity reference radius [m], not altitude radius
J2 = 1.08263e-3  # positive, dimensionless; symmetry axis is inertial +z
INITIAL_ALTITUDE = 400_000.0  # [m]
ORBIT_RADIUS = EARTH_RADIUS + INITIAL_ALTITUDE
CIRCULAR_SPEED = np.sqrt(MU / ORBIT_RADIUS)  # [m / s]
MEAN_MOTION = np.sqrt(MU / ORBIT_RADIUS**3)  # [rad / s]
PERIOD = 2.0 * np.pi / MEAN_MOTION  # [s]
DURATION = 86_400.0  # [s]
SAMPLE_INTERVAL = 10.0  # output spacing [s], not the integration step size
RTOL = 1e-10  # dimensionless
ATOL = np.array([1e-4] * 3 + [1e-7] * 3)  # [m, m, m, m/s, m/s, m/s]


def initial_state(inclination=0.0):
    """Circular Keplerian state on +x; velocity tilted about x by inclination [rad]."""
    return np.array([ORBIT_RADIUS, 0.0, 0.0, 0.0,
                     CIRCULAR_SPEED * np.cos(inclination),
                     CIRCULAR_SPEED * np.sin(inclination)])


def j2_acceleration(position, j2=J2):
    """J2 contribution [m/s^2] for positions [m] of shape (3,) or (3, N).

    Derived as -grad(Phi_J2) from the potential in specific_potential.
    J2 is positive for Earth's oblateness; it equals -C20 when unnormalized.
    """
    position = np.asarray(position)
    radius = np.linalg.norm(position, axis=0)
    z_ratio_squared = (position[2] / radius)**2
    factor = 1.5 * j2 * MU * EQUATORIAL_RADIUS**2 / radius**5
    return factor * position * np.array([
        5.0 * z_ratio_squared - 1.0,
        5.0 * z_ratio_squared - 1.0,
        5.0 * z_ratio_squared - 3.0,
    ])


def specific_potential(position, j2=0.0):
    """Central + J2 gravitational potential [m^2/s^2], zero at infinity.

    Phi = -mu/r * [1 - J2*(R_eq/r)^2*P2(z/r)]. The Hamiltonian/potential
    convention is given in Arnas (2022), Eq. (1): https://arxiv.org/abs/2212.08746
    Shape (3,) gives a scalar; shape (3, N) gives N potential values.
    """
    position = np.asarray(position)
    radius = np.linalg.norm(position, axis=0)
    potential = -MU / radius
    if j2 != 0.0:
        p2 = 0.5 * (3.0 * (position[2] / radius)**2 - 1.0)
        potential += MU / radius * j2 * (EQUATORIAL_RADIUS / radius)**2 * p2
    return potential


def two_body_rhs(t, state, j2=0.0):
    """Return [velocity, central + optional J2 acceleration], independent of t."""
    position = state[:3]
    velocity = state[3:]
    acceleration = -MU * position / np.linalg.norm(position)**3
    if j2 != 0.0:
        acceleration += j2_acceleration(position, j2)
    return np.concatenate((velocity, acceleration))


def analytical_state(times, inclination=0.0):
    """Exact two-body circular states (6, N); elapsed seconds, inclination radians."""
    phase = MEAN_MOTION * np.atleast_1d(times)
    cosine, sine = np.cos(phase), np.sin(phase)
    cos_i, sin_i = np.cos(inclination), np.sin(inclination)
    return np.vstack((
        ORBIT_RADIUS * cosine, ORBIT_RADIUS * sine * cos_i, ORBIT_RADIUS * sine * sin_i,
        -CIRCULAR_SPEED * sine, CIRCULAR_SPEED * cosine * cos_i,
        CIRCULAR_SPEED * cosine * sin_i,
    ))


def propagate(rtol=RTOL, atol=ATOL, *, state0=None, duration=DURATION,
              sample_interval=SAMPLE_INTERVAL, j2=0.0):
    """Integrate Cartesian states in SI; defaults preserve milestone 1.

    J2 = 0 disables the perturbation. Always include the final time, even when
    the duration is not an integer multiple of the output sample interval.
    """
    if not np.isfinite(duration) or duration <= 0.0:
        raise ValueError("Duration must be finite and positive.")
    if not np.isfinite(sample_interval) or sample_interval <= 0.0:
        raise ValueError("Sample interval must be finite and positive.")
    if state0 is None:
        state0 = initial_state()
    times = np.append(np.arange(0.0, duration, sample_interval), duration)
    solution = solve_ivp(
        two_body_rhs, (0.0, duration), state0, args=(j2,),
        method="DOP853", t_eval=times, rtol=rtol, atol=atol,
    )
    if not solution.success:
        raise RuntimeError(f"Orbit integration failed: {solution.message}")
    if solution.t[-1] != duration or not np.all(np.isfinite(solution.y)):
        raise RuntimeError("Integration did not return a finite full-duration trajectory.")
    return solution


def diagnostics(times, states):
    """Milestone 1 equatorial two-body errors; not the invariants for J2 motion."""
    position, velocity = states[:3], states[3:]
    radius = np.linalg.norm(position, axis=0)
    energy = 0.5 * np.sum(velocity**2, axis=0) - MU / radius
    angular_momentum = np.cross(position.T, velocity.T).T
    exact = analytical_state(times)
    exact_energy = -MU / (2.0 * ORBIT_RADIUS)
    exact_h = np.array([0.0, 0.0, ORBIT_RADIUS * CIRCULAR_SPEED])
    energy_error = energy - exact_energy
    # Use the full vector: a magnitude alone could hide a change of orbital plane.
    h_error = np.linalg.norm(angular_momentum - exact_h[:, None], axis=0)
    return {
        "altitude_m": radius - EARTH_RADIUS,
        "altitude_error_m": radius - ORBIT_RADIUS,
        "position_error_m": np.linalg.norm(position - exact[:3], axis=0),
        "velocity_error_m_s": np.linalg.norm(velocity - exact[3:], axis=0),
        "specific_energy_error_m2_s2": energy_error,
        "relative_energy_error": energy_error / abs(exact_energy),
        "angular_momentum_error_m2_s": h_error,
        "relative_angular_momentum_error": h_error / np.linalg.norm(exact_h),
    }


def summarize(errors):
    """Maximum absolute and final errors on the saved output grid."""
    summary = {}
    for name, values in errors.items():
        if name != "altitude_m":
            summary[f"max_abs_{name}"] = float(np.max(np.abs(values)))
            summary[f"final_{name}"] = float(values[-1])
    return summary


def plot_results(times, states, errors, output_dir, show=False):
    """Save figures; use a headless backend unless a plot window is requested."""
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle

    fig, (orbit_ax, altitude_ax) = plt.subplots(1, 2, figsize=(12, 5))
    orbit_ax.add_patch(Circle((0.0, 0.0), EARTH_RADIUS, color="lightsteelblue",
                             label="Spherical Earth"))
    orbit_ax.plot(states[0], states[1], label="Numerical: 24 hours")
    reference = analytical_state(np.linspace(0.0, PERIOD, 600))
    orbit_ax.plot(reference[0], reference[1], "k--", lw=1, label="Analytical circle")
    orbit_ax.plot(states[0, 0], states[1, 0], "ro", ms=4, label="Initial position")
    orbit_ax.set(xlabel="ECI x [m]", ylabel="ECI y [m]", title="Orbit in the xy-plane (z = 0)")
    orbit_ax.set_aspect("equal", adjustable="box")
    orbit_ax.legend(loc="upper left", fontsize=8)

    altitude_ax.plot(times, errors["altitude_m"], label="Numerical altitude")
    altitude_ax.axhline(INITIAL_ALTITUDE, color="black", ls="--", label="Exact altitude")
    altitude_ax.set(xlabel="Elapsed time [s]", ylabel="Altitude [m]",
                    title="Altitude above the spherical surface")
    # A +/-1 m range keeps the physical altitude explicit despite tiny errors.
    altitude_ax.set_ylim(INITIAL_ALTITUDE - 1.0, INITIAL_ALTITUDE + 1.0)
    altitude_ax.ticklabel_format(axis="y", useOffset=False, style="plain")
    max_altitude_error = np.max(np.abs(errors["altitude_error_m"]))
    altitude_ax.text(0.04, 0.08, f"Max altitude error: {max_altitude_error:.3e} m",
                     transform=altitude_ax.transAxes)
    altitude_ax.legend(fontsize=8)
    for ax in (orbit_ax, altitude_ax):
        ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_dir / "orbit_and_altitude.png", dpi=180)

    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    panels = [
        ("position_error_m", "Position error norm [m]"),
        ("velocity_error_m_s", "Velocity error norm [m/s]"),
        ("relative_energy_error", "Signed relative specific energy error"),
        ("relative_angular_momentum_error", "Relative specific angular momentum vector error"),
    ]
    for ax, (name, title) in zip(axes.flat, panels):
        ax.plot(times, errors[name])
        ax.set(title=title, xlabel="Elapsed time [s]")
        ax.grid(alpha=0.25)
    fig.suptitle("Validation against the exact circular orbit")
    fig.tight_layout()
    fig.savefig(output_dir / "validation.png", dpi=180)
    if show:
        plt.show()
    plt.close("all")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--show", action="store_true", help="Also open plot windows.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args()

    solution = propagate()
    errors = diagnostics(solution.t, solution.y)
    nominal = summarize(errors)
    # Same dynamics, initial state, method and output times; change only tolerances.
    tighter = propagate(rtol=RTOL / 100.0, atol=ATOL / 100.0)
    tighter_summary = summarize(diagnostics(tighter.t, tighter.y))
    improvement = (nominal["max_abs_position_error_m"] /
                   tighter_summary["max_abs_position_error_m"])
    report = {
        "versions": {"python": platform.python_version(), "numpy": np.__version__,
                     "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
        "model": {"mu_m3_s2": MU, "earth_radius_m": EARTH_RADIUS,
                  "initial_altitude_m": INITIAL_ALTITUDE, "orbit_radius_m": ORBIT_RADIUS,
                  "initial_speed_m_s": float(CIRCULAR_SPEED),
                  "period_s": float(PERIOD), "duration_s": DURATION,
                  "exact_specific_energy_m2_s2": -MU / (2.0 * ORBIT_RADIUS),
                  "exact_specific_angular_momentum_m2_s": ORBIT_RADIUS * float(CIRCULAR_SPEED)},
        "solver": {"method": "DOP853", "rtol": RTOL, "atol": ATOL.tolist(),
                   "output_interval_s": SAMPLE_INTERVAL, "output_samples": len(solution.t),
                   "rhs_evaluations": solution.nfev, "message": solution.message},
        "nominal": nominal,
        "tighter": {"rtol": RTOL / 100.0, "atol": (ATOL / 100.0).tolist(),
                    "rhs_evaluations": tighter.nfev, **tighter_summary},
        "max_position_error_improvement_factor": improvement,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "validation.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8",
    )
    plot_results(solution.t, solution.y, errors, args.output_dir, args.show)

    print(f"Initial circular speed: {CIRCULAR_SPEED:.9f} m/s")
    print(f"Analytical period: {PERIOD:.9f} s")
    print(f"Revolutions in 24 hours: {DURATION / PERIOD:.9f}")
    print(f"Solver: {solution.message}")
    print(f"Saved samples: {len(solution.t)}; RHS evaluations: {solution.nfev}")
    print("\nMaximum absolute errors over 24 hours:")
    for name, value in nominal.items():
        if name.startswith("max_abs_"):
            print(f"  {name.removeprefix('max_abs_')}: {value:.6e}")
    print(f"\nTighter-run maximum position error: "
          f"{tighter_summary['max_abs_position_error_m']:.6e} m")
    print(f"Position error improvement: {improvement:.2f} times")
    print(f"Figures and report saved in: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
