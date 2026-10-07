"""Milestone 2: compare two-body and J2 motion from one inclined initial state."""

import argparse
import json
from pathlib import Path
import platform

import matplotlib
import numpy as np
import scipy

import orbit
from elements import NODE_TOL, osculating_elements


EXAMPLE_INCLINATION = np.deg2rad(51.6)  # radians internally
COMPARISON_DURATION = 7.0 * 86_400.0  # [s]
OUTPUT_INTERVAL = 60.0  # [s], independent of adaptive integration steps


def secular_raan_rate(a, eccentricity, inclination, j2=orbit.J2):
    """First-order mean-element J2 node rate [rad/s]; a [m], inclination [rad].

    NASA Rickman (2014), slide 103:
    https://ntrs.nasa.gov/api/citations/20200007746/downloads/20200007746.pdf#page=103
    A circular inclined orbit has a node; an equatorial orbit does not.
    """
    if not (np.isfinite(a) and a > 0.0 and 0.0 <= eccentricity < 1.0):
        raise ValueError("This secular prediction requires a bound Keplerian ellipse.")
    if not np.isfinite(inclination) or not 0.0 <= inclination <= np.pi:
        raise ValueError("Inclination must be in [0, pi] radians.")
    if abs(np.sin(inclination)) <= NODE_TOL:
        raise ValueError("RAAN and its drift are undefined for an equatorial orbit.")
    p = a * (1.0 - eccentricity**2)
    n = np.sqrt(orbit.MU / a**3)
    return -1.5 * j2 * n * (orbit.EQUATORIAL_RADIUS / p)**2 * np.cos(inclination)


def fit_raan_drift(times, raan):
    """Unwrap RAAN and fit intercept + slope*time, with a free intercept.

    Sample finely enough that physical node motion between samples is less
    than pi. Reject undefined angles instead of filling gaps or inventing nodes.
    """
    times, raan = np.asarray(times), np.asarray(raan)
    if (times.ndim != 1 or times.shape != raan.shape or len(times) < 2
            or not np.all(np.isfinite(times)) or not np.all(np.isfinite(raan))
            or np.any(np.diff(times) <= 0.0)):
        raise ValueError("RAAN fitting requires defined angles and increasing finite times.")
    unwrapped = np.unwrap(raan)
    centered = times - np.mean(times)
    slope = np.dot(centered, unwrapped - np.mean(unwrapped)) / np.dot(centered, centered)
    intercept = np.mean(unwrapped) - slope * np.mean(times)
    return {
        "slope_rad_s": float(slope),
        "intercept_rad": float(intercept),
        "unwrapped_rad": unwrapped,
        "residual_rad": unwrapped - (intercept + slope * times),
    }


def conservation_diagnostics(states, j2):
    """Errors in matching total specific energy and h_z, relative to t=0.

    Normalize h_z drift by the initial |h|, so the error scale is also defined
    for polar orbits with h_z = 0. Full h is not a J2 conservation diagnostic.
    """
    position, velocity = states[:3], states[3:]
    energy = 0.5 * np.sum(velocity**2, axis=0) + orbit.specific_potential(position, j2)
    h = np.cross(position.T, velocity.T)
    energy_error = energy - energy[0]
    hz_error = h[:, 2] - h[0, 2]
    return {
        "specific_energy_error_m2_s2": energy_error,
        "relative_energy_error": energy_error / abs(energy[0]),
        "hz_error_m2_s": hz_error,
        "relative_hz_error": hz_error / np.linalg.norm(h[0]),
    }


def run_comparison(duration=COMPARISON_DURATION):
    """Run one central case and three J2 tolerance levels; return data and report."""
    state0 = orbit.initial_state(EXAMPLE_INCLINATION)
    settings = {"state0": state0, "duration": duration, "sample_interval": OUTPUT_INTERVAL}
    central = orbit.propagate(**settings, j2=0.0)
    nominal = orbit.propagate(**settings, j2=orbit.J2)
    tighter = orbit.propagate(**settings, j2=orbit.J2,
                              rtol=orbit.RTOL / 100.0, atol=orbit.ATOL / 100.0)
    reference = orbit.propagate(**settings, j2=orbit.J2,
                                rtol=orbit.RTOL / 1000.0, atol=orbit.ATOL / 1000.0)
    solutions = {"two_body": central, "nominal": nominal,
                 "tighter": tighter, "reference": reference}
    elements = {name: osculating_elements(solution.y) for name, solution in solutions.items()}
    fits = {name: fit_raan_drift(solution.t, elements[name]["raan_rad"])
            for name, solution in solutions.items()}
    invariants = {name: conservation_diagnostics(solution.y, 0.0 if name == "two_body" else orbit.J2)
                  for name, solution in solutions.items()}
    initial = elements["nominal"]
    predicted = secular_raan_rate(initial["semimajor_axis_m"][0],
                                  initial["eccentricity"][0], initial["inclination_rad"][0])
    # This sample average is a diagnostic proxy, not a formal mean-element transformation.
    mean_proxy = {name: float(np.mean(initial[name])) for name in
                  ("semimajor_axis_m", "eccentricity", "inclination_rad")}
    proxy_rate = secular_raan_rate(mean_proxy["semimajor_axis_m"], mean_proxy["eccentricity"],
                                   mean_proxy["inclination_rad"])
    measured = fits["nominal"]["slope_rad_s"]
    circular = orbit.analytical_state(central.t, EXAMPLE_INCLINATION)
    convergence = {}
    for name in ("nominal", "tighter"):
        difference = solutions[name].y - reference.y
        convergence[f"{name}_position_difference_m"] = np.linalg.norm(difference[:3], axis=0)
        convergence[f"{name}_velocity_difference_m_s"] = np.linalg.norm(difference[3:], axis=0)

    h = np.cross(nominal.y[:3].T, nominal.y[3:].T)
    keplerian_energy = (0.5 * np.sum(nominal.y[3:]**2, axis=0)
                        - orbit.MU / np.linalg.norm(nominal.y[:3], axis=0))
    report = {
        "versions": {"python": platform.python_version(), "numpy": np.__version__,
                     "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
        "model": {"mu_m3_s2": orbit.MU, "j2": orbit.J2,
                  "j2_reference_radius_m": orbit.EQUATORIAL_RADIUS,
                  "altitude_reference_radius_m": orbit.EARTH_RADIUS,
                  "symmetry_axis": "Fixed inertial z-axis", "duration_s": duration,
                  "inclination_rad": float(EXAMPLE_INCLINATION), "initial_state_si": state0.tolist()},
        "initial_osculating_elements": {
            name: float(values[0]) if np.isfinite(values[0]) else None
            for name, values in initial.items() if name not in ("circular", "equatorial")
        },
        "solver": {"method": "DOP853", "output_interval_s": OUTPUT_INTERVAL,
                   "output_samples": len(nominal.t),
                   "tolerance_levels": {
                       name: {"rtol": orbit.RTOL / scale,
                              "atol": (orbit.ATOL / scale).tolist(),
                              "rhs_evaluations": solutions[name].nfev}
                       for name, scale in (("two_body", 1.0), ("nominal", 1.0),
                                           ("tighter", 100.0), ("reference", 1000.0))}},
        "two_body": {
            "max_position_error_vs_circular_m": float(np.max(np.linalg.norm(central.y[:3] - circular[:3], axis=0))),
            "max_velocity_error_vs_circular_m_s": float(np.max(np.linalg.norm(central.y[3:] - circular[3:], axis=0))),
            "fitted_raan_rate_rad_s": fits["two_body"]["slope_rad_s"],
        },
        "invariants": {"hz_relative_error_normalization": "Initial angular momentum magnitude |h(0)|",
                       **{name: orbit.summarize(errors) for name, errors in invariants.items()}},
        "altitude_bounds_m": {
            name: {"min": float(np.min(np.linalg.norm(solutions[name].y[:3], axis=0) - orbit.EARTH_RADIUS)),
                   "max": float(np.max(np.linalg.norm(solutions[name].y[:3], axis=0) - orbit.EARTH_RADIUS))}
            for name in ("two_body", "nominal")},
        "raan": {"measured_rad_s": measured,
                 "predicted_from_initial_elements_rad_s": float(predicted),
                 "relative_discrepancy_initial_elements": float(abs(measured - predicted) / abs(predicted)),
                 "mean_proxy_elements": mean_proxy,
                 "mean_proxy_note": "Arithmetic sample averages over the run; not a formal osculating-to-mean conversion.",
                 "predicted_from_mean_proxy_rad_s": float(proxy_rate),
                 "relative_discrepancy_mean_proxy": float(abs(measured - proxy_rate) / abs(proxy_rate)),
                 "max_abs_fit_residual_rad": float(np.max(np.abs(fits["nominal"]["residual_rad"]))),
                 "final_unwrapped_rad": float(fits["nominal"]["unwrapped_rad"][-1]),
                 "nominal_tighter_rate_difference_rad_s": abs(measured - fits["tighter"]["slope_rad_s"])},
        "convergence": {f"max_{name}": float(np.max(values)) for name, values in convergence.items()},
        "physical_variations_not_conservation_errors": {
            "max_two_body_specific_energy_change_m2_s2": float(np.max(np.abs(keplerian_energy - keplerian_energy[0]))),
            "max_full_h_vector_change_relative": float(np.max(np.linalg.norm(h - h[0], axis=1)) / np.linalg.norm(h[0])),
        },
    }
    report["convergence"]["position_improvement_factor"] = (
        report["convergence"]["max_nominal_position_difference_m"] /
        report["convergence"]["max_tighter_position_difference_m"])
    return {"solutions": solutions, "fits": fits, "invariants": invariants,
            "convergence": convergence, "report": report}


def plot_comparison(data, output_dir, show=False):
    """Save the trajectories, spherical altitude, RAAN and convergence figures."""
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    central, nominal = data["solutions"]["two_body"], data["solutions"]["nominal"]
    times = nominal.t
    fit = data["fits"]["nominal"]
    fig = plt.figure(figsize=(12, 9))
    trajectory_ax = fig.add_subplot(2, 2, 1, projection="3d")
    altitude_ax = fig.add_subplot(2, 2, 2)
    raan_ax = fig.add_subplot(2, 2, 3)
    residual_ax = fig.add_subplot(2, 2, 4)
    # Draw the fixed central orbit last so it stays visible over the precessing paths.
    for name, solution, color, width, alpha in (
            ("Central + J2", nominal, "tab:orange", 0.5, 0.5),
            ("Two body", central, "tab:blue", 1.3, 1.0)):
        trajectory_ax.plot(*solution.y[:3], color=color, lw=width, alpha=alpha, label=name)
        altitude_ax.plot(times, np.linalg.norm(solution.y[:3], axis=0) - orbit.EARTH_RADIUS,
                         color=color, lw=0.8, label=name)
    limit = orbit.ORBIT_RADIUS * 1.05
    trajectory_ax.set(xlabel="ECI x [m]", ylabel="ECI y [m]", zlabel="ECI z [m]",
                      xlim=(-limit, limit), ylim=(-limit, limit), zlim=(-limit, limit),
                      title="Trajectories over the full comparison")
    trajectory_ax.set_box_aspect((1, 1, 1))
    trajectory_ax.view_init(elev=20, azim=35)
    trajectory_ax.legend(fontsize=8)
    altitude_ax.set(xlabel="Elapsed time [s]", ylabel="Spherical altitude [m]", title="Altitude")
    altitude_ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    altitude_ax.legend(fontsize=8)
    raan_ax.plot(times, data["fits"]["two_body"]["unwrapped_rad"], label="Two body")
    raan_ax.plot(times, fit["unwrapped_rad"], color="tab:orange", label="J2 osculating RAAN")
    raan_ax.plot(times, fit["intercept_rad"] + fit["slope_rad_s"] * times,
                 "k--", lw=1, label="Fitted secular trend")
    predicted = data["report"]["raan"]["predicted_from_initial_elements_rad_s"]
    raan_ax.plot(times, predicted * times, color="tab:red", ls=":", label="First order, initial elements")
    raan_ax.set(xlabel="Elapsed time [s]", ylabel="Unwrapped RAAN [rad]", title="Nodal regression")
    raan_ax.legend(fontsize=8)
    early = times <= 2.0 * orbit.PERIOD
    residual_ax.plot(times[early], fit["residual_rad"][early], color="tab:orange")
    residual_ax.set(xlabel="Elapsed time [s]", ylabel="RAAN minus fitted trend [rad]",
                    title="Short-period variation: first two nominal periods")
    for ax in (altitude_ax, raan_ax, residual_ax):
        ax.grid(alpha=0.25)
    fig.suptitle("Shared initial state: 400,000 m altitude, inclination 51.6 degrees")
    fig.tight_layout()
    fig.savefig(output_dir / "comparison.png", dpi=160)

    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    for name, color in (("nominal", "tab:orange"), ("tighter", "tab:green")):
        errors = data["invariants"][name]
        axes[0, 0].plot(times, errors["relative_energy_error"], color=color, lw=0.8, label=name)
        axes[0, 1].plot(times, errors["relative_hz_error"], color=color, lw=0.8, label=name)
        for ax, quantity in ((axes[1, 0], "position_difference_m"),
                              (axes[1, 1], "velocity_difference_m_s")):
            values = data["convergence"][f"{name}_{quantity}"]
            ax.semilogy(times[1:], values[1:], color=color, label=name)
    titles = ("Relative total specific energy error", "h_z error / initial |h|",
              "Position difference from reference [m]", "Velocity difference from reference [m/s]")
    for ax, title in zip(axes.flat, titles):
        ax.set(title=title, xlabel="Elapsed time [s]")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.25)
    fig.suptitle("J2 verification: matching potential, axial momentum, and three tolerance levels")
    fig.tight_layout()
    fig.savefig(output_dir / "validation.png", dpi=160)
    if show:
        plt.show()
    plt.close("all")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--days", type=float, default=7.0, help="Duration in days; internal time is seconds.")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/j2"))
    parser.add_argument("--show", action="store_true", help="Also open plot windows.")
    args = parser.parse_args()
    if not np.isfinite(args.days) or args.days <= 0.0:
        parser.error("--days must be finite and positive")
    data = run_comparison(args.days * 86_400.0)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plot_comparison(data, args.output_dir, args.show)
    report = data["report"]
    (args.output_dir / "validation.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Duration: {report['model']['duration_s']:.0f} s; samples: {report['solver']['output_samples']}")
    print(f"Measured RAAN rate: {report['raan']['measured_rad_s']:.12e} rad/s")
    print(f"First-order rate, initial elements: {report['raan']['predicted_from_initial_elements_rad_s']:.12e} rad/s")
    print(f"Relative rate discrepancy: {report['raan']['relative_discrepancy_initial_elements']:.6%}")
    print(f"First-order rate, sample mean proxy: {report['raan']['predicted_from_mean_proxy_rad_s']:.12e} rad/s")
    for name in ("nominal", "tighter", "reference"):
        errors = report["invariants"][name]
        print(f"{name}: max relative energy error {errors['max_abs_relative_energy_error']:.6e}; "
              f"max h_z error / |h(0)| {errors['max_abs_relative_hz_error']:.6e}")
    print(f"Position differences from reference: nominal "
          f"{report['convergence']['max_nominal_position_difference_m']:.6e} m; tighter "
          f"{report['convergence']['max_tighter_position_difference_m']:.6e} m")
    print(f"Position improvement: {report['convergence']['position_improvement_factor']:.2f} times")
    print(f"Figures and report saved in: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
