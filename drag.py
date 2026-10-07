"""Milestone 3: seven-day density sensitivity with J2 and illustrative drag."""

import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path
import platform

import matplotlib
import numpy as np
import scipy
from scipy.optimize import brentq

import orbit
from atmosphere import DragModel
from elements import osculating_elements
from j2 import COMPARISON_DURATION, EXAMPLE_INCLINATION, OUTPUT_INTERVAL


DENSITY_MULTIPLIERS = (0.5, 1.0, 2.0)
CHECKPOINTS = (86_400.0, 259_200.0, 604_800.0)
WORK_ATOL = 1e-5  # specific work [m^2/s^2], the seventh integrated component


def rtn_difference(reference_states, states):
    """Project r - r_reference onto reference R,T,N; return shape (3,N) [m].

    R = r/|r|, N = (r cross v)/|r cross v|, T = N cross R. T is the
    transverse/along-track axis; it need not be parallel to v on an ellipse.
    Both inputs must be finite, same-time Cartesian states (6,) or (6,N).
    """
    reference = np.asarray(reference_states, dtype=float)
    states = np.asarray(states, dtype=float)
    if reference.ndim == 1:
        reference = reference[:, None]
    if states.ndim == 1:
        states = states[:, None]
    if (reference.shape != states.shape or reference.ndim != 2
            or reference.shape[0] != 6 or not np.all(np.isfinite(reference))
            or not np.all(np.isfinite(states))):
        raise ValueError("RTN differences require matching finite Cartesian states.")
    radius = np.linalg.norm(reference[:3], axis=0)
    h = np.cross(reference[:3].T, reference[3:].T).T
    h_norm = np.linalg.norm(h, axis=0)
    if np.any(radius == 0.0) or np.any(h_norm == 0.0):
        raise ValueError("The reference state must define an orbital RTN frame.")
    radial, normal = reference[:3] / radius, h / h_norm
    transverse = np.cross(normal.T, radial.T).T
    difference = states[:3] - reference[:3]
    return np.vstack([np.sum(axis * difference, axis=0)
                      for axis in (radial, transverse, normal)])


def nodal_windows(solution):
    """Complete ascending-node to ascending-node time windows [s].

    Bracket with the saved states, refine with dense output; assumes output
    spacing resolves node crossings. Equatorial orbits have undefined nodes.
    The example starts at an ascending node. Partial end orbits are omitted.
    """
    if np.any(osculating_elements(solution.y)["equatorial"]):
        raise ValueError("Nodal averaging is undefined for equatorial orbits.")
    z = solution.y[2]
    crossings = [solution.t[0]] if z[0] == 0.0 and solution.y[5, 0] > 0.0 else []
    brackets = np.flatnonzero((z[:-1] < 0.0) & (z[1:] >= 0.0))
    for index in brackets:
        crossings.append(brentq(lambda t: solution.sol(t)[2],
                                 solution.t[index], solution.t[index + 1], xtol=1e-7))
    return np.array(list(zip(crossings[:-1], crossings[1:])), dtype=float).reshape(-1, 2)


def average_semimajor_axis(solution, windows, quadrature_order=32):
    """Time-average osculating a in each window using dense-output quadrature.

    Gauss-Legendre integration avoids aliasing short-period J2 oscillations on
    a coarse output grid. These averages are not formal mean orbital elements.
    """
    nodes, weights = np.polynomial.legendre.leggauss(quadrature_order)
    averages = []
    for start, end in windows:
        times = (start + end) / 2.0 + (end - start) / 2.0 * nodes
        a = osculating_elements(solution.sol(times))["semimajor_axis_m"]
        averages.append(0.5 * np.dot(weights, a))
    return np.asarray(averages)


def work_energy_diagnostics(solution, model, j2=orbit.J2):
    """Check E(t)-E(0)-W(t), with matching gravity and inertial-velocity work."""
    position, velocity = solution.y[:3], solution.y[3:]
    altitude = np.linalg.norm(position, axis=0) - orbit.EARTH_RADIUS
    acceleration = model.acceleration(position, velocity)
    energy = 0.5 * np.sum(velocity**2, axis=0) + orbit.specific_potential(position, j2)
    energy_change = energy - energy[0]
    residual = energy_change - solution.work_m2_s2
    return {
        "altitude_m": altitude,
        "density_kg_m3": model.density(altitude),
        "drag_acceleration_m_s2": np.linalg.norm(acceleration, axis=0),
        "specific_energy_change_m2_s2": energy_change,
        "specific_work_m2_s2": solution.work_m2_s2,
        "work_energy_residual_m2_s2": residual,
        "relative_work_energy_residual": residual / abs(energy[0]),
    }


def run_experiments(duration=COMPARISON_DURATION, model=None):
    """Return runs, plots' arrays and a report; do not assume guards reach t_end."""
    model = DragModel() if model is None else model
    state0 = orbit.initial_state(EXAMPLE_INCLINATION)
    settings = {"state0": state0, "duration": duration, "sample_interval": OUTPUT_INTERVAL,
                "j2": orbit.J2, "stop_altitude": model.min_altitude_m, "dense_output": True}
    baseline = orbit.propagate(**settings)
    zero = orbit.propagate(**settings, drag=replace(model, density_multiplier=0.0))
    windows = nodal_windows(baseline)
    baseline_means = average_semimajor_axis(baseline, windows)
    cases, case_reports = {}, {}
    for multiplier in DENSITY_MULTIPLIERS:
        case_model = replace(model, density_multiplier=multiplier)
        levels = {}
        for name, scale in (("nominal", 1.0), ("tighter", 100.0), ("reference", 1000.0)):
            solution = orbit.propagate(**settings, drag=case_model,
                                       rtol=orbit.RTOL / scale, atol=orbit.ATOL / scale,
                                       work_atol=WORK_ATOL / scale)
            levels[name] = {"solution": solution,
                            "diagnostics": work_energy_diagnostics(solution, case_model)}
        nominal = levels["nominal"]["solution"]
        # Never extrapolate a stopped trajectory, even for checkpoints or plots.
        shared = nominal.t <= baseline.t[-1]
        rtn = rtn_difference(baseline.sol(nominal.t[shared]), nominal.y[:, shared])
        available_windows = windows[windows[:, 1] <= nominal.t[-1]]
        means = average_semimajor_axis(nominal, available_windows)
        checkpoints = {}
        for checkpoint in CHECKPOINTS:
            if checkpoint > min(nominal.t[-1], baseline.t[-1]):
                checkpoints[str(int(checkpoint))] = {"available": False,
                                                    "reason": "Beyond requested or stopped trajectory."}
                continue
            state, ref_state = nominal.sol(checkpoint), baseline.sol(checkpoint)
            a = osculating_elements(state)["semimajor_axis_m"][0]
            ref_a = osculating_elements(ref_state)["semimajor_axis_m"][0]
            complete = np.flatnonzero(available_windows[:, 1] <= checkpoint)
            index = complete[-1] if len(complete) else None
            checkpoints[str(int(checkpoint))] = {
                "available": True,
                "position_difference_rtn_m": rtn_difference(ref_state, state)[:, 0].tolist(),
                "position_difference_norm_m": float(np.linalg.norm(state[:3] - ref_state[:3])),
                "altitude_difference_m": float(np.linalg.norm(state[:3]) - np.linalg.norm(ref_state[:3])),
                "osculating_semimajor_axis_difference_m": float(a - ref_a),
                "last_complete_reference_orbit_window_s": (
                    available_windows[index].tolist() if index is not None else None),
                "orbit_averaged_semimajor_axis_difference_m": (
                    float(means[index] - baseline_means[index]) if index is not None else None),
            }
        convergence, convergence_arrays = {}, {}
        reference = levels["reference"]["solution"]
        for name in ("nominal", "tighter"):
            solution = levels[name]["solution"]
            common = solution.t <= reference.t[-1]
            difference = solution.y[:, common] - reference.sol(solution.t[common])
            position = np.linalg.norm(difference[:3], axis=0)
            velocity = np.linalg.norm(difference[3:], axis=0)
            convergence[f"max_{name}_position_difference_m"] = float(np.max(position))
            convergence[f"max_{name}_velocity_difference_m_s"] = float(np.max(velocity))
            convergence_arrays[name] = {"time_s": solution.t[common],
                                         "position_m": position, "velocity_m_s": velocity}
        denominator = convergence["max_tighter_position_difference_m"]
        convergence["position_improvement_factor"] = (
            convergence["max_nominal_position_difference_m"] / denominator if denominator else None)
        summaries = {}
        for name, scale in (("nominal", 1.0), ("tighter", 100.0), ("reference", 1000.0)):
            solution = levels[name]["solution"]
            diagnostics = levels[name]["diagnostics"]
            summaries[name] = {
                "rtol": orbit.RTOL / scale, "cartesian_atol_si": (orbit.ATOL / scale).tolist(),
                "work_atol_m2_s2": WORK_ATOL / scale, "rhs_evaluations": solution.nfev,
                "end_time_s": float(solution.t[-1]), "stopped_at_altitude_limit": solution.terminated_early,
                "output_samples": len(solution.t),
                "max_abs_work_energy_residual_m2_s2": float(np.max(np.abs(diagnostics["work_energy_residual_m2_s2"]))),
                "max_abs_relative_work_energy_residual": float(np.max(np.abs(diagnostics["relative_work_energy_residual"]))),
                "final_specific_energy_change_m2_s2": float(diagnostics["specific_energy_change_m2_s2"][-1]),
                "final_specific_work_m2_s2": float(solution.work_m2_s2[-1]),
            }
        diagnostics = levels["nominal"]["diagnostics"]
        case_reports[str(multiplier)] = {
            "density_multiplier": multiplier, "solver_levels": summaries,
            "altitude_bounds_m": [float(np.min(diagnostics["altitude_m"])), float(np.max(diagnostics["altitude_m"]))],
            "density_bounds_kg_m3": [float(np.min(diagnostics["density_kg_m3"])), float(np.max(diagnostics["density_kg_m3"]))],
            "drag_acceleration_bounds_m_s2": [float(np.min(diagnostics["drag_acceleration_m_s2"])),
                                             float(np.max(diagnostics["drag_acceleration_m_s2"]))],
            "checkpoints_elapsed_s": checkpoints, "convergence": convergence,
        }
        cases[multiplier] = {"model": case_model, "levels": levels,
                             "rtn_times_s": nominal.t[shared], "rtn_m": rtn,
                             "windows_s": available_windows, "mean_a_m": means,
                             "convergence": convergence_arrays}
    report = {
        "study": "Illustrative density sensitivity; not a validated space-weather forecast.",
        "versions": {"python": platform.python_version(), "numpy": np.__version__,
                     "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
        "model": {"mu_m3_s2": orbit.MU, "j2": orbit.J2,
                  "j2_reference_radius_m": orbit.EQUATORIAL_RADIUS, "spherical_radius_m": orbit.EARTH_RADIUS,
                  "initial_state_si": state0.tolist(), "inclination_rad": float(EXAMPLE_INCLINATION),
                  "requested_duration_s": duration, "spacecraft_and_atmosphere": asdict(model)},
        "solver": {"method": "DOP853", "output_interval_s": OUTPUT_INTERVAL,
                   "max_internal_step_s": orbit.PERIOD / 16.0,
                   "work_method": "Seventh ODE state, W'=inertial v dot a_drag; W(0)=0.",
                   "relative_work_residual_normalization": "Absolute initial total central + J2 specific energy."},
        "j2_only": {"end_time_s": float(baseline.t[-1]), "output_samples": len(baseline.t),
                    "stopped_at_altitude_limit": baseline.terminated_early,
                    "rtol": orbit.RTOL, "cartesian_atol_si": orbit.ATOL.tolist()},
        "zero_density_check": {"same_output_times": bool(np.array_equal(zero.t, baseline.t)),
                               "max_abs_state_difference_si": float(np.max(np.abs(zero.y - baseline.y)))},
        "averaging": {"description": "Time averages over common complete J2-only ascending-node periods; partial orbits omitted.",
                      "quadrature": "32-point Gauss-Legendre on dense Cartesian output, per period.",
                      "checkpoint_note": "Mean-a differences use the last complete reference orbit before each checkpoint, not the instantaneous time."},
        "cases": case_reports,
    }
    return {"baseline": baseline, "baseline_windows_s": windows, "baseline_mean_a_m": baseline_means,
            "cases": cases, "report": report}


def plot_experiments(data, output_dir, show=False):
    """Save decay, atmosphere, RTN differences and numerical verification."""
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {0.5: "tab:blue", 1.0: "tab:orange", 2.0: "tab:green"}
    baseline = data["baseline"]
    fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
    axes[0].plot(baseline.t, np.linalg.norm(baseline.y[:3], axis=0) - orbit.EARTH_RADIUS,
                 color="black", lw=0.8, alpha=0.65, label="J2 only")
    centers = np.mean(data["baseline_windows_s"], axis=1)
    axes[1].plot(centers, data["baseline_mean_a_m"], color="black", label="J2 only")
    for multiplier, case in data["cases"].items():
        solution = case["levels"]["nominal"]["solution"]
        altitude = case["levels"]["nominal"]["diagnostics"]["altitude_m"]
        label = f"J2 + drag: density x{multiplier:g}"
        axes[0].plot(solution.t, altitude, color=colors[multiplier], lw=0.65, alpha=0.7, label=label)
        axes[1].plot(np.mean(case["windows_s"], axis=1), case["mean_a_m"],
                     color=colors[multiplier], label=label)
    axes[0].set(title="Instantaneous spherical altitude: J2 oscillations and decay", ylabel="Altitude [m]")
    axes[1].set(title="Time-averaged osculating a over complete reference nodal periods",
                ylabel="Orbit-averaged semimajor axis [m]", xlabel="Elapsed time [s]")
    for ax in axes:
        ax.ticklabel_format(axis="y", style="plain", useOffset=False)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.25)
    fig.suptitle("Illustrative density sensitivity: identical spacecraft and initial state")
    fig.tight_layout()
    fig.savefig(output_dir / "decay.png", dpi=150)

    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    for multiplier, case in data["cases"].items():
        solution = case["levels"]["nominal"]["solution"]
        diagnostics = case["levels"]["nominal"]["diagnostics"]
        for ax, key in zip(axes, ("density_kg_m3", "drag_acceleration_m_s2")):
            ax.plot(solution.t, diagnostics[key], color=colors[multiplier], lw=0.7,
                    label=f"Density x{multiplier:g}")
    axes[0].set(title="Exponential density encountered along each trajectory", ylabel="Density [kg/m³]")
    axes[1].set(title="Acceleration opposing atmosphere-relative velocity", ylabel="Drag magnitude [m/s²]",
                xlabel="Elapsed time [s]")
    for ax in axes:
        ax.legend(fontsize=8)
        ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_dir / "atmosphere.png", dpi=150)

    fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    for multiplier, case in data["cases"].items():
        for ax, values in zip(axes, case["rtn_m"]):
            ax.plot(case["rtn_times_s"], values, color=colors[multiplier], lw=0.8,
                    label=f"Density x{multiplier:g}")
    for ax, name in zip(axes, ("Radial R", "Transverse/along-track T", "Cross-track N")):
        ax.set(ylabel=f"{name} difference [m]")
        ax.axhline(0.0, color="black", lw=0.5)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.25)
    axes[0].set_title("Position minus J2-only position, projected in the J2-only RTN frame")
    axes[-1].set_xlabel("Elapsed time [s]")
    fig.tight_layout()
    fig.savefig(output_dir / "rtn.png", dpi=150)

    case = data["cases"][1.0]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    nominal = case["levels"]["nominal"]
    axes[0, 0].plot(nominal["solution"].t, nominal["diagnostics"]["specific_energy_change_m2_s2"],
                    label="E(t) - E(0)")
    axes[0, 0].plot(nominal["solution"].t, nominal["solution"].work_m2_s2, "k--", label="Integrated drag work")
    for name, color in (("nominal", "tab:orange"), ("tighter", "tab:green")):
        level = case["levels"][name]
        axes[0, 1].plot(level["solution"].t, level["diagnostics"]["work_energy_residual_m2_s2"],
                        color=color, lw=0.8, label=name)
        convergence = case["convergence"][name]
        for ax, key in ((axes[1, 0], "position_m"), (axes[1, 1], "velocity_m_s")):
            ax.semilogy(convergence["time_s"][1:], convergence[key][1:], color=color, label=name)
    titles = ("Energy change and work [m²/s²]", "E(t) - E(0) - W(t) [m²/s²]",
              "Position difference from numerical reference [m]", "Velocity difference from reference [m/s]")
    for ax, title in zip(axes.flat, titles):
        ax.set(title=title, xlabel="Elapsed time [s]")
        ax.legend(fontsize=8)
        ax.grid(alpha=0.25)
    fig.suptitle("Numerical verification, density x1: matching J2 potential and inertial work")
    fig.tight_layout()
    fig.savefig(output_dir / "validation.png", dpi=150)
    if show:
        plt.show()
    plt.close("all")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--days", type=float, default=7.0)
    parser.add_argument("--mass-kg", type=float, default=1000.0)
    parser.add_argument("--cd", type=float, default=2.2)
    parser.add_argument("--area-m2", type=float, default=10.0)
    parser.add_argument("--rho-ref-kg-m3", type=float, default=4e-12)
    parser.add_argument("--h-ref-m", type=float, default=400_000.0)
    parser.add_argument("--scale-height-m", type=float, default=60_000.0)
    parser.add_argument("--min-altitude-m", type=float, default=200_000.0)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/drag"))
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()
    if not np.isfinite(args.days) or args.days <= 0.0:
        parser.error("--days must be finite and positive")
    try:
        model = DragModel(mass_kg=args.mass_kg, cd=args.cd, area_m2=args.area_m2,
                          rho_ref_kg_m3=args.rho_ref_kg_m3, h_ref_m=args.h_ref_m,
                          scale_height_m=args.scale_height_m, min_altitude_m=args.min_altitude_m)
        data = run_experiments(args.days * 86_400.0, model)
    except ValueError as error:
        parser.error(str(error))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plot_experiments(data, args.output_dir, args.show)
    report = data["report"]
    (args.output_dir / "validation.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(report["study"])
    print(f"Zero-density maximum state difference: {report['zero_density_check']['max_abs_state_difference_si']:.3e}")
    for multiplier, case in report["cases"].items():
        nominal = case["solver_levels"]["nominal"]
        print(f"\nDensity x{multiplier}: end {nominal['end_time_s']:.1f} s; "
              f"altitude guard triggered: {nominal['stopped_at_altitude_limit']}")
        for time, checkpoint in case["checkpoints_elapsed_s"].items():
            if checkpoint["available"]:
                print(f"  t={time} s: RTN {checkpoint['position_difference_rtn_m']} m; "
                      f"mean-a difference {checkpoint['orbit_averaged_semimajor_axis_difference_m']} m")
            else:
                print(f"  t={time} s: unavailable ({checkpoint['reason']})")
        for name in ("nominal", "tighter", "reference"):
            print(f"  {name} max work-energy residual: "
                  f"{case['solver_levels'][name]['max_abs_work_energy_residual_m2_s2']:.6e} m²/s²")
        print(f"  Position convergence improvement: {case['convergence']['position_improvement_factor']}")
    print(f"Figures and report saved in: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
