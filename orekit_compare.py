"""Milestone 4: force matching and independent SciPy/Orekit comparisons."""

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import platform

import matplotlib
import numpy as np
import scipy

import orbit
from atmosphere import DragModel
from drag import CHECKPOINTS, WORK_ATOL, rtn_difference
from j2 import EXAMPLE_INCLINATION


LEVELS = {"nominal": 1.0, "tighter": 100.0, "reference": 1000.0}


def smoke_test():
    """Import/start the JVM and propagate one central orbit for 60 s first."""
    from orekit_backend import OrekitModel, versions
    state0 = orbit.initial_state(EXAMPLE_INCLINATION)
    solution = OrekitModel().propagate(state0, 60.0, sample_interval=10.0)
    exact = orbit.analytical_state(solution.t, EXAMPLE_INCLINATION)
    position_error = np.linalg.norm(solution.y[:3] - exact[:3], axis=0)
    velocity_error = np.linalg.norm(solution.y[3:] - exact[3:], axis=0)
    if np.max(position_error) > 0.01 or np.max(velocity_error) > 1e-5:
        raise RuntimeError("Orekit smoke test does not match the circular reference.")
    return {"versions": versions(), "duration_s": 60.0, "force_models": solution.force_names,
            "max_position_error_vs_circular_m": float(np.max(position_error)),
            "max_velocity_error_vs_circular_m_s": float(np.max(velocity_error)),
            "final_mass_kg": solution.final_mass_kg}


def check_force_models():
    """Check forces/environment at shared states and dates BEFORE integrating.

    Python functions are used only on the comparison side. Orekit evaluates
    NewtonianAttraction, HolmesFeatherstoneAttractionModel and DragForce.
    """
    from orekit_backend import OrekitModel
    state0 = orbit.initial_state(EXAMPLE_INCLINATION)
    states = {
        "ascending_equator": state0,
        "opposite_equator": -state0,
        "north_oblique": np.r_[(orbit.ORBIT_RADIUS - 10_000.0)
                               * np.array([0.5, 0.4, 0.75]) / np.linalg.norm([0.5, 0.4, 0.75]),
                               [-5000.0, 3500.0, 4000.0]],
        "south_oblique": np.r_[(orbit.ORBIT_RADIUS + 30_000.0)
                               * np.array([0.45, -0.7, -0.55]) / np.linalg.norm([0.45, -0.7, -0.55]),
                               [5000.0, 3000.0, -4000.0]],
        "near_pole": np.array([1000.0, 0.0, np.sqrt(orbit.ORBIT_RADIUS**2 - 1000.0**2),
                               7400.0, 1000.0, 0.0]),
    }
    limits = {"central_acceleration_m_s2": 2e-12, "j2_acceleration_m_s2": 2e-13,
              "drag_acceleration_m_s2": 1e-16, "altitude_m": 1e-6,
              "relative_density": 1e-12, "atmospheric_velocity_m_s": 1e-9}
    maxima = dict.fromkeys(limits, 0.0)
    records = []
    for multiplier in (0.5, 1.0, 2.0):
        config = DragModel(density_multiplier=multiplier)
        model = OrekitModel(orbit.J2, config)
        for label, state in states.items():
            for time in (0.0, 86_400.0, 604_800.0):
                position, velocity = state[:3], state[3:]
                altitude = np.linalg.norm(position) - orbit.EARTH_RADIUS
                local = {"central": orbit.two_body_rhs(time, state)[3:],
                         "j2": orbit.j2_acceleration(position),
                         "drag": config.acceleration(position, velocity)}
                other = model.accelerations(state, time)
                environment = model.environment(position, time)
                density = float(config.density(altitude))
                atmospheric_velocity = velocity - config.relative_velocity(position, velocity)
                errors = {f"{name}_acceleration_m_s2": float(np.linalg.norm(other[name] - values))
                          for name, values in local.items()}
                errors.update({"altitude_m": abs(environment["altitude_m"] - altitude),
                               "relative_density": abs(environment["density_kg_m3"] - density) / density,
                               "atmospheric_velocity_m_s": float(np.linalg.norm(
                                   environment["atmospheric_velocity_m_s"] - atmospheric_velocity))})
                for name, error in errors.items():
                    if not np.isfinite(error) or error > limits[name]:
                        raise RuntimeError(f"Force matching failed at {label}, t={time}, "
                                           f"density x{multiplier}: {name} = {error}. Investigate before integration.")
                    maxima[name] = max(maxima[name], error)
                records.append({"state": label, "elapsed_time_s": time, "density_multiplier": multiplier,
                                "differences": errors, "local_density_kg_m3": density,
                                "orekit_density_kg_m3": environment["density_kg_m3"],
                                "local_atmospheric_velocity_m_s": atmospheric_velocity.tolist(),
                                "orekit_atmospheric_velocity_m_s": environment["atmospheric_velocity_m_s"].tolist(),
                                "local_accelerations_m_s2": {name: value.tolist() for name, value in local.items()},
                                "orekit_accelerations_m_s2": {name: value.tolist() for name, value in other.items()}})
    return {"passed": True, "states_si": {name: state.tolist() for name, state in states.items()},
            "limits": limits, "max_differences": maxima, "records": records}


def differences(reference, other, times):
    """Other minus reference at identical times; RTN axes follow reference."""
    ref_states, states = reference.sol(times), other.sol(times)
    delta = states - ref_states
    return {"time_s": times, "position_m": np.linalg.norm(delta[:3], axis=0),
            "velocity_m_s": np.linalg.norm(delta[3:], axis=0),
            "rtn_m": rtn_difference(ref_states, states)}


def summarize_difference(reference, other, times):
    values = differences(reference, other, times)
    summary = {"max_position_difference_m": float(np.max(values["position_m"])),
               "max_velocity_difference_m_s": float(np.max(values["velocity_m_s"]))}
    checkpoints = {}
    for time in CHECKPOINTS:
        if time > times[-1]:
            checkpoints[str(int(time))] = {"available": False}
        else:
            point = differences(reference, other, np.array([time]))
            checkpoints[str(int(time))] = {"available": True,
                                          "position_difference_m": float(point["position_m"][0]),
                                          "velocity_difference_m_s": float(point["velocity_m_s"][0]),
                                          "position_difference_rtn_m": point["rtn_m"][:, 0].tolist()}
    summary["checkpoints_elapsed_s"] = checkpoints
    return summary, values


def run_case(name, j2, config, duration, interval):
    """Tighten SciPy and Orekit independently, then compare both tight/reference."""
    from orekit_backend import MASS_ATOL, OrekitModel
    state0 = orbit.initial_state(EXAMPLE_INCLINATION)
    minimum = DragModel().min_altitude_m if config is None else config.min_altitude_m
    settings = {"state0": state0, "duration": duration, "sample_interval": interval,
                "j2": j2, "drag": config, "stop_altitude": minimum, "dense_output": True}
    local, other = {}, {}
    for level, scale in LEVELS.items():
        local[level] = orbit.propagate(**settings, rtol=orbit.RTOL / scale,
                                       atol=orbit.ATOL / scale, work_atol=WORK_ATOL / scale)
        other[level] = OrekitModel(j2, config).propagate(
            state0, duration, interval, rtol=orbit.RTOL / scale, atol=orbit.ATOL / scale,
            mass_atol=MASS_ATOL / scale, stop_altitude=minimum)
    common_end = min(solution.t[-1] for solution in (*local.values(), *other.values()))
    times = np.append(np.arange(0.0, common_end, interval), common_end)
    pairs = {
        "both_nominal": (local["nominal"], other["nominal"]),
        "only_scipy_tighter": (local["tighter"], other["nominal"]),
        "only_orekit_tighter": (local["nominal"], other["tighter"]),
        "both_tighter": (local["tighter"], other["tighter"]),
        "both_reference": (local["reference"], other["reference"]),
        "scipy_nominal_vs_reference": (local["reference"], local["nominal"]),
        "scipy_tighter_vs_reference": (local["reference"], local["tighter"]),
        "orekit_nominal_vs_reference": (other["reference"], other["nominal"]),
        "orekit_tighter_vs_reference": (other["reference"], other["tighter"]),
    }
    summaries, arrays = {}, {}
    for label, pair in pairs.items():
        summaries[label], arrays[label] = summarize_difference(*pair, times)
    solver_reports = {}
    for implementation, solutions in (("scipy", local), ("orekit", other)):
        solver_reports[implementation] = {}
        for level, solution in solutions.items():
            scale = LEVELS[level]
            solver_reports[implementation][level] = {
                "rtol": orbit.RTOL / scale, "cartesian_atol_si": (orbit.ATOL / scale).tolist(),
                "end_time_s": float(solution.t[-1]), "output_samples": len(solution.t),
                "stopped_at_altitude_limit": solution.terminated_early,
                "rhs_evaluations": solution.nfev if implementation == "scipy" else solution.rhs_evaluations,
            }
            if implementation == "scipy" and config is not None:
                solver_reports[implementation][level]["work_atol_m2_s2"] = WORK_ATOL / scale
            if implementation == "orekit":
                solver_reports[implementation][level].update({"mass_atol_kg": MASS_ATOL / scale,
                                                              "force_models": solution.force_names,
                                                              "final_mass_kg": solution.final_mass_kg})
    report = {"name": name, "j2": j2, "drag": asdict(config) if config is not None else None,
              "requested_duration_s": duration, "output_interval_s": interval,
              "stop_altitude_m": minimum, "common_comparison_end_s": float(common_end),
              "difference_convention": "Orekit minus SciPy for cross-library pairs; RTN frame follows the stated SciPy level.",
              "solver_runs": solver_reports, "comparisons": summaries}
    if j2 == 0.0 and config is None:
        exact = orbit.analytical_state(times, EXAMPLE_INCLINATION)
        report["circular_reference"] = {
            implementation: {
                level: {"max_position_error_m": float(np.max(np.linalg.norm(solution.sol(times)[:3] - exact[:3], axis=0))),
                        "max_velocity_error_m_s": float(np.max(np.linalg.norm(solution.sol(times)[3:] - exact[3:], axis=0)))}
                for level, solution in solutions.items()}
            for implementation, solutions in (("scipy", local), ("orekit", other))}
    return {"local": local, "orekit": other, "arrays": arrays, "report": report}


def run_comparisons(progress=print):
    """Ordered comparisons; density factors run only after nominal agreement."""
    from orekit_backend import EPOCH_LABEL, EVENT_THRESHOLD, MAX_STEP, MIN_STEP
    smoke = smoke_test()
    progress("Orekit import/60-second central smoke test passed.")
    preflight = check_force_models()
    progress(f"Force/environment preflight passed: {preflight['max_differences']}")
    cases = {}
    definitions = (("central", 0.0, None, 86_400.0, 10.0),
                   ("j2", orbit.J2, None, 604_800.0, 60.0),
                   ("drag_1", orbit.J2, DragModel(), 604_800.0, 60.0),
                   ("drag_0.5", orbit.J2, DragModel(density_multiplier=0.5), 604_800.0, 60.0),
                   ("drag_2", orbit.J2, DragModel(density_multiplier=2.0), 604_800.0, 60.0))
    for definition in definitions:
        name = definition[0]
        progress(f"Running {name}: three independent tolerance levels per propagator ...")
        case = run_case(*definition)
        match = case["report"]["comparisons"]["both_reference"]
        if match["max_position_difference_m"] > 0.05 or match["max_velocity_difference_m_s"] > 5e-5:
            raise RuntimeError(f"{name} reference-level discrepancy needs investigation before subsequent cases: {match}")
        cases[name] = case
        progress(f"{name}: both-reference max position difference {match['max_position_difference_m']:.6e} m")
    report = {
        "purpose": "Independent implementation comparison for matched ideal models; Orekit is not exact truth or real-world validation.",
        "versions": {"python": platform.python_version(), "numpy": np.__version__,
                     "scipy": scipy.__version__, "matplotlib": matplotlib.__version__, **smoke["versions"]},
        "matching": {"inertial_frame": "GCRF, same Cartesian components interpreted in both implementations.",
                     "epoch": EPOCH_LABEL, "time_definition": "Elapsed SI seconds from the common TAI epoch.",
                     "initial_state_si": orbit.initial_state(EXAMPLE_INCLINATION).tolist(),
                     "inclination_rad": float(EXAMPLE_INCLINATION), "mu_m3_s2": orbit.MU,
                     "j2": orbit.J2, "gravity_reference_radius_m": orbit.EQUATORIAL_RADIUS,
                     "normalized_c20": -orbit.J2 / float(np.sqrt(5.0)),
                     "coefficient_convention": "Fully normalized Pbar20=sqrt(5)*P2; Cbar20=-J2/sqrt(5). Degree 2/order 0 only; all other non-central coefficients zero.",
                     "central_gravity": "One explicit NewtonianAttraction; automatic addition disabled. Holmes-Featherstone is non-central only.",
                     "spherical_altitude_radius_m": orbit.EARTH_RADIUS,
                     "orekit_body_shape": "OneAxisEllipsoid with flattening exactly zero (sphere). No ellipsoidal substitution.",
                     "atmosphere_frame": "Uniform passive rotation about GCRF +z; axes coincide at epoch. No ITRF, EOP, precession/nutation or winds.",
                     "spacecraft_and_nominal_atmosphere": asdict(DragModel()),
                     "termination": "Stop on descending spherical-altitude crossing at 200000 m; reject initial state at/below limit. Compare only common available times.",
                     "orekit_event_threshold_s": EVENT_THRESHOLD,
                     "orekit_min_step_s": MIN_STEP, "orekit_max_step_s": MAX_STEP,
                     "scipy_guard_max_step_s": MAX_STEP,
                     "data_files": "None needed: built-in GCRF/TAI, explicitly supplied constants, ideal rotation and toy atmosphere.",
                     "rtol_levels": {name: orbit.RTOL / scale for name, scale in LEVELS.items()},
                     "integrators": "SciPy DOP853 and Hipparchus DormandPrince853Integrator; independent implementations of the same RK family."},
        "smoke_test": smoke, "force_preflight": preflight,
        "agreement_gate": {"both_reference_max_position_m": 0.05, "both_reference_max_velocity_m_s": 5e-5,
                           "meaning": "Consistency gate before subsequent cases; not a real-world accuracy bound."},
        "cases": {name: case["report"] for name, case in cases.items()},
    }
    return {"cases": cases, "report": report}


def plot_comparisons(data, output_dir, show=False):
    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for name, case in data["cases"].items():
        fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
        for label, color in (("both_nominal", "tab:orange"), ("both_tighter", "tab:green"),
                              ("both_reference", "tab:blue")):
            values = case["arrays"][label]
            axes[0, 0].plot(values["time_s"], values["position_m"], color=color, lw=0.8, label=label)
            for ax, component in zip((axes[0, 1], axes[1, 0], axes[1, 1]), values["rtn_m"]):
                ax.plot(values["time_s"], component, color=color, lw=0.8, label=label)
        for ax, title in zip(axes.flat, ("Cartesian position difference magnitude [m]", "Radial R difference [m]",
                                         "Transverse/along-track T difference [m]", "Cross-track N difference [m]")):
            ax.set(title=title, xlabel="Elapsed time [s]")
            ax.legend(fontsize=8)
            ax.grid(alpha=0.25)
        fig.suptitle(f"{name}: Orekit minus SciPy; each RTN frame follows the matching SciPy level")
        fig.tight_layout()
        fig.savefig(output_dir / f"{name}.png", dpi=150)

    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    case = data["cases"]["drag_1"]
    for implementation, color in (("scipy", "tab:blue"), ("orekit", "tab:orange")):
        for level, style in (("nominal", "-"), ("tighter", "--")):
            values = case["arrays"][f"{implementation}_{level}_vs_reference"]
            for ax, key in zip(axes, ("position_m", "velocity_m_s")):
                ax.semilogy(values["time_s"][1:], values[key][1:], color=color, ls=style,
                            label=f"{implementation}: {level} vs own reference")
    axes[0].set(title="Independent convergence for nominal-density drag", ylabel="Position difference [m]")
    axes[1].set(ylabel="Velocity difference [m/s]", xlabel="Elapsed time [s]")
    for ax in axes:
        ax.legend(fontsize=8)
        ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_dir / "convergence.png", dpi=150)
    if show:
        plt.show()
    plt.close("all")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smoke-only", action="store_true")
    parser.add_argument("--forces-only", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/orekit"))
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()
    try:
        if args.smoke_only:
            print(json.dumps(smoke_test(), indent=2, allow_nan=False))
            return
        if args.forces_only:
            result = check_force_models()
            print(json.dumps({key: value for key, value in result.items() if key != "records"}, indent=2))
            return
        data = run_comparisons(lambda message: print(message, flush=True))
    except ModuleNotFoundError as error:
        parser.error(f"Optional Orekit environment is required: {error}. See README setup and requirements-orekit.txt.")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plot_comparisons(data, args.output_dir, args.show)
    (args.output_dir / "validation.json").write_text(
        json.dumps(data["report"], indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Figures and report saved in: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
