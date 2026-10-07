"""Drag physics, stopping events, RTN geometry and seven-day work/convergence."""

from dataclasses import replace
from types import SimpleNamespace
import unittest

import numpy as np
from scipy.integrate import quad

import drag
from atmosphere import DragModel, EARTH_ROTATION_RATE
import orbit


class DragForceTests(unittest.TestCase):
    def setUp(self):
        self.model = DragModel()
        state = orbit.initial_state(drag.EXAMPLE_INCLINATION)
        self.position, self.velocity = state[:3], state[3:]

    def test_exponential_density_and_reference_value(self):
        model = self.model
        heights = model.h_ref_m + model.scale_height_m * np.array([-1.0, 0.0, 1.0])
        np.testing.assert_allclose(model.density(heights), model.rho_ref_kg_m3 * np.exp([1.0, 0.0, -1.0]))
        np.testing.assert_allclose(replace(model, density_multiplier=2.0).density(heights),
                                   2.0 * model.density(heights), rtol=1e-15)

    def test_fixed_state_scaling(self):
        acceleration = self.model.acceleration(self.position, self.velocity)
        changes = ({"density_multiplier": 2.0}, {"rho_ref_kg_m3": 8e-12},
                   {"cd": 4.4}, {"area_m2": 20.0}, {"mass_kg": 500.0})
        for change in changes:
            with self.subTest(change=change):
                scaled = replace(self.model, **change).acceleration(self.position, self.velocity)
                np.testing.assert_allclose(scaled, 2.0 * acceleration, rtol=1e-15, atol=0.0)

    def test_drag_opposes_relative_velocity_and_has_expected_magnitude(self):
        relative = self.model.relative_velocity(self.position, self.velocity)
        acceleration = self.model.acceleration(self.position, self.velocity)
        self.assertLess(np.dot(acceleration, relative), 0.0)
        np.testing.assert_allclose(np.cross(acceleration, relative), 0.0, atol=1e-17)
        expected = (0.5 * self.model.rho_ref_kg_m3 * self.model.cd
                    * self.model.area_m2 / self.model.mass_kg * np.dot(relative, relative))
        self.assertAlmostEqual(np.linalg.norm(acceleration), expected, places=18)
        # Co-rotation moves +x atmosphere toward +y, not toward +z.
        np.testing.assert_allclose(relative, self.velocity - np.array([
            0.0, EARTH_ROTATION_RATE * orbit.ORBIT_RADIUS, 0.0]), rtol=1e-15)

    def test_corotating_spacecraft_has_zero_drag(self):
        atmospheric_velocity = np.cross([0.0, 0.0, EARTH_ROTATION_RATE], self.position)
        np.testing.assert_array_equal(self.model.acceleration(self.position, atmospheric_velocity), 0.0)

    def test_vectorized_force_matches_individual_states(self):
        states = orbit.analytical_state([0.0, 1200.0, 3400.0], drag.EXAMPLE_INCLINATION)
        vectorized = self.model.acceleration(states[:3], states[3:])
        for index in range(states.shape[1]):
            np.testing.assert_allclose(vectorized[:, index],
                                       self.model.acceleration(states[:3, index], states[3:, index]),
                                       rtol=2e-15, atol=1e-20)

    def test_invalid_parameters_and_disabled_forces(self):
        for change in ({"mass_kg": 0.0}, {"scale_height_m": -1.0}, {"cd": -1.0},
                       {"area_m2": -1.0}, {"rho_ref_kg_m3": -1.0},
                       {"density_multiplier": -1.0}, {"min_altitude_m": np.nan}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                replace(self.model, **change)
        for change in ({"rho_ref_kg_m3": 0.0}, {"density_multiplier": 0.0},
                       {"cd": 0.0}, {"area_m2": 0.0}):
            disabled = replace(self.model, **change)
            np.testing.assert_array_equal(disabled.acceleration(self.position, self.velocity), 0.0)

    def test_zero_density_exactly_reproduces_j2(self):
        settings = {"state0": np.concatenate((self.position, self.velocity)), "duration": 86_400.0,
                    "j2": orbit.J2, "stop_altitude": self.model.min_altitude_m}
        baseline = orbit.propagate(**settings)
        for change in ({"rho_ref_kg_m3": 0.0}, {"density_multiplier": 0.0}):
            disabled = orbit.propagate(**settings, drag=replace(self.model, **change))
            np.testing.assert_array_equal(disabled.t, baseline.t)
            np.testing.assert_array_equal(disabled.y, baseline.y)
            np.testing.assert_array_equal(disabled.work_m2_s2, 0.0)


class DragGuardAndGeometryTests(unittest.TestCase):
    def test_stopping_event_returns_the_actual_boundary_state(self):
        model = DragModel(min_altitude_m=395_000.0)
        settings = {"state0": orbit.initial_state(drag.EXAMPLE_INCLINATION), "duration": 20_000.0,
                    "sample_interval": 60.0, "j2": orbit.J2, "drag": model, "dense_output": True}
        solution = orbit.propagate(**settings)
        self.assertTrue(solution.success)
        self.assertTrue(solution.terminated_early)
        self.assertLess(solution.t[-1], settings["duration"])
        self.assertEqual(solution.t[-1], solution.t_events[0][-1])
        self.assertEqual(solution.y.shape, (6, len(solution.t)))
        self.assertEqual(len(solution.work_m2_s2), len(solution.t))
        self.assertLess(abs(np.linalg.norm(solution.y[:3, -1]) - orbit.EARTH_RADIUS
                            - model.min_altitude_m), 1e-6)
        np.testing.assert_array_equal(solution.y[:, -1], solution.y_events[0][-1])
        np.testing.assert_allclose(solution.sol(solution.t[-1]), solution.y[:, -1], atol=1e-8)
        self.assertAlmostEqual(float(solution.work_sol(solution.t[-1])), solution.work_m2_s2[-1])
        tighter = orbit.propagate(**settings, rtol=orbit.RTOL / 100.0, atol=orbit.ATOL / 100.0,
                                  work_atol=drag.WORK_ATOL / 100.0)
        self.assertLess(abs(tighter.t[-1] - solution.t[-1]), 1e-3)

    def test_invalid_initial_altitude_is_rejected(self):
        for minimum in (400_000.0, 410_000.0):
            with self.assertRaises(ValueError):
                orbit.propagate(drag=DragModel(min_altitude_m=minimum))

    def test_rtn_axes_have_correct_orientation_and_preserve_distance(self):
        state = orbit.initial_state(drag.EXAMPLE_INCLINATION)
        inclination = drag.EXAMPLE_INCLINATION
        axes = np.array([[1.0, 0.0, 0.0],
                         [0.0, np.cos(inclination), np.sin(inclination)],
                         [0.0, -np.sin(inclination), np.cos(inclination)]])
        expected = np.array([10.0, -20.0, 30.0])
        displaced = state.copy()
        displaced[:3] += axes.T @ expected
        measured = drag.rtn_difference(state, displaced)[:, 0]
        np.testing.assert_allclose(measured, expected, atol=1e-9)
        self.assertAlmostEqual(np.linalg.norm(measured), np.linalg.norm(displaced[:3] - state[:3]))
        with self.assertRaises(ValueError):
            drag.rtn_difference(np.array([1.0, 0.0, 0.0, 1.0, 0.0, 0.0]), displaced)

    def test_orbit_averaging_recovers_constant_circular_semimajor_axis(self):
        solution = SimpleNamespace(sol=orbit.analytical_state)
        windows = np.array([[0.0, orbit.PERIOD], [orbit.PERIOD, 2.0 * orbit.PERIOD]])
        np.testing.assert_allclose(drag.average_semimajor_axis(solution, windows),
                                   orbit.ORBIT_RADIUS, rtol=0.0, atol=1e-8)
        equatorial = orbit.propagate(duration=2.0 * orbit.PERIOD, dense_output=True)
        with self.assertRaises(ValueError):
            drag.nodal_windows(equatorial)

    def test_work_uses_inertial_velocity_and_agrees_with_independent_quadrature(self):
        model = DragModel()
        solution = orbit.propagate(state0=orbit.initial_state(drag.EXAMPLE_INCLINATION),
                                   duration=orbit.PERIOD, j2=orbit.J2, drag=model, dense_output=True,
                                   rtol=orbit.RTOL / 100.0, atol=orbit.ATOL / 100.0)

        def power(time, relative=False):
            state = solution.sol(time)
            velocity = model.relative_velocity(state[:3], state[3:]) if relative else state[3:]
            return np.dot(velocity, model.acceleration(state[:3], state[3:]))

        independent_work = quad(power, 0.0, orbit.PERIOD, epsabs=1e-7, epsrel=1e-10)[0]
        wrong_frame_work = quad(lambda t: power(t, relative=True), 0.0, orbit.PERIOD,
                                epsabs=1e-7, epsrel=1e-10)[0]
        self.assertLess(abs(solution.work_m2_s2[-1] - independent_work), 1e-6)
        self.assertGreater(abs(solution.work_m2_s2[-1] - wrong_frame_work), abs(independent_work) * 0.01)


class SevenDayDragTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = drag.run_experiments()
        cls.report = cls.data["report"]

    def test_shared_state_full_duration_and_zero_density(self):
        initial = self.data["baseline"].y[:, 0]
        self.assertEqual(self.data["baseline"].t[-1], 604_800.0)
        for case in self.data["cases"].values():
            for level in case["levels"].values():
                solution = level["solution"]
                np.testing.assert_array_equal(solution.y[:, 0], initial)
                self.assertEqual(solution.t[-1], 604_800.0)
                self.assertFalse(solution.terminated_early)
        self.assertTrue(self.report["zero_density_check"]["same_output_times"])
        self.assertEqual(self.report["zero_density_check"]["max_abs_state_difference_si"], 0.0)

    def test_density_sensitivity_and_checkpoint_differences(self):
        for time in ("86400", "259200", "604800"):
            mean_differences = []
            for multiplier in drag.DENSITY_MULTIPLIERS:
                checkpoint = self.report["cases"][str(multiplier)]["checkpoints_elapsed_s"][time]
                self.assertTrue(checkpoint["available"])
                self.assertLess(checkpoint["orbit_averaged_semimajor_axis_difference_m"], 0.0)
                mean_differences.append(checkpoint["orbit_averaged_semimajor_axis_difference_m"])
                rtn = np.array(checkpoint["position_difference_rtn_m"])
                np.testing.assert_allclose(np.linalg.norm(rtn), checkpoint["position_difference_norm_m"], rtol=1e-14)
            self.assertTrue(np.all(np.diff(mean_differences) < 0.0))

    def test_work_energy_balance_converges_for_all_densities(self):
        for case in self.report["cases"].values():
            levels = case["solver_levels"]
            for name in ("nominal", "tighter", "reference"):
                self.assertLess(levels[name]["max_abs_relative_work_energy_residual"], 1e-8)
                self.assertLess(levels[name]["final_specific_work_m2_s2"], -1000.0)
                self.assertLess(levels[name]["final_specific_energy_change_m2_s2"], -1000.0)
            self.assertLess(levels["tighter"]["max_abs_work_energy_residual_m2_s2"],
                            levels["nominal"]["max_abs_work_energy_residual_m2_s2"] / 10.0)

    def test_tighter_trajectories_converge_for_all_densities(self):
        for case in self.report["cases"].values():
            convergence = case["convergence"]
            for quantity in ("position_difference_m", "velocity_difference_m_s"):
                self.assertLess(convergence[f"max_tighter_{quantity}"],
                                convergence[f"max_nominal_{quantity}"] / 10.0)

    def test_averaging_quadrature_is_resolved(self):
        baseline = self.data["baseline"]
        windows = self.data["baseline_windows_s"][[0, -1]]
        coarse = drag.average_semimajor_axis(baseline, windows, quadrature_order=32)
        fine = drag.average_semimajor_axis(baseline, windows, quadrature_order=64)
        np.testing.assert_allclose(coarse, fine, rtol=0.0, atol=0.01)


if __name__ == "__main__":
    unittest.main()
