"""Optional Orekit checks; skipped when the wrapper is absent from the environment."""

from dataclasses import replace
from importlib.util import find_spec
import unittest

import numpy as np

import orbit
from atmosphere import DragModel
from j2 import EXAMPLE_INCLINATION
import orekit_compare


HAS_OREKIT = find_spec("orekit_jpype") is not None


@unittest.skipUnless(HAS_OREKIT, "Install requirements-orekit.txt in the optional Orekit environment.")
class OrekitModelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from orekit_backend import OrekitModel
        cls.model_class = OrekitModel
        cls.state0 = orbit.initial_state(EXAMPLE_INCLINATION)

    def test_smoke_import_and_numerical_central_propagation(self):
        smoke = orekit_compare.smoke_test()
        self.assertEqual(smoke["force_models"], ["NewtonianAttraction"])
        self.assertEqual(smoke["final_mass_kg"], 1000.0)
        self.assertLess(smoke["max_position_error_vs_circular_m"], 0.01)

    def test_forces_density_and_atmospheric_velocity_match_before_integration(self):
        preflight = orekit_compare.check_force_models()
        self.assertTrue(preflight["passed"])
        self.assertEqual(len(preflight["records"]), 45)
        for name, value in preflight["max_differences"].items():
            self.assertLessEqual(value, preflight["limits"][name])

    def test_only_normalized_c20_and_one_central_model_are_enabled(self):
        model = self.model_class(orbit.J2)
        self.assertEqual(model.provider.getMaxDegree(), 2)
        self.assertEqual(model.provider.getMaxOrder(), 0)
        coefficient = model.provider.onDate(model.epoch).getNormalizedCnm(2, 0)
        self.assertEqual(coefficient, -orbit.J2 / np.sqrt(5.0))
        self.assertEqual(model.sphere.getFlattening(), 0.0)
        self.assertEqual(model.sphere.getEquatorialRadius(), orbit.EARTH_RADIUS)
        result = model.propagate(self.state0, 60.0)
        self.assertEqual(result.force_names.count("NewtonianAttraction"), 1)
        self.assertEqual(set(result.force_names), {"NewtonianAttraction", "HolmesFeatherstoneAttractionModel"})

    def test_zero_density_reproduces_orekit_j2_exactly(self):
        baseline = self.model_class(orbit.J2).propagate(self.state0, 3600.0, stop_altitude=200_000.0)
        zero = self.model_class(orbit.J2, DragModel(density_multiplier=0.0)).propagate(
            self.state0, 3600.0, stop_altitude=200_000.0)
        np.testing.assert_array_equal(zero.t, baseline.t)
        np.testing.assert_array_equal(zero.y, baseline.y)
        self.assertNotIn("DragForce", zero.force_names)

    def test_configuration_changes_match_the_independent_drag_force(self):
        from orekit_backend import vector
        for changes in ({"mass_kg": 500.0}, {"area_m2": 20.0}, {"cd": 4.4},
                        {"scale_height_m": 50_000.0}, {"omega_rad_s": 0.0}):
            with self.subTest(changes=changes):
                config = replace(DragModel(), **changes)
                model = self.model_class(orbit.J2, config)
                state = orbit.analytical_state([1000.0], EXAMPLE_INCLINATION)[:, 0]
                measured = model.accelerations(state, 12345.0)["drag"]
                np.testing.assert_allclose(measured, config.acceleration(state[:3], state[3:]),
                                           rtol=1e-12, atol=1e-20)
                if config.omega_rad_s == 0.0:
                    np.testing.assert_array_equal(model.environment(state[:3], 12345.0)["atmospheric_velocity_m_s"], 0.0)
        # Positive physical rotation also has the correct passive coordinate sign.
        model = self.model_class(orbit.J2, DragModel())
        time = 12345.0
        transformed = model.frame.getTransformTo(model.body_frame, model.epoch.shiftedBy(time)).transformPosition(
            vector([orbit.ORBIT_RADIUS, 0.0, 0.0]))
        angle = model.config.omega_rad_s * time
        np.testing.assert_allclose([transformed.getX(), transformed.getY(), transformed.getZ()],
                                   orbit.ORBIT_RADIUS * np.array([np.cos(angle), -np.sin(angle), 0.0]), atol=1e-8)

    def test_uniform_rotation_field_transform(self):
        from org.hipparchus.util import Binary64Field
        from org.orekit.time import FieldAbsoluteDate
        model = self.model_class(orbit.J2, DragModel())
        field_date = FieldAbsoluteDate(Binary64Field.getInstance(), model.epoch.shiftedBy(12345.0))
        transform = model.rotation.getTransform(field_date)
        self.assertAlmostEqual(transform.getRotationRate().getZ().getReal(), model.config.omega_rad_s)

    def test_spherical_altitude_event_matches_scipy_and_prevents_extrapolation(self):
        config = DragModel(min_altitude_m=395_000.0)
        other = self.model_class(orbit.J2, config).propagate(
            self.state0, 20_000.0, rtol=1e-12, atol=orbit.ATOL / 100.0, stop_altitude=config.min_altitude_m)
        local = orbit.propagate(state0=self.state0, duration=20_000.0, j2=orbit.J2, drag=config,
                                rtol=1e-12, atol=orbit.ATOL / 100.0, dense_output=True,
                                work_atol=1e-7)
        self.assertTrue(other.terminated_early)
        self.assertTrue(local.terminated_early)
        self.assertLess(abs(other.t[-1] - local.t[-1]), 1e-5)
        self.assertLess(abs(np.linalg.norm(other.y[:3, -1]) - orbit.EARTH_RADIUS - config.min_altitude_m), 1e-6)
        end = min(other.t[-1], local.t[-1])
        np.testing.assert_allclose(other.sol(end), local.sol(end), rtol=0.0, atol=0.001)
        with self.assertRaises(ValueError):
            other.sol(other.t[-1] + 1.0)
        with self.assertRaises(ValueError):
            self.model_class().propagate(self.state0, 60.0, stop_altitude=400_000.0)


@unittest.skipUnless(HAS_OREKIT, "Install requirements-orekit.txt in the optional Orekit environment.")
class OrekitTrajectoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = orekit_compare.run_comparisons(progress=lambda message: None)

    def test_initial_states_times_and_constant_mass_match(self):
        initial = orbit.initial_state(EXAMPLE_INCLINATION)
        for case in self.data["cases"].values():
            for level in orekit_compare.LEVELS:
                local, other = case["local"][level], case["orekit"][level]
                np.testing.assert_array_equal(local.y[:, 0], initial)
                np.testing.assert_array_equal(other.y[:, 0], initial)
                np.testing.assert_array_equal(local.t, other.t)
                self.assertFalse(other.terminated_early)
                self.assertEqual(other.final_mass_kg, 1000.0)

    def test_central_runs_agree_with_exact_circular_solution(self):
        references = self.data["cases"]["central"]["report"]["circular_reference"]
        for implementation in ("scipy", "orekit"):
            self.assertLess(references[implementation]["nominal"]["max_position_error_m"], 0.1)
            self.assertLess(references[implementation]["reference"]["max_position_error_m"], 1e-4)

    def test_independent_position_and_velocity_convergence(self):
        for case in self.data["cases"].values():
            comparisons = case["report"]["comparisons"]
            for implementation in ("scipy", "orekit"):
                nominal = comparisons[f"{implementation}_nominal_vs_reference"]
                tighter = comparisons[f"{implementation}_tighter_vs_reference"]
                for quantity in ("max_position_difference_m", "max_velocity_difference_m_s"):
                    self.assertLess(tighter[quantity], nominal[quantity] / 10.0)
            self.assertLess(comparisons["both_reference"]["max_position_difference_m"], 0.05)
            self.assertLess(comparisons["both_reference"]["max_velocity_difference_m_s"], 5e-5)

    def test_checkpoints_and_rtn_reference_are_consistent(self):
        for name, case in self.data["cases"].items():
            summary = case["report"]["comparisons"]["both_tighter"]
            for time, point in summary["checkpoints_elapsed_s"].items():
                self.assertEqual(point["available"], float(time) <= case["report"]["requested_duration_s"])
                if point["available"]:
                    np.testing.assert_allclose(np.linalg.norm(point["position_difference_rtn_m"]),
                                               point["position_difference_m"], rtol=1e-12, atol=1e-12)
                    delta = case["orekit"]["tighter"].sol(float(time)) - case["local"]["tighter"].sol(float(time))
                    self.assertAlmostEqual(np.linalg.norm(delta[:3]), point["position_difference_m"])


if __name__ == "__main__":
    unittest.main()
