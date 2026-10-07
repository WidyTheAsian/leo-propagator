"""J2 force/potential consistency, orbital-angle cases and seven-day verification."""

import unittest

import numpy as np

import j2
import orbit
from elements import osculating_elements


class J2ForceTests(unittest.TestCase):
    def test_j2_zero_reproduces_two_body(self):
        settings = {"state0": orbit.initial_state(j2.EXAMPLE_INCLINATION), "duration": 3600.0}
        central = orbit.propagate(**settings)
        disabled = orbit.propagate(**settings, j2=0.0)
        np.testing.assert_array_equal(disabled.t, central.t)
        np.testing.assert_array_equal(disabled.y, central.y)
        np.testing.assert_array_equal(orbit.j2_acceleration(central.y[:3], j2=0.0), 0.0)
        np.testing.assert_array_equal(orbit.specific_potential(central.y[:3], j2=0.0),
                                      -orbit.MU / np.linalg.norm(central.y[:3], axis=0))

    def test_acceleration_is_negative_potential_gradient(self):
        # Finite differences provide a check independent of the acceleration algebra.
        for position in (np.array([6.8e6, -1.2e6, 2.3e6]),
                         np.array([7.0e6, 0.0, 0.0]), np.array([0.0, 0.0, 7.0e6])):
            with self.subTest(position=position):
                gradient = np.zeros(3)
                step = 10.0  # metres: small relative to radius, above cancellation noise
                for axis in range(3):
                    offset = np.zeros(3)
                    offset[axis] = step
                    gradient[axis] = (orbit.specific_potential(position + offset, orbit.J2)
                                      - orbit.specific_potential(position - offset, orbit.J2)) / (2.0 * step)
                state = np.concatenate((position, np.zeros(3)))
                acceleration = orbit.two_body_rhs(0.0, state, orbit.J2)[3:]
                np.testing.assert_allclose(acceleration, -gradient, rtol=2e-8, atol=2e-9)

    def test_j2_direction_on_equator_and_pole(self):
        equator = orbit.j2_acceleration(np.array([orbit.ORBIT_RADIUS, 0.0, 0.0]))
        pole = orbit.j2_acceleration(np.array([0.0, 0.0, orbit.ORBIT_RADIUS]))
        self.assertLess(equator[0], 0.0)  # strengthens the equatorial inward acceleration
        self.assertGreater(pole[2], 0.0)  # weakens the polar inward acceleration
        np.testing.assert_array_equal(equator[1:], 0.0)
        np.testing.assert_array_equal(pole[:2], 0.0)


class OrbitalAngleTests(unittest.TestCase):
    def test_circular_inclined_and_polar_angles(self):
        for inclination in (j2.EXAMPLE_INCLINATION, np.pi / 2.0):
            elements = osculating_elements(orbit.initial_state(inclination))
            self.assertTrue(elements["circular"][0])
            self.assertFalse(elements["equatorial"][0])
            self.assertAlmostEqual(elements["raan_rad"][0], 0.0)
            self.assertAlmostEqual(elements["inclination_rad"][0], inclination)
            self.assertTrue(np.isnan(elements["argument_of_periapsis_rad"][0]))
            self.assertTrue(np.isnan(elements["true_anomaly_rad"][0]))
            self.assertAlmostEqual(elements["argument_of_latitude_rad"][0], 0.0)

    def test_circular_equatorial_angles_include_retrograde_limit(self):
        for inclination in (0.0, 1e-12, np.pi):
            elements = osculating_elements(orbit.initial_state(inclination))
            self.assertTrue(elements["equatorial"][0])
            for name in ("raan_rad", "argument_of_periapsis_rad", "true_anomaly_rad",
                         "argument_of_latitude_rad"):
                self.assertTrue(np.isnan(elements[name][0]), name)
            self.assertAlmostEqual(elements["true_longitude_rad"][0], 0.0)

    def test_eccentric_equatorial_orbit_has_true_anomaly_but_no_node(self):
        state = orbit.initial_state()
        state[4] *= np.sqrt(1.1)  # perigee state with osculating e = 0.1
        elements = osculating_elements(state)
        self.assertFalse(elements["circular"][0])
        self.assertTrue(np.isnan(elements["raan_rad"][0]))
        self.assertTrue(np.isnan(elements["argument_of_periapsis_rad"][0]))
        self.assertAlmostEqual(elements["true_anomaly_rad"][0], 0.0)

    def test_known_eccentric_inclined_orientation(self):
        raan, inclination, periapsis, eccentricity = 0.4, 0.8, 0.3, 0.1

        def rotate_z(angle):
            return np.array([[np.cos(angle), -np.sin(angle), 0.0],
                             [np.sin(angle), np.cos(angle), 0.0], [0.0, 0.0, 1.0]])

        rotate_x = np.array([[1.0, 0.0, 0.0], [0.0, np.cos(inclination), -np.sin(inclination)],
                             [0.0, np.sin(inclination), np.cos(inclination)]])
        rotation = rotate_z(raan) @ rotate_x @ rotate_z(periapsis)
        position = rotation @ np.array([orbit.ORBIT_RADIUS, 0.0, 0.0])
        velocity = rotation @ np.array([0.0, orbit.CIRCULAR_SPEED * np.sqrt(1.0 + eccentricity), 0.0])
        elements = osculating_elements(np.concatenate((position, velocity)))
        self.assertAlmostEqual(elements["eccentricity"][0], eccentricity)
        self.assertAlmostEqual(elements["semimajor_axis_m"][0] / orbit.ORBIT_RADIUS,
                               1.0 / (1.0 - eccentricity))
        for name, expected in (("raan_rad", raan), ("inclination_rad", inclination),
                               ("argument_of_periapsis_rad", periapsis), ("true_anomaly_rad", 0.0)):
            difference = (elements[name][0] - expected + np.pi) % (2.0 * np.pi) - np.pi
            self.assertAlmostEqual(difference, 0.0, places=12)

    def test_raan_unwrapping_and_undefined_angle_rejection(self):
        times = np.arange(4.0)
        for intercept, slope in ((6.2, 0.2), (0.0, -0.2)):
            wrapped = np.mod(intercept + slope * times, 2.0 * np.pi)
            fit = j2.fit_raan_drift(times, wrapped)
            self.assertAlmostEqual(fit["slope_rad_s"], slope)
            np.testing.assert_allclose(fit["unwrapped_rad"], intercept + slope * times, atol=1e-14)
        with self.assertRaises(ValueError):
            j2.fit_raan_drift(times, np.full(4, np.nan))
        with self.assertRaises(ValueError):
            j2.secular_raan_rate(orbit.ORBIT_RADIUS, 0.0, 0.0)
        with self.assertRaises(ValueError):
            osculating_elements(np.array([orbit.ORBIT_RADIUS, 0.0, 0.0, 100.0, 0.0, 0.0]))


class SevenDayJ2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = j2.run_comparison()
        cls.report = cls.data["report"]

    def test_shared_initial_state_and_seven_day_duration(self):
        solutions = self.data["solutions"]
        for solution in solutions.values():
            self.assertEqual(solution.t[-1], 604_800.0)
            np.testing.assert_array_equal(solution.y[:, 0], solutions["two_body"].y[:, 0])
        elements = osculating_elements(solutions["nominal"].y[:, 0])
        self.assertAlmostEqual(elements["inclination_rad"][0], np.deg2rad(51.6))

    def test_disabled_j2_matches_inclined_circular_solution(self):
        self.assertLess(self.report["two_body"]["max_position_error_vs_circular_m"], 5.0)
        self.assertLess(self.report["two_body"]["max_velocity_error_vs_circular_m_s"], 0.01)
        self.assertLess(abs(self.report["two_body"]["fitted_raan_rate_rad_s"]), 1e-15)

    def test_matching_total_energy_and_axial_momentum_are_conserved(self):
        for name in ("nominal", "tighter", "reference"):
            errors = self.report["invariants"][name]
            self.assertLess(errors["max_abs_relative_energy_error"], 1e-8)
            self.assertLess(errors["max_abs_relative_hz_error"], 1e-8)
        for quantity in ("relative_energy_error", "relative_hz_error"):
            self.assertLess(self.report["invariants"]["tighter"][f"max_abs_{quantity}"],
                            self.report["invariants"]["nominal"][f"max_abs_{quantity}"] / 10.0)

    def test_secular_raan_rate_agrees_with_first_order_theory(self):
        raan = self.report["raan"]
        self.assertLess(raan["measured_rad_s"], 0.0)
        # First-order averaging is approximate, unlike the exact milestone 1 reference.
        self.assertLess(raan["relative_discrepancy_initial_elements"], 0.01)
        self.assertLess(raan["relative_discrepancy_mean_proxy"],
                        raan["relative_discrepancy_initial_elements"])

    def test_three_level_trajectory_convergence_and_stable_raan_rate(self):
        convergence = self.report["convergence"]
        for quantity in ("position_difference_m", "velocity_difference_m_s"):
            self.assertLess(convergence[f"max_tighter_{quantity}"],
                            convergence[f"max_nominal_{quantity}"] / 10.0)
        self.assertLess(self.report["raan"]["nominal_tighter_rate_difference_rad_s"], 1e-11)


if __name__ == "__main__":
    unittest.main()
