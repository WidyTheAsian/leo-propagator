"""Physics and full-day accuracy checks; run with python -m unittest -v."""

import unittest

import numpy as np

import orbit


class CircularOrbitTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.solution = orbit.propagate()
        cls.summary = orbit.summarize(
            orbit.diagnostics(cls.solution.t, cls.solution.y)
        )

    def test_initial_state_supports_circular_motion(self):
        state = orbit.initial_state()
        position, velocity = state[:3], state[3:]
        acceleration = orbit.two_body_rhs(0.0, state)[3:]
        self.assertEqual(np.linalg.norm(position) - orbit.EARTH_RADIUS, 400_000.0)
        self.assertEqual(np.dot(position, velocity), 0.0)
        self.assertLess(np.dot(position, acceleration), 0.0)
        np.testing.assert_allclose(
            acceleration,
            -np.dot(velocity, velocity) * position / np.dot(position, position),
            rtol=1e-14, atol=0.0,
        )

    def test_full_day_trajectory_matches_analytical_orbit(self):
        self.assertEqual(self.solution.t[0], 0.0)
        self.assertEqual(self.solution.t[-1], 86_400.0)
        self.assertLess(self.summary["max_abs_position_error_m"], 0.1)
        self.assertLess(self.summary["max_abs_velocity_error_m_s"], 1e-4)
        self.assertLess(self.summary["max_abs_altitude_error_m"], 0.01)
        np.testing.assert_array_equal(self.solution.y[[2, 5]], 0.0)

    def test_energy_and_angular_momentum_are_conserved(self):
        self.assertLess(self.summary["max_abs_relative_energy_error"], 1e-8)
        self.assertLess(self.summary["max_abs_relative_angular_momentum_error"], 1e-8)

    def test_tighter_tolerances_reduce_trajectory_error(self):
        tighter = orbit.propagate(rtol=orbit.RTOL / 100.0, atol=orbit.ATOL / 100.0)
        summary = orbit.summarize(orbit.diagnostics(tighter.t, tighter.y))
        self.assertLess(summary["max_abs_position_error_m"],
                        self.summary["max_abs_position_error_m"] / 10.0)
        self.assertLess(summary["max_abs_velocity_error_m_s"],
                        self.summary["max_abs_velocity_error_m_s"] / 10.0)


if __name__ == "__main__":
    unittest.main()
