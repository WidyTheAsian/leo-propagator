"""Instantaneous Keplerian (osculating) elements from ECI states, in SI units."""

import numpy as np

from orbit import MU


ECCENTRICITY_TOL = 1e-8  # below this, periapsis-based angles are undefined
NODE_TOL = 1e-10  # below this sin(i), the node-based angles are undefined


def osculating_elements(states):
    """Return arrays for states of shape (6,) or (6, N); all angles in radians.

    RAAN is measured from +x toward +y in the equatorial plane. Directed
    in-plane angles follow the angular momentum vector, using atan2.
    Undefined angles are NaN, never an arbitrary zero. The state must have
    finite components, nonzero radius, and a well-defined orbital plane.
    """
    states = np.asarray(states, dtype=float)
    if states.ndim == 1:
        states = states[:, None]
    if states.ndim != 2 or states.shape[0] != 6 or not np.all(np.isfinite(states)):
        raise ValueError("Expected finite states with shape (6,) or (6, N).")
    position, velocity = states[:3].T, states[3:].T
    radius = np.linalg.norm(position, axis=1)
    h = np.cross(position, velocity)
    h_norm = np.linalg.norm(h, axis=1)
    if np.any(radius == 0.0) or np.any(h_norm == 0.0):
        raise ValueError("Orbital elements require nonzero radius and angular momentum.")
    h_unit = h / h_norm[:, None]
    node = np.column_stack((-h[:, 1], h[:, 0], np.zeros(len(radius))))
    node_norm = np.linalg.norm(node, axis=1)
    e_vector = np.cross(velocity, h) / MU - position / radius[:, None]
    eccentricity = np.linalg.norm(e_vector, axis=1)
    # Keplerian energy defines osculating a; it is not the conserved J2 energy.
    keplerian_energy = 0.5 * np.sum(velocity**2, axis=1) - MU / radius
    with np.errstate(divide="ignore"):
        semimajor_axis = -MU / (2.0 * keplerian_energy)
    circular = eccentricity <= ECCENTRICITY_TOL
    equatorial = node_norm / h_norm <= NODE_TOL

    def angle(first, second, defined):
        sine = np.sum(np.cross(first, second) * h_unit, axis=1)
        cosine = np.sum(first * second, axis=1)
        return np.where(defined, np.mod(np.arctan2(sine, cosine), 2.0 * np.pi), np.nan)

    raan = np.where(~equatorial,
                    np.mod(np.arctan2(node[:, 1], node[:, 0]), 2.0 * np.pi), np.nan)
    return {
        "semimajor_axis_m": semimajor_axis,
        "eccentricity": eccentricity,
        "inclination_rad": np.arctan2(node_norm, h[:, 2]),
        "raan_rad": raan,
        "argument_of_periapsis_rad": angle(node, e_vector, ~equatorial & ~circular),
        "true_anomaly_rad": angle(e_vector, position, ~circular),
        "argument_of_latitude_rad": angle(node, position, ~equatorial),
        # Geometric azimuth replaces the undefined angles for a circular equatorial orbit.
        "true_longitude_rad": np.where(equatorial & circular,
                                        np.mod(np.arctan2(position[:, 1], position[:, 0]),
                                               2.0 * np.pi), np.nan),
        "circular": circular,
        "equatorial": equatorial,
    }
