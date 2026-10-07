"""Illustrative exponential density and rigidly co-rotating atmospheric drag."""

from dataclasses import dataclass

import numpy as np

from orbit import EARTH_RADIUS


# IERS numerical standards, mean Earth angular velocity [rad/s]:
# https://hpiers.obspm.fr/eop-pc/index.php?index=constants
EARTH_ROTATION_RATE = 7.2921150e-5


@dataclass(frozen=True)
class DragModel:
    """Constant spacecraft properties and illustrative atmosphere, all in SI.

    No epoch, weather indices, latitude/local-time structure, winds or attitude
    model. min_altitude_m is a propagation guard, not a reentry model. Density
    parameters are teaching choices, not a fitted or validated Earth profile.
    """

    mass_kg: float = 1000.0
    cd: float = 2.2  # dimensionless
    area_m2: float = 10.0  # constant effective area normal to relative flow
    rho_ref_kg_m3: float = 4e-12  # illustrative density at h_ref
    h_ref_m: float = 400_000.0
    scale_height_m: float = 60_000.0
    density_multiplier: float = 1.0
    omega_rad_s: float = EARTH_ROTATION_RATE  # fixed inertial +z axis
    min_altitude_m: float = 200_000.0

    def __post_init__(self):
        for name in self.__dataclass_fields__:
            if not np.isfinite(getattr(self, name)):
                raise ValueError(f"{name} must be finite.")
        for name in ("mass_kg", "scale_height_m"):
            if getattr(self, name) <= 0.0:
                raise ValueError(f"{name} must be positive.")
        for name in ("cd", "area_m2", "rho_ref_kg_m3", "h_ref_m",
                     "density_multiplier", "min_altitude_m"):
            if getattr(self, name) < 0.0:
                raise ValueError(f"{name} must be nonnegative.")

    @property
    def enabled(self):
        return (self.rho_ref_kg_m3 > 0.0 and self.density_multiplier > 0.0
                and self.cd > 0.0 and self.area_m2 > 0.0)

    def density(self, altitude_m):
        """rho [kg/m^3] = multiplier * rho_ref * exp(-(h - h_ref)/H)."""
        altitude_m = np.asarray(altitude_m)
        if self.rho_ref_kg_m3 == 0.0 or self.density_multiplier == 0.0:
            return np.zeros_like(altitude_m, dtype=float)
        # Do not silently clip a physically different atmosphere on overflow.
        with np.errstate(over="raise", invalid="raise"):
            return (self.density_multiplier * self.rho_ref_kg_m3
                    * np.exp(-(altitude_m - self.h_ref_m) / self.scale_height_m))

    def relative_velocity(self, position, velocity):
        """v_rel [m/s] for (3,) or (3,N) arrays; rigid rotation without winds."""
        position, velocity = np.asarray(position), np.asarray(velocity)
        atmospheric_velocity = self.omega_rad_s * np.array([
            -position[1], position[0], np.zeros_like(position[0])])
        return velocity - atmospheric_velocity

    def acceleration(self, position, velocity):
        """Drag [m/s^2], opposite relative flow; accepts (3,) or (3,N).

        NASA atmospheric rotation/drag equations, slide 22:
        https://ntrs.nasa.gov/api/citations/20160012488/downloads/20160012488.pdf#page=22
        """
        position, velocity = np.asarray(position), np.asarray(velocity)
        if not self.enabled:
            return np.zeros_like(velocity, dtype=float)
        relative = self.relative_velocity(position, velocity)
        altitude = np.linalg.norm(position, axis=0) - EARTH_RADIUS
        return (-0.5 * self.density(altitude) * self.cd * self.area_m2 / self.mass_kg
                * np.linalg.norm(relative, axis=0) * relative)
