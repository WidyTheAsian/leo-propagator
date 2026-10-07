"""Optional Orekit force models and Cartesian numerical propagation.

Only this module needs the Java wrapper. No local acceleration or density
function is called by an Orekit force model. The Python callback defines
only the prescribed uniform frame rotation.
"""

from importlib.metadata import version
from types import SimpleNamespace

import numpy as np
import orekit_jpype
import jpype

if not jpype.isJVMStarted():
    orekit_jpype.initVM()

# Java imports must follow JVM startup, as documented by the JPype wrapper.
from java.lang import System
from java.util import Properties
from org.hipparchus.geometry.euclidean.threed import (
    FieldRotation, FieldVector3D, Rotation, RotationConvention, Vector3D)
from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
from org.orekit.bodies import OneAxisEllipsoid
from org.orekit.forces.drag import DragForce, IsotropicDrag
from org.orekit.forces.gravity import HolmesFeatherstoneAttractionModel, NewtonianAttraction
from org.orekit.forces.gravity.potential import GravityFieldFactory, TideSystem
from org.orekit.frames import FieldTransform, Frame, FramesFactory, Transform, TransformProvider
from org.orekit.models.earth.atmosphere import SimpleExponentialAtmosphere
from org.orekit.orbits import CartesianOrbit, OrbitType
from org.orekit.propagation import SpacecraftState
from org.orekit.propagation.events import AltitudeDetector
from org.orekit.propagation.numerical import NumericalPropagator
from org.orekit.time import AbsoluteDate, TimeScalesFactory
from org.orekit.utils import PVCoordinates

import orbit
from atmosphere import DragModel


EPOCH_LABEL = "2026-01-01T00:00:00 TAI"
MAX_STEP = orbit.PERIOD / 16.0  # [s], also the altitude event max-check interval
MIN_STEP = 0.001  # [s]
MASS_ATOL = 1e-9  # [kg]; mass is a constant seventh Orekit state
EVENT_THRESHOLD = 1e-9  # Orekit root-finding time tolerance [s]


def vector(values):
    return Vector3D(*[float(value) for value in values])


def array(value):
    return np.array([value.getX(), value.getY(), value.getZ()])


@jpype.JImplements(TransformProvider)
class UniformRotation:
    """GCRF to atmosphere: passive rotation by omega*t about fixed +z.

    Angular rate is expressed in the new frame; z is common to both frames.
    Field transforms use the same formula with field-valued elapsed time.
    """

    def __init__(self, epoch, omega):
        self.epoch, self.omega = epoch, omega

    @jpype.JOverride
    def getTransform(self, date):
        if isinstance(date, AbsoluteDate):
            rotation = Rotation(Vector3D.PLUS_K, self.omega * date.durationFrom(self.epoch),
                                RotationConvention.FRAME_TRANSFORM)
            return Transform(date, rotation, vector([0.0, 0.0, self.omega]))
        field = date.getField()
        rotation = FieldRotation(FieldVector3D(field, Vector3D.PLUS_K),
                                 date.durationFrom(self.epoch).multiply(self.omega),
                                 RotationConvention.FRAME_TRANSFORM)
        rate = FieldVector3D(field, vector([0.0, 0.0, self.omega]))
        return FieldTransform(date, rotation, rate)

    # Explicit overrides avoid JPype's overloaded-default-method dispatch error.
    @jpype.JOverride
    def getStaticTransform(self, date):
        return self.getTransform(date).toStaticTransform()

    @jpype.JOverride
    def getKinematicTransform(self, date):
        return self.getTransform(date)


class OrekitModel:
    """Explicitly matched gravity, sphere and built-in exponential atmosphere."""

    def __init__(self, j2=0.0, drag=None):
        self.j2 = j2
        self.config = DragModel() if drag is None else drag
        self.frame = FramesFactory.getGCRF()
        self.epoch = AbsoluteDate(2026, 1, 1, 0, 0, 0.0, TimeScalesFactory.getTAI())
        self.rotation = UniformRotation(self.epoch, self.config.omega_rad_s)
        self.body_frame = Frame(self.frame, self.rotation, "Uniformly rotating toy atmosphere", False)
        # Zero flattening makes geodetic height exactly the spherical |r|-R_E.
        self.sphere = OneAxisEllipsoid(orbit.EARTH_RADIUS, 0.0, self.body_frame)
        self.atmosphere = SimpleExponentialAtmosphere(
            self.sphere, self.config.density_multiplier * self.config.rho_ref_kg_m3,
            self.config.h_ref_m, self.config.scale_height_m)
        self.forces = {"central": NewtonianAttraction(orbit.MU)}
        if j2 != 0.0:
            # Fully normalized Pbar_20 = sqrt(5)*P2. Only degree 2, order 0.
            coefficients = jpype.JArray(jpype.JDouble, 2)([[1.0], [0.0], [-j2 / np.sqrt(5.0)]])
            sine = jpype.JArray(jpype.JDouble, 2)([[0.0], [0.0], [0.0]])
            self.provider = GravityFieldFactory.getNormalizedProvider(
                orbit.EQUATORIAL_RADIUS, orbit.MU, TideSystem.UNKNOWN, coefficients, sine)
            # This force supplies the NON-CENTRAL part; central is added once above.
            self.forces["j2"] = HolmesFeatherstoneAttractionModel(self.frame, self.provider)
        if drag is not None and drag.enabled:
            self.forces["drag"] = DragForce(
                self.atmosphere, IsotropicDrag(drag.area_m2, drag.cd))

    def state(self, state_si, time=0.0):
        date = self.epoch.shiftedBy(float(time))
        state_si = np.asarray(state_si, dtype=float)
        pv = PVCoordinates(vector(state_si[:3]), vector(state_si[3:]))
        initial = CartesianOrbit(pv, self.frame, date, orbit.MU)
        return SpacecraftState(initial).withMass(self.config.mass_kg)

    def accelerations(self, state_si, time=0.0):
        state = self.state(state_si, time)
        return {name: array(force.acceleration(state, force.getParameters(state.getDate())))
                for name, force in self.forces.items()}

    def environment(self, position_si, time=0.0):
        date, position = self.epoch.shiftedBy(float(time)), vector(position_si)
        height = self.sphere.transform(position, self.frame, date).getAltitude()
        return {"altitude_m": float(height),
                "density_kg_m3": float(self.atmosphere.getDensity(date, position, self.frame)),
                "atmospheric_velocity_m_s": array(self.atmosphere.getVelocity(date, position, self.frame))}

    def propagate(self, state0, duration, sample_interval=60.0, rtol=orbit.RTOL,
                  atol=orbit.ATOL, mass_atol=MASS_ATOL, stop_altitude=None):
        """Hipparchus DOP853, numerical Cartesian integration plus dense ephemeris.

        Output grid and downward spherical-altitude stop match orbit.propagate.
        No force evaluations are delegated to the Python propagator.
        """
        if not np.isfinite(duration) or duration <= 0.0:
            raise ValueError("Duration must be finite and positive.")
        if not np.isfinite(sample_interval) or sample_interval <= 0.0:
            raise ValueError("Sample interval must be finite and positive.")
        state0 = np.asarray(state0, dtype=float)
        if state0.shape != (6,) or not np.all(np.isfinite(state0)):
            raise ValueError("Expected six finite initial Cartesian components.")
        if stop_altitude is not None:
            if not np.isfinite(stop_altitude) or stop_altitude < 0.0:
                raise ValueError("Stopping altitude must be finite and nonnegative.")
            if np.linalg.norm(state0[:3]) - orbit.EARTH_RADIUS <= stop_altitude:
                raise ValueError("Initial altitude must be above the stopping altitude.")
        integrator = DormandPrince853Integrator(
            MIN_STEP, MAX_STEP, np.append(np.broadcast_to(atol, (6,)), mass_atol).tolist(),
            [float(rtol)] * 7)
        propagator = NumericalPropagator(integrator)
        propagator.setOrbitType(OrbitType.CARTESIAN)
        # Disable automatic central attraction, then add our explicit built-in model once.
        propagator.setIgnoreCentralAttraction(True)
        propagator.setInitialState(self.state(state0))
        for force in self.forces.values():
            propagator.addForceModel(force)
        if stop_altitude is not None:
            propagator.addEventDetector(AltitudeDetector(
                float(stop_altitude), self.sphere).withMaxCheck(MAX_STEP).withThreshold(EVENT_THRESHOLD))
        generator = propagator.getEphemerisGenerator()
        final = propagator.propagate(self.epoch.shiftedBy(float(duration)))
        ephemeris = generator.getGeneratedEphemeris()
        end = float(final.getDate().durationFrom(self.epoch))
        force_names = [str(force.getClass().getSimpleName()) for force in propagator.getAllForceModels()]

        def sol(times):
            times = np.asarray(times, dtype=float)
            if np.any(times < 0.0) or np.any(times > end):
                raise ValueError("Cannot extrapolate an Orekit trajectory beyond its available interval.")
            columns = []
            for time in np.atleast_1d(times):
                state = final if time == end else ephemeris.propagate(self.epoch.shiftedBy(float(time)))
                pv = state.getPVCoordinates(self.frame)
                columns.append(np.concatenate((array(pv.getPosition()), array(pv.getVelocity()))))
            values = np.column_stack(columns)
            return values[:, 0] if times.ndim == 0 else values

        times = np.append(np.arange(0.0, end, sample_interval), end)
        return SimpleNamespace(t=times, y=sol(times), sol=sol, end_time_s=end,
                               terminated_early=end < duration, rhs_evaluations=integrator.getEvaluations(),
                               force_names=force_names, final_mass_kg=float(final.getMass()))


def versions():
    def java_artifact_version(java_class, resource):
        properties = Properties()
        stream = java_class.class_.getResourceAsStream(resource)
        try:
            properties.load(stream)
            return str(properties.getProperty("version"))
        finally:
            stream.close()

    return {"wrapper": version("orekit-jpype"), "jpype": version("JPype1"),
            "orekit": java_artifact_version(SpacecraftState, "/META-INF/maven/org.orekit/orekit/pom.properties"),
            "hipparchus": java_artifact_version(DormandPrince853Integrator, "/META-INF/maven/org.hipparchus/hipparchus-ode/pom.properties"),
            "java": str(System.getProperty("java.runtime.version")),
            "java_vendor": str(System.getProperty("java.vendor")),
            "java_vm": str(System.getProperty("java.vm.name"))}
