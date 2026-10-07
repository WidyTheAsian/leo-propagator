# leo-propagator

A learning project for building and understanding satellite orbit propagation
in Python, with an eventual comparison against Orekit. The implementation is
kept small so that the physical model, numerical integration and verification
can be read and explained directly.

**Milestone 1 remains the default:** a satellite initially in a circular orbit
400,000 m above a spherical Earth, subject only to Earth's central gravity.
All calculations, saved values and plot axes use SI units; angles use radians.

**Milestone 2 adds optional J2 gravity:** a seven-day comparison of two-body
and J2 motion from the same 400,000 m, 51.6-degree initial Cartesian state.
The original tests and example outputs are preserved. J2 results are saved
separately under `outputs/j2/`.

## Implemented capabilities

- Integrate all six Cartesian position and velocity components for 86,400 s
  in a clearly defined Earth-centred inertial frame using `solve_ivp` and DOP853.
- Calculate the circular speed, analytical orbital period and exact trajectory.
- Compare numerical position and velocity with the analytical solution at
  matching times, and quantify altitude, specific energy and angular momentum errors.
- Check convergence by repeating the propagation with tighter tolerances.
- Save orbit, altitude and validation figures, plus a JSON report, and run four
  original milestone 1 physics and numerical verification tests.
- Enable an axisymmetric J2 acceleration and verify the matching total energy
  and axial angular momentum; compare secular nodal regression with theory.
- Extract osculating elements with explicit undefined-angle handling, and run
  13 additional tests for J2, orbital geometry and seven-day convergence.

## Setup and run

Run the commands from the repository root. Use Python 3.9 through 3.12 with
the versions pinned in `requirements.txt`. Local environments are ignored by
Git; the environment definitions remain in the repository.

To create an environment with Conda (Python 3.9, matching the validation run's
Python version family):

```powershell
conda env create --prefix ./.conda --file environment.yml
conda activate ./.conda
python orbit.py
python j2.py --days 7
python -m unittest -v
python -m pip check
```

Alternatively, on Windows with Python 3.12 installed:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe orbit.py
.\.venv\Scripts\python.exe j2.py --days 7
.\.venv\Scripts\python.exe -m unittest -v
.\.venv\Scripts\python.exe -m pip check
```

On Linux/macOS with Python 3.12 installed:

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python orbit.py
.venv/bin/python j2.py --days 7
.venv/bin/python -m unittest -v
.venv/bin/python -m pip check
```

If the project's `.conda` environment already exists, activate it instead of
creating it again. For a direct invocation in PowerShell without Conda activation:

```powershell
# Supply the DLL paths normally set by conda activation.
$env:PATH = "$PWD\.conda;$PWD\.conda\Library\bin;$PWD\.conda\Scripts;$env:PATH"
.\.conda\python.exe -m pip install -r requirements.txt
.\.conda\python.exe orbit.py
.\.conda\python.exe j2.py --days 7
.\.conda\python.exe -m unittest -v
```

Both scripts save figures without opening windows. Add `--show` for plot
windows or `--output-dir path` to select the output directory. `orbit.py`
runs the original 24-hour equatorial two-body example. `j2.py` defaults to
seven days; `--days` selects another positive duration. Each example performs
its own tolerance comparison. Use `python -m unittest -v test` to run only
the four original tests, or `python -m unittest -v` to run all 17.

## Files and outputs

- `orbit.py`: shared constants, central/J2 gravity, initial state and solver,
  plus the original milestone 1 reference, diagnostics and plotting.
- `j2.py`: inclined comparison, correct J2 invariants, secular RAAN analysis
  and plots, including three integration tolerance levels.
- `elements.py`: instantaneous Keplerian elements and undefined-angle masks.
- `test.py`: four original physics and full-day checks, unchanged.
- `test_j2.py`: 13 J2, element geometry and seven-day verification checks.
- `requirements.txt`: NumPy, SciPy and Matplotlib versions used for this milestone.
- `outputs/orbit_and_altitude.png`: the numerical orbit overlaid on the exact
  circle, and altitude versus elapsed time.
- `outputs/validation.png`: position, velocity, energy and angular momentum errors.
- `outputs/validation.json`: model constants, solver settings, package versions,
  maximum and final errors, and the tighter-tolerance comparison.
- `outputs/j2/comparison.png`, `outputs/j2/validation.png` and
  `outputs/j2/validation.json`: separate milestone 2 figures and numerical report.

![Orbit and altitude](outputs/orbit_and_altitude.png)

![Numerical validation](outputs/validation.png)

## Frame, units and assumptions

The frame is an idealised, right-handed Earth-centred inertial (ECI) frame:

- Origin: the centre of a stationary, spherical Earth.
- `+z`: a fixed reference north-pole direction.
- `+x`: an arbitrary fixed direction in the reference equatorial plane.
- `+y`: completes the right-handed system, with unit vector `z cross x`.

The axes stay fixed throughout the integration; they do not rotate with Earth.
The initial orbit lies in the `xy` plane and travels counterclockwise when
viewed from `+z`. This defines an idealised frame, not a full implementation of
GCRF, EME2000 or Earth orientation transformations. Time is elapsed seconds
from an arbitrary initial instant, so no calendar epoch or time-scale conversion
is needed for this autonomous model.

| Quantity | Adopted value | SI unit |
| --- | ---: | --- |
| Earth gravitational parameter, `mu = GM` | 3.986004418e14 | m³/s² |
| Spherical Earth radius, `R_E` | 6,371,000 | m |
| Initial altitude, `h_0` | 400,000 | m |
| Initial orbital radius, `r_0 = R_E + h_0` | 6,771,000 | m |
| Propagation duration | 86,400 | s |

The radius is the volumetric mean radius listed in
[NASA's Earth fact sheet](https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html).
`mu` is an adopted model parameter, used identically in the force law and reference
solution. Altitude means `norm(position) - R_E`; it is not ellipsoidal/geodetic
altitude. The satellite is a test particle: its mass cancels from the equations.
Milestone 1 uses only Earth's central gravity. Milestone 2 adds the J2 gravity
term, with the symmetry axis fixed along inertial `z`, while keeping spherical
altitude as the reporting convention. Both modes omit atmosphere and drag,
other bodies, radiation pressure, relativity and manoeuvres. Higher gravity
harmonics, time-varying coefficients and pole motion are not included.

## Equations and initial conditions

Define the state as

$$
\mathbf y = [x,y,z,v_x,v_y,v_z]^T.
$$

Newton's law for a central gravitational field gives

$$
\frac{d\mathbf r}{dt}=\mathbf v,\qquad
\frac{d\mathbf v}{dt}=-\frac{\mu}{\|\mathbf r\|^3}\mathbf r.
$$

The vector acceleration points inward. Its magnitude is `mu / r²`: the third
power in the vector expression includes the extra factor needed to turn the
position vector into a unit vector. `two_body_rhs` returns these six derivatives,
not the next state. Its `t` argument is unused because this force law has no
explicit dependence on time.

### Why is that the initial velocity?

For a circle, the required inward acceleration is `v²/r_0`. Equating this to
`mu/r_0²` yields

$$
v_c=\sqrt{\frac{\mu}{r_0}}=7672.598648385\ \mathrm{m/s}.
$$

The velocity must also be tangent to the circle. Starting on the positive
`x` axis, choose motion along positive `y`:

$$
\mathbf r(0)=(6771000,0,0)\ \mathrm m,\qquad
\mathbf v(0)=(0,7672.598648385,0)\ \mathrm{m/s}.
$$

Tangency alone does not make an orbit circular; the speed must have this value
as well. The inertial orbital velocity is used directly: no Earth-rotation speed
is added. With central gravity and both `z` and `v_z` initially zero, the orbit
remains in this plane, although the code integrates all three spatial dimensions.

### How does the integrator advance the orbit?

`solve_ivp` calls `two_body_rhs` at trial states. DOP853 is an explicit
eighth-order Runge–Kutta method: each step combines derivative estimates from
several intermediate states to advance position and velocity together. In
general Runge–Kutta notation, with step length `dt`,

$$
\mathbf k_i=f\left(t+c_i\,dt,
  \mathbf y+dt\sum_{j<i}a_{ij}\mathbf k_j\right),\qquad
\mathbf y_{\rm next}=\mathbf y+dt\sum_i b_i\mathbf k_i.
$$

The method's coefficients specify its trial stages and weights. An error
estimate controls whether a step is accepted and how large the next step is.
It is scaled using `atol_i + rtol * abs(y_i)` for each state component. See the
[SciPy solver documentation](https://docs.scipy.org/doc/scipy-1.13.1/reference/generated/scipy.integrate.solve_ivp.html).

The nominal run uses `rtol = 1e-10`, position `atol = 1e-4 m`, and velocity
`atol = 1e-7 m/s`. Relative tolerance is dimensionless; each absolute tolerance
has the units of its own component. These control estimated local error, and
do not guarantee a particular final position error after a full day.

`t_eval` requests saved states every 10 s, including both endpoints. It does
not impose 10 s integration steps. The solver adapts its internal steps and
interpolates to the requested times. The report's RHS evaluation count measures
derivative calls, including trial stages and interpolation work; it is not a
step count. The script checks successful completion before using the result.

### What evidence shows the result is correct?

The angular rate and period follow from the circular speed:

$$
n=\sqrt{\frac{\mu}{r_0^3}},\qquad
T=\frac{2\pi}{n}=5544.855095981\ \mathrm s.
$$

The satellite completes about 15.582012245 revolutions in 86,400 s. A full day
is not an integer number of periods, so the final point should not coincide
with the initial point.

The exact trajectory at each output time is

$$
\mathbf r_{\rm exact}(t)=r_0(\cos nt,\sin nt,0),\qquad
\mathbf v_{\rm exact}(t)=v_c(-\sin nt,\cos nt,0).
$$

The script compares full position and velocity vectors at matching times:
`norm(r_num - r_exact)` and `norm(v_num - v_exact)`. This detects radial,
orbital-plane and phase errors; matching the shape of a circle alone would
not test whether the satellite moves around it at the correct rate.

Two independent physical invariants are also checked:

$$
\epsilon=\tfrac12\|\mathbf v\|^2-\frac{\mu}{\|\mathbf r\|},\qquad
\mathbf h=\mathbf r\times\mathbf v.
$$

Their exact values for this initial condition are

$$
\epsilon_0=-\frac{\mu}{2r_0}=-29434385.009599764\ \mathrm{m^2/s^2},
\qquad
\mathbf h_0=(0,0,51951165448.21493)\ \mathrm{m^2/s}.
$$

Energy is constant because gravity is conservative; angular momentum is constant
because a central force has zero torque, `r cross a = 0`. Diagnostics use the
exact invariant values as references, which agree with the initial numerical
values to floating-point precision. Relative energy error is
`(epsilon - epsilon_0) / abs(epsilon_0)`; relative angular momentum error is
`norm(h - h_0) / norm(h_0)`. The latter compares the whole vector, so it can
detect an incorrect orbital-plane direction as well as a magnitude change.
The reported maxima take absolute values; final altitude and energy errors
retain their signs.

DOP853 does not enforce these invariants exactly. Small conservation errors,
small time-matched trajectory errors and improvement when tolerances tighten
together validate this implementation within the stated model.

## Measured validation results

Run with Python 3.9.13, NumPy 2.0.2, SciPy 1.13.1 and Matplotlib 3.9.4.
The solver reached 86,400 s successfully and returned 8,641 saved states, using
6,191 RHS evaluations. These maxima are measured on the 10 s output grid,
including interpolated states; they are not rigorous bounds between samples.

| Error quantity | Maximum absolute value | Final value at 86,400 s |
| --- | ---: | ---: |
| Altitude error [m] | 1.420762e-3 | -2.268264e-4 |
| Position error norm [m] | 7.625973e-3 | 6.848396e-3 |
| Velocity error norm [m/s] | 8.388148e-6 | 7.838255e-6 |
| Specific energy error [m²/s²] | 8.693080e-3 | -2.459142e-3 |
| Relative specific energy error | 2.953376e-10 | -8.354659e-11 |
| Specific angular momentum vector error norm [m²/s] | 7.671570 | 2.170174 |
| Relative specific angular momentum vector error | 1.476689e-10 | 4.177334e-11 |

Reducing every tolerance by a factor of 100 gives `rtol = 1e-12`, position
`atol = 1e-6 m`, and velocity `atol = 1e-9 m/s`. That run uses 10,991 RHS
evaluations. Its maximum position error is **1.199673e-4 m**, an improvement
of **63.57 times**. Its maximum velocity error is 1.335773e-7 m/s, relative
energy error is 3.111666e-12, and relative angular momentum error is
1.555950e-12. The dynamics, initial state and saved times are identical in both
runs. Error reduction supports convergence; adaptive local tolerances do not
imply an identical reduction factor in accumulated trajectory error.

The orbit curves overlap at the figure's scale. The altitude panel uses a
400,000 +/- 1 m range, so the small altitude error is shown in its annotation.
The validation figure makes the trajectory and conservation errors visible.
Errors can oscillate because both integration and output interpolation
contribute; they need not grow monotonically with time.

All **4 original milestone 1 tests pass**. The trajectory checks require
maximum position error below 0.1 m, velocity error below 1e-4 m/s and altitude
error below 0.01 m.
Both relative conservation errors must stay below 1e-8, and the tighter run
must reduce both position and velocity errors by more than a factor of 10.
The initial-state check verifies tangency, inward gravity and the centripetal
acceleration condition; the trajectory check also verifies the full duration
and preservation of `z = v_z = 0`. `pip check` reports no broken requirements.

Full-precision values are saved in `outputs/validation.json`. Small differences
in the last digits may occur across platforms. These are numerical errors
relative to the specified ideal model.

## Numerical verification and prediction accuracy

The errors above measure agreement with the exact solution of the same
idealised two-body model. The trajectory comparisons, conservation checks and
tolerance experiment verify the implementation and its numerical accuracy
for this case. The analytical reference uses the same adopted Earth constants,
so these checks do not establish their suitability for a particular real satellite.

Real-world orbit-prediction accuracy also depends on the initial-state
uncertainty, force models, atmospheric density and satellite properties, and
the physical reference frame and epoch. Earth's oblateness and atmospheric
drag matter for an actual low Earth orbit. Milestone 1 omits both; milestone 2
includes J2 but still omits drag and higher gravity harmonics.
This project has not been compared with tracking observations or Orekit yet;
the small numerical errors are not a claim of millimetre accuracy for real
satellite position predictions.

## Reading and explaining the code

Follow `initial_state` → `two_body_rhs` → `propagate` first. These contain the
physical model and numerical integration. Then read `analytical_state` and
`diagnostics` to see how the result is checked. `summarize`, `plot_results` and
`main` turn those arrays into reports and figures.

For an interactive inspection:

```python
import orbit

solution = orbit.propagate()
print(solution.y.shape)  # (6, 8641): state components by saved times
print(solution.t[-1])    # 86400.0 s
position = solution.y[:3]  # metres
velocity = solution.y[3:]  # metres per second
errors = orbit.diagnostics(solution.t, solution.y)
print(orbit.summarize(errors))
```

When later comparing with Orekit, match `mu`, the inertial axes, initial
Cartesian state, elapsed times and force model. Match the spherical radius
when comparing altitudes. Using a different Earth radius changes the initial
orbital radius even if both cases are described as a 400,000 m orbit.

## Milestone 2: J2 gravitational perturbation

### Force model, potential and constants

J2 is the leading zonal term describing Earth's oblateness. The model is
axisymmetric about the fixed inertial `z` axis. Rotation about that axis does
not change this zonal potential, so no Earth-rotation transformation is needed
for this idealised model.

| Quantity | Value | Unit |
| --- | ---: | --- |
| `mu` (unchanged) | 3.986004418e14 | m³/s² |
| Positive `J2` | 1.08263e-3 | dimensionless |
| J2 reference equatorial radius `R_eq` | 6,378,137 | m |
| Spherical altitude reference `R_E` | 6,371,000 | m |
| Initial inclination | 0.9005898940290741 (51.6 degrees) | rad |
| Comparison duration | 604,800 | s |

`J2` and `R_eq` come from [NASA's Earth fact sheet](https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html).
The two radii have different roles: `R_eq` sets the gravity harmonic's scale;
`R_E` defines the existing spherical altitude. The satellite still starts at
radius 6,771,000 m. A positive J2 corresponds to negative **unnormalized** C20.

Using `r = norm(position)` and `s = z/r`, the specific potential, zero at infinity, is

$$
\Phi(\mathbf r)=-\frac{\mu}{r}
\left[1-J_2\left(\frac{R_{\rm eq}}r\right)^2P_2(s)\right],
\qquad P_2(s)=\tfrac12(3s^2-1).
$$

This sign convention agrees with the Hamiltonian in
[Arnas (2022), equation (1)](https://arxiv.org/pdf/2212.08746#page=3).
Differentiating the potential gives `a = -grad(Phi)`: central gravity plus

$$
\mathbf a_{J2}=\frac{3J_2\mu R_{\rm eq}^2}{2r^5}
\begin{bmatrix}
x(5s^2-1)\\
y(5s^2-1)\\
z(5s^2-3)
\end{bmatrix}.
$$

Acceleration has units m/s² and potential has units m²/s². For positive J2,
the additional acceleration is inward on the equator and outward on the pole.
`two_body_rhs(t, state, j2=0.0)` and `propagate(..., j2=0.0)` retain pure
central gravity. Passing `j2=orbit.J2` adds the perturbation; zero skips it
exactly. `specific_potential` uses the same coefficient and reference radius.

### Shared initial conditions

At the ascending node, set initial RAAN to zero and rotate the original
circular velocity about `x` by the inclination:

$$
\mathbf r_0=(6771000,0,0)\ \mathrm m,
\qquad
\mathbf v_0=(0,4765.817619411,6012.965361426)\ \mathrm{m/s}.
$$

The speed is still `sqrt(mu/r_0)`, and velocity remains perpendicular to
position. Both simulations use **this identical Cartesian state**. It is
circular in the two-body model, and initially has zero osculating Keplerian
eccentricity in the J2 case. It is not an exact circular solution of the J2
field: periodic radius and element changes are physical responses to the
additional force. The velocity is not retuned for J2.

```python
import numpy as np
import orbit

state0 = orbit.initial_state(np.deg2rad(51.6))
central = orbit.propagate(state0=state0, duration=604800.0, j2=0.0)
perturbed = orbit.propagate(state0=state0, duration=604800.0, j2=orbit.J2)
```

### Osculating elements, RAAN and undefined angles

`elements.osculating_elements` treats each Cartesian state as an instantaneous
Keplerian orbit. It computes `h = r cross v`, node vector `N = z_hat cross h`,
eccentricity vector `e = (v cross h)/mu - r/r`, and `a` from the instantaneous
two-body energy. That energy is used to define an **osculating element**, not
as the conserved J2 energy.

RAAN is `atan2(N_y, N_x)`, measured from `+x` toward `+y` in the equatorial
plane. Returned angles are radians in `[0, 2*pi)`. The comparison applies
`np.unwrap` before fitting RAAN as `intercept + rate*time`. The free intercept
and the fit over all 10,081 output states reduce the influence of periodic
oscillations compared with a difference of two endpoints. Saving every 60 s
is sufficiently fine to avoid ambiguous RAAN jumps in this example.

| Case | Undefined angles returned as `NaN` | Useful defined angle |
| --- | --- | --- |
| Circular, inclined | Argument of periapsis, true anomaly | RAAN and argument of latitude |
| Eccentric, equatorial | RAAN, argument of periapsis, argument of latitude | True anomaly |
| Circular, equatorial | RAAN, argument of periapsis, true anomaly, argument of latitude | Geometric true longitude `atan2(y, x)` |

Circular means `e <= 1e-8`; equatorial means `sin(i) <= 1e-10`, including
the retrograde `i = pi` limit. Longitude is measured geometrically from `+x`
toward `+y`. Radial states with zero angular momentum have no orbital plane
and raise `ValueError`. RAAN fitting and the secular node prediction reject
undefined equatorial nodes explicitly. The JSON report represents undefined
initial angles as `null` rather than nonstandard JSON NaNs.

### First-order secular node prediction

For mean elements `a`, `e`, `i`, the standard first-order J2 prediction is

$$
\dot\Omega=-\frac32J_2 n\left(\frac{R_{\rm eq}}p\right)^2\cos i,
\qquad n=\sqrt{\frac{\mu}{a^3}},\quad p=a(1-e^2).
$$

See [Rickman, NASA (2014), slide 103](https://ntrs.nasa.gov/api/citations/20200007746/downloads/20200007746.pdf#page=103).
The negative rate for this prograde orbit describes westward nodal regression.
The main comparison first inserts the initial osculating values
`a = 6,771,000 m`, `e = 0`, `i = 0.900589894 rad`. It also reports a second
estimate using arithmetic averages of the sampled osculating `a`, `e`, `i`.
These averages are a diagnostic **mean-element proxy**, not a formal
osculating-to-mean transformation.

### Actual seven-day results

Run with Python 3.9.13, NumPy 2.0.2, SciPy 1.13.1 and Matplotlib 3.9.4.
Both cases reach 604,800 s and share 10,081 output times, spaced by 60 s.
All maxima below refer to those saved samples, including interpolated states.

| RAAN quantity | Measured value |
| --- | ---: |
| Fitted two-body drift [rad/s] | 5.631485e-21 |
| Fitted J2 secular drift [rad/s] | -1.018645894336e-6 |
| First order, initial elements [rad/s] | -1.014233568650e-6 |
| Relative discrepancy using initial elements | 0.435040% |
| First order, sample mean proxy [rad/s] | -1.017827427685e-6 |
| Relative discrepancy using mean proxy | 0.080413% |
| Maximum absolute RAAN fit residual [rad] | 4.512172e-4 |
| Final unwrapped J2 RAAN [rad] | -0.616358533 |

The measured RAAN contains short-period oscillations as well as secular drift.
The residual panel shows four oscillations in roughly two nominal orbital
periods; a finite-duration line fit does not remove them perfectly. More
fundamentally, first-order averaged theory expects mean elements, while the
specified initial state supplies osculating elements. The sample means are
`a = 6,765,019.335 m`, `e = 0.000945600`, and `i = 0.900239016 rad`; they differ
from the initial values. Their use improves agreement, supporting the
mean-versus-osculating explanation. The remaining discrepancy also includes
averaging choices and higher-order effects omitted by first-order theory.
This interpretation is consistent with the distinction and transformations
discussed in [Arnas (2022)](https://arxiv.org/abs/2212.08746).

Tightening tolerances changes the fitted rate by only **7.311555e-16 rad/s**,
far below the 4.412326e-9 rad/s initial-element theory discrepancy. The latter
is therefore not explained by the integration error at these settings.

The J2 spherical altitude ranges from **387,786.713 m to 400,000.000 m**.
The two-body altitude stays within 0.004 m of 400,000 m. The substantial J2
radius variation follows from the shared Keplerian initial state; total
energy conservation rules out interpreting it as dissipative orbital decay.

![J2 trajectories, altitude and RAAN](outputs/j2/comparison.png)

### J2 conservation and convergence

The time-independent potential conserves

$$
E=\tfrac12\|\mathbf v\|^2+\Phi(\mathbf r),
\qquad h_z=(\mathbf r\times\mathbf v)_z.
$$

Axial symmetry gives zero torque about `z`, preserving `h_z`. The full
angular momentum vector rotates as the orbital plane precesses. The two-body
energy `v²/2 - mu/r` exchanges with the J2 potential contribution and also
varies physically. For this run, its maximum change is **52,148.049 m²/s²**,
and the full angular momentum vector changes by **0.475124** of its initial
magnitude. Neither is required to be constant in the J2 tests.

Energy error is referenced to the initial **matching total energy** and
normalized by its absolute value. Axial momentum error is `h_z(t)-h_z(0)`,
normalized by initial `norm(h)` so that the scale remains defined for polar
orbits with zero `h_z`. Absolute errors are also recorded in the JSON report.

| Solver level | `rtol` | Position `atol` [m] | Velocity `atol` [m/s] |
| --- | ---: | ---: | ---: |
| Nominal | 1e-10 | 1e-4 | 1e-7 |
| Tighter | 1e-12 | 1e-6 | 1e-9 |
| Numerical reference | 1e-13 | 1e-7 | 1e-10 |

| Seven-day J2 check | Nominal | Tighter |
| --- | ---: | ---: |
| Maximum absolute total specific energy error [m²/s²] | 1.719590e-2 | 1.114383e-4 |
| Maximum relative total specific energy error | 5.836507e-10 | 3.782358e-12 |
| Maximum absolute `h_z` error [m²/s] | 9.391926 | 0.061234 |
| Maximum `abs(h_z error) / norm(h_0)` | 1.807837e-10 | 1.178675e-12 |
| Maximum position difference from reference [m] | 1.420498 | 0.006375425 |
| Maximum velocity difference from reference [m/s] | 0.001610305 | 7.225751e-6 |

The tighter position difference is **222.81 times smaller**. The reference
run's relative energy and axial momentum errors are 3.975307e-13 and
1.257096e-13. J2 trajectory differences are measured against this finite-
accuracy numerical reference, not an exact J2 trajectory or observations.
The inclined two-body run independently matches the analytical circular
trajectory within 1.544617 m in position over seven days.

All **17 tests pass**, including all four unchanged milestone 1 tests.
New checks cover zero-J2 equivalence, force versus finite-difference potential
gradient, equatorial/polar force signs, known element orientations, circular
and equatorial singularities, RAAN unwrapping, correct invariant conservation,
first-order rate agreement within 1%, and tolerance convergence. The nominal
and tighter J2 position and velocity differences must improve by more than
10 times, and both correct relative invariant errors must remain below 1e-8.
`orbit.py` also reproduces the committed milestone 1 plots and report unchanged.

![J2 conservation and convergence](outputs/j2/validation.png)

Read `orbit.j2_acceleration` and `orbit.specific_potential` together, then
`elements.osculating_elements`, `j2.fit_raan_drift` and `j2.run_comparison`.
The force integrates Cartesian states; orbital elements are extracted for
analysis, not used as singular integration coordinates. Full-precision
results, initial state, constants and all solver settings are in
[`outputs/j2/validation.json`](outputs/j2/validation.json).

## Remaining milestones

J2 is implemented; the following remain future work.

3. **Atmospheric drag:** introduce an explicit density model and satellite
   drag properties, then examine orbital decay and sensitivity to those inputs.
4. **Orekit comparison:** match the initial Cartesian state, inertial frame,
   epoch, Earth constants, force models and output times; compare trajectories
   while accounting for numerical tolerances and model differences.
