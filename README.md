# initial orbit determination
Double-R method for preliminary (initial) orbit determination using three sets of right-ascension and declination angles.

This Fortran program is designed to solve the problem of orbit determination from three optical sightings using the Double-R iteration method. The solution provides the state vectors (position and velocity) of the object at the second observation time (t2).

## Input
- ra1, ra2, ra3: Right ascension at observation times t1, t2, and t3 (in radians).
- dec1, dec2, dec3: Declination at observation times t1, t2, and t3 (in radians).
- Mjd1, Mjd2, Mjd3: Modified Julian Date of observation times t1, t2, and t3.
- rSite1, rSite2, rSite3: Position vectors of observation sites in the Earth-Centered Inertial (ECI) frame (in kilometers).

## Output
- r: Position vector at observation time t2 (in kilometers).
- v: Velocity vector at observation time t2 (in kilometers per second).

## Execution
To run the program, execute the compiled binary in a Fortran environment. The program uses several helper functions from the **doubler_functions** module.

### Helping Module "doubler_functions"
This Fortran module provides subroutines and functions for performing the Double-R iteration. The main subroutine `doubler_itr` executes the core Double-R iteration, while `rv2ceo` calculates classical orbital elements from position and velocity vectors. Additionally, there are utility functions such as `site_vector`, `sidereal_time`, and `universal_time` for handling observation site coordinates and time conversions.

- **Subroutines:**
    - **doubler_itr:** Executes the main Double-R iteration. It takes various inputs related to observations, station positions, and times and computes values like true anomalies, eccentricities, and differences between calculated and observed times.
    - **degrad:** Transforms angles between degrees and radians based on the specified type.

- **Functions:**
    - **rv2ceo:** Calculates classical orbital elements (semi-major axis, eccentricity, inclination, RA of Ascending Node, argument of perigee, true anomaly) from position and velocity vectors.
    - **site_vector:** Computes the observation site position vector in Earth-Centered Inertial (ECI) coordinates from geodetic coordinates.
    - **sidereal_time:** Calculates the local sidereal time from the Julian Date, UT, and East Longitude.
    - **universal_time:** Computes Universal Time (UT) from a given Julian Date.
    - **CROSS:** Calculates the cross product between two vectors.
    - **angl2v:** Calculates the angle between two 3D vectors.

- **Additional Information:**
The module includes declarations for constants like Earth's radius, gravitational constant, and mathematical constants.
Implicit none is used for variable declaration to ensure explicit variable typing.
The module uses a selected real kind for numerical precision (ikind).

## Iteration
The program uses the Double-r iteration method to refine the estimates of the state vectors. The iteration process is repeated until convergence or a maximum iteration count is reached.

## Results
The final estimates of the state vectors (position and velocity) and classical orbital elements (semi-major axis, eccentricity, inclination, right ascension of ascending node, argument of perigee, and true anomaly) are provided.

## Note
- The program is initialized with initial guesses for the position vector at t1 and t3.
- The iteration process refines these guesses to find the state vectors at t2.

**Caution:** This script is intended for educational purposes and may require adaptations for specific use cases or real-world applications.
