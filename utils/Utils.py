from dataclasses import dataclass, field
import numpy as np
from scipy.optimize import newton
import spiceypy as spy
from typing import Dict, Set
import multimin as mm

def Geo2Eclip(lon, lat, alt, date=None, et=None, frame='ITRF93'):
    """
    Converts geodetic coordinates (latitude, longitude, altitude) of an impact 
    event on Earth to ecliptic J2000 coordinates.

    Parameters:
    ----------
    lon : float
        Geodetic longitude of the impact site [degrees].
    lat : float
        Geodetic latitude of the impact site [degrees].
    alt : float
        Altitude of the impact site above Earth's reference spheroid [km].
    date : str
        UTC date and time of the impact event in format 'YYYY-MM-DD HH:MM:SS'.

    Returns:
    -------
    r_earth_ecl : ndarray
        Cartesian coordinates [km] of the impact site in the Ecliptic J2000 frame.

    Notes:
    ------
    1. The function first converts geodetic coordinates to Earth-centered 
       Cartesian coordinates in the ITRF93 frame.
    2. Then, it applies a transformation matrix to convert from ITRF93 
       (Earth-fixed) to ECLIPJ2000 (inertial, ecliptic-based) frame.
    3. This transformation is necessary for orbital calculations, ensuring 
       the position is in an inertial reference frame.
    """
    deg = np.pi/180

    lon = lon*deg
    lat = lat*deg

    n, props = spy.bodvrd('399','RADII',3)
    RE_spice = props[0]  #Equatorial radius of the reference spheroid.
    RP_spice = props[2]  #Polar radius of the reference spheroid.
    f_spice = (RE_spice-RP_spice)/RE_spice # Flattening coefficient.
    r_earth_fixed = spy.georec(lon, lat, alt, RE_spice, f_spice)
    #print("Equatorial and Polar Radios: ", RE_spice, RP_spice)

    if date: 
        et = spy.utc2et(date)  #Convert from UTC to ephemerides time
        #print("ET", et)
    
        r_earth_fixed = spy.georec(lon, lat, alt, RE_spice, f_spice)  #Convert geodetic coordinates to rectangular coordinates in the ITRF93 frame (rotante)
        #print("GeoRec", r_earth_fixed)
        M_ecl = spy.pxform(frame, 'ECLIPJ2000', et) 
        r_earth_ecl = spy.mxv(M_ecl, r_earth_fixed)  #from ITRF93 (rotante) frame to inertial frame ECLIPJ2000
    else: 
        print('et: ephemeris time is not provided')    
    return r_earth_ecl


def Geo2Rec(lon, lat, alt):
    """
    lon: (float) [°]
    lat: (float) [°]
    alt: (float) km
    date: (str) '2000-08-16 00:00:00'
    """
    deg = np.pi/180

    lon = lon*deg
    lat = lat*deg

    n, props = spy.bodvrd('399','RADII',3)
    RE_spice = props[0]
    RP_spice = props[2]
    #print("Equatorial and Polar Radios: ", RE_spice, RP_spice)
    f_spice = (RE_spice-RP_spice)/RE_spice
    #print(RE_spice, RP_spice)

    r_earth_fixed = spy.georec(lon, lat, alt, RE_spice, f_spice)

    return r_earth_fixed

def Geo2Eclip2(lon, lat, alt, date):
    """
    lon: (float) [°]
    lat: (float) [°]
    alt: (float) km
    date: (str) '2000-08-16 00:00:00'
    """

    et = spy.utc2et(date)
    r_earth_fixed = Geo2Rec(lon, lat, alt)
    mx = spy.pxform('IAU_EARTH', 'ECLIPJ2000', et)
    r_earth_ecl = spy.mxv(mx, r_earth_fixed)

    return r_earth_ecl

def z_axis_rotation(x: float) -> np.ndarray:
    return np.array([[np.cos(x), -np.sin(x), 0],[np.sin(x), np.cos(x), 0],[0,0,1]])

def mag(x) -> float:
    """
    Compute the magnitude of a vector.
    
    Parameters:
    -----------
    x : array-like
        Input vector (can be list, tuple, numpy array, etc.)
        
    Returns:
    --------
    float
        Magnitude of the vector
    """
    return (x@x)**0.5

def get_velocity_ecliptic(vx: float, vy: float, vz: float, lon: float, lat: float, alt: float, 
                         date: str | None = None, et: float | None = None) -> np.ndarray:
    """
    Convert velocity vector from Earth-fixed to ecliptic J2000 coordinates.
    
    This function takes a velocity vector in Earth-fixed coordinates and converts it
    to ecliptic J2000 coordinates, accounting for Earth's rotation.
    
    Parameters:
    -----------
    vx, vy, vz : float
        Velocity components in Earth-fixed coordinates [km/s]
    lon, lat : float
        Geodetic longitude and latitude of the observation point [degrees]
    alt : float
        Altitude above Earth's reference spheroid [km]
    date : str, optional
        UTC date and time in format 'YYYY-MM-DD HH:MM:SS'
        Must be provided if et is not provided
    et : float, optional
        Ephemeris time (ET) in seconds past J2000
        Must be provided if date is not provided
        
    Returns:
    --------
    v_eclip : np.ndarray
        Velocity vector in ecliptic J2000 coordinates [km/s]
        
    Raises:
    -------
    ValueError
        If neither date nor et is provided, or if both are provided
        
    Notes:
    ------
    The function accounts for Earth's rotation by adding the cross product
    of Earth's angular velocity and the position vector to the input velocity.
    """
    # Input validation
    if date is not None and et is not None:
        raise ValueError("Provide either 'date' or 'et', not both")
    if date is None and et is None:
        raise ValueError("Either 'date' or 'et' must be provided")
    
    # Convert velocity to numpy array
    v = np.array([vx, vy, vz])
    
    # Get position vector in Earth-fixed coordinates
    r = Geo2Rec(lon, lat, alt)
    
    # Earth's rotation parameters
    t_sidereal = 86164.09053083288  # Sidereal day in seconds
    w_earth = 2 * np.pi / t_sidereal  # Earth's angular velocity [rad/s]
    omega = np.array([0, 0, w_earth])
    
    # Add Earth's rotation contribution to velocity
    v_E = v + spy.vcrss(omega, r)  # Velocity in Earth-fixed frame [km/s]
    
    # Convert ephemeris time if date is provided
    if date is not None:
        et = spy.utc2et(date)
    
    # Transform from Earth-fixed to ecliptic J2000 coordinates
    # At this point, et is guaranteed to be a float due to input validation
    mx = spy.pxform('ITRF93', 'ECLIPJ2000', et)  # type: ignore
    v_eclip = spy.mxv(mx, v_E)
    
    return v_eclip


def Kepler(E, M, e):
    return E - e * np.sin(E) - M

@dataclass
class CanonicalUnits:
    """
    Data class to provide canonical units and compute the gravitational parameter mu in canonical units.
    """
    deg: float = field(default=np.pi / 180, init=False)
    AU_m: float = field(default=1.496e11, init=False)  # meters
    M_sun: float = field(default=1.9891e30, init=False)  # kg
    G: float = field(default=6.67430e-11, init=False)  # m^3 / (kg s^2)
    year: float = field(default=365.25 * 24 * 3600, init=False)  # seconds

    @property
    def G_cu(self) -> float:
        # Gravitational constant in canonical units
        return self.G * (1 / self.AU_m) ** 3 * self.M_sun * (self.year) ** 2

    @property
    def mu(self) -> float:
        # Gravitational parameter mu in canonical units
        return self.G_cu

@dataclass
class OrbitalElements:
    a: float
    e: float
    i: float
    Omega: float
    w: float
    E: float | None = None
    M: float | None = None

    def __post_init__(self):
        if self.E is None and self.M is None:
            raise ValueError("Either 'E' (Eccentric Anomaly) or 'M' (Mean Anomaly) must be provided.")
        if self.M is not None and self.E is None:
            self.E = newton(Kepler, 1, args=(self.M, self.e))
        if self.E is not None and self.M is None:
            # Optionally, you could allow both, or raise an error if both are provided
            pass

@dataclass
class GravitationalParameters:
    G: float | None = None
    M: float | None = None
    m: float | None = None
    mu: float | None = None
    
    def __post_init__(self):
        """Validate and calculate missing parameters after initialization"""
        # Check if we have enough information to calculate mu
        if self.mu is None:
            if self.G is not None and self.M is not None:
                self.mu = self.G * self.M
            else:
                raise ValueError("Either 'mu' must be provided, or both 'G' and 'M' must be provided")
        
        # If we have mu but not G or M, we can't calculate them uniquely
        # This is fine - the user can access mu directly
    
    @classmethod
    def from_mu(cls, mu: float):
        """Create instance with only mu parameter"""
        return cls(mu=mu)
    
    @classmethod
    def from_GM(cls, G: float, M: float, m: float | None = None):
        """Create instance with G and M parameters, optionally m"""
        return cls(G=G, M=M, m=m)
    
    @classmethod
    def from_standard(cls, m: float | None = None):
        """Create instance with standard gravitational constant and solar mass"""
        return cls(G=6.67430e-11, M=1.9891e30, m=m)

class OrbitTrasformations:

    @staticmethod
    def sqrt_e(e: float) -> float:
        return (1-e**2)**0.5

    @staticmethod
    def nu(a: float, mu: float) -> float:
        return (mu*a)**0.5
    
    @staticmethod
    def r(a: float, e: float, E: float) -> float:
        return a*(1 - e*np.cos(E))

    @staticmethod
    def h(a: float, e: float, mu: float) -> float:
        return (mu*a*(1-e**2))**0.5

    @staticmethod
    def Matrix_trasformation(i: float, w: float, Omega: float) -> np.ndarray:
        Matrix = np.zeros((3, 3))

        Matrix[0, 0] = np.cos(Omega)*np.cos(w) - np.sin(Omega)*np.cos(i)*np.sin(w)
        Matrix[0, 1] = -np.cos(Omega)*np.sin(w) - np.sin(Omega)*np.cos(i)*np.cos(w)
        Matrix[0, 2] = np.sin(Omega)*np.sin(i)

        Matrix[1, 0] = np.sin(Omega)*np.cos(w) + np.cos(Omega)*np.cos(i)*np.sin(w)
        Matrix[1, 1] = -np.sin(Omega)*np.sin(w) + np.cos(Omega)*np.cos(i)*np.cos(w)
        Matrix[1, 2] = -np.cos(Omega)*np.sin(i)

        Matrix[2, 0] = np.sin(i)*np.sin(w)
        Matrix[2, 1] = np.sin(i)*np.cos(w)
        Matrix[2, 2] = np.cos(i)

        return Matrix


def compute_functions(i: float, w: float, Omega: float, 
                    which: Set[str]) -> Dict[str, float]:
    """Compute transformation functions A, B, C, D, F, G"""
    results = {}
    
    if 'A' in which:
        results['A'] = (np.cos(Omega)*np.cos(w) - np.sin(Omega)*np.cos(i)*np.sin(w))
    if 'B' in which:
        results['B'] = (-np.cos(Omega)*np.sin(w) - np.sin(Omega)*np.cos(i)*np.cos(w))
    if 'C' in which:
        results['C'] = (np.sin(Omega)*np.cos(w) + np.cos(Omega)*np.cos(i)*np.sin(w))
    if 'D' in which:
        results['D'] = (-np.sin(Omega)*np.sin(w) + np.cos(Omega)*np.cos(i)*np.cos(w))
    if 'F' in which:
        results['F'] = np.sin(w) * np.sin(i)
    if 'G' in which:
        results['G'] = np.cos(w) * np.sin(i)
    
    return results

def compute_state_vector(E: np.array, mu: float) -> np.array:
    state_vector = spy.conics(E+[0, mu], 0)
    return state_vector

def trasformation_X_to_E(x: float, y: float, z: float, vx: float, vy: float, vz: float, mu: float) -> tuple[float, float, float, float, float, float]:

    elements = spy.oscelt([x, y, z, vx, vy, vz], et=0, mu=mu)
    q = elements[0]
    e = elements[1]
    i = elements[2]
    Omega = elements[3]
    w = elements[4]
    M = elements[5]
    a = q/(1-e)

    return q, e, i, Omega, w, M, a

def trasformation_E_to_X(q: float, e: float, i: float, Omega: float, w: float, M: float, mu: float) -> tuple[float, float, float, float, float, float]:

    state_vec = spy.conics([q, e, i, Omega, w, M]+[0, mu], 0)
    x = state_vec[0]
    y = state_vec[1]
    z = state_vec[2]
    vx = state_vec[3]
    vy = state_vec[4]
    vz = state_vec[5]

    return x, y, z, vx, vy, vz

def compute_jacobian_XoE(a: float, e: float, i: float, Omega: float, w: float, M: float, mu: float) -> np.array:
    E = newton(Kepler, M, args=(M, e))

    functions = compute_functions(i, w, Omega, {'A', 'B', 'C', 'D', 'F', 'G'})
    A = functions['A']
    B = functions['B']
    C = functions['C']
    D = functions['D']
    F = functions['F']
    G = functions['G']
    
    ot = OrbitTrasformations()
    r = ot.r(a, e, E)
    eps = ot.sqrt_e(e)
    nu = ot.nu(a, mu)
    nur = nu/r

    q = a*(1-e)
    state_vector = compute_state_vector([q, e, i, Omega, w, M], mu)

    #partial a

    partial_a_x = state_vector[0]/a
    partial_a_y = state_vector[1]/a
    partial_a_z = state_vector[2]/a
    partial_a_vx = -state_vector[3]/(2*a)
    partial_a_vy = -state_vector[4]/(2*a)
    partial_a_vz = -state_vector[5]/(2*a)

    partial_a = [partial_a_x, partial_a_y, partial_a_z, partial_a_vx, partial_a_vy, partial_a_vz]

    #partial e

    #dX/de
    dcosEde=-a*np.sin(E)**2/r       
    dsinEde=a*np.cos(E)*np.sin(E)/r
    dnurde=(nu*a/r**2)*(np.cos(E)-(a/r)*e*np.sin(E)**2)
    depsde=-e/eps

    drAde=a*(dcosEde-1)
    drBde=a*(depsde*np.sin(E)+eps*dsinEde)

    dvAde=-(dnurde*np.sin(E)+nur*dsinEde)
    dvBde=(dnurde*eps*np.cos(E)+nur*depsde*np.cos(E)+nur*eps*dcosEde)

    partial_e = np.array([
        drAde*A+drBde*B,
        drAde*C+drBde*D,
        drAde*F+drBde*G,
        dvAde*A+dvBde*B,
        dvAde*C+dvBde*D,
        dvAde*F+dvBde*G
    ])

    #partial i
    partial_i_x = state_vector[2]*np.sin(Omega)
    partial_i_y = -state_vector[2]*np.cos(Omega)
    partial_i_z = -state_vector[0]*np.sin(Omega) + state_vector[1]*np.cos(Omega)

    partial_i_vx = state_vector[5]*np.sin(Omega)
    partial_i_vy = -state_vector[5]*np.cos(Omega)
    partial_i_vz = -state_vector[3]*np.sin(Omega) + state_vector[4]*np.cos(Omega)

    partial_i = [partial_i_x, partial_i_y, partial_i_z, partial_i_vx, partial_i_vy, partial_i_vz]

    #partial Omega
    partial_Omega_x = -state_vector[1]
    partial_Omega_y = state_vector[0]
    partial_Omega_z = 0

    partial_Omega_vx = -state_vector[4]
    partial_Omega_vy = state_vector[3]
    partial_Omega_vz = 0

    partial_Omega = [partial_Omega_x, partial_Omega_y, partial_Omega_z, partial_Omega_vx, partial_Omega_vy, partial_Omega_vz]

    #partial w

    partial_w_x = -state_vector[1]*np.cos(i) - state_vector[2]*np.sin(i)*np.cos(Omega)
    partial_w_y = state_vector[0]*np.cos(i) - state_vector[2]*np.sin(i)*np.sin(Omega)
    partial_w_z = state_vector[0]*np.sin(i)*np.cos(Omega) + state_vector[1]*np.sin(i)*np.sin(Omega)
    partial_w_vx = -state_vector[4]*np.cos(i) - state_vector[5]*np.sin(i)*np.cos(Omega)
    partial_w_vy = state_vector[3]*np.cos(i) - state_vector[5]*np.sin(i)*np.sin(Omega)
    partial_w_vz = state_vector[3]*np.sin(i)*np.cos(Omega) + state_vector[4]*np.sin(i)*np.sin(Omega)

    partial_w = [partial_w_x, partial_w_y, partial_w_z, partial_w_vx, partial_w_vy, partial_w_vz]

    #partial M
    n = (mu/a**3)**0.5
    factor = -(mu*a**3)**0.5/r**3

    partial_M_x = (1/n) * state_vector[3]
    partial_M_y = (1/n) * state_vector[4]
    partial_M_z = (1/n) * state_vector[5]
    partial_M_vx = factor * state_vector[0]
    partial_M_vy = factor * state_vector[1]
    partial_M_vz = factor * state_vector[2]

    partial_M = [partial_M_x, partial_M_y, partial_M_z, partial_M_vx, partial_M_vy, partial_M_vz]

    J = np.zeros((6,6))
    J[:,0] = partial_a
    J[:,1] = partial_e
    J[:,2] = partial_i
    J[:,3] = partial_Omega
    J[:,4] = partial_w
    J[:,5] = partial_M

    Je2c=np.eye(6)
    Je2c[0,0]=1/(1-e)
    Je2c[0,1]=q/(1-e)**2
    JX2c=np.matmul(J,Je2c)
    
    return JX2c

def Jacobian_qei_to_QEI(q: float, e: float, i: float, q_max: float, e_max: float, i_max: float) -> np.array:
    partialQ_q = q_max/(q*(q_max - q))
    partialQ_e = 0
    partialQ_i = 0
    partialE_q = 0
    partialE_e = e_max/(e*(e_max - e))
    partialE_i = 0
    partialI_q = 0
    partialI_e = 0
    partialI_i = i_max/(i*(i_max - i))

    Jacobian = np.array([[partialQ_q, partialQ_e, partialQ_i], [partialE_q, partialE_e, partialE_i], [partialI_q, partialI_e, partialI_i]])
    return Jacobian

def P_E_CMND(q: float, e: float, i: float, F: mm.FitCMND) -> float:
    max = 2*np.pi; min = 0

    element = np.array([q, e, i])
    scales=[1.35,1.00,np.pi]
    u_element = mm.Util.tIF(element, scales, mm.Util.f2u)

    P_qei = F.cmnd.pdf(u_element)
    P_WwM = 1/(max - min)**3

    return P_qei * P_WwM

def P_X_CMND(x: np.ndarray, y: np.ndarray, z: np.ndarray, vx: np.ndarray, vy: np.ndarray, vz: np.ndarray, 
              q_max: float, e_max: float, i_max: float, mu: float, F: mm.FitCMND) -> np.ndarray:
    """
    Compute probability density function of NEAs in phase space (X).
    
    Vectorized version: x, y, z, vx, vy, vz are arrays (or scalars).
    Returns array of probability values.
    
    Parameters:
    -----------
    x, y, z : np.ndarray
        Position coordinates (arrays or scalars)
    vx, vy, vz : np.ndarray
        Velocity coordinates (arrays or scalars)
    q_max : float
        Maximum perihelion distance for NEA classification (AU)
    e_max : float
        Maximum eccentricity for valid orbits
    i_max : float
        Maximum inclination
    mu : float
        Gravitational parameter
    F : mm.FitCMND
        CMND fit model for orbital elements PDF
        
    Returns:
    --------
    np.ndarray
        Probability density values, same shape as broadcasted input arrays
    """
    # Constants for NEA orbit validation
    Q_MAX_NEA = 1.3  # Maximum perihelion distance for NEA classification (AU)
    E_MAX_VALID = 1.0  # Maximum eccentricity for valid elliptical orbits
    
    # Convert inputs to arrays and get broadcast shape
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)
    vx = np.asarray(vx)
    vy = np.asarray(vy)
    vz = np.asarray(vz)
    
    shape = np.broadcast(x, y, z, vx, vy, vz).shape
    P = np.empty(shape, dtype=float)

    # Flatten arrays for iteration
    x_flat = x.ravel()
    y_flat = y.ravel()
    z_flat = z.ravel()
    vx_flat = vx.ravel()
    vy_flat = vy.ravel()
    vz_flat = vz.ravel()

    # Process each point in the flattened arrays
    for idx in range(x_flat.size):
        # Transform from Cartesian to orbital elements
        q, e, i, Omega, w, M, a = trasformation_X_to_E(
            x_flat[idx], y_flat[idx], z_flat[idx], 
            vx_flat[idx], vy_flat[idx], vz_flat[idx], mu
        )
        
        # Validate orbit: must be NEA (q <= Q_MAX_NEA) and elliptical (e < E_MAX_VALID)
        if q > Q_MAX_NEA or e >= E_MAX_VALID:
            P.flat[idx] = 0.0
            continue

        # Compute Jacobian matrices
        jacobian_X_to_E = compute_jacobian_XoE(a, e, i, Omega, w, M, mu)
        jacobian_qei_to_QEI = Jacobian_qei_to_QEI(q, e, i, q_max, e_max, i_max)

        #compute determinants
        det_X_to_E = np.linalg.det(jacobian_X_to_E)
        det_QEI = jacobian_qei_to_QEI[0,0] * jacobian_qei_to_QEI[1,1] * jacobian_qei_to_QEI[2,2]
        
        # Check for computational inconsistencies
        if det_X_to_E == 0 or not np.isfinite(det_X_to_E):
            P.flat[idx] = np.nan
            continue
        
        # Compute probability density
        inv_det_X_to_E = 1.0 / det_X_to_E
        pdf_orbital_elements = P_E_CMND(q, e, i, F)
        P.flat[idx] = pdf_orbital_elements * abs(inv_det_X_to_E) * abs(det_QEI)
                    
    return P.reshape(shape)

def hypercube_surface_integral_P_X_CMND(center, widths, max_elements: list, n_points=8, mu=1, F=mm.FitCMND):
    """
    Calculate the surface integral of P_xyvxvy in a hypercube centered at (x, y, vx, vy)
    with dimensions (dx, dy, dvx, dvy) using Gauss-Legendre quadrature.

    Parameters:
        center: tuple/list/array of (x, y, vx, vy) center
        widths: tuple/list/array of (dx, dy, dvx, dvy) side lengths
        max_elemetns: tuple/list/array of (a_max, e_max, i_max) max elements
        n_points: number of quadrature points per dimension

    Returns:
        Integral (float)
    """
    from numpy.polynomial.legendre import leggauss

    x0, y0, z0, vx0, vy0, vz0 = center
    dx, dy, dz, dvx, dvy, dvz = widths
    q_max, e_max, i_max = max_elements

    # Get Gauss-Legendre points and weights for [-1, 1]
    pts, wts = leggauss(n_points)

    # Map points from [-1, 1] to [center-width/2, center+width/2] for each dimension
    x_pts = x0 + 0.5*dx*pts
    y_pts = y0 + 0.5*dy*pts
    z_pts = z0 + 0.5*dz*pts
    vx_pts = vx0 + 0.5*dvx*pts
    vy_pts = vy0 + 0.5*dvy*pts
    vz_pts = vz0 + 0.5*dvz*pts

    # Create meshgrid of all quadrature points
    X, Y, Z, VX, VY, VZ = np.meshgrid(x_pts, y_pts, z_pts, vx_pts, vy_pts, vz_pts, indexing='ij')
    WX, WY, WZ, WVX, WVY, WVZ = np.meshgrid(wts, wts, wts, wts, wts, wts, indexing='ij')

    # Flatten for vectorized evaluation
    Xf = X.ravel()
    Yf = Y.ravel()
    Zf = Z.ravel()
    VXf = VX.ravel()
    VYf = VY.ravel()
    VZf = VZ.ravel()
    WF = (WX * WY * WZ * WVX * WVY * WVZ).ravel()
    # Evaluate P at all points
    Pf = P_X_CMND(Xf, Yf, Zf, VXf, VYf, VZf, q_max, e_max, i_max, mu, F=F)

    # Integral is sum(P * weight) * volume factor
    integral = np.sum(Pf * WF) * (0.5*dx) * (0.5*dy) * (0.5*dz) * (0.5*dvx) * (0.5*dvy) * (0.5*dvz)
    #return integral, y_points
    return integral

def box_surface_integral_P_X_CMND(center, widths, fixed_values, max_elements: list, integrate_over='velocity', n_points=8, mu=1, F=mm.FitCMND):
    """
    Generalized function to calculate the integral of P_X_CMND over a box in either position or velocity space,
    while keeping the other fixed.
    
    Parameters:
        center: tuple/list/array of (c1, c2, c3) center of the integration box
        widths: tuple/list/array of (d1, d2, d3) side lengths of the integration box
        fixed_values: tuple/list/array of (f1, f2, f3) fixed values for the non-integrated dimensions
        max_elements: tuple/list/array of (q_max, e_max, i_max) max elements
        integrate_over: 'position' or 'velocity' - specifies which dimensions to integrate over
        n_points: number of quadrature points per dimension
        mu: gravitational parameter
        F: FitCMND object
    
    Returns:
        Integral (float)
    
    Examples:
        # Integrate over velocity box with fixed position (equivalent to v_box_surface_integral_P_X_CMND)
        integral = box_surface_integral_P_X_CMND(
            center=(vx_c, vy_c, vz_c), 
            widths=(dvx, dvy, dvz),
            fixed_values=(x, y, z),  # fixed position
            max_elements=[q_max, e_max, i_max],
            integrate_over='velocity'
        )
        
        # Integrate over position box with fixed velocity (equivalent to r_box_surface_integral_P_X_CMND)
        integral = box_surface_integral_P_X_CMND(
            center=(x_c, y_c, z_c), 
            widths=(dx, dy, dz),
            fixed_values=(vx, vy, vz),  # fixed velocity
            max_elements=[q_max, e_max, i_max],
            integrate_over='position'
        )
    """
    from numpy.polynomial.legendre import leggauss
    
    if integrate_over not in ['position', 'velocity']:
        raise ValueError("integrate_over must be 'position' or 'velocity'")
    
    c1, c2, c3 = center
    d1, d2, d3 = widths
    f1, f2, f3 = fixed_values
    q_max, e_max, i_max = max_elements
    
    # Get Gauss-Legendre points and weights for [-1, 1]
    pts, wts = leggauss(n_points)
    
    # Map points from [-1, 1] to [center-width/2, center+width/2] for integration dimensions
    dim1_pts = c1 + 0.5*d1*pts
    dim2_pts = c2 + 0.5*d2*pts
    dim3_pts = c3 + 0.5*d3*pts
    
    # Create meshgrid for integration dimensions
    DIM1, DIM2, DIM3 = np.meshgrid(dim1_pts, dim2_pts, dim3_pts, indexing='ij')
    W1, W2, W3 = np.meshgrid(wts, wts, wts, indexing='ij')
    
    # Flatten for vectorized evaluation
    DIM1f = DIM1.ravel()
    DIM2f = DIM2.ravel()
    DIM3f = DIM3.ravel()
    WF = (W1 * W2 * W3).ravel()
    
    # Create constant arrays for fixed dimensions (same size as integration arrays)
    F1f = np.full_like(DIM1f, f1)
    F2f = np.full_like(DIM2f, f2)
    F3f = np.full_like(DIM3f, f3)
    
    # Assign to position and velocity based on integrate_over parameter
    if integrate_over == 'velocity':
        # Integrating over velocities, positions are fixed
        Xf, Yf, Zf = F1f, F2f, F3f  # fixed positions
        VXf, VYf, VZf = DIM1f, DIM2f, DIM3f  # integrated velocities
    else:  # integrate_over == 'position'
        # Integrating over positions, velocities are fixed
        Xf, Yf, Zf = DIM1f, DIM2f, DIM3f  # integrated positions
        VXf, VYf, VZf = F1f, F2f, F3f  # fixed velocities
    
    # Evaluate P at all points
    Pf = P_X_CMND(Xf, Yf, Zf, VXf, VYf, VZf, q_max, e_max, i_max, mu, F=F)
    
    # Integral is sum(P * weight) * volume factor
    integral = np.sum(Pf * WF) * (0.5*d1) * (0.5*d2) * (0.5*d3)
    return integral