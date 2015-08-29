from __future__ import division
from . constants.Angle import REV2RAD
from . constants.Coordinates import EARTH_ANGULAR_SPEED
from . constants.Orbit import DAY2SEC
from . constants.Orbit import ECCENTRIC_ANOMALY_TOLERANCE
from . constants.Orbit import SECOND_ZONAL_TERM_J2
from . constants.Orbit import STANDARD_GRAVITATIONAL_PARAMETER
from . decorators import accepts
from . decorators import returns
from . Coordinates import Coordinates
from . Ephemeris import Ephemeris
from . Error import OrbitError
import numpy as np


class Orbit(object):

    __slots__ = [
        "_argument_of_perigee",
        "_datetime",
        "_eccentric_anomaly",
        "_ephemeris",
        "_longitude_of_the_ascending_node",
        "_mean_anomaly",
        "_mean_motion",
        "_position_ecf",
        "_position_eci",
        "_position_geo",
        "_semimajor_axis",
        "_timedelta",
        "_true_anomaly",
        "_velocity_ecf",
        "_velocity_eci",
        ]

    @accepts(Ephemeris)
    def __init__(self, ephemeris):
        """Constructor of a generic Orbit instance.

        Parameters:

        ephemeris
            Ephemeris instance containing the Kepler parameters
        """

        for item in self.__slots__:
            self.__setattr__(item, None)

        if ephemeris:
            self._ephemeris = ephemeris.to_copy()
        else:
            msg = "Ephemeris instance is not complete"
            raise OrbitError(msg)

    @accepts(np.ndarray)
    def compute(self, datetime):
        """Compute elliptical orbit for specified datetime.

        Parameters:

        datetime
            datetime for predicted elliptical orbit parameters
        """

        sup = self.__class__
        sup._add_datetime(self, datetime)
        sup._calc_timedelta(self)
        sup._calc_mean_motion(self)
        sup._calc_mean_anomaly(self)
        sup._calc_eccentric_anomaly(self)
        sup._calc_true_anomaly(self)
        sup._calc_semimajor_axis(self)
        sup._calc_argument_of_perigee_and_longitude_of_the_ascending_node(self)
        sup._calc_coordinates_from_orb_to_eci(self)
        sup._calc_coordinates_from_eci_to_ecf(self)
        sup._calc_coordinates_from_ecf_to_geo(self)

    @classmethod
    def _add_datetime(cls, obj, date):
        """Safe setter for datetime attribute."""

        try:
            obj._datetime = np.asarray([
                x.strftime("%Y-%m-%dT%H:%M:%S.%fZ")
                for x in date.ravel()], dtype="datetime64")[:, None]
        except (TypeError, ValueError):
            msg = "Invalid array for input time"
            err = OrbitError(msg)
            err.__cause__ = None
            raise err

    @classmethod
    def _calc_timedelta(cls, obj):
        """Compute time difference with respect to epoch."""

        try:
            fmt = "%Y-%m-%dT%H:%M:%S.%fZ"
            epoch = np.asarray([
                obj.ephemeris.epoch_datetime.strftime(fmt)
                ]).astype("datetime64")
            delta = (obj.datetime - epoch)
            obj._timedelta = (delta / np.timedelta64(1, "D")).astype(float)
        except TypeError:
            msg = "Ephemeris attribute from Orbit instance is not complete"
            err = OrbitError(msg)
            err.__cause__ = None
            raise err

    @classmethod
    def _calc_mean_motion(cls, obj):
        """Compute mean motion for a specific datetime."""

        # Call necessary properties.
        dt = obj.timedelta
        n0 = obj.ephemeris.mean_motion
        n1 = obj.ephemeris.mean_motion_first_dif
        # Set mean motion property.
        obj._mean_motion = n0 + n1 * dt

    @classmethod
    def _calc_mean_anomaly(cls, obj):
        """Compute mean anomaly for a specific datetime."""

        # Call necessary properties.
        dt = obj.timedelta
        m0 = obj.ephemeris.mean_anomaly.rad
        n = obj.mean_motion * REV2RAD
        # Set mean anomaly property.
        obj._mean_anomaly = m0 + n * dt

    @classmethod
    def _calc_eccentric_anomaly(cls, obj):
        """Compute eccentric anomaly for a specific datetime."""

        # Call necessary properties.
        eccent = obj.ephemeris.eccentricity
        m = obj.mean_anomaly
        ean_old = m
        dif = np.asarray([1])
        # Solve Kepler's equation by recursive methods until a tolerance
        # is reached.
        while (dif > 0).all():
            ean_new = np.where(dif > 0, m + eccent*np.sin(ean_old), ean_old)
            dif = np.abs(ean_new - ean_old) - ECCENTRIC_ANOMALY_TOLERANCE
            ean_old = ean_new
        # Set eccentric anomaly property.
        obj._eccentric_anomaly = ean_old

    @classmethod
    def _calc_true_anomaly(cls, obj):
        """Compute true anomaly for a specific datetime."""

        # Call necessary properties.
        eccent = obj.ephemeris.eccentricity
        factor = np.sqrt((1+eccent)/(1-eccent))
        ean = obj.eccentric_anomaly
        # Set true anomaly property.
        obj._true_anomaly = 2 * np.arctan(factor * np.tan(ean/2))

    @classmethod
    def _calc_semimajor_axis(cls, obj):
        """Compute semimajor axis for a specific datetime."""

        # Call necessary constants and properties.
        mu = STANDARD_GRAVITATIONAL_PARAMETER
        n = obj.mean_motion * REV2RAD / DAY2SEC
        # Set semimajor axis property using Kepler's third law.
        obj._semimajor_axis = (mu / n**2)**(1/3)

    @classmethod
    def _calc_argument_of_perigee_and_longitude_of_the_ascending_node(cls, obj):
        """Compute argument of perigee and longitude of the ascending
        node for a specific datetime."""

        # Call necessary constants.
        j2 = SECOND_ZONAL_TERM_J2
        mu = STANDARD_GRAVITATIONAL_PARAMETER
        # Call necessary properties.
        dt = obj.timedelta
        i = obj.ephemeris.inclination.rad
        n = obj.mean_motion
        eccent = obj.ephemeris.eccentricity
        a = obj.semimajor_axis
        p = a * (1-eccent**2)
        # Compute first derivatives of argument of perigee and longitude
        # of the ascending node.
        factor = -REV2RAD * j2 / (mu*p**2) * 3 * n
        w0 = obj.ephemeris.argument_of_perigee.rad
        w1 = factor * ((5/4)*np.sin(i)**2 - 1)
        o0 = obj.ephemeris.longitude_of_the_ascending_node.rad
        o1 = (factor/2) * np.cos(i)
        # Set argument of perigee property.
        obj._argument_of_perigee = w0 + w1 * dt
        # Set longitude of the ascending node property.
        obj._longitude_of_the_ascending_node = o0 + o1 * dt

    @classmethod
    def _calc_coordinates_from_orb_to_eci(cls, obj):
        """Compute ECI coordinates for a specific datetime."""

        # Call necessary constants.
        mu = STANDARD_GRAVITATIONAL_PARAMETER
        # Call necessary properties.
        i = obj.ephemeris.inclination.rad
        o = obj.longitude_of_the_ascending_node
        w = obj.argument_of_perigee
        v = obj.true_anomaly
        eccent = obj.ephemeris.eccentricity
        a = obj.semimajor_axis
        # Compute intermediate factors.
        p = a * (1 - eccent**2)
        r = p / (1 + eccent * np.cos(v))
        l = np.sqrt(mu * p)
        # Compute (X,Y,Z) coordinates in ECI reference system.
        rx = r * (np.cos(o)*np.cos(w+v) - np.sin(o)*np.sin(w+v)*np.cos(i))
        ry = r * (np.sin(o)*np.cos(w+v) + np.cos(o)*np.sin(w+v)*np.cos(i))
        rz = r * (np.sin(w+v)*np.sin(i))
        # Compute (Vx,Vy,Vz) velocities in ECI reference system.
        vx = (rx * l * eccent)/(r * p) * np.sin(v) - l/r * (
            np.cos(o)*np.sin(w+v) + np.sin(o)*np.cos(w+v)*np.cos(i))
        vy = (ry * l * eccent)/(r * p) * np.sin(v) - l/r * (
            np.sin(o)*np.sin(w+v) - np.cos(o)*np.cos(w+v)*np.cos(i))
        vz = (rz * l * eccent)/(r * p) * np.sin(v) + l/r * (
            np.cos(w+v)*np.sin(i))
        # Set ECI satellite position and velocity properties.
        obj._position_eci = np.hstack([rx, ry, rz])
        obj._velocity_eci = np.hstack([vx, vy, vz])

    @classmethod
    def _calc_coordinates_from_eci_to_ecf(cls, obj):
        """Compute ECF coordinates for a specific datetime."""

        rx, ry, rz = [x[:, None] for x in obj.position_eci.T]
        vx, vy, vz = [x[:, None] for x in obj.velocity_eci.T]
        vx += EARTH_ANGULAR_SPEED * ry
        vy -= EARTH_ANGULAR_SPEED * rx

        obj._position_ecf = Coordinates.from_eci_to_ecf(
            np.hstack([rx, ry, rz]), obj.datetime)
        obj._velocity_ecf = Coordinates.from_eci_to_ecf(
            np.hstack([vx, vy, vz]), obj.datetime)

    @classmethod
    def _calc_coordinates_from_ecf_to_geo(cls, obj):
        """Compute geodetic coordinates for a specific datetime."""
        obj._position_geo = Coordinates.from_ecf_to_geo(obj.position_ecf)

    @classmethod
    def _calc_coordinates_from_geo_to_ecf(cls, obj):
        """Compute ECF coordinates using geodetic coordinates."""
        obj._position_ecf = Coordinates.from_geo_to_ecf(obj.position_geo)

    @property
    @returns(np.ndarray)
    def argument_of_perigee(self):
        """argument of perigee in radians"""
        return self._argument_of_perigee

    @property
    @returns(np.ndarray)
    def datetime(self):
        """datetime for predicted elliptical orbit parameters"""
        return self._datetime

    @property
    @returns(np.ndarray)
    def eccentric_anomaly(self):
        """eccentric anomaly in radians"""
        return self._eccentric_anomaly

    @property
    @returns(Ephemeris)
    def ephemeris(self):
        """Ephemeris object containing the reference orbital parameters"""
        return self._ephemeris

    @property
    @returns(np.ndarray)
    def longitude_of_the_ascending_node(self):
        """longitude of the ascending node in radians"""
        return self._longitude_of_the_ascending_node

    @property
    @returns(np.ndarray)
    def mean_anomaly(self):
        """mean anomaly in radians"""
        return self._mean_anomaly

    @property
    @returns(np.ndarray)
    def mean_motion(self):
        """mean motion in revs per day"""
        return self._mean_motion

    @property
    @returns(np.ndarray)
    def position_ecf(self):
        """satellite position in ECF reference system"""
        return self._position_ecf

    @property
    @returns(np.ndarray)
    def position_eci(self):
        """satellite position in ECI reference system"""
        return self._position_eci

    @property
    @returns(np.ndarray)
    def position_geo(self):
        """satellite position in geodetic reference system"""
        return self._position_geo

    @property
    @returns(np.ndarray)
    def velocity_ecf(self):
        """satellite velocity in ECF reference system"""
        return self._velocity_ecf

    @property
    @returns(np.ndarray)
    def velocity_eci(self):
        """satellite velocity in ECI reference system"""
        return self._velocity_eci

    @property
    @returns(np.ndarray)
    def semimajor_axis(self):
        """semimajor axis of the elliptic orbit in meters"""
        return self._semimajor_axis

    @property
    @returns(np.ndarray)
    def timedelta(self):
        """increment of time since the epoch, in days"""
        return self._timedelta

    @property
    @returns(np.ndarray)
    def true_anomaly(self):
        """true anomaly in radians"""
        return self._true_anomaly
