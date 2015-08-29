from __future__ import division
from . constants.Angle import REV2RAD
from . constants.Coordinates import EARTH_FLATTENING_FACTOR
from . constants.Coordinates import EARTH_SEMIMAJOR_AXIS
from . constants.Coordinates import GEODETIC_COORDINATES_TOLERANCE
from . constants.Orbit import DAY2SEC
from . constants.Orbit import EARTH_ANGULAR_SPEED
from . constants.Orbit import ECCENTRIC_ANOMALY_TOLERANCE
from . constants.Orbit import SECOND_ZONAL_TERM_J2
from . constants.Orbit import STANDARD_GRAVITATIONAL_PARAMETER
from . decorators import accepts
from . decorators import returns
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
                for x in date.ravel()
                ]).astype("datetime64")
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
            delta = (obj.datetime - epoch)[:, None]
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
    def _calc_greenwich_mean_sidereal_time(cls, obj):
        """Compute Greenwich mean sidereal time for an Orbit instance."""

        # Call necessary constants.
        we = EARTH_ANGULAR_SPEED
        t1 = np.datetime64("2000-01-01T12:00:00Z")
        t2 = np.datetime64("2000-01-01T00:00:00Z")
        # Compute Julian centuries since epoch 2000/01/01 12:00:00.
        dtm1 = ((obj.datetime - t1) / np.timedelta64(1, "D")).astype(float)
        dtm1 = dtm1[:, None] / 36525
        # Compute UTC time for the datetime of interest.
        dtm2 = ((obj.datetime - t2) / np.timedelta64(1, "D")).astype(float)
        dtm2 = ((dtm2 % dtm2.astype(int)) * 86400)[:, None]
        # Compute the Greenwich Sidereal Time (GST), that is, the angle
        # between the vernal point and the Greenwich meridian (which is
        # also the angle between the ECI and ECEF reference systems).
        c0, c1, c2, c3 = [24110.54841, 8640184.812866, 0.093104, -6.2e-6]
        gst0 = c0 + c1 * dtm1 + c2 * dtm1**2 + c3 * dtm1**3
        gst = (gst0/86400 * 2*np.pi) + we * dtm2 / 1.00278

        return gst

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

        # Call necessary properties.
        rx, ry, rz = [x[:, None] for x in obj.position_eci.T]
        vx, vy, vz = [x[:, None] for x in obj.velocity_eci.T]
        # Compute the Greenwich Sidereal Time (GST), that is, the angle
        # between the vernal point and the Greenwich meridian (which is
        # also the angle between the ECI and ECF reference systems).
        gst = cls._calc_greenwich_mean_sidereal_time(obj)
        # Compute ECF coordinates as a rotation of ECI coordinates.
        rx_s = + np.cos(gst)*rx + np.sin(gst)*ry
        ry_s = - np.sin(gst)*rx + np.cos(gst)*ry
        rz_s = rz
        # Compute (Vx,Vy,Vz) velocities in ECF reference system.
        vx += EARTH_ANGULAR_SPEED*ry
        vy -= EARTH_ANGULAR_SPEED*rx
        vx_s = + np.cos(gst)*vx + np.sin(gst)*vy
        vy_s = - np.sin(gst)*vx + np.cos(gst)*vy
        vz_s = vz
        # Set ECF satellite position and velocity properties.
        obj._position_ecf = np.hstack([rx_s, ry_s, rz_s])
        obj._velocity_ecf = np.hstack([vx_s, vy_s, vz_s])

    @classmethod
    def _calc_coordinates_from_ecf_to_geo(cls, obj):
        """Compute geodetic coordinates for a specific datetime."""

        # Call necessary constants.
        ae = EARTH_SEMIMAJOR_AXIS
        e2 = EARTH_FLATTENING_FACTOR
        # Call necessary properties.
        rx_s, ry_s, rz_s = [x[:, None] for x in obj.position_ecf.T]
        p_factor = np.sqrt(rx_s**2 + ry_s**2)

        def _estimate_n_factor(latitude):
            return ae / np.sqrt(1 - e2 * np.sin(latitude)**2)

        def _estimate_height(n_factor, latitude):
            return p_factor / np.cos(latitude) - n_factor

        def _estimate_latitude(n_factor, altitude):
            sin_lat = rz_s / (n_factor*(1-e2) + altitude)
            return np.arctan((rz_s + e2*n_factor*sin_lat)/p_factor)

        def _estimate_longitude():
            return np.arctan2(ry_s, rx_s)

        # Estimate recursively the values of altitude and latitude.
        n_factor_old = ae
        alt_old = 0
        lat_old = _estimate_latitude(n_factor_old, alt_old)
        dif = np.asarray([1])
        # Repeat recursive method until a tolerance is reached.
        while (dif > 0).all():
            n_factor_new = np.where(
                dif > 0, _estimate_n_factor(lat_old), n_factor_old)
            alt_new = np.where(
                dif > 0, _estimate_height(n_factor_new, lat_old), alt_old)
            lat_new = np.where(
                dif > 0, _estimate_latitude(n_factor_new, alt_new), lat_old)
            dif1 = np.abs(alt_new - alt_old)
            dif2 = np.abs(lat_new - lat_old)
            dif = dif1 + dif2 - GEODETIC_COORDINATES_TOLERANCE
            # Overwrite old temporary variables with new values.
            n_factor_old = n_factor_new
            alt_old = alt_new
            lat_old = lat_new
        # Compute longitude for a specific time.
        alt = alt_old
        lat = lat_old
        lon = _estimate_longitude()
        # Set geodetic satellite position property.
        obj._position_geo = np.hstack([lat, lon, alt])

    @classmethod
    def _calc_coordinates_from_geo_to_ecf(cls, obj):
        """Compute ECF coordinates using geodetic coordinates."""

        # Call necessary constants.
        ae = EARTH_SEMIMAJOR_AXIS
        e2 = EARTH_FLATTENING_FACTOR
        # Call necessary properties.
        lat, lon, alt = [x[:, None] for x in obj.position_geo.T]

        # Compute ECF coordinates.
        n = ae / np.sqrt(1 - e2 * np.sin(lat)**2)
        rx_s = (n + alt) * np.cos(lat) * np.cos(lon)
        ry_s = (n + alt) * np.cos(lat) * np.sin(lon)
        rz_s = (n * (1 - e2) + alt) * np.sin(lat)
        # Set ECF satellite position property
        obj._position_ecf = np.hstack([rx_s, ry_s, rz_s])

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
