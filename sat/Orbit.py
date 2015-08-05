from __future__ import division
from . _constants.Angle import REV2RAD
from . _constants.Orbit import DAY2SEC
from . _constants.Orbit import EARTH_ANGULAR_SPEED
from . _constants.Orbit import EARTH_FLATTENING_FACTOR
from . _constants.Orbit import EARTH_SEMIMAJOR_AXIS
from . _constants.Orbit import ECCENTRIC_ANOMALY_TOLERANCE
from . _constants.Orbit import GEODETIC_COORDINATES_TOLERANCE
from . _constants.Orbit import SECOND_ZONAL_TERM_J2
from . _constants.Orbit import STANDARD_GRAVITATIONAL_PARAMETER
from . _decorators import accepts
from . _decorators import returns
from . _errors.Orbit import OrbitError
from . Ephemeris import Ephemeris
import numpy as np


class Orbit(object):

    _properties = [
        "_argument_of_perigee",
        "_eccentric_anomaly",
        "_ephemeris",
        "_longitude_of_the_ascending_node",
        "_mean_anomaly",
        "_mean_motion",
        "_satellite_position_ecef",
        "_satellite_position_geod",
        "_semimajor_axis",
        "_time",
        "_true_anomaly",
        ]

    def __init__(self, ephemeris):

        for item in self._properties:
            self.__setattr__("_{}".format(item), None)

        if ephemeris:
            self._ephemeris = ephemeris.to_copy()
        else:
            raise OrbitError("Ephemeris instance is not complete")
        pass

    @property
    @returns(np.ndarray)
    def argument_of_perigee(self):
        """argument of perigee in radians"""
        return self._argument_of_perigee

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
    def satellite_position_ecef(self):
        """satellite position in ECEF reference system"""
        return self._satellite_position_ecef

    @property
    @returns(np.ndarray)
    def satellite_position_geod(self):
        """satellite position in geodetic reference system"""
        return self._satellite_position_geod

    @property
    @returns(np.ndarray)
    def semimajor_axis(self):
        """semimajor axis of the elliptic orbit in meters"""
        return self._semimajor_axis

    @property
    @returns(np.ndarray)
    def time(self):
        """instant of time with respect to the epoch, in days"""
        return self._time

    @time.setter
    @accepts(np.ndarray)
    def time(self, val):
        try:
            epoch = np.datetime64(self.ephemeris.epoch_datetime)
            delta = np.datetime64(val) - epoch
            self._time = np.asarray(delta / np.timedelta(1, "D")).astype(float)
        except TypeError:
            msg = "Ephemeris attribute from Orbit instance is not set"
            err = OrbitError(msg)
            err.__cause__ = None
            raise err

    @property
    @returns(np.ndarray)
    def true_anomaly(self):
        """true anomaly in radians"""
        return self._true_anomaly

    def _calc_mean_motion(self):

        dt = self.time
        n0 = self.ephemeris.mean_motion
        n1 = self.ephemeris.mean_motion_first_dif
        n2 = self.ephemeris.mean_motion_second_dif

        self._mean_motion = np.asarray(n0 + n1 * dt + 0.5 * n2 * dt**2)

    def _calc_mean_anomaly(self):

        dt = self.time
        m0 = self.ephemeris.mean_anomaly.rad
        n = self.mean_motion * REV2RAD

        self._mean_anomaly = np.asarray(m0 + n * dt)

    def _calc_eccentric_anomaly(self):

        eccent = self.ephemeris.eccentricity
        m = self.mean_anomaly
        ean_old = m
        dif = 1

        while (dif > 0).all():
            ean_new = np.where(dif > 0, m + eccent*np.sin(ean_old), ean_old)
            dif = np.abs(ean_new - ean_old) - ECCENTRIC_ANOMALY_TOLERANCE
            ean_old = ean_new

        self._eccentric_anomaly = ean_old

    def _calc_true_anomaly(self):

        eccent = self.ephemeris.eccentricity
        factor = np.sqrt((1+eccent)/(1-eccent))
        ean = self.eccentric_anomaly

        self._true_anomaly = 2 * np.arctan(factor * np.tan(ean/2))

    def _calc_semimajor_axis(self):

        mu = STANDARD_GRAVITATIONAL_PARAMETER
        n = self.mean_motion * REV2RAD / DAY2SEC

        self._semimajor_axis = (mu / n**2)**(1/3)

    def _calc_argument_of_perigee_and_longitude_of_the_ascending_node(self):

        j2 = SECOND_ZONAL_TERM_J2
        mu = STANDARD_GRAVITATIONAL_PARAMETER
        we = EARTH_ANGULAR_SPEED

        dt = self.time
        i = self.ephemeris.inclination
        n = self.mean_motion * REV2RAD

        eccent = self.ephemeris.eccentricity
        a = self.semimajor_axis
        p = a * (1-eccent**2)

        factor = -REV2RAD * j2 / (mu*p**2) * 3 * n
        w0 = self.ephemeris.argument_of_perigee
        w1 = factor * ((5/4)*np.sin(i)**2 - 1)
        omega0 = self.ephemeris.longitude_of_the_ascending_node
        omega1 = (factor/2) * np.cos(i) - we

        self._argument_of_perigee = np.asarray(w0 + w1 * dt)
        self._longitude_of_the_ascending_node = np.asarray(omega0 + omega1 * dt)

    def _calc_cartesian_and_geodetic_coordinates(self):

        i = self.ephemeris.inclination
        o = self.longitude_of_the_ascending_node
        w = self.argument_of_perigee
        v = self.true_anomaly

        eccent = self.ephemeris.eccentricity
        a = self.semimajor_axis
        p = a * (1 - eccent**2)

        r_s = p / (1 + eccent * np.cos(v))

        rx_s = r_s * (np.cos(o)*np.cos(w+v) - np.sin(o)*np.sin(w+v)*np.cos(i))
        ry_s = r_s * (np.sin(o)*np.cos(w+v) + np.cos(o)*np.sin(w+v)*np.cos(i))
        rz_s = r_s * (np.sin(w+v)*np.sin(i))

        self._satellite_position_ecef = np.asarray([rx_s, ry_s, rz_s]).T

        ae = EARTH_SEMIMAJOR_AXIS
        e2 = EARTH_FLATTENING_FACTOR
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

        n_factor_old = ae
        alt_old = 0
        lat_old = _estimate_latitude(n_factor_old, alt_old)
        dif = 1

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

            n_factor_old = n_factor_new
            alt_old = alt_new
            lat_old = lat_new

        alt = alt_old
        lat = lat_old
        lon = _estimate_longitude()

        self._satellite_position_geod = np.asarray([lat, lon, alt]).T
