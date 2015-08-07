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
        "argument_of_perigee",
        "datetime",
        "eccentric_anomaly",
        "ephemeris",
        "longitude_of_the_ascending_node",
        "mean_anomaly",
        "mean_motion",
        "satellite_position_ecef",
        "satellite_position_geod",
        "semimajor_axis",
        "timedelta",
        "true_anomaly",
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
    def timedelta(self):
        """increment of time since the epoch, in days"""
        return self._timedelta

    @property
    @returns(np.ndarray)
    def true_anomaly(self):
        """true anomaly in radians"""
        return self._true_anomaly

    @accepts(np.ndarray)
    def compute(self, datetime):
        """Compute elliptical orbit for specified datetime.

        Parameters:

        datetime
            datetime for predicted elliptical orbit parameters
        """

        def _add_datetime(self, datetime):
            try:
                self._datetime = np.asarray([
                    x.strftime("%Y-%m-%dT%H:%M:%S.%fZ")
                    for x in datetime.ravel()
                    ]).astype("datetime64")
            except (TypeError, ValueError):
                msg = "Invalid array for input time"
                err = OrbitError(msg)
                print(err)
                exit()

        def _calc_timedelta(self):
            try:
                epoch = np.asarray([
                    self.ephemeris.epoch_datetime.strftime("%Y-%m-%dT%H:%M:%S.%fZ")
                    ]).astype("datetime64")
                delta = (self.datetime - epoch)[:, np.newaxis]
                self._timedelta = (delta / np.timedelta64(1, "D")).astype(float)
            except TypeError:
                msg = "Ephemeris attribute from Orbit instance is not complete"
                err = OrbitError(msg)
                print(err)
                exit()

        def _calc_mean_motion(self):

            dt = self.timedelta
            n0 = self.ephemeris.mean_motion
            n1 = self.ephemeris.mean_motion_first_dif
            n2 = self.ephemeris.mean_motion_second_dif

            self._mean_motion = np.asarray(n0 + n1 * dt + 0.5 * n2 * dt**2)

        def _calc_mean_anomaly(self):

            dt = self.timedelta
            m0 = self.ephemeris.mean_anomaly.rad
            n = self.mean_motion * REV2RAD
            self._mean_anomaly = np.asarray(m0 + n * dt)

        def _calc_eccentric_anomaly(self):

            eccent = self.ephemeris.eccentricity
            m = self.mean_anomaly
            ean_old = m
            dif = np.asarray([1])

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

            dt = self.timedelta
            i = self.ephemeris.inclination.rad
            n = self.mean_motion

            eccent = self.ephemeris.eccentricity
            a = self.semimajor_axis
            p = a * (1-eccent**2)

            factor = -REV2RAD * j2 / (mu*p**2) * 3 * n
            w0 = self.ephemeris.argument_of_perigee.rad
            w1 = factor * ((5/4)*np.sin(i)**2 - 1)
            omega0 = self.ephemeris.longitude_of_the_ascending_node.rad
            omega1 = (factor/2) * np.cos(i)

            self._argument_of_perigee = np.asarray(w0 + w1 * dt)
            self._longitude_of_the_ascending_node = np.asarray(omega0 + omega1 * dt)

        def _calc_cartesian_coordinates(self):

            we = EARTH_ANGULAR_SPEED

            i = self.ephemeris.inclination.rad
            o = self.longitude_of_the_ascending_node
            w = self.argument_of_perigee
            v = self.true_anomaly

            eccent = self.ephemeris.eccentricity
            a = self.semimajor_axis
            p = a * (1 - eccent**2)

            r_s = p / (1 + eccent * np.cos(v))

            rx_s1 = r_s * (np.cos(o)*np.cos(w+v) - np.sin(o)*np.sin(w+v)*np.cos(i))
            ry_s1 = r_s * (np.sin(o)*np.cos(w+v) + np.cos(o)*np.sin(w+v)*np.cos(i))
            rz_s1 = r_s * (np.sin(w+v)*np.sin(i))

            ref1 = np.datetime64("2000-01-01T12:00:00Z")
            ref2 = np.datetime64("2000-01-01T00:00:00Z")
            dtm1 = ((self.datetime - ref1) / np.timedelta64(1, "D")).astype(float)[:, np.newaxis] / 36525
            dtm2 = ((self.datetime - ref2) / np.timedelta64(1, "D")).astype(float)[:, np.newaxis]
            dtm2 = ((dtm2 % dtm2.astype(int)) * 86400)


            gst0 = 24110.54841 + 8640184.812866 * dtm1 + 0.093104 * dtm1**2 - 6.2e-6 * dtm1**3
            gst2 = (gst0/86400 * 2*np.pi) + we * dtm2

            rx_s = + np.cos(gst2)*rx_s1 + np.sin(gst2)*ry_s1
            ry_s = - np.sin(gst2)*rx_s1 + np.cos(gst2)*ry_s1
            rz_s = rz_s1

            #rx_s = rx_s1
            #ry_s = ry_s1

            self._satellite_position_ecef = np.asarray([rx_s, ry_s, rz_s]).T

        def _calc_geodetic_coordinates(self):

            ae = EARTH_SEMIMAJOR_AXIS
            e2 = EARTH_FLATTENING_FACTOR

            rx_s, ry_s, rz_s = self._satellite_position_ecef.T
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
            dif = np.asarray([1])

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

        _add_datetime(self, datetime)
        _calc_timedelta(self)

        _calc_mean_motion(self)
        _calc_mean_anomaly(self)
        _calc_eccentric_anomaly(self)
        _calc_true_anomaly(self)
        _calc_semimajor_axis(self)
        _calc_argument_of_perigee_and_longitude_of_the_ascending_node(self)
        _calc_cartesian_coordinates(self)
        _calc_geodetic_coordinates(self)
