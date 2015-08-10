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
        "satellite_position_ecf",
        "satellite_position_eci",
        "satellite_position_geo",
        "satellite_velocity_ecf",
        "satellite_velocity_eci",
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
            msg = "Ephemeris instance is not complete"
            err = OrbitError(msg)
            print(err)
            exit()

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
    def satellite_position_ecf(self):
        """satellite position in ECEF reference system"""
        return self._satellite_position_ecf

    @property
    @returns(np.ndarray)
    def satellite_position_eci(self):
        """satellite position in ECI reference system"""
        return self._satellite_position_eci

    @property
    @returns(np.ndarray)
    def satellite_velocity_ecf(self):
        """satellite velocity in ECEF reference system"""
        return self._satellite_velocity_ecf

    @property
    @returns(np.ndarray)
    def satellite_velocity_eci(self):
        """satellite velocity in ECI reference system"""
        return self._satellite_velocity_eci

    @property
    @returns(np.ndarray)
    def satellite_position_geo(self):
        """satellite position in geodetic reference system"""
        return self._satellite_position_geo

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
            """Safe setter for datetime attribute."""
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
            """Compute time difference with respect to epoch."""
            try:
                fmt = "%Y-%m-%dT%H:%M:%S.%fZ"
                epoch = np.asarray([
                    self.ephemeris.epoch_datetime.strftime(fmt)
                    ]).astype("datetime64")
                delta = (self.datetime - epoch)[:, np.newaxis]
                self._timedelta = (delta / np.timedelta64(1, "D")).astype(float)
            except TypeError:
                msg = "Ephemeris attribute from Orbit instance is not complete"
                err = OrbitError(msg)
                print(err)
                exit()

        def _calc_mean_motion(self):
            """Compute mean motion for a specific datetime."""

            # Call necessary properties.
            dt = self.timedelta
            n0 = self.ephemeris.mean_motion
            n1 = self.ephemeris.mean_motion_first_dif
            n2 = self.ephemeris.mean_motion_second_dif
            # Set mean motion property.
            self._mean_motion = n0 + n1 * dt + 0.5 * n2 * dt**2

        def _calc_mean_anomaly(self):
            """Compute mean anomaly for a specific datetime."""

            # Call necessary properties.
            dt = self.timedelta
            m0 = self.ephemeris.mean_anomaly.rad
            n = self.mean_motion * REV2RAD
            # Set mean anomaly property.
            self._mean_anomaly = m0 + n * dt

        def _calc_eccentric_anomaly(self):
            """Compute eccentric anomaly for a specific datetime."""

            # Call necessary properties.
            eccent = self.ephemeris.eccentricity
            m = self.mean_anomaly
            ean_old = m
            dif = np.asarray([1])
            # Solve Kepler's equation by recursive methods until a tolerance
            # is reached.
            while (dif > 0).all():
                ean_new = np.where(dif > 0, m + eccent*np.sin(ean_old), ean_old)
                dif = np.abs(ean_new - ean_old) - ECCENTRIC_ANOMALY_TOLERANCE
                ean_old = ean_new
            # Set eccentric anomaly property.
            self._eccentric_anomaly = ean_old

        def _calc_true_anomaly(self):
            """Compute true anomaly for a specific datetime."""

            # Call necessary properties.
            eccent = self.ephemeris.eccentricity
            factor = np.sqrt((1+eccent)/(1-eccent))
            ean = self.eccentric_anomaly
            # Set true anomaly property.
            self._true_anomaly = 2 * np.arctan(factor * np.tan(ean/2))

        def _calc_semimajor_axis(self):
            """Compute semimajor axis for a specific datetime."""

            # Call necessary constants and properties.
            mu = STANDARD_GRAVITATIONAL_PARAMETER
            n = self.mean_motion * REV2RAD / DAY2SEC
            # Set semimajor axis property using Kepler's third law.
            self._semimajor_axis = (mu / n**2)**(1/3)

        def _calc_argument_of_perigee_and_longitude_of_the_ascending_node(self):
            """Compute argument of perigee and longitude of the ascending
            node for a specific datetime."""

            # Call necessary constants.
            j2 = SECOND_ZONAL_TERM_J2
            mu = STANDARD_GRAVITATIONAL_PARAMETER
            # Call necessary properties.
            dt = self.timedelta
            i = self.ephemeris.inclination.rad
            n = self.mean_motion
            eccent = self.ephemeris.eccentricity
            a = self.semimajor_axis
            p = a * (1-eccent**2)
            # Compute first derivatives of argument of perigee and longitude
            # of the ascending node.
            factor = -REV2RAD * j2 / (mu*p**2) * 3 * n
            w0 = self.ephemeris.argument_of_perigee.rad
            w1 = factor * ((5/4)*np.sin(i)**2 - 1)
            o0 = self.ephemeris.longitude_of_the_ascending_node.rad
            o1 = (factor/2) * np.cos(i)
            # Set argument of perigee property.
            self._argument_of_perigee = w0 + w1 * dt
            # Set longitude of the ascending node property.
            self._longitude_of_the_ascending_node = o0 + o1 * dt

        def _calc_eci_coordinates(self):
            """Compute ECI coordinates for a specific datetime."""

            # Call necessary constants.
            mu = STANDARD_GRAVITATIONAL_PARAMETER
            # Call necessary properties.
            i = self.ephemeris.inclination.rad
            o = self.longitude_of_the_ascending_node
            w = self.argument_of_perigee
            v = self.true_anomaly
            eccent = self.ephemeris.eccentricity
            a = self.semimajor_axis
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
            self._satellite_position_eci = np.asarray([rx, ry, rz]).T
            self._satellite_velocity_eci = np.asarray([vx, vy, vz]).T

        def _calc_ecf_coordinates(self):
            """Compute ECEF coordinates for a specific datetime."""

            # Call necessary constants.
            we = EARTH_ANGULAR_SPEED
            t1 = np.datetime64("2000-01-01T12:00:00Z")
            t2 = np.datetime64("2000-01-01T00:00:00Z")
            # Call necessary properties.
            rx, ry, rz = self.satellite_position_eci.T
            vx, vy, vz = self.satellite_velocity_eci.T
            # Compute Julian centuries since epoch 2000/01/01 12:00:00.
            dtm1 = ((self.datetime - t1) / np.timedelta64(1, "D")).astype(float)
            dtm1 = dtm1[:, np.newaxis] / 36525
            # Compute UTC time for the datetime of interest.
            dtm2 = ((self.datetime - t2) / np.timedelta64(1, "D")).astype(float)
            dtm2 = dtm2[:, np.newaxis]
            dtm2 = ((dtm2 % dtm2.astype(int)) * 86400)
            # Compute the Greenwich Sidereal Time (GST), that is, the angle
            # between the vernal point and the Greenwich meridian (which is
            # also the angle between the ECI and ECEF reference systems).
            c0, c1, c2, c3 = [24110.54841, 8640184.812866, 0.093104, -6.2e-6]
            gst0 = c0 + c1 * dtm1 + c2 * dtm1**2 + c3 * dtm1**3
            gst2 = (gst0/86400 * 2*np.pi) + we * dtm2
            # Compute ECEF coordinates as a rotation of ECI coordinates.
            rx_s = + np.cos(gst2)*rx + np.sin(gst2)*ry
            ry_s = - np.sin(gst2)*rx + np.cos(gst2)*ry
            rz_s = rz
            # Compute (Vx,Vy,Vz) velocities in ECF reference system.
            vx_s = + np.cos(gst2)*vx + np.sin(gst2)*vy
            vy_s = - np.sin(gst2)*vx + np.cos(gst2)*vy
            vz_s = vz
            # Set ECEF satellite position and velocity properties.
            self._satellite_position_ecf = np.hstack([rx_s, ry_s, rz_s])
            self._satellite_velocity_ecf = np.hstack([vx_s, vy_s, vz_s])

        def _calc_geo_coordinates(self):
            """Compute geodetic coordinates for a specific datetime."""

            # Call necessary constants.
            ae = EARTH_SEMIMAJOR_AXIS
            e2 = EARTH_FLATTENING_FACTOR
            # Call necessary properties.
            rx_s, ry_s, rz_s = self.satellite_position_ecf.T
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
            self._satellite_position_geo = np.asarray([lat, lon, alt]).T

        _add_datetime(self, datetime)
        _calc_timedelta(self)
        _calc_mean_motion(self)
        _calc_mean_anomaly(self)
        _calc_eccentric_anomaly(self)
        _calc_true_anomaly(self)
        _calc_semimajor_axis(self)
        _calc_argument_of_perigee_and_longitude_of_the_ascending_node(self)
        _calc_eci_coordinates(self)
        _calc_ecf_coordinates(self)
        _calc_geo_coordinates(self)
