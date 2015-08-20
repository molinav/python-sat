from __future__ import division
from . constants.Scene import DATASET_COORDINATE_TOLERANCE
from . constants.Scene import DATASET_DRIVER
from . constants.Scene import DATE_PATTERN
from . constants.Scene import DEFAULT_SCAN_STEP_ANGLE
from . constants.Orbit import EARTH_FLATTENING_FACTOR
from . constants.Orbit import EARTH_SEMIMAJOR_AXIS
from . decorators import accepts
from . decorators import returns
from . Angle import Angle
from . Ephemeris import Ephemeris
from . Error import SceneError
from . Orbit import Orbit
from datetime import datetime
from datetime import timedelta
from gdal import Dataset
from gdal import Open
from six import text_type
import numpy as np
import re


class Scene(object):

    __slots__ = [
        "_dataset",
        "_end_datetime",
        "_orbit",
        "_position_ecf",
        "_position_geo",
        "_position_pix",
        "_scan_step_angle",
        "_scan_timedelta",
        "_start_datetime",
        ]

    @accepts(text_type, Ephemeris)
    def __init__(self, path, ephemeris):
        """Constructor of a generic Scene instance.

        Parameters:

        path
            string path where the image is stored
        ephemeris
            Ephemeris instance containing the Kepler parameters
        """

        # Check that the dataset is valid.
        dataset = Open(path)
        if not dataset:
            msg = "Invalid path for image file"
            raise SceneError(msg)
        if dataset.GetDriver().LongName != DATASET_DRIVER:
            msg = "Image file must be a {}".format(DATASET_DRIVER)
            raise SceneError(msg)
        for item in self.__slots__:
            self.__setattr__(item, None)
        # Build Orbit instance from list of TLE lines.
        orbit = Orbit(ephemeris.to_copy())
        # Set main attributes.
        self._dataset = dataset
        self._orbit = orbit

    @accepts(int, int)
    def compute(self, x_density, y_density):
        """Compute scene coordinates for corresponding datetime."""

        # Proceed only if the Scene instance's attributes are well-defined.
        try:
            # Call necessary properties.
            xsize = self.dataset.RasterXSize
            ysize = self.dataset.RasterYSize
            sdate = self.start_datetime
            dtime = self.scan_timedelta
            # Check that point density is valid both along and across track.
            if x_density > (xsize // 10):
                msg = "{} (minimum {}, found {})".format(
                    "Insufficient across-track point density",
                    xsize // 4,
                    x_density)
                raise SceneError(msg)
            if y_density > (ysize // 10):
                msg = "{} (minimum {}, found {})".format(
                    "Insufficient along-track point density",
                    ysize // 4,
                    y_density)
                raise SceneError(msg)
            # Set point spacing.
            xran = np.arange((xsize // 2) % x_density, xsize, x_density)
            yran = np.arange((ysize // 2) % y_density, ysize, y_density)
            xnum = len(xran)
            ynum = len(yran)
            # Compute expected datetimes for every scan line and propagate
            # attributes from Orbit property.
            times = np.asarray([sdate + n*dtime for n in range(ysize)])[:, None]
            self.orbit.compute(times)
        except AttributeError:
            parent = self.__class__.__name__
            msg = "{} attributes are not complete".format(parent)
            err = SceneError(msg)
            err.__cause__ = None
            raise err
        # Call necessary constants.
        ae = EARTH_SEMIMAJOR_AXIS
        be = ae * np.sqrt(1 - EARTH_FLATTENING_FACTOR)
        # Call necessary properties.
        r_vec = self.orbit.position_ecf
        u_orb = self._compute_spacecraft_fixed_unit_vectors()
        w_xyz = self._compute_spacecraft_geodetic_vectors()
        # Set attitude angles to 0 radians.
        roll, pitch, yaw = np.radians([0, 0, 0])
        rot_a = np.array([[1, pitch, -roll], [-pitch, 1, yaw], [roll, -yaw, 1]])
        # Create temporary arrays where ECF and geodetic coordinates will be
        # stored during the computation.
        tmp_ecf = np.empty((xnum*ynum, 3))
        tmp_geo = np.empty((xnum*ynum, 3))
        tmp_pix = np.empty((xnum*ynum, 2))
        for i, j in enumerate(yran):
            rot_t = w_xyz[:, :, j]
            u_xyz = np.dot(np.dot(rot_t.T, rot_a.T), u_orb[:, xran])
            # Estimate recursively the ECF coordinates for every image's pixel.
            r_vec_ii = r_vec[j][:, None]
            pix_norm = ae
            sat_norm = np.linalg.norm(r_vec_ii, axis=0)[None, :]
            dotru = np.sum(r_vec_ii * u_xyz, axis=0)[None, :]
            dif = np.array([1])
            # Repeat recursive method until a tolerance is reached.
            while (dif > DATASET_COORDINATE_TOLERANCE).all():
                # Calculate new ECF vectors for every pixel.
                d_norm = - dotru - np.sqrt(pix_norm**2 - sat_norm**2 + dotru**2)
                pix_vec_ii = d_norm * u_xyz + r_vec_ii
                # Update image position coordinates.
                self._position_ecf = pix_vec_ii.T
                self._calc_coordinates_from_ecf_to_geo()
                # Compute new norm for pixel vectors.
                lat = self.position_geo[:, 0]
                pix_norm2 = np.sqrt(
                    ((ae**2 * np.cos(lat))**2 + (be**2 * np.sin(lat))**2) /
                    ((ae * np.cos(lat))**2 + (be * np.sin(lat))**2))
                # Overwrite old temporary variables with new values.
                dif = (pix_norm2 - pix_norm)**2
                pix_norm = pix_norm2
            # Grow temporary variables which store ECF and geodetic coordinates.
            tmp_ecf[i*xnum:(i+1)*xnum, :] = self.position_ecf
            tmp_geo[i*xnum:(i+1)*xnum, :] = self.position_geo
            tmp_pix[i*xnum:(i+1)*xnum, :] = np.asarray([
                j * np.ones((xsize,))[xran] - 0.5,
                xsize - xran - 0.5,
                #(sign==-1)*xsize + sign*xran - 0.5,
                ]).T
        # Update ECF and geodetic coordinates.
        self._position_ecf = tmp_ecf
        self._position_geo = tmp_geo
        self._position_pix = tmp_pix

    def _calc_coordinates_from_ecf_to_geo(self):
        """Compute geodetic coordinates for a specific datetime."""
        Orbit._calc_coordinates_from_ecf_to_geo(self)

    def _compute_spacecraft_fixed_unit_vectors(self):
        """Return unit vectors in spacecraft-fixed coordinate system."""

        delta_step = self.scan_step_angle.rad
        xsize = self.dataset.RasterXSize
        delta = (xsize//2 - np.arange(xsize)[:, None]) * delta_step

        ub = np.hstack([np.cos(delta), np.zeros((xsize, 1)), np.sin(delta)]).T
        return ub

    def _compute_spacecraft_geodetic_vectors(self):
        """Return unit vectors in spacecraft-geodetic coordinate system."""

        # Call necessary properties.
        lat, lon, alt = [x[:, None] for x in self.orbit.position_geo.T]
        vs = self.orbit.velocity_ecf

        # Compute unit vectors (wx, wy, wz).
        wx = np.hstack([
            -np.cos(lat)*np.cos(lon), -np.cos(lat)*np.sin(lon), -np.sin(lat)])
        wz = np.cross(vs, wx)
        wz /= np.linalg.norm(wz, axis=1)[:, None]
        wy = np.cross(wz, wx)

        # Store all the unit vectors as 3x3 arrays with [wx, wy, wz].T
        # structure concatenated along the third axis.
        w_xyz = np.vstack([wx.T[None, :], wy.T[None, :], wz.T[None, :]])
        return w_xyz

    @classmethod
    def _parse_datetime(cls, text_date):
        """Return a datetime object from a valid datetime string."""

        match = re.match(DATE_PATTERN, text_date)
        if match:
            yr = int(match.group(1))
            dy = int(match.group(2)) - 1
            ms = int("".join([match.group(3), "000"]))
            return datetime(yr, 1, 1) + timedelta(days=dy, microseconds=ms)
        else:
            msg = "Input text does not match datetime pattern"
            raise ValueError(msg)

    @property
    @returns(Dataset)
    def dataset(self):
        """dataset containing the image"""
        return self._dataset

    @property
    @returns(datetime)
    def end_datetime(self):
        """datetime for last scanline"""

        if self._end_datetime is None:
            try:
                text_date = self.dataset.GetMetadata()["STOP"]
                self._end_datetime = self._parse_datetime(text_date)
            except KeyError:
                msg = "Scene dataset does not contain end datetime"
                err = SceneError(msg)
                err.__cause__ = None
                raise err
        return self._end_datetime

    @property
    def orbit(self):
        """Orbit instance for the associated Earth satellite"""
        return self._orbit

    @property
    @returns(np.ndarray)
    def position_ecf(self):
        """dataset position in ECF reference system"""
        return self._position_ecf

    @property
    @returns(np.ndarray)
    def position_geo(self):
        """dataset position in geodetic reference system"""
        return self._position_geo

    @property
    @returns(np.ndarray)
    def position_pix(self):
        """dataset position in pixel (row, col) reference system"""
        return self._position_pix

    @property
    @returns(timedelta)
    def scan_timedelta(self):
        """period of time which is necessary to measure a scanline"""

        if self._scan_timedelta is None:
            sec = (self.end_datetime - self.start_datetime).total_seconds()             
            sec /= (self.dataset.RasterYSize - 1)
            self._scan_timedelta = timedelta(seconds=sec)
        return self._scan_timedelta

    @property
    def scan_step_angle(self):
        """angle between across-track pixels within the image"""

        if self._scan_step_angle is None:
            self._scan_step_angle = Angle(deg=DEFAULT_SCAN_STEP_ANGLE)
        return self._scan_step_angle

    @property
    @returns(datetime)
    def start_datetime(self):
        """datetime for first scanline"""

        if self._start_datetime is None:
            try:
                text_date = self.dataset.GetMetadata()["START"]
                self._start_datetime = self._parse_datetime(text_date)
            except KeyError:
                msg = "Scene dataset does not contain start datetime"
                err = SceneError(msg)
                err.__cause__ = None
                raise err
        return self._start_datetime
