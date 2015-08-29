from __future__ import division
from __future__ import print_function
from . constants.Angle import DEG2RAD
from . constants.Coordinates import EARTH_FLATTENING_FACTOR
from . constants.Coordinates import EARTH_SEMIMAJOR_AXIS
from . constants.Scene import DATASET_COORDINATE_TOLERANCE
from . constants.Scene import DATASET_DRIVER
from . constants.Scene import DATE_PATTERN
from . constants.Scene import DEFAULT_SCAN_STEP_LOOK_ANGLE
from . constants.Scene import HPRT_BORDER_MASK_MARGIN
from . constants.Scene import HRPT_LATITUDE_STEP
from . constants.Scene import HRPT_LONGITUDE_STEP
from . constants.Scene import HRPT_CLOUD_THRESHOLD_B5
from . constants.Scene import SPACECRAFT_ASCENDING
from . constants.Scene import SPACECRAFT_DESCENDING
from . decorators import accepts
from . decorators import returns
from . Angle import Angle
from . Coordinates import Coordinates
from . Ephemeris import Ephemeris
from . Error import SceneError
from . Orbit import Orbit
from datetime import datetime
from datetime import timedelta
from scipy.interpolate import griddata
from scipy.interpolate import NearestNDInterpolator
from scipy.interpolate import RectBivariateSpline
from scipy.ndimage import binary_dilation
from scipy.ndimage import binary_erosion
from scipy.ndimage import distance_transform_edt
from scipy.sparse import coo_matrix
from six import text_type
import gdal
import numpy as np
import re


class Scene(object):

    __slots__ = [
        "_dataset",
        "_clock_drift",
        "_end_datetime",
        "_gcp_position_ecf",
        "_gcp_position_geo",
        "_gcp_position_pix",
        "_gcp_xrange",
        "_gcp_yrange",
        "_hrpt_geo",
        "_hrpt_img",
        "_hrpt_msk",
        "_hrpt_gtr",
        "_hrpt_pix",
        "_orbit",
        "_scan_step_look_angle",
        "_scan_timedelta",
        "_scan_xsize",
        "_scan_ysize",
        "_spacecraft_attitude",
        "_spacecraft_direction",
        "_spacecraft_fixed_look_angle",
        "_spacecraft_fixed_unit_look_vector",
        "_spacecraft_geodetic_unit_vector",
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
        gdal.PushErrorHandler("CPLQuietErrorHandler")
        dataset = gdal.Open(path)
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
        """Compute scene coordinates for corresponding datetime.

        Parameters:

        x_density
            x spacing between ground control points
        y_density
            y spacing between ground control points
        """

        # Proceed only if the Scene instance's attributes are well-defined.
        try:
            # Check that point density is valid both along and across track.
            if x_density > (self.scan_xsize // 10):
                msg = "{} (minimum {}, found {})".format(
                    "Insufficient across-track point density",
                    self.scan_xsize // 4,
                    x_density)
                raise SceneError(msg)
            if y_density > (self.scan_ysize // 10):
                msg = "{} (minimum {}, found {})".format(
                    "Insufficient along-track point density",
                    self.scan_ysize // 4,
                    y_density)
                raise SceneError(msg)
            # Set control point indices.
            self._gcp_xrange = np.arange(
                0.5 + (self.scan_xsize // 2) % x_density,
                self.scan_xsize,
                x_density)
            self._gcp_yrange = np.arange(
                0.5 + (self.scan_ysize // 2) % y_density,
                self.scan_ysize,
                y_density)
            # Compute expected datetimes for every scan line and propagate
            # attributes from Orbit property.
            times = np.asarray([
                self.start_datetime + self.clock_drift + n * self.scan_timedelta
                for n in range(self.scan_ysize)])[:, None]
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
        u_orb = self._compute_spacecraft_fixed_unit_look_vectors()
        w_xyz = self._compute_spacecraft_geodetic_unit_vectors()
        # Read attitude angles and build rotation array.
        roll, pitch, yaw = self.spacecraft_attitude
        rot_a = np.array([[1, pitch, -roll], [-pitch, 1, yaw], [roll, -yaw, 1]])
        # Create temporary arrays where ECF and geodetic coordinates will be
        # stored during the computation.
        xnum = len(self.gcp_xrange)
        ynum = len(self.gcp_yrange)
        xran = self.gcp_xrange.astype(int)
        tmp_ecf = np.empty((xnum*ynum, 3))
        tmp_geo = np.empty((xnum*ynum, 3))
        tmp_pix = np.empty((xnum*ynum, 2))
        for i, j in enumerate(self.gcp_yrange):
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
                self._gcp_position_ecf = pix_vec_ii.T
                self._gcp_position_geo = Coordinates.from_ecf_to_geo(
                    self.gcp_position_ecf)
                # Compute new norm for pixel vectors.
                lat = self.gcp_position_geo[:, 0]
                pix_norm2 = np.sqrt(
                    ((ae**2 * np.cos(lat))**2 + (be**2 * np.sin(lat))**2) /
                    ((ae * np.cos(lat))**2 + (be * np.sin(lat))**2))
                # Overwrite old temporary variables with new values.
                dif = (pix_norm2 - pix_norm)**2
                pix_norm = pix_norm2
            # Grow temporary variables which store ECF and geodetic coordinates.
            tmp_ecf[i*xnum:(i+1)*xnum, :] = self.gcp_position_ecf
            tmp_geo[i*xnum:(i+1)*xnum, :] = self.gcp_position_geo
            tmp_pix[i*xnum:(i+1)*xnum, :] = np.asarray([
                j * np.ones((xnum,)), self.gcp_xrange]).T
        # Update ECF and geodetic coordinates.
        self._gcp_position_ecf = tmp_ecf
        self._gcp_position_geo = tmp_geo
        self._gcp_position_pix = tmp_pix

    def build_geographic_window(self, slices=400, buffer=100, fill_value=-999):
        """Calc geodetic coordinates for all the hrpt image and build
        geodetic grid for pixel indices and its geotransform.

        Optional parameters:

        slices
            number of rows in which the method is split (default 400)
        buffer
            number of overlapping rows between two slices (default 100)
        fill_value
            no data value (default -999)
        """

        # Verify the types and ranges of all the input arguments.
        if not isinstance(slices, int):
            msg = "slices must be a positive int"
            raise SceneError(msg)
        if slices < 0 or slices >= self.scan_ysize:
            msg = "slices must be positive and not greater than image length"
            raise SceneError(msg)
        if not isinstance(buffer, int):
            msg = "buffer must be a positive int"
            raise SceneError(msg)
        if buffer < 0 or buffer >= self.scan_ysize:
            msg = "buffer must be positive and not greater than image length"
            raise SceneError(msg)
        if not isinstance(fill_value, int):
            msg = "fill value must be an int"
            raise SceneError(msg)

        # Show start message.
        msg = "\nComputing HRPT attributes for {} instance with id {}"
        print(msg.format(self.__class__.__name__, id(self)))

        # Call necessary properties.
        xstep = 0.5 + np.arange(self.scan_xsize)
        ystep = 0.5 + np.arange(self.scan_ysize)

        # Show info message.
        msg = "\n  Calculating geodetic coordinates for all the HRPT pixels"
        print(msg)

        def calc_hrpt_geo():
            """Inner function where self._hrpt_geo is computed."""

            # Call necessary properties.
            yran = self.gcp_yrange
            xran = self.gcp_xrange
            gcps = self.gcp_position_geo[:, :2] / DEG2RAD
            gcps = np.rollaxis(gcps.reshape((len(yran), len(xran), 2)), 2)
            # Compute latitude and longitude interpolators.
            lat_spl = RectBivariateSpline(x=yran, y=xran, z=gcps[0])
            lon_spl = RectBivariateSpline(x=yran, y=xran, z=gcps[1])
            # Compute latitude and longitude for all the pixels.
            latlon = np.asarray([f(ystep, xstep) for f in [lat_spl, lon_spl]])
            height = np.zeros((1,) + latlon.shape[1:])
            self._hrpt_geo = np.vstack([latlon * DEG2RAD, height])

        # Compute latitudes and longitudes for all the HRPT pixels.
        calc_hrpt_geo()
        raw_ymin, raw_xmin = self.hrpt_geo.min(axis=1).min(axis=1)[:2] / DEG2RAD
        raw_ymax, raw_xmax = self.hrpt_geo.max(axis=1).max(axis=1)[:2] / DEG2RAD

        # Show info message.
        msg = "\n  Calculating geotransform for gridded HRPT indices image"
        print(msg)

        def calc_hrpt_gtr():
            """Inner function where self._hrpt_gtr is computed."""

            # Compute geotransform parameters.
            prj_st = np.asarray([HRPT_LATITUDE_STEP, HRPT_LONGITUDE_STEP])
            prj_nw = np.asarray([raw_ymax, raw_xmin]) - 0.5 * prj_st
            prj_se = np.asarray([raw_ymin, raw_xmax]) + 0.5 * prj_st
            prj_sz = np.round((prj_se - prj_nw) / prj_st).astype(int)

            # Set geotransform property.
            self._hrpt_gtr = (prj_nw[1], prj_st[1], 0, prj_nw[0], 0, prj_st[0])

            # Compute pixel indices in the new grid.
            prj_st3 = prj_st[:, None, None]
            prj_nw3 = prj_nw[:, None, None]
            prj_rc = (self.hrpt_geo[:2]/DEG2RAD - prj_nw3) // prj_st3

            return prj_sz, prj_rc

        # Compute geotransform for the new grid.
        prj_size, prj_rowcol = calc_hrpt_gtr()

        # Show info message.
        msg = "\n  Calculating gridded HRPT border mask"
        print(msg)

        def calc_border_mask():
            """Inner function where gridded border mask is obtained."""

            margin = HPRT_BORDER_MASK_MARGIN
            iter_n = int(margin//2)
            msk_sz = prj_size + 2 * margin
            msk_rowcol = prj_rowcol + margin
            row, col = msk_rowcol.reshape((2, np.prod(msk_rowcol.shape[1:])))
            # Arrange test array.
            ones = np.ones(row.shape)
            test = coo_matrix((ones, (row, col)), shape=msk_sz).toarray()
            test = test > 0
            kern = np.ones((3, 3), dtype=bool)
            # Dilate and erode test array to fill holes.
            test = binary_dilation(test, structure=kern, iterations=iter_n)
            test = binary_erosion(test, structure=kern, iterations=iter_n+5)
            return np.logical_not(test[None, margin:-margin, margin:-margin])

        # Set pixels within border mask to fill value.
        border_mask = calc_border_mask()

        # Show info message.
        msg = "\n  Calculating gridded HRPT indices"
        print(msg)

        def calc_hrpt_pts():
            """Inner function that returns xx, yy and vv values."""

            # Compute intermediate arrays.
            sz = (-1, 2)
            tp = np.float32
            raw_rc = np.mgrid[:self.scan_ysize, :self.scan_xsize] + 0.5
            qry_rc = np.mgrid[:prj_size[0], :prj_size[1]] + 0.5

            # Compute outputs.
            rc = [prj_rowcol, raw_rc, qry_rc]
            return [np.rollaxis(x, 0, 3).reshape(sz).astype(tp) for x in rc]

        # Set pixels from gridded image as xx and pixels from hrpt as yy. Set
        # also pixels to be evaluated as vv.
        xx, yy, vv = calc_hrpt_pts()

        def calc_hrpt_pix_legacy1():
            """Legacy inner function where self._hrpt_pix is computed."""

            func = NearestNDInterpolator(xx, yy)
            temp = func(vv).reshape((prj_size[0], prj_size[1], 2))
            self._hrpt_pix = np.rollaxis(temp, 2)
            self._hrpt_pix[border_mask] = fill_value

        def calc_hrpt_pix_legacy2():
            """Legacy inner function where self._hrpt_pix is computed."""

            # Prepare info message for loop.
            nfmt = len(str(prj_size[0]))
            msg2 = "\n    rows {1:{0}d} - {2:{0}d} out of {3}"

            # For every slice project the pixel indices.
            temp = np.empty((2, prj_size[0], prj_size[1]), dtype=np.int16)

            for lim1 in range(0, prj_size[0], slices):

                lim2 = min(prj_size[0], lim1 + slices)

                # Show info message.
                print(msg2.format(nfmt, lim1, lim2, prj_size[0]))

                # Get array slices for the current loop.
                flag = (xx[:, 0] > lim1-buffer) & (xx[:, 0] < lim2+buffer)
                xx_n = xx[flag]
                yy_n = yy[flag]
                vv_n = vv[lim1*prj_size[1]:lim2*prj_size[1]]

                # Convert to gridded data.
                size = (lim2-lim1, prj_size[1], 2)
                temp[:, lim1:lim2] = np.rollaxis(
                    griddata(xx_n, yy_n, vv_n, method="nearest",
                             fill_value=fill_value,).reshape(size), 2)

            # Set array into self._hrpt_pix and mask border values.
            self._hrpt_pix = temp
            self._hrpt_pix[border_mask] = fill_value

        def calc_hrpt_pix():
            """Inner function where self._hrpt_pix is computed."""

            # Arrange test array.
            row, col = xx.T
            ones_array = np.ones((yy.shape[0],), dtype=int)
            test_array = np.asarray([
                coo_matrix((val, (row, col)), shape=prj_size).toarray()
                for val in [yy[:, 0], yy[:, 1], ones_array]])

            # Divide array of indices by the array of weights.
            nums = test_array[None, 2]
            mass = np.where(nums, nums, 1)
            data = (test_array[:2] // mass).astype(np.int16)

            # Build necessary masks.
            nodata_mask = np.equal(nums, 0)
            values_mask = np.logical_not(border_mask)
            distance_mask = np.repeat(nodata_mask, 2, axis=0) & values_mask

            # Prepare info message for loop.
            nfmt = len(str(prj_size[0]))
            msg2 = "\n    rows {1:{0}d} - {2:{0}d} out of {3}"

            n = data.shape[1]
            temp = np.empty(data.shape, dtype=np.int16)

            for lim1 in range(0, n, slices):

                lim2 = min(n, lim1 + slices)
                lim1_b = max(0, lim1 - buffer)
                lim2_b = min(n, lim2 + buffer)

                # Show info message.
                print(msg2.format(nfmt, lim1, lim2, data.shape[1]))

                # Store array of indices in a temporary array.
                indices = distance_transform_edt(
                    distance_mask[:, lim1_b:lim2_b, :],
                    return_distances=False,
                    return_indices=True)
                temp_b = data[:, lim1_b:lim2_b, :][tuple(indices)]
                temp[:, lim1:lim2, :] = temp_b[:, lim1-lim1_b:lim2-lim1_b, :]

            # Store temporary array of indices in self._hrpt_pix.
            self._hrpt_pix = temp
            self._hrpt_pix[:, border_mask[0]] = fill_value

        # Calculate HRPT indices image.
        calc_hrpt_pix()

        # Show info message.
        msg = "\n  Calculating gridded HRPT bands"
        print(msg)

        def calc_hrpt_img():
            """Inner function where self._hrpt_img is computed"""

            # Extract HRPT bands from GDAL Dataset instance.
            bnd_raw = np.vstack([
                self.dataset.GetRasterBand(i+1).ReadAsArray(
                    0, 0, self.scan_xsize, self.scan_ysize)[None, :]
                for i in range(self.dataset.RasterCount)])

            # Generate the array of absolute indices that connect the final
            # geodetic dataset with the original raw dataset.
            index = self.hrpt_pix.astype(int)
            index = index[0] * bnd_raw.shape[2] + index[1]
            bnd_raw = bnd_raw.reshape((bnd_raw.shape[0], -1))

            # Assign digital-numbers using the indices array.
            bnd_geo_shape = (bnd_raw.shape[0],) + index.shape
            bnd_geo = bnd_raw.take(index.ravel(), axis=1, mode="warp")
            bnd_geo = bnd_geo.reshape(bnd_geo_shape).astype(np.int16)
            bnd_geo[:, border_mask[0]] = 0

            self._hrpt_img = bnd_geo

        # Calculate HRPT 5-band digital-number image.
        calc_hrpt_img()

        # Show info message.
        msg = "\n  Calculating gridded HRPT cloud masks"
        print(msg)

        def calc_hrpt_msk():
            """Inner function where self._hrpt_msk is computed"""

            # Set thermal evaluator (channel 5).
            cld = self.hrpt_img[4] > HRPT_CLOUD_THRESHOLD_B5

            self._hrpt_msk = 255 * cld[None, :].astype(np.uint8)

        # Calculate HRPT cloud mask.
        calc_hrpt_msk()

        # Show info message.
        msg = "\n  HRPT channels were successfully gridded with shape {}"
        print(msg.format(self.hrpt_img.shape))

    def _compute_spacecraft_fixed_unit_look_vectors(self):
        """Return unit look vectors in spacecraft-fixed coordinates."""

        self._spacecraft_fixed_unit_look_vector = np.vstack([
            np.cos(self.spacecraft_fixed_look_angle),
            np.zeros((1, self.scan_xsize)),
            np.sin(self.spacecraft_fixed_look_angle)])
        return self._spacecraft_fixed_unit_look_vector

    def _compute_spacecraft_geodetic_unit_vectors(self):
        """Return unit vectors in spacecraft-geodetic coordinates."""

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
        self._spacecraft_geodetic_unit_vector = w_xyz
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
    @returns(timedelta)
    def clock_drift(self):
        """on board clock variation"""
        if self._clock_drift is None:
            self._clock_drift = timedelta(seconds=0)
        return self._clock_drift

    @property
    @returns(gdal.Dataset)
    def dataset(self):
        """GDAL dataset containing the image"""
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
    @returns(np.ndarray)
    def gcp_xrange(self):
        """column indices where control points are defined"""
        return self._gcp_xrange

    @property
    @returns(np.ndarray)
    def gcp_yrange(self):
        """row indices where control points are defined"""
        return self._gcp_yrange

    @property
    @returns(np.ndarray)
    def hrpt_geo(self):
        """interpolated latitude and longitude values for all the pixels"""
        return self._hrpt_geo

    @property
    @returns(np.ndarray)
    def hrpt_img(self):
        """HRPT image in gridded window"""
        return self._hrpt_img

    @property
    @returns(np.ndarray)
    def hrpt_msk(self):
        """HRPT cloud mask in gridded window"""
        return self._hrpt_msk

    @property
    @returns(tuple)
    def hrpt_gtr(self):
        """geotransform for HRPT gridded window"""
        return self._hrpt_gtr

    @property
    @returns(np.ndarray)
    def hrpt_pix(self):
        """gridded image of HRPT indices"""
        return self._hrpt_pix

    @property
    @returns(Orbit)
    def orbit(self):
        """Orbit instance for the associated Earth satellite"""
        return self._orbit

    @property
    @returns(np.ndarray)
    def gcp_position_ecf(self):
        """dataset gcp position in ECF reference system"""
        return self._gcp_position_ecf

    @property
    @returns(np.ndarray)
    def gcp_position_geo(self):
        """dataset gcp position in geodetic reference system"""
        return self._gcp_position_geo

    @property
    @returns(np.ndarray)
    def gcp_position_pix(self):
        """dataset gcp position in pixel (row, col) reference system"""
        return self._gcp_position_pix

    @property
    @returns(timedelta)
    def scan_timedelta(self):
        """period of time which is necessary to measure a scanline"""

        if self._scan_timedelta is None:
            sec = (self.end_datetime - self.start_datetime).total_seconds()             
            sec /= (self.scan_ysize - 1)
            self._scan_timedelta = timedelta(seconds=sec)
        return self._scan_timedelta

    @property
    @returns(Angle)
    def scan_step_look_angle(self):
        """angle between across-track pixels within the image"""

        if self._scan_step_look_angle is None:
            self._scan_step_look_angle = Angle(deg=DEFAULT_SCAN_STEP_LOOK_ANGLE)
        return self._scan_step_look_angle

    @property
    @returns(int)
    def scan_xsize(self):
        """number of columns in the Scene instance"""
        if self._scan_xsize is None:
            self._scan_xsize = self.dataset.RasterXSize
        return self._scan_xsize

    @property
    @returns(int)
    def scan_ysize(self):
        """number of rows in the Scene instance"""
        if self._scan_ysize is None:
            self._scan_ysize = self.dataset.RasterYSize
        return self._scan_ysize

    @property
    @returns(np.ndarray)
    def spacecraft_attitude(self):
        """spacecraft (roll, pitch, yaw) angles"""
        if self._spacecraft_attitude is None:
            self._spacecraft_attitude = np.zeros((3, 1))
        return self._spacecraft_attitude

    @property
    @returns(int)
    def spacecraft_direction(self):
        """spacecraft direction (+1 ascending, -1 descending)"""
        if self._spacecraft_direction is None:
            direction = self.dataset.GetMetadata()["LOCATION"].lower()
            if direction == "ascending":
                self._spacecraft_direction = SPACECRAFT_ASCENDING
            else:
                self._spacecraft_direction = SPACECRAFT_DESCENDING
        return self._spacecraft_direction

    @property
    @returns(np.ndarray)
    def spacecraft_fixed_look_angle(self):
        """look angles in spacecraft-fixed coordinate system"""

        if self._spacecraft_fixed_look_angle is None:
            # Call necessary properties.
            dstep = self.scan_step_look_angle.rad
            xsize = self.scan_xsize
            # Compute look angles.
            dlook = -((xsize//2 - 0.5 - np.arange(xsize)[None, :]) * dstep)
            self._spacecraft_fixed_look_angle = dlook
        return self._spacecraft_fixed_look_angle

    @property
    @returns(np.ndarray)
    def spacecraft_fixed_unit_look_vector(self):
        """unit look vectors in spacecraft-fixed coordinate system"""
        return self._spacecraft_fixed_unit_look_vector

    @property
    @returns(np.ndarray)
    def spacecraft_geodetic_unit_vector(self):
        """unit vectors in spacecraft-geodetic coordinate system"""
        return self._spacecraft_geodetic_unit_vector

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
