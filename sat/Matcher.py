from __future__ import division
from . constants.Angle import DEG2RAD
from . constants.Matcher import MAX_ATTITUDE_ANGLE
from . constants.Matcher import MAX_DISTANCE
from . constants.Matcher import MAX_HAMMING
from . constants.Scene import HRPT_CLOUD_THRESHOLD_B5
from . constants.Scene import SPACECRAFT_ASCENDING
from . decorators import returns
from . Coordinates import Coordinates
from . Error import MatcherError
from . Scene import Scene
from datetime import timedelta
from scipy.ndimage import binary_dilation
import cv2
import gdal
import os.path
import numpy as np


class Matcher(object):

    __slots__ = [
        "_bfmatcher",
        "_descriptor",
        "_detector",
        "_match_table",
        "_source",
        "_source_cloudmask",
        "_source_descriptors",
        "_source_gtr",
        "_source_keypoints",
        "_source_matchband",
        "_target",
        "_target_cloudmask",
        "_target_descriptors",
        "_target_gtr",
        "_target_keypoints",
        "_target_matchband",
        ]

    def __init__(self, source, fast_lim=40, brisk_lim=50):
        """Constructor of a generic Matcher instance.

        Parameters:

        source
            code for reference dataset

        Optional parameters:

        fast_lim
            limit value introduced in the FAST constructor
        brisk_lim
            limit value introduced in the BRISK constructor
        """

        for item in self.__slots__:
            self.__setattr__(item, None)

        # Create detector, descriptor and Brute-Force matcher.
        try:
            self._detector = cv2.FastFeatureDetector(fast_lim)
        except AttributeError:
            self._detector = cv2.FastFeatureDetector_create(fast_lim)
        try:
            self._descriptor = cv2.BRISK(thresh=brisk_lim)
        except AttributeError:
            self._descriptor = cv2.BRISK_create(thresh=brisk_lim)
        self._bfmatcher = cv2.BFMatcher(cv2.NORM_HAMMING, crossCheck=True)

        # Build source path.
        dirname = os.path.dirname(__file__)
        srcpath = os.path.join(dirname, "constants", "_matcher", source)

        # Open the source and build its attributes.
        gdal.PushErrorHandler("CPLQuietErrorHandler")
        self._source = gdal.Open(srcpath)
        self._build_source_geotransform()

    def draw_matches(self, path, color=(0, 255, 255), radius=5, thickness=1):
        """Export matches and the matches bands to an image file.

        Parameters:

        path
            output path where the image will be stored

        Optional parameters:

        color
            BGR tuple defining the color of circles and lines
            (default (0, 255, 255) = yellow)
        radius
            circle radius (default 5)
        thickness
            circle and line thickness (default 1)
        """

        # Call necessary properties.
        srcrows, srccols = self.source_matchband.shape
        tgtrows, tgtcols = self.target_matchband.shape

        # Create the output image with all its values set to 0.
        shp = (max(srcrows, tgtrows), srccols+tgtcols, 3)
        img = np.zeros(shp, dtype=np.uint8)

        # Place the source and target match bands in the output image.
        img[:srcrows, :srccols] = np.dstack(3*[self.source_matchband])
        img[:tgtrows, srccols:] = np.dstack(3*[self.target_matchband])

        # For each match, draw yellow circles over the matching pixels and a
        # line between them.
        for elem in self.match_table:
            # Read points.
            y1, x1 = elem[0:2].astype(int)
            y2, x2 = elem[2:4].astype(int)
            # Draw a circle representing every point.
            cv2.circle(img, (x1, y1), radius, color, thickness)
            cv2.circle(img, (x2+srccols, y2), radius, color, thickness)
            # Draw a line between the two circles.
            cv2.line(img, (x1, y1), (x2+srccols, y2), color, thickness)

        # Save the image to a file.
        cv2.imwrite(path, img)

    def draw_source_keypoints(self, path):
        """Export source keypoints and the match band to an image file.

        Parameters:

        path
            output path where image will be stored
        """
        self.__class__._draw_keypoints(
            path, self.source_matchband, self.source_keypoints)

    def draw_target_keypoints(self, path):
        """Export target keypoints and the match band to an image file.

        Parameters:

        path
            output path where the image will be stored
        """
        self.__class__._draw_keypoints(
            path, self.target_matchband, self.target_keypoints)

    def fix_attitude(self, max_dist=MAX_DISTANCE, max_hamm=MAX_HAMMING):
        """Evaluate optimal values for target's attitude angles.

        Parameters:

        max_dist
            maximum allowed geographic distance to validate a matching
        max_hamm
            maximum allowed descriptor distance to validate a matching
        """

        # Verify that the target's attitude angles have not been evaluated yet.
        if not np.allclose(self.target.spacecraft_attitude, 0):
            msg = "Target spacecraft already has attitude angles defined"
            raise MatcherError(msg)

        # Call finder.
        new_r_pix = self._find_matches(max_dist, max_hamm)[0]

        # Calculate new look unit vector.
        num = len(new_r_pix)
        new_r_sat = np.empty(new_r_pix.shape)
        new_d = np.empty(new_r_pix.shape)
        new_ub = np.empty(new_r_pix.shape)
        new_ug = np.empty(new_r_pix.shape)
        for i in range(num):
            # Get the row and column on the spacecraft original imagery.
            row1, col1 = self.match_table[i, [2, 3]].astype(int)
            row2, col2 = self.target.hrpt_pix[:, row1, col1]
            # Get the new look vector in ECF coordinate system.
            new_r_sat[i] = self.target.orbit.position_ecf[row2, :]
            new_d[i] = new_r_pix[i] - new_r_sat[i]
            # Transform look vector into spacecraft geodetic coordinate system.
            rot_T = self.target.spacecraft_geodetic_unit_vector[:, :, row2].T
            new_d[i] = np.dot(new_d[i], rot_T)
            # Calculate new spacecraft-geodetic and spacecraft-fixed unit look
            # vectors
            new_ug[i] = new_d[i] / np.linalg.norm(new_d[i])
            new_ub[i] = self.target.spacecraft_fixed_unit_look_vector[:, col2]

        # Calculate optimal attitude values which relates spacecraft-geodetic
        # and spacecraft-fixed unit look vectors by LSM.
        if num is not 1:
            cnt = -1
            flag = True
            while flag:
                values = list(range(-(cnt < 0), cnt)) + list(range(cnt+1, num))
                ux, uy, uz = new_ug.T[:, values]
                cd, nd, sd = new_ub.T[:, values]
                # Create temporary values for LSM solution.
                g1 = + np.sum(sd * sd)
                g2 = - np.sum(sd * cd)
                g3 = + np.sum(cd * cd)
                h1 = - np.sum(uy * sd)
                h2 = + np.sum(ux * sd - uz * cd)
                h3 = + np.sum(uy * cd)
                # Calculate optimal roll, pitch and yaw.
                r_angle = h2 / len(cd)
                p_angle = (g1*h3 - g2*h1) / (g1*g3 - g2*g2)
                y_angle = (g3*h1 - g2*h3) / (g1*g3 - g2*g2)
                attitude = np.asarray([r_angle, p_angle, y_angle])
                # Check the values are within the valid range.
                flag = np.any(np.abs(attitude) > MAX_ATTITUDE_ANGLE)
                cnt += 1
                # Alternative if the common LSM method does not give a solution.
                if cnt == num + 1:
                    num = 1
                    break
        if num is 1:
            ux, uy, uz = new_ug.T
            cd, nd, sd = new_ub.T
            # Estimate roll, pitch and yaw minimizing pitch and yaw.
            r_angle = + ux * sd - uz * cd
            p_angle = + uy * cd
            y_angle = - uy * sd
            # Choose the median values.
            attitude = np.asarray([r_angle, p_angle, y_angle])
            attitude = np.min(attitude, axis=1).flatten()
            flag = np.any(np.abs(attitude) > MAX_ATTITUDE_ANGLE)
            if flag:
                msg = "attitude angles outside range"
                raise MatcherError(msg)

        # Store the attitude angles in the Scene target and return their values.
        self.target._spacecraft_attitude = attitude
        return attitude

    def fix_clock_drift(self, max_dist=MAX_DISTANCE, max_hamm=MAX_HAMMING):
        """Evaluate optimal value for the target's clock drift's optimal.

        Parameters:

        max_dist
            maximum allowed geographic distance to validate a matching
        max_hamm
            maximum allowed descriptor distance to validate a matching
        """

        # Call finder.
        dr_xyz_alongtrack = self._find_matches(max_dist, max_hamm)[1]

        # Calculate velocity vectors for every ECF distance vector.
        num = len(dr_xyz_alongtrack)
        v_sat_alongtrack = np.empty(dr_xyz_alongtrack.shape)
        for i in range(num):
            # Get the row and column on the spacecraft original imagery.
            row1, col1 = self.match_table[i, [2, 3]].astype(int)
            row2, col2 = self.target.hrpt_pix[:, row1, col1]
            # Get the velocity vector.
            v_sat_i = self.target.orbit.velocity_ecf[row2, :]
            old_wy_i = self.target.spacecraft_geodetic_unit_vector[1, :, row2]
            v_sat_alongtrack[i] = np.sum(v_sat_i * old_wy_i)

        # The list of vectors dr_xyz is related to the list of velocities v_sat
        # so that v_sat * dt = - dr_xyz, where dt is the elapsed time. Solve it
        # by using LSM method: A*X = B.
        aa = + v_sat_alongtrack[:, None]
        bb = - dr_xyz_alongtrack[:, None]
        xx = np.linalg.lstsq(aa, bb)[0]
        rs = np.abs(np.dot(aa, xx) - bb)

        # Repeat the LSM method with the most reliable matches.
        try:
            aa = aa[rs < 500]
            bb = bb[rs < 500]
            xx = np.linalg.lstsq(aa, bb)[0]
        except np.linalg.linalg.LinAlgError:
            pass

        # Calculate clock drift as a timedelta instance.
        dtime = xx[0, 0]
        drift = timedelta(seconds=dtime)

        # Update target clock drift and return drift value.
        self.target._clock_drift += drift
        return drift

    def set_target(self, target):
        """Add a Scene instance to be matched and compute its attributes.

        Parameters:

        target
            a Scene instance that needs matching correction
        """

        # Open the target and build its attributes.
        self._target = target
        self._build_target_cloudmask()
        self._build_target_geotransform()
        self._build_source_cloudmask()

        # Build the match bands for both the source and the target.
        self._build_matchbands()
        self._build_keypoints()
        self._build_descriptors()

    def _build_source_cloudmask(self):
        """Create cloud mask based on channel 5 brightness temperature."""

        # Call necessary properties.
        k = np.ones((3, 3))
        xsize = self.source.RasterXSize
        ysize = self.source.RasterYSize
        b5min = HRPT_CLOUD_THRESHOLD_B5

        # Create cloud mask for source dataset and dilate.
        src_cld = self.source.GetRasterBand(5).ReadAsArray(
            0, 0, xsize, ysize).astype(np.int16)
        src_cld = binary_dilation((src_cld > b5min), k, iterations=20)

        self._source_cloudmask = 255 * src_cld.astype(np.uint8)[None, :]

    def _build_target_cloudmask(self):
        """Create cloud mask based on channel 5 brightness temperature."""

        # Call necessary properties.
        k = np.ones((3, 3))

        # Dilate the already existing cloud mask from target dataset.
        tgt_cld = self.target.hrpt_msk[0]
        tgt_cld = binary_dilation(tgt_cld, k, iterations=20)

        self._target_cloudmask = 255 * tgt_cld.astype(np.uint8)[None, :]

    def _build_source_geotransform(self):
        """Create geotransform array for source dataset."""

        srcgtr = np.asarray(self.source.GetGeoTransform()).reshape((2, 3))
        self._source_gtr = srcgtr

    def _build_target_geotransform(self):
        """Create geotransform array for target dataset."""
        self._target_gtr = np.asarray(self.target.hrpt_gtr).reshape((2, 3))

    def _build_matchbands(self):
        """Create match bands based on channels 2/3 for day/night."""

        # Get index of the dataset band used in matching (it starts in 1).
        idx = 2 + int(self.target.spacecraft_direction is SPACECRAFT_ASCENDING)
        xsize = self.source.RasterXSize
        ysize = self.source.RasterYSize

        def enclose(band, qmin, qmax):
            cmin = np.percentile(band[band > 0], qmin)
            cmax = np.percentile(band[band > 0], qmax)
            return np.maximum(0, np.minimum(255, 255/(cmax-cmin) * (band-cmin)))

        # Create match band for source dataset.
        srcband = self.source.GetRasterBand(idx).ReadAsArray(
            0, 0, xsize, ysize).astype(np.int16)
        srcband[srcband < 0] = 0
        srcband = enclose(srcband, 1, 99)
        srcband = np.dstack(3 * [srcband]).astype(np.uint8)
        srcband = cv2.cvtColor(srcband, cv2.COLOR_BGR2GRAY)
        self._source_matchband = srcband

        # Create match band for target dataset.
        tgtband = self.target.hrpt_img[idx-1]
        tgtband[tgtband < 0] = 0
        tgtband = enclose(tgtband, 1, 99)
        tgtband = np.dstack(3 * [tgtband]).astype(np.uint8)
        tgtband = cv2.cvtColor(tgtband, cv2.COLOR_BGR2GRAY)
        self._target_matchband = tgtband

    def _build_keypoints(self):
        """Create lists of keypoints using the FAST detector."""

        k = np.ones((3, 3))
        srcclds = self.source_cloudmask[0]
        tgtclds = self.target_cloudmask[0]
        srcmask = binary_dilation((self.source_matchband == 0), k, iterations=5)
        tgtmask = binary_dilation((self.target_matchband == 0), k, iterations=5)

        # Create keypoints for source dataset.
        srckps = self.detector.detect(self.source_matchband, None)
        self._source_keypoints = np.asarray([
            x for x in srckps if
            not srcclds[x.pt[1], x.pt[0]] and not srcmask[x.pt[1], x.pt[0]]])

        # Create keypoints for target dataset.
        tgtkps = self.detector.detect(self.target_matchband, None)
        self._target_keypoints = np.asarray([
            x for x in tgtkps if
            not tgtclds[x.pt[1], x.pt[0]] and not tgtmask[x.pt[1], x.pt[0]]])

    def _build_descriptors(self):
        """Create lists of descriptors using the BRISK algorithm."""

        # Create source descriptors.
        srckps, srcdes = self.descriptor.compute(
            self.source_matchband, self.source_keypoints)
        self._source_keypoints = np.asarray(srckps)
        self._source_descriptors = srcdes

        # Create target descriptors.
        tgtkps, tgtdes = self.descriptor.compute(
            self.target_matchband, self.target_keypoints)
        self._target_keypoints = np.asarray(tgtkps)
        self._target_descriptors = tgtdes

    @staticmethod
    def _draw_keypoints(path, img, kps):
        """Export keypoints over the dataset to an image file."""

        fmt = cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS
        try:
            img2 = cv2.drawKeypoints(img, kps, flags=fmt)
        except TypeError:
            img2 = np.empty(img.shape, dtype=np.uint8)
            img2 = cv2.drawKeypoints(img, kps, img2, flags=fmt)
        cv2.imwrite(path, img2)

    def _find_matches(self, max_dist, max_hamm):
        """Find feature matches between the source and the target."""

        # Call necessary properties.
        srckp1 = self.source_keypoints
        tgtkp1 = self.target_keypoints
        srcdes = self.source_descriptors
        tgtdes = self.target_descriptors
        srcgtr = self.source_gtr
        tgtgtr = self.target_gtr

        # Call brute-force matcher.
        m = self.bfmatcher.match(srcdes, tgtdes)
        num = len(m)
        if num is 0:
            msg = "No matches source-target were found"
            raise MatcherError(msg)
        m = sorted(m, key=lambda x: x.distance)
        match_difference = np.asarray([x.distance for x in m])
        srckp2 = np.asarray([[y+0.5 for y in srckp1[x.queryIdx].pt] for x in m])
        tgtkp2 = np.asarray([[y+0.5 for y in tgtkp1[x.trainIdx].pt] for x in m])

        # Convert pixel coordinates to geodetic and ECF coordinates.
        alt = np.zeros((len(m), 1))
        new_lonlat = srcgtr[:, 0][None, :] + np.dot(srckp2, srcgtr[:, 1:].T)
        old_lonlat = tgtgtr[:, 0][None, :] + np.dot(tgtkp2, tgtgtr[:, 1:].T)
        new_xyz = Coordinates.from_geo_to_ecf(
            np.hstack([new_lonlat[:, [1, 0]], alt]) * DEG2RAD)
        old_xyz = Coordinates.from_geo_to_ecf(
            np.hstack([old_lonlat[:, [1, 0]], alt]) * DEG2RAD)
        distance_xyz = old_xyz - new_xyz

        # Compute geographic and temporal distances between matches.
        distance_xyz_alongtrack = np.empty((num,))
        time_difference = np.empty((num,))
        for i in range(num):
            # Get the row and column on the spacecraft original imagery.
            row1, col1 = tgtkp2[i, [1, 0]].astype(int)
            row2, col2 = self.target.hrpt_pix[:, row1, col1]
            # Compute the along-track geographic distance.
            old_wy_i = self.target.spacecraft_geodetic_unit_vector[1, :, row2]
            distance_xyz_alongtrack[i] = np.sum(distance_xyz[i] * old_wy_i)
            # Compute the time difference between matches positions.
            dstep = self.target.scan_timedelta.total_seconds() / 1100
            time_difference[i] = distance_xyz_alongtrack[i] * dstep

        # Create match table with information about feature matches.
        match_table = np.hstack([
            srckp2[:, [1, 0]],
            tgtkp2[:, [1, 0]],
            new_lonlat[:, [1, 0]],
            old_lonlat[:, [1, 0]],
            match_difference[:, None],
            np.linalg.norm(distance_xyz, axis=1)[:, None],
            distance_xyz_alongtrack[:, None],
            time_difference[:, None]
            ])
        # Validate matches by using their geographic distances.
        mask1 = np.abs(match_table[:, 9]) < max_dist
        match_table = match_table[mask1]
        if len(match_table) is 0:
            msg = "No matches source-target after distance restraint"
            raise MatcherError(msg)
        # Validate matches by using their descriptor distances.
        mask2 = match_table[:, 8] < max_hamm
        match_table = match_table[mask2]
        if len(match_table) is 0:
            msg = "No matches source-target after distance+descriptor restraint"
            raise MatcherError(msg)

        self._match_table = match_table

        return new_xyz[mask1][mask2], distance_xyz_alongtrack[mask1][mask2]

    @property
    def bfmatcher(self):
        """brute-force matcher"""
        return self._bfmatcher

    @property
    def descriptor(self):
        """BRISK descriptor"""
        return self._descriptor

    @property
    def detector(self):
        """FAST detector"""
        return self._detector

    @property
    @returns(np.ndarray)
    def match_table(self):
        """table containing all the relevant matches attributes"""
        return self._match_table

    @property
    @returns(gdal.Dataset)
    def source(self):
        """source dataset used as reference"""
        return self._source

    @property
    @returns(np.ndarray)
    def source_cloudmask(self):
        """cloud mask for source dataset"""
        return self._source_cloudmask

    @property
    @returns(np.ndarray)
    def source_descriptors(self):
        """descriptors from source dataset keypoints"""
        return self._source_descriptors

    @property
    @returns(np.ndarray)
    def source_gtr(self):
        """tuple geotransform for source dataset"""
        return self._source_gtr

    @property
    @returns(np.ndarray)
    def source_keypoints(self):
        """keypoints from source dataset"""
        return self._source_keypoints

    @property
    @returns(np.ndarray)
    def source_matchband(self):
        """source dataset's band used in matching process"""
        return self._source_matchband

    @property
    @returns(Scene)
    def target(self):
        """target dataset to be geometrically corrected"""
        return self._target

    @property
    @returns(np.ndarray)
    def target_cloudmask(self):
        """cloud mask for target dataset"""
        return self._target_cloudmask

    @property
    @returns(np.ndarray)
    def target_descriptors(self):
        """descriptors from target dataset keypoints"""
        return self._target_descriptors

    @property
    @returns(np.ndarray)
    def target_gtr(self):
        """tuple geotransform for target dataset"""
        return self._target_gtr

    @property
    @returns(np.ndarray)
    def target_keypoints(self):
        """keypoints from target dataset"""
        return self._target_keypoints

    @property
    @returns(np.ndarray)
    def target_matchband(self):
        """target dataset's band used in matching process"""
        return self._target_matchband
