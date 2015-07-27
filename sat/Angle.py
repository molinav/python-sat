from __future__ import division
from . _constants.Angle import DEG
from . _constants.Angle import GRD
from . _constants.Angle import RAD
from . _constants.Angle import REV
from . _constants.Angle import DEG2RAD
from . _constants.Angle import GRD2RAD
from . _constants.Angle import RAD2RAD
from . _constants.Angle import REV2RAD
from . _constants.Angle import DEFAULT_ANGLE_VALUE
from . _constants.Angle import INFINITY
from . _constants.Angle import PI
from . _decorators import accepts
from . _decorators import limits
from . _decorators import returns
from numbers import Real


class Angle(object):
    """Base class which handles angles in different notations."""

    _angle_min = -INFINITY
    _angle_max = +INFINITY
    _limit_min = "["
    _limit_max = "]"

    def __init__(self, **kwargs):
        self._angle = DEFAULT_ANGLE_VALUE
        # Set angle value.
        num_args = len(kwargs)
        if num_args is 0:
            pass
        elif num_args is 1:
            key, val = list(kwargs.items())[0]
            try:
                factor = {
                    DEG: DEG2RAD, GRD: GRD2RAD, RAD: RAD2RAD, REV: REV2RAD}
                self.rad = val * factor[key]
            except KeyError:
                msg = "invalid key '{}' (expected '{}', '{}', '{}', '{}')"\
                    .format(key, DEG, GRD, RAD, REV)
                raise AttributeError(msg)
        else:
            msg = "too many arguments (expected 1)"
            raise AttributeError(msg)

    def norm(self):
        """return an Angle instance within the interval [0, 2*pi) rad"""
        val = self.rad % (2*PI)
        return Angle(rad=val)

    def wrap(self):
        """return an Angle instance within the interval (-pi, pi] rad"""
        val = self.rad % (2*PI)
        val = val - (2*PI) if val > PI else val
        return Angle(rad=val)

    def normed(self):
        """normalize an angle within the interval [0, 2*pi) rad"""
        val = self.rad % (2*PI)
        self.rad = val

    def wrapped(self):
        """wrap the instance value within the interval (-pi, pi] rad"""
        val = self.rad % (2*PI)
        self.rad = val - (2*PI) if val > PI else val

    def _return_angle(self, factor):
        """generic getter that returns the angle value"""
        return self._angle / factor

    def _update_angle(self, factor, val):
        """generic setter that updates the angle value"""

        value_min = self._angle_min / factor
        value_max = self._angle_max / factor

        @limits(self._limit_min, value_min, value_max, self._limit_max)
        def _calc_safe(s, v):
            s._angle = v * factor

        _calc_safe(self, val)

    @property
    @returns(Real)
    def deg(self):
        """angle value in degrees"""
        return self._return_angle(DEG2RAD)

    @deg.setter
    @accepts(Real)
    def deg(self, val):
        self._update_angle(DEG2RAD, val)

    @property
    @returns(Real)
    def grd(self):
        """angle value in gradians"""
        return self._return_angle(GRD2RAD)

    @grd.setter
    @accepts(Real)
    def grd(self, val):
        self._update_angle(GRD2RAD, val)

    @property
    @returns(Real)
    def rad(self):
        """angle value in radians"""
        return self._return_angle(RAD2RAD)

    @rad.setter
    @accepts(Real)
    def rad(self, val):
        self._update_angle(RAD2RAD, val)

    @property
    @returns(Real)
    def rev(self):
        """angle value in revolutions"""
        return self._return_angle(REV2RAD)

    @rev.setter
    @accepts(Real)
    def rev(self, val):
        self._update_angle(REV2RAD, val)


class AziAngle(Angle):
    """Angle subclass whose angle is restricted to [0, 2*pi] rad."""

    _angle_min = 0
    _angle_max = 2*PI
    _limit_min = "["
    _limit_max = ")"


class ZenAngle(Angle):
    """Angle subclass whose angle is restricted to [0, pi] rad."""

    _angle_min = 0
    _angle_max = PI
    _limit_min = "["
    _limit_max = "]"


class LatAngle(Angle):
    """Angle subclass whose angle is restricted to [-pi/2, +pi/2] rad."""

    _angle_min = -PI/2
    _angle_max = +PI/2
    _limit_min = "["
    _limit_max = "]"


class LonAngle(Angle):
    """Angle subclass whose angle is restricted to [-pi, +pi] rad."""

    _angle_min = -PI
    _angle_max = +PI
    _limit_min = "("
    _limit_max = "]"
