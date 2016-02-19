"""
Routines for astronomical related calculation.

"""
import datetime

import numpy as np

import astropy.constants as const
import astropy.units as u

def beam_area(*args):
    """
    Calculate the Gaussian beam area.

    Parameters
    ----------
    args: float
        FWHM of the beam.
        If args is a single argument, a symmetrical beam is assumed.
        If args has two arguments, the two arguments are bmaj and bmin,
        the width of the major and minor axes of the beam in that order.

    Return
    ------
    out: float
        Beam area. No unit conversion is performed, i.e. the unit will depend
        on the input arguments. For example, beam width in degree wil return
        the beam area in square degree. Likewise, beam width in pixel will
        return the beam area in pixel.

    """
    if len(args) > 2:
        raise ValueError('Input argument must be a single beam width for a '
                         'symmetrical beam, or widths of the major and minor '
                         'axes of the beam.')
    if len(args) == 2:
        bmaj, bmin = args
    else:
        bmaj = args[0]
        bmin = bmaj
    return np.pi * bmaj * bmin / (4 * np.log(2))


def h2hms24(h):
    """
    Convert decimal hours to hh:mm:ss, rounding value to +0h to +24h.

    """
    hm, seconds = divmod((h % 24) * 3600, 60)
    hours, minutes = divmod(hm, 60)
    return '{:02.0f}:{:02.0f}:{:02.3f}'.format(hours, minutes, seconds)


def h2hms_signed(h):
    """
    Convert decimal hours to hh:mm:ss, rounding value to -24h to +24h.

    """
    inverse = ''
    if h < 0:
        inverse = '-'
    hms = str(datetime.timedelta(hours=abs(h))).rsplit(', ')[-1]
    return inverse + hms


def d2dms(d, delimiter=':', precision=6):
    """
    Convert decimal degrees to dd:mm:ss

    """
    inverse = ''
    if d < 0:
        inverse = '-'
    minutes, seconds = divmod(abs(d) * 3600, 60)
    degrees, minutes = divmod(minutes, 60)
    return '{0:s}{1:.0f}{4}{2:02.0f}{4}{3:0{5:d}.{6:d}f}'\
        .format(inverse, degrees, minutes, seconds, delimiter, 3 + precision,
                precision)


def lst2gha(lst, site_long=116.670456):
    """
    Convert LST in decimal degree to GHA in decimal hours.

    Longitude of the observer can be specify with site_long = +/- decimal
    degree, where + is east and - is west of the prime meridian. Else assume
    MWA 128T location, site_long = 116.670456.

    """
    gha = lst/15. - site_long/15.
    if gha > 24.0:
        gha -= 24.0
    elif gha < 0.0:
        gha += 24.0
    return gha


def jysr2k(intensity, freq):
    """
    Convert Jy/sr to K.

    Parameters
    ----------
    intensity: array-like
        Intensity (brightness) in Jy/sr
    freq: float
        Frequency of the map in MHz

    Return
    ------
    out: array-like or float
        Brightness temperature in Kelvin


    """
    ba = 1 * u.sr
    equiv = u.brightness_temperature(ba, freq * u.MHz)
    return (intensity * u.Jy).to(u.K, equivalencies=equiv).value


def k2jysr(temp, freq):
    """
    Convert K to Jy/sr.

    Parameters
    ----------
    temp: array-like
        Brightness temperature in Kelvin
    freq: float
        Frequency of the map in MHz

    Return
    ------
    out: array-like or float
        Intensity (brightness) in Jy/sr

    """
    ba = 1 * u.sr
    equiv = u.brightness_temperature(ba, freq * u.MHz)
    return (temp * u.K).to(u.Jy, equivalencies=equiv).value



def jybeam2k(intensity, freq, beam_width):
    """
    Convert Jy/beam to K.

    Parameters
    ----------
    intensity: array-like
        Intensity (brightness) in Jy/beam
    freq: float
        Frequency of the map in MHz
    beam_width: float
        The Gaussian FWHM width in degree

    Return
    ------
    out: array-like or float
        Brightness temperature in Kelvin

    """
    ba = beam_area(beam_width) * u.Unit('deg2')
    equiv = u.brightness_temperature(ba, freq * u.MHz)
    return (intensity * u.Jy).to(u.K, equivalencies=equiv).value


def k2jybeam(temp, freq, beam_width):
    """
    Convert K to Jy/beam.

    Parameters
    ----------
    temp: array-like
        Brightness temperature in Kelvin
    freq: float
        Frequency of the map in MHz
    beam_width: float
        The Gaussian FWHM width in degree

    Return
    ------
    out: array-like or float
        Intensity (brightness) in Jy/beam

    """
    ba = beam_area(beam_width) * u.Unit('deg2')
    equiv = u.brightness_temperature(ba, freq * u.MHz)
    return (temp * u.K).to(u.Jy, equivalencies=equiv).value
