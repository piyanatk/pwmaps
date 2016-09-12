"""
Store system and configuration information for running MAPS.

"""
from __future__ import print_function, division

import os

import numpy as np


H21CM = 1420.40575177


class MAPS(object):
    """
    MAPS configurations.

    """
    try:
        ROOT_DIR = os.environ['SIM']
    except KeyError:
        print("'SIM' environment variable is not in your path.\nMAPS_DIR is "
              "temporary set to {:s}. You can set this path manually by "
              "calling `set_MPAS_DIR()` (not implemented yet)"
              .format(os.getcwd()))
        ROOT_DIR = os.getcwd()
    ARRAY_DIR = ROOT_DIR + '/array'
    ARRAY_CONFIG = {
        'mwa_128': ARRAY_DIR + '/mwa_128_crossdipole_gp_20110225.txt',
        'vla_d': ARRAY_DIR + '/VLA_D.txt'}
    ARRAY_LOC = {
        'mwa_128': ('-26.7033', '116.671', '377.830'),
        'vla_d': ('34.025778', '252.3210278', '2125.3704')}
    MAPS_GHA = {
        'mwa_128': -7.778066666666667,
        'vla_d': 16.821401853333334}


class MWA(object):
    EOR0 = (0.0, -30.0)
    EOR1 = (4.0 * 15.0, -30.0)
    EOR2 = (10.33 * 15.0, -10.0)
    ZENITH_DEC = -26.7033
    FREQ_EOR_LOW_40KHZ = np.arange(138.895, 167.055, 0.04)
    FREQ_EOR_HI_40KHZ = np.arange(167.055, 195.255, 0.04)
    FREQ_EOR_ALL_40KHZ = np.arange(138.895, 195.255, 0.04)
    FREQ_EOR_LOW_80KHZ = np.arange(138.915, 167.075, 0.08)
    FREQ_EOR_HI_80KHZ = np.arange(167.075, 195.275, 0.08)
    FREQ_EOR_ALL_80KHZ = np.arange(138.915, 195.275, 0.08)
