"""
This module provides and automated objects to run drift scan simulations
of the MWA using MAPS. It requires a python wrapper for MAPS (pwmaps.py).

"""
import os
import multiprocessing
from datetime import datetime

import numpy as np
import astropy.constants as const

from . import astro, pwmaps
from . import settings as s


# TODO: Get rid of this class.
class _InputError(Exception):
    """Class for input error exceptions.

    """
    def __init__(self, message, func_name, object_name):
        self.err = message
        self.func_name = func_name
        self.object_name = object_name

    def __str__(self):
        return '{0}: {1}.{2}:\n>>>> {3}'\
               .format(str(datetime.now()), self.object_name, self.func_name,
                       self.err)


class Drift:
    """
    This class provides an easy setup object for a drift scan simulation
    of the MWA EoR observation.


    """
    def __init__(self, target_ra, target_ha, sky_img=None, oobs=None,
                 pointing_center='zenith', fov_size=(412530.0, 412530.0),
                 duration=2.0, frequency=140.0, corr_int_time=1.0,
                 corr_chan_bw=0.04, scan_start='gha', site='MWA_128',
                 name=None, convert_k2jysr=False):
        """
        Initialize a drift scan.

        Parameters
        ----------
        target_ra: float
            Right ascension at the center of a star/target to be observed
            [decimal hours]
        target_ha: float
            Hour angle of the star/target from the zenith meridian
            at the star of an observation [decimal hours]
        sky_img: string
            Name and path of the input sky image. The image must be a SIN
            projection of equal size in FITS format
        oobs: string
            Out-of-bound source list file.
            For example, a file containing,
                # RA(decimal hour) Dec(decimal degree) I Q U V
                0.0 -26.0 20.0 0.0 0.0 0.0
            will add an unpolarized source of 20.0 Jy at (0.0h -26.0d) to the
            visibility.
        pointing_center: 'zenith' or ('hh:mm:ss', 'dd:mm:ss')
            Pointing center of the drift scan. Default mode 'zenith' uses
            the zenith calculated from target_ra and target_ha as a pointing.
        fov_size: array-like of float, 2 elements, optional
            Angular size of the field of view (RA[arcseconds, Dec[arcseconds])
            This is equivalent to the angular size of the sky image.
            The default value is 2 radians which produce a whole sky image
            in SIN projection.
        duration: float
            Duration of an observation [second]
        frequency: float
            Observing frequency at the center of the channel [MHz]
        corr_int_time: float
            Correlator integration time [second]
        corr_chan_bw: float
            Correlator channel bandwidth [MHz]
        scan_start: 'gha' or 'year:day-of-year:hour:minute:second'
            Start time of a scan. If 'gha', the program will determine a
            Greenwich hour angle from the site keyword and use GHA convention
            to use relative time hard-coded in visgen.
            Else, user can give an absolute by giving a string of
            'year:day-of-year:hour:minute:second'. When specifying the absolute
            time, visgen will calculate the HA from this using the array
            location, but precision effects such as precession and nutation
            are not included.
        site: 'MWA_128', 'VLA_D'
            Defining array to use.
        name: string, optional
            Name of the observation. The name of sky_img - ".fits" with out path
            will be use if None
        convert_k2jysr: {True, False}
            Perform conversion on sky_img from Kelvin to Jy/sr.

        """
        # TODO: assert that eitehr sky_img or oobs exist
        # TODO: Use astropy.coordinates
        self._ha = target_ha
        self._ra = target_ra
        self.site = site
        if pointing_center == 'zenith':
            self.fov_center_ra = astro.h2hms24(self._ra + self._ha)
            self.fov_center_dec = astro.d2dms(float(s.MAPS.ARRAY_LOC[self.site.lower()][0]))
        else:
            self.fov_center_ra = pointing_center[0]
            self.fov_center_dec = pointing_center[1]
        self.fov_size_ra = str(fov_size[0])
        self.fov_size_dec = str(fov_size[1])
        self._center_frequency = frequency
        self.frequency = str(self._center_frequency - corr_chan_bw / 2.)
        self.corr_int_time = str(corr_int_time)
        self.corr_chan_bw = str(corr_chan_bw)
        self.channel = self.frequency + ':' + self.corr_chan_bw
        # TODO: recalculate and check GHA
        if scan_start == 'gha':
            self.scan_start = 'GHA {0:f}'.format(s.MAPS.MAPS_GHA[site.lower()])
        else:
            self.scan_start = scan_start
        self.scan_duration = str(duration)
        # We do not need time and frequency average.
        self.time_cells = '0'
        self.freq_cells = '0'
        if name is None:
            if sky_img is not None:
                self.name = sky_img.rsplit('/', 1)[-1][0:-5]
            else:
                self.name = 'visgen_{:.2f}h_{:.2f}ha_{:.3f}MHz_{:.3f}kHz_{:.2f}sec'\
                    .format(target_ra, target_ha, frequency, corr_chan_bw, duration)
        else:
            self.name = name
        self.sky_img = sky_img
        self.oobs = oobs
        self.spec_file = None
        self.vis_in = None
        self.vis_out = None
        self.vislog = None
        self.uvfits = None
        self.convert_k2jysr = convert_k2jysr
        self.__spec = ''
        self.update_spec()
        self.__log = ''

    def __str__(self):
        return self.__spec + self.__log

    def update_spec(self):
        speckeys = np.array(['FOV_center_RA', 'FOV_center_Dec', 'FOV_size_RA',
                             'FOV_size_Dec', 'Corr_int_time', 'Corr_chan_bw',
                             'Scan_start', 'Scan_duration', 'Channel'])
        header = ('# {0}\n'
                  '# MAPS drift scan simulation\n'
                  '# name: {1}\n'
                  '# sky image: {2}\n'
                  '# OOB sources: {3}\n'
                  '# sky uvgrid: {4}\n'
                  '# visibility: {5}\n'
                  '# uvfits: {6}\n'
                  '# visgen log: {7}\n'
                  '# visgen specification file: {8}\n'
                  .format(str(datetime.now()), self.name, self.sky_img,
                          self.oobs, self.vis_in, self.vis_out, self.uvfits,
                          self.vislog, self.spec_file))
        spec = ''
        for k in speckeys:
            spec += '{0} = {1}\n'.format(k, self.__dict__.get(k.lower()))
        self.__spec = header + spec + 'Endscan\n\n'

    def print_spec(self):
        print self.__spec

    def write_spec(self):
        """
        Write observation specification file (*.ospec)

        Parameters
        ----------
        filename: str, optional
            Name of the ospec file
        verbose: boolean, optional
            Print the ospec file if True

        """
        outfile = self.name + '.ospec'
        self.spec_file = outfile
        self.update_spec()
        self.append_log('# $> write_spec()\n'
                        '# >>> visgen spec: {0}\n'
                        .format(self.spec_file))
        with open(outfile, 'w') as f:
            f.write(self.__spec)

    def append_log(self, string):
        self.__log += '# {0}\n{1}\n'.format(str(datetime.now()), string)

    def print_log(self):
        print self.__str__()

    def write_log(self):
        self.update_spec()
        self.append_log('# $> write_log()\n'
                        '# >>> log file: {0}\n'
                        .format(self.name + '.log'))
        with open(self.name + '.log', 'w') as f:
            f.write(self.__str__())

    def im2uv(self):
        if self.sky_img is None:
            raise _InputError('No imagae file', self.im2uv.__name__,
                              self.__name__)
        else:
            print '# im2uv: ' + self.name
            if self.convert_k2jysr:
                normalizer = 2 * (self._center_frequency * 1e6) ** 2 \
                    * const.k_B.si.value / (const.c.si.value ** 2)
            else:
                normalizer = None
            pwmaps.im2uv(self.sky_img, verbose=False, normalizer=normalizer)
            self.vis_in = self.sky_img.rsplit('/', 1)[-1][0:-5] + '.dat'
            self.update_spec()
            self.append_log('# $> im2uv({0})\n'
                            '# >>> sky uvgrid: {1}\n'
                            .format(self.sky_img, self.vis_in))

    def visgen(self, mpi=1):
        if self.spec_file is None:
            raise _InputError('No oobs file')
        if self.vis_in is None and self.oobs is None:
            raise _InputError('Neither uvgrid file nor oob source list exist.',
                              self.visgen.__name__, self.name)
        else:
            print '# visgen: ' + self.name
            pwmaps.visgen(self.name, self.spec_file, oobs=self.oobs,
                        uvgrid=self.vis_in, mpi=mpi, site=self.site)
            self.vis_out = self.name + '.vis'
            self.vislog = self.name + '.vislog'
            self.update_spec()
            self.append_log('# $> visgen()\n'
                            '# >>>> visgen visibility: {0}\n'
                            '# >>>> visgen log file: {1}\n'
                            .format(self.vis_out, self.vislog))

    def maps2uvfits(self):
        if self.vis_out is None:
            raise _InputError('visibility from visgen is not present',
                              self.maps2uvfits.__name__, self.name)
        else:
            print '# maps2uvfits: ' + self.name
            pwmaps.maps2uvfits(self.vis_out, site=self.site, verbose=False)
            self.uvfits = self.name + '.uvfits'
            self.update_spec()
            self.append_log('# $> maps2uvfits({0})\n'
                            '# >>>> uvfits: {1}\n'
                            .format(self.vis_out, self.uvfits))

    def run(self):
        # TODO: Need to check if input exist
        if self.sky_img is not None:
            self.im2uv()
        self.write_spec()
        self.visgen()
        if self.vis_in is not None:
            os.remove(self.vis_in)
            self.append_log('# remove ' + self.vis_in)
        self.maps2uvfits()
        self.write_spec()
        self.write_log()


def __call_go(instance):
    """
    Wrapper to make Drift.run() pickle-able.

    """
    return instance.run()


def batch_drift(instance, nprocs=4):
    pool = multiprocessing.Pool(nprocs)
    pool.map(__call_go, instance)
    pool.close()
    pool.join()
