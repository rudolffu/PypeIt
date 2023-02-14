"""
Module for the XLT/BFOSC instrument

.. include:: ../include/links.rst
"""
import numpy as np

from astropy.time import Time

from pypeit import msgs
from pypeit import telescopes
from pypeit import io
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.core import parse
from pypeit.images import detector_container


class XLTBfosc(spectrograph.Spectrograph):
    """
    Child to handle BFOSC specific code for each camera
    """
    ndet = 1
    telescope = telescopes.XLTTelescopePar()
    name = 'xlt_bfosc'
    camera = 'BFOSC'
    # supported = True
    comment = 'Grisms 4, 7'
    pypeline = 'MultiSlit'

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        Detector data from `here
        <http://www.not.iac.es/instruments/detectors/CCD14/>`__.

        .. warning::

            Many of the necessary detector parameters are read from the file
            header, meaning the ``hdu`` argument is effectively **required** for
            NOT/ALFOSC.  The optional use of ``hdu`` is only viable for
            automatically generated documentation.

        Args:
            det (:obj:`int`):
                1-indexed detector number.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent parameters are set to a default.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """
        # http://www.not.iac.es/instruments/detectors/CCD14/

        if hdu is None:
            binning = '1,1'
            gain = None
            ronoise = None
        else:
            # binning = self.get_meta_value(self.get_headarr(hdu), 'XBINNING')
            binning = '1,1'
            # gain = np.atleast_1d(hdu[0].header['GAIN'])  # e-/ADU
            # ronoise = np.atleast_1d(hdu[0].header['RDNOISE'])  # e-
            gain = np.atleast_1d(2.2)
            ronoise = np.atleast_1d(7.8)

        # Detector 1
        detector_dict = dict(
            binning         = binning,
            det             = 1,
            dataext         = 1,
            specaxis        = 0,
            specflip        = True,
            spatflip        = False,
            xgap            = 0.,
            ygap            = 0.,
            ysize           = 1.,
            platescale      = 0.2138,
            mincounts       = -1e10,
            darkcurr        = 1.3,      # e-/pix/hr
            saturation      = 700000.,  # ADU
            nonlinear       = 0.86,
            # datasec         = np.atleast_1d('[:,{}:{}]'.format(1, 2062)),  # Unbinned
            oscansec        = None,
            numamplifiers   = 1,
            gain            = gain,     # e-/ADU
            ronoise         = ronoise   # e-
        )

#        # Parse datasec, oscancsec from the header
#        head1 = hdu[1].header
#        detector_dict['gain'] = np.atleast_1d(head1['GAIN'])  # e-/ADU
#        detector_dict['ronoise'] = np.atleast_1d(head1['RDNOISE'])  # e-

        # Return
        return detector_container.DetectorContainer(**detector_dict)

    def configuration_keys(self):
        """
        Return the metadata keys that define a unique instrument
        configuration.

        This list is used by :class:`~pypeit.metadata.PypeItMetaData` to
        identify the unique configurations among the list of frames read
        for a given reduction.

        Returns:
            :obj:`list`: List of keywords of data pulled from file headers
            and used to constuct the :class:`~pypeit.metadata.PypeItMetaData`
            object.
        """
        return ['dispname', 'decker', 'binning']

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the ``PypeIt``-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='RA')
        self.meta['dec'] = dict(ext=0, card='DEC')
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['decker'] = dict(ext=0, card='FILTER')
        self.meta['binning'] = dict(card=None, compound=True)
        self.meta['mjd'] = dict(ext=0, card=None, compound=True)
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='WHEEL3')
        self.meta['idname'] = dict(ext=0, card='OBSTYPE')
        # self.meta['obstype'] = dict(ext=0, card='OBSTYPE')
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')

    def compound_meta(self, headarr, meta_key):
        """
        Methods to generate metadata requiring interpretation of the header
        data, instead of simply reading the value of a header card.

        Args:
            headarr (:obj:`list`):
                List of `astropy.io.fits.Header`_ objects.
            meta_key (:obj:`str`):
                Metadata keyword to construct.

        Returns:
            object: Metadata value read from the header(s).
        """
        if meta_key == 'binning':
            binspec, binspatial = [headarr[0]['XBINNING'], headarr[0]['YBINNING']]
            return parse.binning2string(binspec, binspatial)
        elif meta_key == 'mjd':
            ttime = Time(headarr[0]['DATE-OBS'], format='isot')
            return ttime.mjd
        else:
            msgs.error("Not ready for this compound meta")

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.

        Args:
            ftype (:obj:`str`):
                Type of frame to check. Must be a valid frame type; see
                frame-type :ref:`frame_type_defs`.
            fitstbl (`astropy.table.Table`_):
                The table with the metadata for one or more frames to check.
            exprng (:obj:`list`, optional):
                Range in the allowed exposure time for a frame of type
                ``ftype``. See
                :func:`pypeit.core.framematch.check_frame_exptime`.

        Returns:
            `numpy.ndarray`_: Boolean array with the flags selecting the
            exposures in ``fitstbl`` that are ``ftype`` type frames.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype == 'science':
            return good_exp & (fitstbl['idname'] == 'SPECLTARGET')
        # if ftype == 'standard':
        #     return good_exp & ((fitstbl['target'] == 'STD')
        #                         | (fitstbl['target'] == 'STD,SLIT'))
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'BIAS')
        if ftype in ['pixelflat', 'illumflat', 'trace']:
            return good_exp & (fitstbl['idname'] == 'SPECLFLAT')
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'SPECLLAMP')
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)
#    def pypeit_file_keys(self):
#        """
#        Define the list of keys to be output into a standard ``PypeIt`` file.
#
#        Returns:
#            :obj:`list`: The list of keywords in the relevant
#            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
#            :ref:`pypeit_file`.
#        """
#        return super().pypeit_file_keys() + ['slitwid']



