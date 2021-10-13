"""
Module for Keck/MOSFIRE specific methods.

.. include:: ../include/links.rst
"""
import os
from pkg_resources import resource_filename

from IPython import embed

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch, meta
from pypeit import utils
from pypeit import io
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container
from scipy import special
from pypeit.spectrographs.slitmask import SlitMask

from pypeit.utils import index_of_x_eq_y

class KeckMOSFIRESpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Keck/MOSFIRE specific code
    """
    ndet = 1
    name = 'keck_mosfire'
    telescope = telescopes.KeckTelescopePar()
    camera = 'MOSFIRE'
    supported = True
    comment = 'Gratings tested: Y, J, K'

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

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
        # Detector 1
        detector_dict = dict(
            binning         = '1,1',
            det             = 1,
            dataext         = 0,
            specaxis        = 1,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.1798,
            darkcurr        = 0.8,
            saturation      = 1e9, # ADU, this is hacked for now
            nonlinear       = 1.00,  # docs say linear to 90,000 but our flats are usually higher
            numamplifiers   = 1,
            mincounts       = -1e10,
            gain            = np.atleast_1d(2.15),  # Taken from MOSFIRE detector webpage
            ronoise         = np.atleast_1d(5.8), # This is for 16 non-destructuve reads, the default readout mode
            datasec         = np.atleast_1d('[5:2044,5:2044]'),
            #oscansec        = np.atleast_1d('[:,:]')
        )
        return detector_container.DetectorContainer(**detector_dict)

    # def get_rawimage(self, raw_file, det):
    #     """
    #     Read raw images and generate a few other bits and pieces
    #     that are key for image processing.
    #
    #     Over-ride standard get_rawimage() for MOSFIRE to deal
    #     with long2pos slitmask
    #
    #     Parameters
    #     ----------
    #     raw_file : :obj:`str`
    #         File to read
    #     det : :obj:`int`
    #         1-indexed detector to read
    #
    #     Returns
    #     -------
    #     detector_par : :class:`pypeit.images.detector_container.DetectorContainer`
    #         Detector metadata parameters.
    #     raw_img : `numpy.ndarray`_
    #         Raw image for this detector.
    #     hdu : `astropy.io.fits.HDUList`_
    #         Opened fits file
    #     exptime : :obj:`float`
    #         Exposure time read from the file header
    #     rawdatasec_img : `numpy.ndarray`_
    #         Data (Science) section of the detector as provided by setting the
    #         (1-indexed) number of the amplifier used to read each detector
    #         pixel. Pixels unassociated with any amplifier are set to 0.
    #     oscansec_img : `numpy.ndarray`_
    #         Overscan section of the detector as provided by setting the
    #         (1-indexed) number of the amplifier used to read each detector
    #         pixel. Pixels unassociated with any amplifier are set to 0.
    #     """
    #
    #     detector, raw_img, hdu, exptime, rawdatasec_img, oscansec_img = super().get_rawimage(raw_file, det)
    #
    #     headarr = self.get_headarr(raw_file)
    #
    #     if 'long2pos' in self.get_meta_value(headarr, 'decker'):
    #         # Select only the 3 slits used for long2pos and neglect the others that are random
    #         rawdatasec_img[:883, :] = 0
    #         rawdatasec_img[1195:, :] = 0
    #
    #     # Return
    #     return detector, raw_img, hdu, exptime, rawdatasec_img, oscansec_img

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of ``PypeIt`` methods.
        """
        par = super().default_pypeit_par()

        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.30 #0.20  # Might be grating dependent..
        par['calibrations']['wavelengths']['sigdetect']=5.0
        par['calibrations']['wavelengths']['fwhm']= 5.0
        par['calibrations']['wavelengths']['n_final']= 4
        par['calibrations']['wavelengths']['lamps'] = ['OH_NIRES']
        #par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['method'] = 'holy-grail'
        # Reidentification parameters
        #par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_nires.fits'
        par['calibrations']['slitedges']['edge_thresh'] = 50.
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Flats
        # Do not illumination correct. We should also not be flat fielding given the bars.
        # TODO Implement imaging flats for MOSFIRE. Do test with/without illumination flats.
        # Turn of illumflat
        turn_off = dict(use_biasimage=False, use_overscan=False, use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)

        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['extraction']['sn_gauss'] = 4.0

        # Flexure
        par['flexure']['spec_method'] = 'skip'

        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] ='nothing'

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 20]
        par['calibrations']['arcframe']['exprng'] = [1, None]
        par['calibrations']['darkframe']['exprng'] = [20, None]
        par['scienceframe']['exprng'] = [20, None]

        # Sensitivity function parameters
        par['sensfunc']['extrap_blu'] = 0.0  # Y-band contaminated by higher order so don't extrap much
        par['sensfunc']['extrap_red'] = 0.0
        par['fluxcalib']['extrap_sens'] = True
        par['sensfunc']['extrap_red'] = 0.0
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 13
        par['sensfunc']['IR']['maxiter'] = 2
        par['sensfunc']['IR']['telgridfile'] \
                = os.path.join(par['sensfunc']['IR'].default_root,
                               'TelFit_MaunaKea_3100_26100_R20000.fits')
        return par

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the ``PypeIt`` parameters to hard-wired values used for
        specific instrument configurations.

        Args:
            scifile (:obj:`str`):
                File to use when determining the configuration and how
                to adjust the input parameters.
            inp_par (:class:`~pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`~pypeit.par.parset.ParSet`: The PypeIt parameter set
            adjusted for configuration specific parameter values.
        """
        par = super().config_specific_par(scifile, inp_par=inp_par)

        headarr = self.get_headarr(scifile)

        # Turn on the use of mask design
        if 'LONGSLIT' not in self.get_meta_value(headarr, 'decker'):
            par['calibrations']['slitedges']['use_maskdesign'] = True
            # Assign RA, DEC, OBJNAME to detected objects
            par['reduce']['slitmask']['assign_obj'] = True
            # force extraction of undetected objects
            par['reduce']['slitmask']['extract_missing_objs'] = True
            # needed for better slitmask design matching
            par['calibrations']['flatfield']['tweak_slits'] = False
            if 'long2pos' in self.get_meta_value(headarr, 'decker'):
                # use dither info in the header to find objects
                par['reduce']['slitmask']['use_dither_offset'] = True
                # assume that the main target is always detected, i.e., skipping force extraction
                par['reduce']['slitmask']['extract_missing_objs'] = False


        # Return
        return par

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
        self.meta['target'] = dict(ext=0, card='TARGNAME')
        self.meta['decker'] = dict(ext=0, card='MASKNAME')
        self.meta['binning'] = dict(ext=0, card=None, default='1,1')

        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
        self.meta['exptime'] = dict(ext=0, card='TRUITIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='OBSMODE')
        self.meta['idname'] = dict(card=None, compound=True)
        self.meta['frameno'] = dict(ext=0, card='FRAMENUM')
        self.meta['object'] = dict(ext=0, card='OBJECT')
        # Filter
        self.meta['filter1'] = dict(ext=0, card='FILTER')
        # Lamps
        # flat lamps on/off
        self.meta['lampstat01'] = dict(ext=0, card='FLATSPEC')
        # lamp_names = ['FLATSPEC']
        # for kk,lamp_name in enumerate(lamp_names):
        #     self.meta['lampstat{:02d}'.format(kk+1)] = dict(ext=0, card=lamp_name)

        # Dithering
        self.meta['dithpat'] = dict(ext=0, card='PATTERN')
        self.meta['dithpos'] = dict(ext=0, card='FRAMEID')
        self.meta['dithoff'] = dict(ext=0, card='YOFFSET')

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
        if meta_key == 'idname':
            if headarr[0].get('KOAIMTYP', None) is not None:
                return headarr[0].get('KOAIMTYP')
            else:
                try:
                    # TODO: This should be changed to except on a specific error.
                    FLATSPEC = int(headarr[0].get('FLATSPEC'))
                    PWSTATA7 = int(headarr[0].get('PWSTATA7'))
                    PWSTATA8 = int(headarr[0].get('PWSTATA8'))
                    if FLATSPEC == 0 and PWSTATA7 == 0 and PWSTATA8 == 0:
                        if 'Flat:Off' in headarr[0].get('OBJECT') or "lamps off" in headarr[0].get('OBJECT'):
                            return 'flatlamp'
                        else:
                            return 'object'
                    elif FLATSPEC == 1:
                        return 'flatlamp'
                    elif PWSTATA7 == 1 or PWSTATA8 == 1:
                        return 'arclamp'
                except:
                    return 'unknown'
        else:
            msgs.error("Not ready for this compound meta")

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
        return ['decker', 'dispname', 'filter1']

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard ``PypeIt`` file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
#        pypeit_keys = super().pypeit_file_keys()
#        # TODO: Why are these added here? See
#        # pypeit.metadata.PypeItMetaData.set_pypeit_cols
#        pypeit_keys += [calib', 'comb_id', 'bkg_id']
#        return pypeit_keys
        return super().pypeit_file_keys() + [ 'lampstat01', 'dithpat', 'dithpos', 'dithoff', 'frameno']

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
        if ftype in ['science', 'standard']:
            return good_exp & (fitstbl['idname'] == 'object')
        if ftype in ['bias', 'dark']:
            return good_exp & self.lamps(fitstbl, 'off') & (fitstbl['idname'] == 'dark')
        if ftype in ['pixelflat', 'trace']:
            return good_exp & (fitstbl['idname'] == 'flatlamp')
        if ftype in ['illumflat']:
            # Flats and trace frames are typed together
            return good_exp & self.lamps(fitstbl, 'dome') & (fitstbl['idname'] == 'flatlamp')
        if ftype == 'pinhole':
            # Don't type pinhole frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            # TODO: This is a kludge.  Allow science frames to also be
            # classified as arcs
            is_arc = fitstbl['idname'] == 'arclamp'
            is_obj = self.lamps(fitstbl, 'off') & (fitstbl['idname'] == 'object')
            return good_exp & (is_arc | is_obj)
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    def lamps(self, fitstbl, status):
        """
        Check the lamp status.

        Args:
            fitstbl (`astropy.table.Table`_):
                The table with the fits header meta data.
            status (:obj:`str`):
                The status to check. Can be ``'off'``, ``'arcs'``, or
                ``'dome'``.

        Returns:
            `numpy.ndarray`_: A boolean array selecting fits files that meet
            the selected lamp status.

        Raises:
            ValueError:
                Raised if the status is not one of the valid options.
        """
        if status == 'off':
            # Check if all are off
            return fitstbl['lampstat01'] == '0'
            # return np.all(np.array([fitstbl[k] == 0 for k in fitstbl.keys() if 'lampstat' in k]),
            #               axis=0)
        # DP: this is not checking the lamp status of the arcs, but of the flats
        # if status == 'arcs':
        #     # Check if any arc lamps are on
        #     arc_lamp_stat = [ 'lampstat{0:02d}'.format(i) for i in range(1,6) ]
        #     return np.any(np.array([ fitstbl[k] == '1' for k in fitstbl.keys()
        #                                     if k in arc_lamp_stat]), axis=0)
        if status == 'dome':
            return fitstbl['lampstat01'] == '1'

        raise ValueError('No implementation for status = {0}'.format(status))



    def parse_dither_pattern(self, file_list, ext=None):
        """
        Parse headers from a file list to determine the dither pattern.

        Parameters
        ----------
        file_list (list of strings):
            List of files for which dither pattern is desired
        ext (int, optional):
            Extension containing the relevant header for these files. Default=None. If None, code uses
            self.primary_hdrext

        Returns
        -------
        dither_pattern, dither_id, offset_arcsec

        dither_pattern (str `numpy.ndarray`_):
            Array of dither pattern names
        dither_id (str `numpy.ndarray`_):
            Array of dither pattern IDs
        offset_arc (float `numpy.ndarray`_):
            Array of dither pattern offsets
        """
        nfiles = len(file_list)
        offset_arcsec = np.zeros(nfiles)
        dither_pattern = []
        dither_id = []
        for ifile, file in enumerate(file_list):
            hdr = fits.getheader(file, self.primary_hdrext if ext is None else ext)
            dither_pattern.append(hdr['PATTERN'])
            dither_id.append(hdr['FRAMEID'])
            offset_arcsec[ifile] = hdr['YOFFSET']
        return np.array(dither_pattern), np.array(dither_id), np.array(offset_arcsec)

    def tweak_standard(self, wave_in, counts_in, counts_ivar_in, gpm_in, meta_table, debug=False):
        """

        This routine is for performing instrument/disperser specific tweaks to standard stars so that sensitivity
        function fits will be well behaved. For example, masking second order light. For instruments that don't
        require such tweaks it will just return the inputs, but for isntruments that do this function is overloaded
        with a method that performs the tweaks.

        Parameters
        ----------
        wave_in: (float np.ndarray) shape = (nspec,)
            Input standard star wavelenghts
        counts_in: (float np.ndarray) shape = (nspec,)
            Input standard star counts
        counts_ivar_in: (float np.ndarray) shape = (nspec,)
            Input inverse variance of standard star counts
        gpm_in: (bool np.ndarray) shape = (nspec,)
            Input good pixel mask for standard
        meta_table: (astropy.table)
            Table containing meta data that is slupred from the specobjs object. See unpack_object routine in specobjs.py
            for the contents of this table.

        Returns
        -------
        wave_out: (float np.ndarray) shape = (nspec,)
            Output standard star wavelenghts
        counts_out: (float np.ndarray) shape = (nspec,)
            Output standard star counts
        counts_ivar_out: (float np.ndarray) shape = (nspec,)
            Output inverse variance of standard star counts
        gpm_out: (bool np.ndarray) shape = (nspec,)
            Output good pixel mask for standard

        """


        # Could check the wavelenghts here to do something more robust to header/meta data issues
        if 'Y-spectroscopy' in meta_table['DISPNAME']:
            #wave_out = np.copy(wave_in)
            #counts_out = np.copy(counts_in)
            #counts_ivar_out = np.copy(counts_ivar_in)
            #gpm_out = np.copy(gpm_in)
            # The blue edge and red edge of the detector are contaiminated by higher order light. These are masked
            # by hand.
            #second_order_region= (wave_in < 9520.0) | (wave_in > 11256.0)
            #gpm_out = gpm_in & np.logical_not(second_order_region)

            # Use a sigmoid function to apodize the spectrum smoothly in the regions where it is bad.
            #dlam = 10.0 # width that determines how shaprly  apodization occurs
            #sigmoid_blue_arg = (wave_in - wave_blue)/dlam
            #sigmoid_red_arg = (wave_red - wave_in)/dlam

            #sigmoid_apodize = special.expit(sigmoid_blue_arg) * special.expit(sigmoid_red_arg)
            #counts = counts_in*sigmoid_apodize
            # No we apodize only the flux. Since there are more counts in the in the unapodized spectrum, there is also
            # more variance in the original counts_ivar_in, and so in this way the S/N ratio is naturally reduced
            # in this region. There is not an obvious way to tweak the error vector here, and we don't to just mask
            # since then the polynomial fits go crazy at the boundaries. This is a reasonable compromise to not mask the
            # counts_ivar_in. The flux is bogus in the regions we are apodizing, so it does not really matter what we do,
            # it is better than operating on the original bogus flux.

            # Inflat the errors in the apodized region so that they don't inform the fits much
            #apo_pix = (sigmoid_apodize < 0.95)
            #mean_counts, med_counts, sigma_counts = sigma_clipped_stats(counts_in[np.logical_not(apo_pix)], sigma=3.0)
            #sigma_apo = med_counts/20.0 # roughly S/N ratio 20 in the apodized region
            #counts_ivar[apo_pix] = 1.0/sigma_apo**2
            #sigma_apo = np.sqrt(np.abs(counts[apo_pix]))
            #counts_ivar[apo_pix] = utils.inverse(sigma_apo**2)
            #counts_ivar[apo_pix] = utils.clip_ivar(counts[apo_pix], counts_ivar_in[apo_pix], 10.0, mask=gpm_in[apo_pix])
            wave_blue = 9520.0  # blue wavelength below which there is contamination
            wave_red = 11256.0  # red wavelength above which the spectrum is containated
            second_order_region= (wave_in < wave_blue) | (wave_in > wave_red)
            wave = wave_in.copy()
            counts = counts_in.copy()
            gpm = gpm_in.copy()
            counts_ivar = counts_ivar_in.copy()
            # By setting the wavelengths to zero, we guarantee that the sensitvity function will only be computed
            # over the valid wavelength region. While we could mask, this would still produce a wave_min and wave_max
            # for the zeropoint that includes the bad regions, and the polynomial fits will extrapolate crazily there
            wave[second_order_region] = 0.0
            counts[second_order_region] = 0.0
            counts_ivar[second_order_region] = 0.0
            gpm[second_order_region] = False
            #if debug:
            #    from matplotlib import pyplot as plt
            #    counts_sigma = np.sqrt(utils.inverse(counts_ivar_in))
            #    plt.plot(wave_in, counts, color='red', alpha=0.7, label='apodized flux')
            #    plt.plot(wave_in, counts_in, color='black', alpha=0.7, label='flux')
            #    plt.plot(wave_in, counts_sigma, color='blue', alpha=0.7, label='flux')
            #    plt.axvline(wave_blue, color='blue')
            #    plt.axvline(wave_red, color='red')
            #    plt.legend()
            #    plt.show()
            return wave, counts, counts_ivar, gpm
        else:
            return wave_in, counts_in, counts_ivar_in, gpm_in

    def get_slitmask(self, filename):
        """
        Parse the slitmask data from a MOSFIRE file into :attr:`slitmask`, a
        :class:`~pypeit.spectrographs.slitmask.SlitMask` object.

        This can be used for multi-object slitmask, but it it's not good
        for "LONGSLIT" nor "long2pos". Both "LONGSLIT" and "long2pos" have emtpy/incomplete
        binTable where the slitmask data are stored.


        Args:
            filename (:obj:`str`):
                Name of the file to read.

        Returns:
            :class:`~pypeit.spectrographs.slitmask.SlitMask`: The slitmask
            data read from the file. The returned object is the same as
            :attr:`slitmask`.
        """
        # Open the file
        hdu = io.fits_open(filename)

        # load slitmask info
        targs = hdu['Target_List'].data
        ssl = hdu['Science_Slit_List'].data
        msl = hdu['Mechanical_Slit_List'].data
        # some needed cleanup
        ssl = ssl[ssl['Slit_Number'] != ' ']
        msl = msl[msl['Slit_Number'] != ' ']

        # Book keeping: Count and check that the # of objects in the SSL matches that of the MSL
        # and if we recover the total number of CSUs
        CSUnumslits = 46
        numslits = np.zeros(len(ssl))
        for i in range(len(ssl)):
            slit = ssl[i]
            numslits[i] = np.where(slit['Target_Name'] == msl['Target_in_Slit'])[0].size

        if (numslits.sum() != CSUnumslits) and ('LONGSLIT' not in self.get_meta_value(filename, 'decker')) \
                and ('long2pos' not in self.get_meta_value(filename, 'decker')):
            msgs.error('The number of allocated CSU slits does not match the number of possible slits. '
                       'Slitmask design matching not possible. Turn parameter `use_maskdesign` off')

        targ_dist_center = np.array(ssl['Target_to_center_of_slit_distance'], dtype=float)

        slit_centers = np.array(ssl['Slit_length'], dtype=float) / 2.
        # Some # adjustment for long2pos
        if 'long2pos' in self.get_meta_value(filename, 'decker'):
            slit_centers_wrong = np.array(ssl['Slit_length'], dtype=float) / 2.
            # correct slit length
            ssl['Slit_length'] = '22.970', '7.010', '22.970'
            slit_centers = np.array(ssl['Slit_length'], dtype=float) / 2.
            centers_diff = slit_centers - slit_centers_wrong
            targ_dist_center[0] -= centers_diff[0]
            targ_dist_center[2] += centers_diff[2]

        # Projected distance (in arcsec) of the object from the left and right (top and bot) edges of the slit
        topdist = np.round(slit_centers + targ_dist_center, 3)
        botdist = np.round(slit_centers - targ_dist_center, 3)

        # Find the index to map the objects in the Science Slit List and the Target list
        indx = index_of_x_eq_y(targs['Target_Name'], ssl['Target_Name'])
        targs_mtch = targs[indx]
        obj_ra = targs_mtch['RA_Hours']+' '+targs_mtch['RA_Minutes']+' '+targs_mtch['RA_Seconds']
        obj_dec = targs_mtch['Dec_Degrees']+' '+targs_mtch['Dec_Minutes']+' '+targs_mtch['Dec_Seconds']
        obj_ra, obj_dec = meta.convert_radec(obj_ra, obj_dec)
        #   - Pull out the slit ID, object ID, name, object coordinates, top and bottom distance
        objects = np.array([np.array(ssl['Slit_Number'], dtype=int),
                           np.zeros(ssl['Slit_Number'].size, dtype=int),   # no object ID
                           obj_ra,
                           obj_dec,
                           np.array(ssl['Target_Name']),
                           np.array(targs_mtch['Magnitude'], dtype=float),
                           ['None']*ssl['Slit_Number'].size,       # no magnitude band
                           topdist,
                           botdist]).T

        # PA corresponding to positive x on detector (spatial)
        posx_pa = hdu[0].header['SKYPA0']
        if posx_pa < 0.:
            posx_pa += 360.

        slit_ra = ssl['Slit_RA_Hours']+' '+ssl['Slit_RA_Minutes']+' '+ssl['Slit_RA_Seconds']
        slit_dec = ssl['Slit_Dec_Degrees']+' '+ssl['Slit_Dec_Minutes']+' '+ssl['Slit_Dec_Seconds']
        slit_ra, slit_dec = meta.convert_radec(slit_ra, slit_dec)

        # Instantiate the slit mask object and return it
        self.slitmask = SlitMask(np.array([np.zeros(ssl['Slit_Number'].size),   # mosfire maskdef has not slit corners
                                           np.zeros(ssl['Slit_Number'].size),
                                           np.zeros(ssl['Slit_Number'].size),
                                           np.zeros(ssl['Slit_Number'].size),
                                           np.zeros(ssl['Slit_Number'].size),
                                           np.zeros(ssl['Slit_Number'].size),
                                           np.zeros(ssl['Slit_Number'].size),
                                           np.zeros(ssl['Slit_Number'].size)]).T.reshape(-1,4,2),
                                 slitid=np.array(ssl['Slit_Number'], dtype=int),
                                 align=ssl['Target_Name'] == 'posB',
                                 science=ssl['Target_Name'] != 'posB',
                                 onsky=np.array([slit_ra,
                                                 slit_dec,
                                                 np.array(ssl['Slit_length'], dtype=float),
                                                 np.array(ssl['Slit_width'], dtype=float),
                                                 np.array([round(hdu[0].header['SKYPA0'],2)]*ssl['Slit_Number'].size)]).T,
                                 objects=objects,
                                 posx_pa=posx_pa)
        return self.slitmask

    def get_maskdef_slitedges(self, ccdnum=None, filename=None, debug=None):
        """
         Provides the slit edges positions predicted by the slitmask design using
         the mask coordinates already converted from mm to pixels by the method
        `mask_to_pixel_coordinates`.

        If not already instantiated, the :attr:`slitmask`, :attr:`amap`,
        and :attr:`bmap` attributes are instantiated.  If so, a file must be provided.

        Args:
        ccdnum (:obj:`int`):
            Detector number
        filename (:obj:`str`, optional):
            The filename to use to (re)instantiate the :attr:`slitmask` and :attr:`grating`.
            Default is None, i.e., to use previously instantiated attributes.
        debug (:obj:`bool`, optional):
            Run in debug mode
        Returns:
            Three `numpy.ndarray`_ and a :class:`~pypeit.spectrographs.slitmask.SlitMask`.
            Two arrays are the predictions of the slit edges from the slitmask design and
            one contains the indices to order the slits from left to right in the PypeIt orientation

        """
        # Re-initiate slitmask and amap and bmap
        if filename is not None:
            # Reset the slitmask
            self.get_slitmask(filename)

        if self.slitmask is None:
            msgs.error('Unable to read slitmask design info. Provide a file.')

        # Hard-Coded for MOSFIRE
        platescale = 0.1798
        slit_gap = round(0.96/platescale)

        # build an array of values containing the bottom (right) edge of the slits
        # starting edge
        edge = 2034 if 'long2pos' not in self.get_meta_value(filename, 'decker') else 305
        bot_edges = np.array([edge], dtype=np.int)
        for i in range(self.slitmask.nslits - 1):
            # target is the slit number
            edge -= (self.slitmask.onsky[:,2][i]/platescale + slit_gap)
            bot_edges = np.append(bot_edges, np.round(edge))

        # build an array of values containing the top (left) edge of the slits
        top_edges = bot_edges - np.round(self.slitmask.onsky[:,2]/platescale)

        # Sort slits from left to right
        sortindx = np.argsort(self.slitmask.slitid[::-1])

        # This print a QA table with info on the slits sorted from left to right.
        if not debug:
            num = 0
            msgs.info('Expected slits')
            msgs.info('*' * 18)
            msgs.info('{0:^6s} {1:^12s}'.format('N.', 'Slit_Number'))
            msgs.info('{0:^6s} {1:^12s}'.format('-' * 5, '-' * 13))
            for i in range(sortindx.shape[0]):
                msgs.info('{0:^6d} {1:^12d}'.format(num, self.slitmask.slitid[sortindx][i]))
                num += 1
            msgs.info('*' * 18)

        # If instead we run this method in debug mode, we print more info
        if debug:
            num = 0
            msgs.info('Expected slits')
            msgs.info('*' * 92)
            msgs.info('{0:^5s} {1:^10s} {2:^12s} {3:^12s} {4:^16s} {5:^16s}'.format('N.', 'Slit_Number',
                                                                                    'slitLen(arcsec)',
                                                                                    'slitWid(arcsec)',
                                                                                    'top_edges(pix)',
                                                                                    'bot_edges(pix)'))
            msgs.info('{0:^5s} {1:^10s} {2:^12s} {3:^12s} {4:^16s} {5:^14s}'.format('-' * 4, '-' * 13, '-' * 11,
                                                                                    '-' * 11, '-' * 18, '-' * 15))
            for i in range(sortindx.size):
                msgs.info('{0:^5d}{1:^14d} {2:^9.3f} {3:^12.3f}    {4:^16.2f} {5:^14.2f}'.format(num,
                            self.slitmask.slitid[sortindx][i], self.slitmask.onsky[:,2][sortindx][i],
                            self.slitmask.onsky[:,3][sortindx][i], top_edges[sortindx][i], bot_edges[sortindx][i]))
                num += 1
            msgs.info('*' * 92)

        return top_edges, bot_edges, sortindx, self.slitmask


def find_longslit_pos(raw_file):
    """
    Given a MOSFIRE science raw file, find the position of the slit
    in the LONGSLIT slitmask

    Args:
        raw_file: (:obj:`str`):
            Name of the file to read.

    Returns:
        :obj:`tuple`: Two integer number indicating the x position of the
        beginning and the end of the slit.

    """

    hdu = io.fits_open(raw_file)
    decker = hdu[0].header['MASKNAME']
    pixelscale = hdu[0].header['PSCALE']
    CSUnum=int(decker.split("x")[0].split('-')[1])
    # Hard-coded for MOSFIRE
    slit_gap = 0.96
    slit_length = CSUnum * 7.01/pixelscale + (CSUnum-1)*slit_gap/pixelscale
    pix_start = hdu[0].header['CRPIX2'] - slit_length/2
    pix_end = hdu[0].header['CRPIX2'] + slit_length/2

    return int(pix_start), int(pix_end)
