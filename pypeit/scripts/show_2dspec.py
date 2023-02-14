"""
This script enables the viewing of a processed FITS file
with extras.  Run above the Science/ folder.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import os

import numpy as np

from IPython import embed

from astropy.io import fits
from astropy.stats import sigma_clipped_stats

from pypeit import msgs
from pypeit import slittrace
from pypeit import specobjs
from pypeit import io

from pypeit.display import display
from pypeit.images.imagebitmask import ImageBitMask
from pypeit.images.detector_container import DetectorContainer
from pypeit import masterframe
from pypeit import spec2dobj
from pypeit.scripts import scriptbase


# TODO: We should not be calling objects by the same names as modules we've
# loaded.  I.e., we've loaded specobjs, and here we're passing specobjs to this
# function...
def show_trace(specobjs, det, viewer, ch):
    """
    Overplot the extracted object traces for this detector in the provided ginga
    channel.

    Args:
        specobjs (:class:`~pypeit.specobjs.SpecObjs`):
            Object holding the 1D spectral extractions.  If None, the function
            doesn't do anything.
        det (:obj:`str`):
            The string identifier for the detector or mosaic used to select the
            extractions to show.
        viewer (?):
        ch (?):
    """
    if specobjs is None:
        return
    in_det = np.where(specobjs.DET == det)[0]
    for kk in in_det:
        trace = specobjs[kk]['TRACE_SPAT']
        obj_id = specobjs[kk].NAME
        maskdef_objname = specobjs[kk].MASKDEF_OBJNAME
        maskdef_extr_flag = specobjs[kk].MASKDEF_EXTRACT
        manual_extr_flag = specobjs[kk].hand_extract_flag
        if maskdef_objname is not None:
            trc_name = '{}     OBJNAME:{}'.format(obj_id, maskdef_objname)
        else:
            trc_name = obj_id
        if maskdef_extr_flag is not None and maskdef_extr_flag is True:
            display.show_trace(viewer, ch, trace, trc_name, color='#f0e442') #hdu.name)
        elif manual_extr_flag is True:
            display.show_trace(viewer, ch, trace, trc_name, color='#33ccff') #hdu.name)
        else:
            display.show_trace(viewer, ch, trace, trc_name, color='orange') #hdu.name)


class Show2DSpec(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Display sky subtracted, spec2d image in a '
                                                 'Ginga viewer.  Run above the Science/ folder',
                                    width=width)

        parser.add_argument('file', type = str, default = None, help = 'PypeIt spec2d file')
        parser.add_argument('--list', default=False, help='List the extensions only?',
                            action='store_true')
        # User can provide only '--det' or '--detname', not both
        parser.add_argument('--det', default='1', type=str,
                            help='Detector name or number.  If a number, the name is constructed '
                                 'assuming the reduction is for a single detector.  If a string, '
                                 'it must match the name of the detector object (e.g., DET01 for '
                                 'a detector, MSC01 for a mosaic).')
        parser.add_argument("--spat_id", type=int, default=None,
                            help='Restrict plotting to this slit (PypeIt ID notation)')
        parser.add_argument("--maskID", type=int, default=None,
                            help='Restrict plotting to this maskID')
        parser.add_argument('--showmask', default=False, help='Overplot masked pixels',
                            action='store_true')
        parser.add_argument('--removetrace', default=False, action="store_true",
                            help='Do not overplot traces in the skysub, sky_resid, and resid '
                                 'channels')
        parser.add_argument('--embed', default=False, action='store_true',
                            help='Upon completion embed in ipython shell')
        parser.add_argument('--ignore_extract_mask', default=False, action='store_true',
                            help='Ignore the extraction mask')
        parser.add_argument("--sensfunc", type=str, default=None,
                            help='Pass in a sensfunc to display the sky-subtracted image with a '
                                 'flux calibration')
        parser.add_argument('--channels', type=str,
                            help='Only show a subset of the channels (0-indexed), e.g. 1,3')
        parser.add_argument('--prefix', type=str, default='',
                            help="Channel name prefix [lets you display more than one set]")
        parser.add_argument('--no_clear', dest='clear', default=True, 
                            action="store_false",
                            help='Do *not* clear all existing tabs')
        parser.add_argument('-v', '--verbosity', type=int, default=2,
                            help='Verbosity level between 0 [none] and 2 [all]')
        return parser

    @staticmethod
    def main(args):

        # List only?
        if args.list:
            io.fits_open(args.file).info()
            return

        # Parse the detector name
        try:
            det = int(args.det)
        except:
            detname = args.det
        else:
            detname = DetectorContainer.get_name(det)

        # Load it up -- NOTE WE ALLOW *OLD* VERSIONS TO GO FORTH
        spec2DObj = spec2dobj.Spec2DObj.from_file(args.file, detname, chk_version=False)
        # Use the appropriate class to get the "detector" number
        det = spec2DObj.detector.parse_name(detname)

        # Setup for PypeIt imports
        msgs.reset(verbosity=args.verbosity)

        # Find the set of channels to show
        if args.channels is not None:
            show_channels = [int(item) for item in args.channels.split(',')]
        else:
            show_channels = [0,1,2,3]

        # Grab the slit edges
        slits = spec2DObj.slits
        if spec2DObj.sci_spat_flexure is not None:
            msgs.info("Offseting slits by {}".format(spec2DObj.sci_spat_flexure))
        all_left, all_right, mask = slits.select_edges(flexure=spec2DObj.sci_spat_flexure)

        # TODO -- This may be too restrictive, i.e. ignore BADFLTCALIB??
        gpm = mask == 0

        # Restrict on spat_id or maskid? 
        if args.spat_id is not None:
            gpm &= slits.spat_id == args.spat_id
        elif args.maskID is not None:
            gpm &= slits.maskdef_id == args.maskID
        
        left = all_left[:, gpm]
        right = all_right[:, gpm]
        slid_IDs = spec2DObj.slits.slitord_id[gpm]
        maskdef_id = None if spec2DObj.slits.maskdef_id is None \
                            else spec2DObj.slits.maskdef_id[gpm]

        # Object traces from spec1d file
        spec1d_file = args.file.replace('spec2d', 'spec1d')
        if args.file[-2:] == 'gz':
            spec1d_file = spec1d_file[:-3]
        if os.path.isfile(spec1d_file):
            sobjs = specobjs.SpecObjs.from_fitsfile(spec1d_file, chk_version=False)
        else:
            sobjs = None
            msgs.warn('Could not find spec1d file: {:s}'.format(spec1d_file) + msgs.newline() +
                      '                          No objects were extracted.')

        display.connect_to_ginga(raise_err=True, allow_new=True)

        # Now show each image to a separate channel

        # Show the bitmask?
        if args.showmask:
            viewer, ch_mask = display.show_image(spec2DObj.bpmmask, chname='MASK',
                                                 waveimg=spec2DObj.waveimg, 
                                                 clear=args.clear)

        channel_names = []
        # SCIIMG
        if 0 in show_channels:
            image = spec2DObj.sciimg  # Processed science image
            gpm = spec2DObj.select_flag(invert=True)
            mean, med, sigma = sigma_clipped_stats(image[gpm], sigma_lower=5.0, sigma_upper=5.0)
            cut_min = mean - 1.0 * sigma
            cut_max = mean + 4.0 * sigma
            chname_sci = args.prefix+f'sciimg-{detname}'
            # Clear all channels at the beginning
            viewer, ch_sci = display.show_image(image, chname=chname_sci,
                                                waveimg=spec2DObj.waveimg, 
                                                clear=args.clear,
                                                cuts=(cut_min, cut_max))

            if sobjs is not None:
                show_trace(sobjs, detname, viewer, ch_sci)
            display.show_slits(viewer, ch_sci, left, right, slit_ids=slid_IDs,
                               maskdef_ids=maskdef_id)
            channel_names.append(chname_sci)

        # SKYSUB
        if 1 in show_channels:
            gpm = spec2DObj.select_flag(invert=True)
            if args.ignore_extract_mask:
                gpm |= spec2DObj.select_flag(flag='EXTRACT')

            image = (spec2DObj.sciimg - spec2DObj.skymodel) * gpm.astype(float)
            mean, med, sigma = sigma_clipped_stats(image[gpm], sigma_lower=5.0, sigma_upper=5.0)
            cut_min = mean - 1.0 * sigma
            cut_max = mean + 4.0 * sigma
            chname_skysub = args.prefix+f'skysub-{detname}'
            viewer, ch_skysub = display.show_image(image, chname=chname_skysub,
                                                   waveimg=spec2DObj.waveimg,
                                                   cuts=(cut_min, cut_max), wcs_match=True)
            if not args.removetrace and sobjs is not None:
                    show_trace(sobjs, detname, viewer, ch_skysub)
            display.show_slits(viewer, ch_skysub, left, right, slit_ids=slid_IDs,
                               maskdef_ids=maskdef_id)
            channel_names.append(chname_skysub)

        # TODO Place holder for putting in sensfunc
        #if args.sensfunc:
        #    # Load the sensitivity function
        #    wave_sens, sfunc, _, _, _ = sensfunc.SensFunc.load(sensfunc_masterframe_name)
        #    # Interpolate the sensitivity function onto the wavelength grid of the data. Since the image is rectified
        #    # this is trivial and we don't need to do a 2d interpolation
        #    sens_factor = flux_calib.get_sensfunc_factor(
        #        pseudo_dict['wave_mid'][:,islit], wave_sens, sfunc, fits.getheader(files[0])['TRUITIME'],
        #        extrap_sens=parset['fluxcalib']['extrap_sens'])
        #    # Compute the median sensitivity and set the sensitivity to zero at locations 100 times the median. This
        #    # prevents the 2d image from blowing up where the sens_factor explodes because there is no throughput
        #    sens_gpm = sens_factor < 100.0*np.median(sens_factor)
        #    sens_factor_masked = sens_factor*sens_gpm
        #    sens_factor_img = np.repeat(sens_factor_masked[:, np.newaxis], pseudo_dict['nspat'], axis=1)
        #    imgminsky = sens_factor_img*pseudo_dict['imgminsky']
        #    imgminsky_gpm = sens_gpm[:, np.newaxis] & pseudo_dict['inmask']
        #else:
        #    imgminsky= pseudo_dict['imgminsky']

        # SKRESIDS
        if 2 in show_channels:
            # the block below is repeated because if showing this channel but
            # not channel 1 it will crash
            gpm = spec2DObj.select_flag(invert=True)
            if args.ignore_extract_mask:
                gpm |= spec2DObj.select_flag(flag='EXTRACT')
            chname_skyresids = args.prefix+f'sky_resid-{detname}'
            image = (spec2DObj.sciimg - spec2DObj.skymodel) * np.sqrt(spec2DObj.ivarmodel) \
                        * gpm.astype(float)
            viewer, ch_sky_resids = display.show_image(image, chname_skyresids,
                                                       waveimg=spec2DObj.waveimg, cuts=(-5.0, 5.0))
            if not args.removetrace and sobjs is not None:
                    show_trace(sobjs, detname, viewer, ch_sky_resids)
            display.show_slits(viewer, ch_sky_resids, left, right, slit_ids=slid_IDs,
                               maskdef_ids=maskdef_id)
            channel_names.append(chname_skyresids)

        # RESIDS
        if 3 in show_channels:
            chname_resids = args.prefix+f'resid-{detname}'
            # full model residual map
            gpm = spec2DObj.select_flag(invert=True)
            image = (spec2DObj.sciimg - spec2DObj.skymodel - spec2DObj.objmodel) \
                        * np.sqrt(spec2DObj.ivarmodel) * gpm.astype(float)
            viewer, ch_resids = display.show_image(image, chname=chname_resids,
                                                   waveimg=spec2DObj.waveimg, cuts=(-5.0, 5.0),
                                                   wcs_match=True)
            if not args.removetrace and sobjs is not None:
                    show_trace(sobjs, detname, viewer, ch_resids)
            display.show_slits(viewer, ch_resids, left, right, slit_ids=slid_IDs,
                               maskdef_ids=maskdef_id)
            channel_names.append(chname_resids)


        # After displaying all the images sync up the images with WCS_MATCH
        shell = viewer.shell()
        shell.start_global_plugin('WCSMatch')
        shell.call_global_plugin_method('WCSMatch', 'set_reference_channel', [channel_names[-1]],
                                        {})

        if args.embed:
            embed()

        # Playing with some mask stuff
        #out = shell.start_operation('TVMask')
        #maskfile = '/Users/joe/python/PypeIt-development-suite/REDUX_OUT/Shane_Kast_blue/600_4310_d55/shane_kast_blue_setup_A/crmask.fits'
        #out = shell.call_local_plugin_method(chname_resids, 'TVMask', 'load_file', [maskfile], {})


