#!/usr/bin/env python

import numpy as np
import healpy
import astropy.io.fits as pyfits
import pylab
import argparse

import ugali.utils.projector
import ugali.utils.healpix
import ugali.candidate.associate

p = argparse.ArgumentParser()
p.add_argument('-l', '--load', default=None)
args = p.parse_args()
load = args.load

NEST=True
NSIDE=4096


# Initiate_mask
if load is not None:
    healpix_mask = ugali.utils.healpix.read_map(load, nest=NEST)
else:
    healpix_mask = np.tile(0, healpy.nside2npix(NSIDE))


infile_dust = '/Users/mcnanna/Research/DES_luminosity/ebv_sfd98_fullres_nside_4096_nest_equatorial.fits.gz'
ebv_map = ugali.utils.healpix.read_map(infile_dust, nest=True)

cut_ebv = (ebv_map > 0.2)
healpix_mask[cut_ebv] |= 0b0001


def cut_from_hotspots(ras, decs, radius=0.1):
    # Blocks pixels within 'radius' degrees of (ra, dec)
    cut = np.tile(False, healpy.nside2npix(NSIDE))

    for ra,dec in zip(ras, decs):
        bad_pixels = ugali.utils.healpix.angToDisc(NSIDE, ra, dec, radius, nest=NEST)
        cut[bad_pixels] = True
    
    return cut

#external_cat_list = ['Nilson73'] #'Harris96', 'ExtraClusters', 'WEBDA14']
external_cat_list = ['McConnachie15', 'Harris96', 'Corwen04', 'Nilson73', 'Webbink85', 'Kharchenko13', 'Bica08', 'WEBDA14', 'ExtraDwarfs','ExtraClusters']
for external_cat in external_cat_list:
    catalog = ugali.candidate.associate.catalogFactory(external_cat)
    external_cut = cut_from_hotspots(catalog['ra'], catalog['dec'], radius=0.1)
    healpix_mask[external_cut] |= 0b0010


reader_bsc = pyfits.open('/Users/mcnanna/Research/DES_luminosity/bsc5.fits')
d_bsc = reader_bsc[1].data
bsc_cut = cut_from_hotspots(d_bsc['RA'], d_bsc['DEC'], radius=0.1)
healpix_mask[bsc_cut] |= 0b0100

healpy.write_map('healpix_mask.fits', healpix_mask, dtype=np.int32, nest=NEST, coord='C', overwrite=True)

pylab.figure()
healpy.mollview(healpix_mask, nest=True, coord='C', cmap='binary')
pylab.savefig('healpix_mask.png')
