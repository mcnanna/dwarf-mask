#!/usr/bin/env python

import numpy as np
import healpy
import astropy.io.fits as pyfits
import pylab
import argparse

import ugali.utils.projector
import ugali.utils.healpix
import ugali.candidate.associate
import associate

p = argparse.ArgumentParser()
p.add_argument('-l', '--load', default=None)
p.add_argument('-w', '--write', action='store_true')
p.add_argument('-m', '--mode', type=int, default=2, help="Mode 1 uses ugali.candidate.associate.py and masks a radius of 6' near known galaxies. Mode 2 uses my associate.py and masks dynamically based on object size.")
args = p.parse_args()

NEST=True
NSIDE=4096


# Initiate_mask
if args.load is not None:
    healpix_mask = ugali.utils.healpix.read_map(args.load, nest=NEST)
else:
    healpix_mask = np.tile(0, healpy.nside2npix(NSIDE))


infile_dust = '/Users/mcnanna/Research/DES_luminosity/ebv_sfd98_fullres_nside_4096_nest_equatorial.fits.gz'
ebv_map = ugali.utils.healpix.read_map(infile_dust, nest=True)

cut_ebv = (ebv_map > 0.2)
healpix_mask[cut_ebv] |= 0b0001


def cut_circles(ras, decs, radii=None, default_radius=0.1, factor=1.0):
    """Blocks pixels near given (ras, decs). 
    If radius (in degrees) are given, it blocks all pixels within disc of factor*radius. 
    If radius is not given or is nan, it defaults to default_radius (in degrees)
    """
    cut = np.tile(False, healpy.nside2npix(NSIDE))

    if radii is None:
        radii = np.tile(default_radius/factor, len(ras))

    for ra,dec,radius in zip(ras, decs, radii):
        if np.isnan(radius) or radius == 0:
            radius = default_radius/factor
        bad_pixels = ugali.utils.healpix.angToDisc(NSIDE, ra, dec, factor*radius, nest=NEST, inclusive=True)
        cut[bad_pixels] = True
    
    return cut


"""
McConnachie12: Half-light radius along major axis. Could incorporate ellipticity info
Harris96: Half-light radius
Corwen04: No radius data
Nilson73: Major axis (no ellipticity info)
Webbink85: 10**log(radius)
Kharchenko13: Radius of the central part
Bica08: Major axis. Could incorporate ellipticity info
WEBDA: Diameter/2

Vizier data on: McConnachie12, Nilson73, Webbink85 (should check log or ln), Kharchenko13 (3 radii, used middle), Bica08
ugali data on: Harris 96, WEBDA14
"""
external_cat_list = ['McConnachie12', 'Harris96', 'Corwen04', 'Nilson73', 'Webbink85', 'Kharchenko13', 'Bica08', 'WEBDA14', 'ExtraDwarfs','ExtraClusters']
for external_cat in external_cat_list:
    if args.mode == 1:
        catalog = ugali.candidate.associate.catalogFactory(external_cat)
        external_cut = cut_circles(catalog['ra'], catalog['dec'], default_radius=0.1)
    elif args.mode == 2:
        catalog = associate.catalogFactory(external_cat)
        external_cut = cut_circles(catalog['ra'], catalog['dec'], catalog['radius'], default_radius=0.1)
    healpix_mask[external_cut] |= 0b0010


reader_bsc = pyfits.open('/Users/mcnanna/Research/DES_luminosity/bsc5.fits')
d_bsc = reader_bsc[1].data

bsc_cut = cut_circles(d_bsc['RA'], d_bsc['DEC'], default_radius=0.1)
healpix_mask[bsc_cut] |= 0b0100


if args.write:
    print('Writing map...')
    healpy.write_map('healpix_mask_{}.fits.gz'.format(args.mode), healpix_mask, dtype=np.int32, nest=NEST, coord='C', overwrite=True)

pylab.figure()
healpy.mollview(healpix_mask, nest=True, coord='C', cmap='binary')
if args.write:
    pylab.savefig('healpix_mask_{}.png'.format(args.mode))


