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
p.add_argument('-w', '--write', action='store_true', help="Save the mask in a .fits.gz file and save the plot as a .png")
p.add_argument('-m', '--mode', type=int, default=2, help="Mode 1 uses ugali.candidate.associate.py and masks a radius of 6' near known galaxies. Mode 2 uses my associate.py and masks dynamically based on object size.")
p.add_argument('-s', '--survey', default=None, help="'des' or 'ps1'")
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
healpix_mask[cut_ebv] |= 0b00001


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
ugali data on: Harris96, WEBDA14
"""
#external_cat_list = ['McConnachie15', 'Harris96', 'Corwen04', 'Nilson73', 'Webbink85', 'Kharchenko13', 'Bica08', 'WEBDA14', 'ExtraDwarfs','ExtraClusters']
external_cat_list = ['Harris96', 'Corwen04', 'Nilson73', 'Webbink85', 'Kharchenko13', 'Bica08', 'WEBDA14', 'ExtraClusters']
for external_cat in external_cat_list:
    if args.mode == 1:
        catalog = ugali.candidate.associate.catalogFactory(external_cat)
        external_cut = cut_circles(catalog['ra'], catalog['dec'], default_radius=0.1)
    elif args.mode == 2:
        catalog = associate.catalogFactory(external_cat)
        external_cut = cut_circles(catalog['ra'], catalog['dec'], catalog['radius'], default_radius=0.1)
    healpix_mask[external_cut] |= 0b00010

known_dwarfs = ['McConnachie15', 'ExtraDwarfs']
for external_cat in known_dwarfs:
    if args.mode == 1:
        catalog = ugali.candidate.associate.catalogFactory(external_cat)
        external_cut = cut_circles(catalog['ra'], catalog['dec'], default_radius=0.1)
    elif args.mode == 2:
        catalog = associate.catalogFactory(external_cat)
        external_cut = cut_circles(catalog['ra'], catalog['dec'], catalog['radius'], default_radius=0.1)
    healpix_mask[external_cut] |= 0b00100


reader_bsc = pyfits.open('/Users/mcnanna/Research/DES_luminosity/bsc5.fits')
d_bsc = reader_bsc[1].data

bsc_cut = cut_circles(d_bsc['RA'], d_bsc['DEC'], default_radius=0.1)
healpix_mask[bsc_cut] |= 0b01000

if args.survey is not None:
    if 'des' in args.survey:
        des_footprint = ugali.utils.healpix.read_map('/Users/mcnanna/Research/DES_luminosity/y3a2_footprint_griz_1exp_v2.0.fits.gz', nest=True) 
        des_cut = des_footprint < 1
        healpix_mask[des_cut] |= 0b10000

    if 'ps1' in args.survey:
        ps1_cut = healpy.query_strip(NSIDE, np.radians(90.0 - -25.0), np.radians(90.0 - -90.0), nest=False) # This function apparently isn't implemented for nest=True
        ps1_cut = healpy.ring2nest(NSIDE, ps1_cut)
        healpix_mask[ps1_cut] |= 0b10000


pylab.figure()
healpy.mollview(healpix_mask, nest=True, coord='C', cmap='binary')
if args.write:
    pylab.savefig('healpix_mask_{}.png'.format(args.survey))

if args.write:
    print('Writing map...')
    healpy.write_map('healpix_mask_{}.fits.gz'.format(args.survey), healpix_mask, dtype=np.int32, nest=NEST, coord='C', overwrite=True)



