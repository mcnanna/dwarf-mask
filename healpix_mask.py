#!/usr/bin/env python

import subprocess
import numpy as np
import healpy as hp
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
plt.ion()
import argparse

import ugali.utils.projector
import ugali.utils.healpix
import ugali.candidate.associate
import associate
from skymap import Skymap

NEST=True
NSIDE=4096

# Custom color map for plotting, based on 'ocean'
n = 4
base = plt.cm.get_cmap('ocean_r')
cmap = base(np.linspace(0,1,n))
cl = base(np.linspace(0,1,n))
color_list = [cl[0], cl[3], cl[1], cl[2]] # Swap green before the blues
color_list[3] = np.array([0., 0, 0.25, 1.]) # Make the dark blue darker
new_cmap = base.from_list('new_cmap', color_list, n)

def cut_circles(ras, decs, radii=None, default_radius=0.1, min_radius=0.05):
    """Blocks pixels near given (ras, decs). 
    If radius (in degrees) are given, it blocks all pixels within disc of radius, with a minimum radius of min_radius.  
    If radius is not given or is nan, it defaults to default_radius (in degrees)
    """
    cut = np.tile(False, hp.nside2npix(NSIDE))

    if radii is None:
        radii = np.tile(default_radius, len(ras))

    for ra,dec,radius in zip(ras, decs, radii):
        if np.isnan(radius) or radius == 0:
            radius = default_radius
        elif radius < min_radius:
            radius = min_radius
        bad_pixels = ugali.utils.healpix.angToDisc(NSIDE, ra, dec, radius, nest=NEST, inclusive=True)
        cut[bad_pixels] = True
    
    return cut


def main(survey, plot=True, write=True):
    print("Creating mask with NSIDE={} and NEST={}...".format(NSIDE, NEST))

    # Initiate_mask
    if args.load is not None:
        healpix_mask = ugali.utils.healpix.read_map(args.load, nest=NEST)
    else:
        healpix_mask = np.tile(0, hp.nside2npix(NSIDE))


    infile_dust = 'ebv_sfd98_fullres_nside_4096_nest_equatorial.fits.gz'
    ebv_map = ugali.utils.healpix.read_map(infile_dust, nest=True)

    cut_ebv = (ebv_map > 0.2)
    healpix_mask[cut_ebv] |= 0b00001

    """
    McConnachie15: Half-light radius along major axis. Could incorporate ellipticity info
    Harris96: Half-light radius
    Corwen04: No radius data
    Nilson73: Major axis (no ellipticity info)
    Webbink85: 10**log(radius)
    Kharchenko13: Radius of the central part
    Bica08: Major axis. Could incorporate ellipticity info
    WEBDA: Diameter/2
    ExtraDwarfs: Half-light radius
    ExtraClusters: Half-light radius

    Vizier data on: Nilson73, Webbink85 (should check log or ln), Kharchenko13 (3 radii, used middle), Bica08
    ugali data on: McConnachi15, Harris96, WEBDA14
    """
    external_cat_list = ['Harris96', 'Corwen04', 'Nilson73', 'Webbink85', 'Kharchenko13', 'Bica08', 'WEBDA14', 'ExtraClusters', 'ExtraStructures']
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


    reader_bsc = pyfits.open('bsc5.fits')
    d_bsc = reader_bsc[1].data

    bsc_cut = cut_circles(d_bsc['RA'], d_bsc['DEC'], default_radius=0.1)
    healpix_mask[bsc_cut] |= 0b01000

    if 'des' == survey:
        des_footprint = ugali.utils.healpix.read_map('y3a2_footprint_griz_1exp_v2.0.fits.gz', nest=True) 
        #print('Masking footprint...')
        #npix = hp.nside2npix(NSIDE)
        #nbad = 3
        #des_cut = np.fromiter( (des_footprint[i] < 1 or sum(des_footprint[hp.get_all_neighbours(NSIDE, i)] < 1) >= nbad for i in range(npix)), dtype=bool, count=npix ) # Mask if at least nbad neighbors outside the footprint
        des_cut = des_footprint < 1
        healpix_mask[des_cut] |= 0b10000

    if 'ps1' == survey:
        ps1_cut = hp.query_strip(NSIDE, np.radians(90.0 - -25.0), np.radians(90.0 - -90.0), nest=False) # This function apparently isn't implemented for nest=True
        ps1_cut = hp.ring2nest(NSIDE, ps1_cut)
        healpix_mask[ps1_cut] |= 0b10000

        failures = pyfits.open('ugali_failures.fits')[1].data # NSIDE = 256
        fail_pix_256 = [ugali.utils.healpix.angToPix(256, fail['ra'], fail['dec'], nest=NEST) for fail in failures]
        fail_pix = ugali.utils.healpix.ud_grade_ipix(fail_pix_256, 256, NSIDE, nest=NEST)
        healpix_mask[fail_pix] |= 0b100000

    if plot:
        print("Simplifying mask for plotting...")
        simplified_mask = np.copy(healpix_mask)
        cut_catalog = np.where((simplified_mask & 0b00010) | (simplified_mask & 0b00100) | (simplified_mask & 0b01000))
        cut_ebv = np.where(simplified_mask & 0b00001)
        cut_footprint = np.where((simplified_mask & 0b10000) | (simplified_mask & 0b100000))
        simplified_mask[cut_catalog] = 1
        simplified_mask[cut_ebv] = 2
        simplified_mask[cut_footprint] = 3

        print("Plotting...")
        title = ''
        if survey == 'des':
            title = 'DES ' + title
        elif survey == 'ps1':
            title = 'Pan-STARRS ' + title

        """
        # Using mollview
        hp.mollview(simplified_mask, nest=True, coord='C', cmap=new_cmap, title=title, xsize=1600)
        ax = plt.gca()
        cbar = ax.images[-1].colorbar
        cbar.set_ticks( (np.arange(n) + 0.5)*(n-1)/n )
        cbar.set_ticklabels(['Unmasked', 'Association', r'$E(B-V)$', 'Footprint'])
        hp.graticule()
        plt.savefig('healpix_mask_{}_v1.png'.format(survey), bbox_inches='tight')
        """

        # Using skymap
        simplified_mask_ring = hp.reorder(simplified_mask, n2r=True)
        fig, ax = plt.subplots(figsize=(12,8))                     
        smap = Skymap(projection = 'mbtfpq',lon_0=0)
        im,lon,lat,values = smap.draw_hpxmap(simplified_mask_ring, xsize=1600, cmap=new_cmap)
        cbar = plt.colorbar(ticks = (np.arange(n) + 0.5)*(n-1)/n, fraction=0.02)
        cbar.set_ticklabels(['Unmasked', 'Association', r'$E(B-V) > 0.2$', 'Footprint'])
        plt.title(title)
        plt.savefig('healpix_mask_{}.png'.format(survey), bbox_inches='tight')

    if write:
        print('Writing mask to healpix_mask_{}.fits.gz ...'.format(survey))
        hp.write_map('healpix_mask_{}.fits.gz'.format(survey), healpix_mask, dtype=np.int32, nest=NEST, coord='C', overwrite=True)

    return healpix_mask, simplified_mask_ring


p = argparse.ArgumentParser()
p.add_argument('-s', '--survey', help="None(=both des and ps1), des, ps1, or gen, short for general, which doesn't apply a footprint mask")
p.add_argument('--no_plot', action='store_true')
p.add_argument('--no_write', action='store_true')
p.add_argument('-l', '--load', default=None, help="Load a map to write on top of. Default: None")
p.add_argument('-m', '--mode', type=int, default=2, help="Mode 1 uses ugali.candidate.associate.py and masks a radius of 6' near known galaxies. Mode 2 uses my associate.py and masks dynamically based on object size. Default: 2")
args = p.parse_args()

if args.survey is None:

    des_mask, des_mask_simplified = main('des', plot=(not args.no_plot), write=(not args.no_write))
    ps1_mask, ps1_mask_simplified = main('ps1', plot=(not args.no_plot), write=(not args.no_write))

    # Put the two healpix masks on the same plot
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(12, 10))
    plt.sca(axes[0])
    smap = Skymap(projection='mbtfpq', lon_0=0)
    im,lon,lat,values = smap.draw_hpxmap(des_mask_simplified, xsize=1600, cmap=new_cmap)
    plt.title('DES', fontsize='large')
    plt.sca(axes[1])
    smap = Skymap(projection='mbtfpq', lon_0=0)
    im,lon,lat,values = smap.draw_hpxmap(ps1_mask_simplified, xsize=1600, cmap=new_cmap)
    plt.title('Pan-STARRS', fontsize='large')
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.175, 0.020, 0.65])
    cbar = fig.colorbar(im, cax=cbar_ax, ticks = (np.arange(n) + 0.5)*(n-1)/n )
    cbar.set_ticklabels(['Unmasked', 'Association', r'$E(B-V) > 0.2$', 'Footprint'])
    plt.savefig('healpix_masks.png', bbox_inches='tight')

else:
    main(args.survey, plot=(not args.no_plot), write=(not args.no_write))
