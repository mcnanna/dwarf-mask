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


# I'm gonna start using a bitmask now

def cut_from_hotspots(ras, decs, radius=0.1):
    # Blocks pixels within 'radius' degrees of (ra, dec)
    cut = np.tile(True, healpy.nside2npix(NSIDE))

    for ra,dec in zip(ras, decs):
        bad_pixels = ugali.utils.healpix.angToDisc(NSIDE, ra, dec, radius, nest=NEST)
        cut[bad_pixels] = False
    
    return cut



# Initiate_mask
if load is not None:
    healpix_mask = ugali.utils.healpix.read_map(load, nest=NEST)
else:
    healpix_mask = np.tile(True, healpy.nside2npix(NSIDE))


# Cuts to apply:
    apply_ebv = True
    apply_Nilson = False
    apply_external_cats = True
    apply_bsc = True

if apply_ebv:
    infile_dust = '/Users/mcnanna/Research/DES_luminosity/ebv_sfd98_fullres_nside_4096_nest_equatorial.fits.gz'
    ebv_map = ugali.utils.healpix.read_map(infile_dust, nest=True)

    cut_ebv = (ebv_map < 0.2)
    healpix_mask &= cut_ebv    

if apply_Nilson:
    external_cat_list = ['Nilson73'] #'Harris96', 'ExtraClusters', 'WEBDA14']
    for external_cat in external_cat_list:
        catalog = ugali.candidate.associate.catalogFactory(external_cat)
        external_cut = cut_from_hotspots(catalog['ra'], catalog['dec'], radius=0.1)
        healpix_mask &= external_cut

if apply_external_cats:
    external_cat_list = ['McConnachie15', 'Harris96', 'Corwen04', 'Nilson73', 'Webbink85', 'Kharchenko13', 'Bica08', 'WEBDA14', 'ExtraDwarfs','ExtraClusters']
    for external_cat in external_cat_list:
        catalog = ugali.candidate.associate.catalogFactory(external_cat)
        external_cut = cut_from_hotspots(catalog['ra'], catalog['dec'], radius=0.1)
        healpix_mask &= external_cut

if apply_bsc:
    reader_bsc = pyfits.open('bsc5.fits')
    d_bsc = reader_bsc[1].data
    bsc_cut = cut_from_hotspots(d_bsc['RA'], d_bsc['DEC'], radius=0.1)
    healpix_mask &= bsc_cut


names = np.array(['ebv', 'Nilson', 'ext', 'bsc'])
bools = [apply_ebv, apply_Nilson, apply_external_cats, apply_bsc]
trailer = '_{}'.format('_'.join(names[bools]))

print('writing healpix_mask'+trailer+' ...')
healpy.write_map('masks/healpix_mask'+trailer+'.fits', healpix_mask, nest=NEST, dtype=bool, coord='C', overwrite=True)

pylab.figure()
healpy.mollview(healpix_mask, nest=True)
pylab.savefig('mask_plots/healpix_mask'+trailer)
