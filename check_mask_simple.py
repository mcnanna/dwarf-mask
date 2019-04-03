#!/usr/bin/env python

import numpy as np
import healpy
import astropy.io.fits as pyfits
import pylab
import argparse

import ugali.utils.projector
import ugali.utils.healpix
#import ugali.candidate.associate as associate
import associate

p = argparse.ArgumentParser()
p.add_argument('survey', help="'ps1' or 'des'")
args = p.parse_args()

pylab.ion()

infile = '/Users/mcnanna/Research/DES_luminosity/candidate_list_{}_simple.fits'.format(args.survey)
r = pyfits.open(infile)
d_original = r[1].data
r.close()

# Minimum significance threshold
d = d_original#[d_original['SIG'] > 5.]


# Iterative approach
for ii in range(0, 6):
    # Consolidate nearby peaks
    match_1, match_2, angsep = ugali.utils.projector.match(d['ra'], d['dec'], d['ra'], d['dec'], tol=0.5, nnearest=2)
    index_exclude = np.where(d['sig'][match_1] > d['sig'][match_2], match_2, match_1)
    print len(index_exclude)
    cut_consolidate = np.tile(True, len(d))
    cut_consolidate[index_exclude] = False
    d = d[cut_consolidate]

pix = ugali.utils.healpix.angToPix(4096, d['ra'], d['dec'], nest=True)
mask = ugali.utils.healpix.read_map('healpix_mask_{}.fits'.format(args.survey), nest=True)

#cut_footprint = (d['DEC'] > -25.)
cut_footprint = np.where(mask[pix] & 0b10000, False, True)
if args.survey == 'ps1':
    cut_modulus = (d['MODULUS'] < 21.75)
elif args.survey == 'des':
    cut_modulus = (d['MODULUS'] < 23.5)
cut_ebv = np.where(mask[pix] & 0b0001, False, True)
cut_associate = np.where(mask[pix] & 0b0010, False, True)
cut_dwarfs = np.where(mask[pix] & 0b0100, False, True)
cut_bsc = np.where(mask[pix] & 0b1000, False, True)
cut_size = (d['r'] >= 0.020)
model = d['N_MODEL'] / (np.pi * d['r']**2)
core = d['N_OBS_HALF'] / (np.pi * (0.5 * d['r'])**2)
full = d['N_OBS'] / (np.pi * d['r']**2)
ratio = (core - model) / (full - model)
cut_core = (ratio > 1.)


# Significance histogram 
bins = np.arange(5., 40.25, 0.5)
pylab.figure()
pylab.yscale('log')
pylab.hist(d['SIG'], bins=bins, color='red', histtype='step', cumulative=-1, label='All')
pylab.hist(d['SIG'][cut_ebv & cut_footprint], bins=bins, color='blue', histtype='step', cumulative=-1, label='E(B-V) < 0.2 mag\n& in footprint')
pylab.hist(d['SIG'][cut_ebv & cut_footprint & cut_modulus], bins=bins, color='green', histtype='step', cumulative=-1, label='above & (m - M) < 21.75') # 23.25
pylab.hist(d['SIG'][cut_ebv & cut_footprint & cut_modulus & cut_bsc], bins=bins, color='magenta', histtype='step', cumulative=-1, label='above & no bsc association') 
pylab.hist(d['SIG'][cut_ebv & cut_footprint & cut_modulus & cut_associate & cut_bsc], bins=bins, color='black', histtype='step', cumulative=-1, label='above & no external association') 
pylab.hist(d['SIG'][cut_ebv & cut_footprint & cut_modulus & cut_associate & cut_bsc & cut_core], bins=bins, color='orange', histtype='step', cumulative=-1, label='above & ratio > 1.') 
pylab.hist(d['SIG'][cut_ebv & cut_footprint & cut_modulus & cut_associate & cut_bsc & cut_core & cut_size], bins=bins, color='cyan', histtype='step', cumulative=-1, label='above & r > 0.025') 
pylab.hist(d['SIG'][cut_ebv & cut_footprint & cut_modulus & cut_associate & cut_bsc & cut_core & cut_size & cut_dwarfs], bins=bins, color='yellow', histtype='step', cumulative=-1, label='above & no dwarf assocation')
pylab.legend(loc='upper right')
pylab.xlabel('SIG')
pylab.ylabel('Cumulative Counts')
pylab.savefig('significance_distribution_{}_simple.png'.format(args.survey))

"""
# Skymap of candidates
infile_dust = '/Users/mcnanna/Research/DES_luminosity/ebv_sfd98_fullres_nside_4096_nest_equatorial.fits.gz'
ebv_map = ugali.utils.healpix.read_map(infile_dust, nest=True)

pylab.figure()
healpy.mollview(ebv_map, xsize=1600, min=0., max=0.5, cmap='binary', title='SIG > 10 Hotspots', unit='E(B-V)', nest=True)
healpy.graticule()
cut = (d['SIG'] > 10.)
healpy.projscatter(d['RA'][cut], d['DEC'][cut], lonlat=True, c='red', marker='o', edgecolor='none', s=2, vmax=0.5)
healpy.projscatter(d['RA'][cut & cut_ebv & cut_footprint], d['DEC'][cut & cut_ebv & cut_footprint], lonlat=True, c='blue', marker='o', edgecolor='none', s=2, vmax=0.5)
healpy.projscatter(d['RA'][cut & cut_ebv & cut_footprint & cut_modulus], d['DEC'][cut & cut_ebv & cut_footprint & cut_modulus], lonlat=True, c='green', marker='o', edgecolor='none', s=2, vmax=0.5)
healpy.projscatter(d['RA'][cut & cut_ebv & cut_footprint & cut_modulus & cut_associate], d['DEC'][cut & cut_ebv & cut_footprint & cut_modulus & cut_associate], lonlat=True, c='magenta', marker='o', edgecolor='none', s=10, vmax=0.5)
pylab.savefig('significance_map.png')


pylab.figure()
pylab.hist(d['MODULUS'][cut & cut_ebv & cut_footprint], bins=np.arange(14.25, 25.25, 0.5), cumulative=False)
pylab.xlabel('m-M')
pylab.ylabel('Counts')
pylab.title('E(B-V) < 0.2 mag & Dec > -25 deg & SIG > 10')
pylab.savefig('modulus_distribution.png')


#cut_sig = (d['SIG'] > 6.)
cut_final = cut_ebv & cut_footprint & cut_modulus & cut_associate & cut_bsc & cut_core & cut_size

known_dwarfs = associate.McConnachie15() 
#known_dwarfs = associate.ExtraDwarfs()
match_candidate, match_known_dwarfs, angsep = known_dwarfs.match(d['RA'], d['DEC'], coord='cel', tol=0.2)

print '%30s%15s%15s%15s%15s%15s%10s%10s'%('name', 'sig', 'modulus', 'r', 'ratio', 'angsep', 'bit', 'cut')
for ii in range(0, len(known_dwarfs.data)):
    if ii in match_known_dwarfs:
        index = np.argmin(np.fabs(ii - match_known_dwarfs))
        sig = d['SIG'][match_candidate[index]]
        a = angsep[index]
        ratio_candidate = ratio[match_candidate[index]]
        r = d['r'][match_candidate[index]]
        modulus = d['MODULUS'][match_candidate[index]]
        bit = mask[pix[match_candidate[index]]]
        wascut = '' if cut_final[match_candidate[index]] else 'cut'
    else:
        sig = 0.
        a = 0.
        ratio_candidate = -9.
        r = 0.
        modulus= 0.
        bit = mask[ugali.utils.healpix.angToPix(4096, known_dwarfs.data['ra'][ii], known_dwarfs.data['dec'][ii], nest=True)]
        wascut = ''
    if bit < 16:
        print '%30s%15.2f%15.2f%15.3f%15.3f%15.3f%10i%10s'%(known_dwarfs.data['name'][ii], sig, modulus, r, ratio_candidate, a, bit, wascut)


def what_matches(name):
    obj = known_dwarfs[known_dwarfs['name'] == name]
    ra, dec = obj['ra'], obj['dec']
    external_cat_list = ['Harris96', 'Corwen04', 'Nilson73', 'Webbink85', 'Kharchenko13', 'Bica08', 'WEBDA14', 'ExtraClusters']
    for cat in external_cat_list:
        catalog = associate.catalogFactory(cat)
        whocares, matches, angseps = catalog.match(ra, dec, coord='cel', tol=6.0)
        for i in range(len(matches)):
            index = matches[i]
            angsep = angseps[i]
            radius = catalog[index]['radius']
            print '%15s%15s%10.3f%10.3f'%(cat, catalog[index]['name'], radius, angsep) 

"""

###############################################################
# Below is all stuff that was included in Keith's original code. 
# I'm not using it, but I didn't want to erase it either 
"""
d_final = d[cut_final]


pylab.figure()
#pylab.scatter(d_final['SIG'], d_final['r'], c=d_final['N_OBS'], s=2, vmax=50)
pylab.hexbin(d_final['SIG'], d_final['r'], mincnt=1)
pylab.colorbar()


pylab.figure()
#pylab.scatter(d_final['SIG'], core / full, s=2)
pylab.hexbin(d_final['SIG'], ratio, mincnt=1)

pylab.figure()
pylab.hist(ratio, bins=51)


def uniqueIdentifier(data):
    #return np.round(np.random.random(10)).astype(int)(data['ra'] * 1.e6) + data['dec']
    #return np.round((data['ra'] * 1.e3) + (data['dec'] * 1.e9)).astype(int)
    #return '%.3f%.3f'%(data['ra'], data['dec'])
    #return 
    return ugali.utils.healpix.angToPix(2**14, data['ra'], data['dec'])


unique_identifier_original = uniqueIdentifier(d_original)

unique_identifier_final = uniqueIdentifier(d_final)
cut_subset = np.in1d(unique_identifier_original, unique_identifier_final)

print int(np.sum(cut_subset))
print len(d_final)
assert int(np.sum(cut_subset)) == len(d_final)

hdu = pyfits.BinTableHDU(d_original[cut_subset])
#hdu.writeto('candidate_list_subset.fits', clobber=True)

#matcR_1, match_2, angsep = ugali.utils.projector.match(d_final['ra'], d_final['dec'], d_final['ra'], d_final['dec'], tol=0.01, nnearest=2)
"""
