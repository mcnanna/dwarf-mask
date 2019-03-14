#!/usr/bin/env python

import numpy as np
import healpy
import astropy.io.fits as pyfits
import pylab

import ugali.utils.projector
import ugali.utils.healpix
import ugali.candidate.associate

pylab.ion()

infile = '/Users/mcnanna/Research/DES_luminosity/candidate_list.fits'
r = pyfits.open(infile)
d_original = r[1].data
r.close()

# Minimum significance threshold
d = d_original[d_original['SIG'] > 5.]


# Iterative approach
for ii in range(0, 5):
    # Consolidate nearby peaks
    match_1, match_2, angsep = ugali.utils.projector.match(d['ra'], d['dec'], d['ra'], d['dec'], tol=0.5, nnearest=2)
    index_exclude = np.where(d['sig'][match_1] > d['sig'][match_2], match_2, match_1)
    print len(index_exclude)
    cut_consolidate = np.tile(True, len(d))
    cut_consolidate[index_exclude] = False
    d = d[cut_consolidate]


pix = ugali.utils.healpix.angToPix(4096, d['ra'], d['dec'], nest=True)
mask = ugali.utils.healpix.read_map('healpix_mask.fits.gz', nest=True)

cut_dec = (d['DEC'] > -25.)
cut_modulus = (d['MODULUS'] < 21.75)
cut_ebv = np.where(mask[pix] & 0b0001, False, True)
cut_associate = np.where(mask[pix] & 0b0010, False, True)
cut_bsc = np.where(mask[pix] & 0b0100, False, True)
cut_size = (d['r'] > 0.025)
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
pylab.hist(d['SIG'][cut_ebv & cut_dec], bins=bins, color='blue', histtype='step', cumulative=-1, label='E(B-V) < 0.2 mag\n& Dec > -25 deg')
pylab.hist(d['SIG'][cut_ebv & cut_dec & cut_modulus], bins=bins, color='green', histtype='step', cumulative=-1, label='above & (m - M) < 21.75') # 23.25
#pylab.hist(d['SIG'][cut_ebv & cut_dec & cut_modulus & cut_Nilson], bins=bins, color='purple', histtype='step', cumulative=-1, label='above & no Nilson73 association') # 23.25
pylab.hist(d['SIG'][cut_ebv & cut_dec & cut_modulus & cut_associate], bins=bins, color='magenta', histtype='step', cumulative=-1, label='above & no galaxy association') # 23.25
pylab.hist(d['SIG'][cut_ebv & cut_dec & cut_modulus & cut_associate & cut_bsc], bins=bins, color='black', histtype='step', cumulative=-1, label='above & no BSC association') # 23.25
pylab.hist(d['SIG'][cut_ebv & cut_dec & cut_modulus & cut_associate & cut_bsc & cut_core], bins=bins, color='orange', histtype='step', cumulative=-1, label='above & ratio > 1.') # 23.25
pylab.hist(d['SIG'][cut_ebv & cut_dec & cut_modulus & cut_associate & cut_bsc & cut_core & cut_size], bins=bins, color='cyan', histtype='step', cumulative=-1, label='above & r > 0.025') # 23.25
pylab.legend(loc='upper right')
pylab.xlabel('SIG')
pylab.ylabel('Cumulative Counts')
pylab.savefig('significance_distribution.png')


# Skymap of candidates
infile_dust = '/Users/mcnanna/Research/DES_luminosity/ebv_sfd98_fullres_nside_4096_nest_equatorial.fits.gz'
ebv_map = ugali.utils.healpix.read_map(infile_dust, nest=True)

pylab.figure()
healpy.mollview(ebv_map, xsize=1600, min=0., max=0.5, cmap='binary', title='SIG > 10 Hotspots', unit='E(B-V)', nest=True)
healpy.graticule()
cut = (d['SIG'] > 10.)
healpy.projscatter(d['RA'][cut], d['DEC'][cut], lonlat=True, c='red', marker='o', edgecolor='none', s=2, vmax=0.5)
healpy.projscatter(d['RA'][cut & cut_ebv & cut_dec], d['DEC'][cut & cut_ebv & cut_dec], lonlat=True, c='blue', marker='o', edgecolor='none', s=2, vmax=0.5)
healpy.projscatter(d['RA'][cut & cut_ebv & cut_dec & cut_modulus], d['DEC'][cut & cut_ebv & cut_dec & cut_modulus], lonlat=True, c='green', marker='o', edgecolor='none', s=2, vmax=0.5)
healpy.projscatter(d['RA'][cut & cut_ebv & cut_dec & cut_modulus & cut_associate], d['DEC'][cut & cut_ebv & cut_dec & cut_modulus & cut_associate], lonlat=True, c='magenta', marker='o', edgecolor='none', s=10, vmax=0.5)
pylab.savefig('significance_map.png')


pylab.figure()
pylab.hist(d['MODULUS'][cut & cut_ebv & cut_dec], bins=np.arange(14.25, 25.25, 0.5), cumulative=False)
pylab.xlabel('m-M')
pylab.ylabel('Counts')
pylab.title('E(B-V) < 0.2 mag & Dec > -25 deg & SIG > 10')
pylab.savefig('modulus_distribution.png')


cut_sig = (d['SIG'] > 6.)
cut_final = cut_ebv & cut_dec & cut_modulus & cut_associate & cut_bsc & cut_core & cut_size & cut_sig


###############################################################
# Below is all stuff that was included in Keith's original code. 
# I'm not using it, but I didn't want to erase it either 
"""
for ii in np.random.choice(np.nonzero(cut_final)[0], 50, replace=False):
    print 'candidate_%.2f_%.2f.png   %.2f'%(d['RA'][ii], d['DEC'][ii], d['SIG'][ii])


d_final = d[cut_final]

model = d_final['N_MODEL'] / (np.pi * d_final['r']**2)
core = d_final['N_OBS_HALF'] / (np.pi * (0.5 * d_final['r'])**2)
full = d_final['N_OBS'] / (np.pi * d_final['r']**2)
ratio = (core - model) / (full - model)

mcconnachie15 = ugali.candidate.associate.McConnachie15() 
match_candidate, match_mcconnachie15, angsep = mcconnachie15.match(d_final['RA'], d_final['DEC'], coord='cel', tol=0.2)

for ii in range(0, len(mcconnachie15.data)):
    if ii in match_mcconnachie15:
        index = np.argmin(np.fabs(ii - match_mcconnachie15))
        sig = d_final['SIG'][match_candidate[index]]
        a = angsep[index]
        ratio_candidate = ratio[match_candidate[index]]
        r = d_final['r'][match_candidate[index]]
    else:
        sig = 0.
        a = 0.
        ratio_candidate = -9.
        r = 0.
    print '%30s%15.2f%15.3f%15.3f%15.3f'%(mcconnachie15.data['name'][ii], sig, r, ratio_candidate, a)


#match = ugali.utils.projector.match(151.77, 16.08, d_final['ra'], d_final['dec'], tol=0.1)[1]


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
