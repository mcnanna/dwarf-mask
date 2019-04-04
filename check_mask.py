#!/usr/bin/env python

import numpy as np
import healpy
import astropy.io.fits as pyfits
import pylab
pylab.ion()
import argparse

import ugali.utils.projector
import ugali.utils.healpix
#import ugali.candidate.associate as associate
import associate

p = argparse.ArgumentParser()
p.add_argument('survey', help="'ps1' or 'des'")
p.add_argument('alg', help = "'u'/'ugali' or 's'/'simple'")
args = p.parse_args()
if 'u' in args.alg:
    args.alg = 'ugali'
elif 's' in args.alg:
    args.alg = 'simple'


infile = '/Users/mcnanna/Research/DES_luminosity/candidate_list_{0}_{1}.fits'.format(args.survey, args.alg)
r = pyfits.open(infile)
d_original = r[1].data
r.close()

# Minimum significance threshold
d = d_original#[d_original['SIG'] > 5.] [d_original['TS'] > 25?]

if args.alg == 'simple':
    SIG = 'SIG'
elif args.alg == 'ugali':
    SIG = 'TS'


### Define and apply cuts

# Consolidate nearby peaks, iterave approach
done = False
while not done:
    match_1, match_2, angsep = ugali.utils.projector.match(d['ra'], d['dec'], d['ra'], d['dec'], tol=0.5, nnearest=2)
    index_exclude = np.where(d[SIG][match_1] > d[SIG][match_2], match_2, match_1)
    if len(index_exclude) == 0:
        done = True
    cut_consolidate = np.tile(True, len(d))
    cut_consolidate[index_exclude] = False
    d = d[cut_consolidate]

# Geometric cuts
pix = ugali.utils.healpix.angToPix(4096, d['ra'], d['dec'], nest=True)
mask = ugali.utils.healpix.read_map('healpix_mask_{}.fits.gz'.format(args.survey), nest=True)

cut_footprint = np.where(mask[pix] & 0b10000, False, True)
cut_ebv = np.where(mask[pix] & 0b0001, False, True)
cut_associate = np.where(mask[pix] & 0b0010, False, True)
cut_dwarfs = np.where(mask[pix] & 0b0100, False, True)
cut_bsc = np.where(mask[pix] & 0b1000, False, True)

# Other cuts (modulous, size, shape)
if args.survey == 'ps1':
    cut_modulus = (d['MODULUS'] < 21.75)
elif args.survey == 'des':
    cut_modulus = (d['MODULUS'] < 23.5)

if args.alg == 'simple':
    cut_size = (d['r'] >= 0.020)
    model = d['N_MODEL'] / (np.pi * d['r']**2)
    core = d['N_OBS_HALF'] / (np.pi * (0.5 * d['r'])**2)
    full = d['N_OBS'] / (np.pi * d['r']**2)
    ratio = (core - model) / (full - model)
    cut_core = (ratio > 1.)

# Significance cut
if args.survey == 'ps1' and args.alg == 'simple':
    min_sig = 6.
elif args.survey == 'ps1' and args.alg == 'ugali':
    min_sig = 80 #TS
elif args.survey == 'des' and args.alg == 'simple':
    min_sig = 7.
elif args.survey == 'des' and args.alg == 'ugali':
    min_sig = 50 #TS
cut_sig = d[SIG] > min_sig

cut_bulk = cut_ebv & cut_footprint & cut_modulus & cut_associate & cut_bsc
if args.alg == 'simple':
    cut_bulk = cut_bulk & cut_core & cut_size

cut_final = cut_bulk & cut_dwarfs & cut_sig

# List of unassociated hotspots
f = open('textfiles/remainers_{}_{}.txt'.format(args.survey, args.alg), 'w')
f.write('%10s%10s%10s%10s\n'%(SIG, 'ra', 'dec', 'modulus'))
for remainer in d[cut_final]:
    f.write('%10.2f%10.2f%10.2f%10.2f\n'%(remainer[SIG], remainer['ra'], remainer['dec'], remainer['modulus']))
f.close()


# Cross-check with other algorithm
# Note: this will cause an error if the text other alg's text file hasn't been written yet (see code above), so comment it out if you haven't made it yet
if args.alg == 'simple':
    d2 = np.genfromtxt('textfiles/remainers_{}_ugali.txt'.format(args.survey), names=True)
if args.alg == 'ugali':
    d2 = np.genfromtxt('textfiles/remainers_{}_simple.txt'.format(args.survey), names=True)

match1, match2, angseps = ugali.utils.projector.match(d['ra'][cut_final], d['dec'][cut_final], d2['ra'], d2['dec'], tol=0.5)
matches = d[cut_final][match1]
cut_cross = np.array([d[i]['ra'] in matches['ra'] for i in range(len(d))])


### Signal Detection


# Detections of known satellites
def print_detections(known_dwarf_catalog, write='w'):
    known_dwarfs = associate.catalogFactory(known_dwarf_catalog)
    match_candidate, match_known_dwarfs, angsep = known_dwarfs.match(d['RA'], d['DEC'], coord='cel', tol=0.2)

    f = open('textfiles/signal_{}_{}.txt'.format(args.survey, args.alg), write)
    if args.alg == 'simple':
        f.write('%25s%15s%15s%15s%15s%15s%10s%10s\n'%('name', SIG, 'modulus', 'r', 'ratio', 'angsep', 'bit', 'cut'))
    elif args.alg == 'ugali':
        f.write('%25s%15s%15s%15s%10s%10s\n'%('name', SIG, 'modulus', 'angsep', 'bit', 'cut'))

    for ii in range(0, len(known_dwarfs.data)):
        if ii in match_known_dwarfs:
            index = np.argmin(np.fabs(ii - match_known_dwarfs))
            sig = d[SIG][match_candidate[index]]
            a = angsep[index]
            if args.alg == 'simple':
                ratio_candidate = ratio[match_candidate[index]]
                r = d['r'][match_candidate[index]]
            modulus = d['MODULUS'][match_candidate[index]]
            bit = mask[pix[match_candidate[index]]]
            wascut = '' if cut_bulk[match_candidate[index]] else 'cut'
        else:
            sig = 0.
            a = 0.
            if args.alg == 'simple':
                ratio_candidate = -9.
                r = 0.
            modulus= 0.
            bit = mask[ugali.utils.healpix.angToPix(4096, known_dwarfs.data['ra'][ii], known_dwarfs.data['dec'][ii], nest=True)]
            wascut = ''
        if bit < 16:
            if args.alg == 'simple':
                f.write('%25s%15.2f%15.2f%15.3f%15.3f%15.3f%10i%10s\n'%(known_dwarfs.data['name'][ii], sig, modulus, r, ratio_candidate, a, bit, wascut))
            elif args.alg == 'ugali':
                f.write('%25s%15.2f%15.2f%15.3f%10i%10s\n'%(known_dwarfs.data['name'][ii], sig, modulus, a, bit, wascut))
    f.close()

print_detections('McConnachie15', 'w')
print_detections('ExtraDwarfs', 'a')

# If a known satellite matches something in another catalog (you can tell from the bit), use this to find out what it's matching
def what_matches(name, known_dwarf_catalog='McConnachie15'):
    known_dwarfs = associate.catalogFactory(known_dwarf_catalog)
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

# Significance histogram 
pylab.figure()
if args.alg == 'simple':
    bins = np.arange(5., 40.25, 0.5)
elif args.alg == 'ugali':
    bins = np.logspace(np.log10(25.), np.log10(200000.), num=70)
    pylab.xscale('log')
pylab.yscale('log')
pylab.hist(d[SIG], bins=bins, color='red', histtype='step', cumulative=-1, label='All')
pylab.hist(d[SIG][cut_ebv & cut_footprint], bins=bins, color='blue', histtype='step', cumulative=-1, label='E(B-V) < 0.2 mag\n& in footprint')
pylab.hist(d[SIG][cut_ebv & cut_footprint & cut_modulus], bins=bins, color='green', histtype='step', cumulative=-1, label='above & (m - M) < {}'.format(21.75 if args.survey == 'ps1' else 23.5))
pylab.hist(d[SIG][cut_ebv & cut_footprint & cut_modulus & cut_bsc], bins=bins, color='orange', histtype='step', cumulative=-1, label='above & no bsc association') 
pylab.hist(d[SIG][cut_ebv & cut_footprint & cut_modulus & cut_associate & cut_bsc], bins=bins, color='black', histtype='step', cumulative=-1, label='above & no external association') # = cut_bulk for ugali
if args.alg == 'simple':
    pylab.hist(d[SIG][cut_ebv & cut_footprint & cut_modulus & cut_associate & cut_bsc & cut_core], bins=bins, color='brown', histtype='step', cumulative=-1, label='above & ratio > 1.') 
    pylab.hist(d[SIG][cut_ebv & cut_footprint & cut_modulus & cut_associate & cut_bsc & cut_core & cut_size], bins=bins, color='olive', histtype='step', cumulative=-1, label='above & r > 0.02') # = cut_bulk for simple
    pylab.hist(d[SIG][cut_bulk & cut_dwarfs], bins=bins, color='purple', histtype='step', cumulative=-1, label='above & no dwarf assocation')
elif args.alg == 'ugali':
    pylab.hist(d[SIG][cut_bulk & cut_dwarfs], bins=bins, color='purple', histtype='step', cumulative=-1, label='above & no dwarf assocation')
pylab.hist(d[SIG][cut_bulk & cut_dwarfs & cut_cross], bins=bins, color='darkturquoise', histtype='step', cumulative=-1, label='above & found by both algorithms') 
pylab.legend(loc='upper right')
pylab.xlabel(SIG)
pylab.ylabel('Cumulative Counts')
pylab.savefig('plots/significance_distribution_{}_{}.png'.format(args.survey, args.alg), bbox_inches='tight')

# Skymap of candidates
infile_dust = '/Users/mcnanna/Research/DES_luminosity/ebv_sfd98_fullres_nside_4096_nest_equatorial.fits.gz'
ebv_map = ugali.utils.healpix.read_map(infile_dust, nest=True)

pylab.figure()
healpy.mollview(ebv_map, xsize=1600, min=0., max=0.5, cmap='binary', title='{0} > {1} Hotspots'.format(SIG, min_sig), unit='E(B-V)', nest=True)
healpy.graticule()
healpy.projscatter(d['RA'][cut_sig], d['DEC'][cut_sig], lonlat=True, c='red', marker='o', edgecolor='none', s=2, vmax=0.5)
healpy.projscatter(d['RA'][cut_sig & cut_ebv & cut_footprint], d['DEC'][cut_sig & cut_ebv & cut_footprint], lonlat=True, c='blue', marker='o', edgecolor='none', s=2, vmax=0.5)
healpy.projscatter(d['RA'][cut_sig & cut_ebv & cut_footprint & cut_modulus], d['DEC'][cut_sig & cut_ebv & cut_footprint & cut_modulus], lonlat=True, c='green', marker='o', edgecolor='none', s=2, vmax=0.5)
healpy.projscatter(d['RA'][cut_final & ~cut_cross], d['DEC'][cut_final & ~cut_cross], lonlat=True, c='purple', marker='o', edgecolor='none', s=10, vmax=0.5)
healpy.projscatter(d['RA'][cut_final & cut_cross], d['DEC'][cut_final & cut_cross], lonlat=True, c='darkturquoise', marker='*', edgecolor='none', s=40, vmax=0.5)
pylab.savefig('plots/significance_map_{}_{}.png'.format(args.survey, args.alg), bbox_inches='tight')

# Modulus distribution of remainers
pylab.figure()
pylab.hist(d['MODULUS'][cut_final], bins=np.arange(14.25, 25.25, 0.5), cumulative=False)
pylab.xlabel('m-M')
pylab.ylabel('Counts')
pylab.title('Unassociated Hotspots')
pylab.savefig('plots/modulus_distribution_{}_{}.png'.format(args.survey, args.alg), bbox_inches='tight')


###############################################################
# Below is all stuff that was included in Keith's original code. 
# I'm not using it, but I didn't want to erase it either 
"""
d_final = d[cut_final]


pylab.figure()
#pylab.scatter(d_final[SIG], d_final['r'], c=d_final['N_OBS'], s=2, vmax=50)
pylab.hexbin(d_final[SIG], d_final['r'], mincnt=1)
pylab.colorbar()


pylab.figure()
#pylab.scatter(d_final[SIG], core / full, s=2)
pylab.hexbin(d_final[SIG], ratio, mincnt=1)

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

#match_1, match_2, angsep = ugali.utils.projector.match(d_final['ra'], d_final['dec'], d_final['ra'], d_final['dec'], tol=0.01, nnearest=2)
"""
