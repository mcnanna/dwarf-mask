#!/usr/bin/env python

import subprocess
import numpy as np
import healpy
import astropy.io.fits as pyfits
import pylab
pylab.ion()
import argparse
from textable import TexTable

import ugali.utils.projector
import ugali.utils.healpix
#import ugali.candidate.associate as associate
import associate

p = argparse.ArgumentParser()
p.add_argument('survey', help="'ps1' or 'des'")
p.add_argument('alg', help = "'u'/'ugali' or 's'/'simple'")
p.add_argument('--sig', type=float, help="Minimum significance/TS threshold.")
p.add_argument('--no_cross', action='store_true', help="Don't cross-check between the two algorithms. Use this is you only have the candidates for one algorithm")
p.add_argument('--no_textfile', action='store_true', help=argparse.SUPPRESS)
args = p.parse_args()
if 'u' in args.alg:
    args.alg = 'ugali'
elif 's' in args.alg:
    args.alg = 'simple'

subprocess.call('mkdir -p textfiles tables diagnostic_plots'.split())

infile = 'candidates/candidate_list_{0}_{1}.fits'.format(args.survey, args.alg)
r = pyfits.open(infile)
d_original = r[1].data
r.close()

# Minimum significance threshold (not used right now)
d = d_original#[d_original['SIG'] > 5.] [d_original['TS'] > 25?]

if args.alg == 'simple':
    SIG = 'SIG'
elif args.alg == 'ugali':
    SIG = 'TS'


### Define and apply cuts

# Consolidate nearby peaks, iterative approach
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

cut_ebv = np.where(mask[pix] & 0b00001, False, True)
cut_associate = np.where(mask[pix] & 0b00010, False, True)
cut_dwarfs = np.where(mask[pix] & 0b00100, False, True)
cut_bsc = np.where(mask[pix] & 0b01000, False, True)
cut_footprint = np.where(mask[pix] & 0b10000, False, True)

# Other cuts (modulus, size, shape)
if args.survey == 'ps1':
    cut_modulus = (d['MODULUS'] < 21.75)
elif args.survey == 'des':
    cut_modulus = (d['MODULUS'] < 23.5)

# Significance cut
if args.sig is not None:
    min_sig = args.sig
elif args.survey == 'ps1' and args.alg == 'simple':
    min_sig = 6.
elif args.survey == 'ps1' and args.alg == 'ugali':
    min_sig = 80. #TS
elif args.survey == 'des' and args.alg == 'simple':
    min_sig = 7.
elif args.survey == 'des' and args.alg == 'ugali':
    min_sig = 50. #TS
cut_sig = d[SIG] > min_sig

# Combine cuts which should filter everything but dwarfs
cut_bulk = cut_ebv & cut_footprint & cut_modulus & cut_associate & cut_bsc

# Final cut also filters known satellites and low sig detections
cut_final = cut_bulk & cut_dwarfs & cut_sig


# List of unassociated hotspots
# As a .txt file
f = open('textfiles/remains_{}_{}.txt'.format(args.survey, args.alg), 'w')
if args.alg == 'ugali':
    f.write('%20s'%('name'))
f.write('%10s%10s%10s%10s\n'%(SIG, 'ra', 'dec', 'modulus'))
for remainer in sorted(d[cut_final], key = lambda x: x[SIG], reverse=True):
    if args.alg == 'ugali':
        f.write('%20s'%(''.join(remainer['NAME'].split())))
    f.write('%10.2f%10.2f%10.2f%10.2f\n'%(remainer[SIG], remainer['ra'], remainer['dec'], remainer['modulus']))
f.close()

# As a LaTeX deluxetable
justs = 'cccc'
header_row1 = [SIG, r"$\alpha_{2000}$", r"$\delta_{2000}$", r"$m - M$"]
header_row2 = ['', '(deg)', '(deg)', '']
if args.alg == 'ugali':
    justs = 'l' + justs
    header_row1 = ['Name'] + header_row1
    header_row2 = [''] + header_row2
t = TexTable(len(justs), justs=justs, comments='', caption='')
t.add_header_row(header_row1)
t.add_header_row(header_row2)
t.add_data(np.genfromtxt('textfiles/remains_{}_{}.txt'.format(args.survey, args.alg), skip_header=1, dtype=object).T)
t.print_table('tables/remains_{}_{}.tex'.format(args.survey, args.alg))
subprocess.call("pdflatex -output-directory tables tables/remains_{}_{}.tex".format(args.survey, args.alg).split())


# Cross-check with other algorithm
if not args.no_cross:
    # Load other alg's results from textfile created above
    if args.alg == 'simple':
        try:
            d2 = np.genfromtxt('textfiles/remains_{}_ugali.txt'.format(args.survey), names=True)
        except IOError:
            subprocess.call('python filter_candidates.py {} ugali --no_textfile'.format(args.survey).split())
            d2 = np.genfromtxt('textfiles/remains_{}_ugali.txt'.format(args.survey), names=True)
    if args.alg == 'ugali':
        try:
            d2 = np.genfromtxt('textfiles/remains_{}_simple.txt'.format(args.survey), names=True)
        except IOError:
            subprocess.call('python filter_candidates.py {} simple --no_textfile'.format(args.survey).split())
            d2 = np.genfromtxt('textfiles/remains_{}_simple.txt'.format(args.survey), names=True)

    match1, match2, angseps = ugali.utils.projector.match(d['ra'][cut_final], d['dec'][cut_final], d2['ra'], d2['dec'], tol=0.2)
    matches = d[cut_final][match1]
    cut_cross = np.array([d[i]['ra'] in matches['ra'] for i in range(len(d))])

    if args.alg == 'ugali': # This 'if' is just so that it only writes this file once, although writing it twice isn't a big deal
        f = open('textfiles/remains_{}_both.txt'.format(args.survey), 'w')
        f.write('%20s%10s%10s%10s%10s%12s%12s%10s\n'%('name', 'TS', 'SIG', 'ra', 'dec', 'mod_ugali', 'mod_simple', 'angsep'))
        for i in range(len(angseps)):
            uga = d[cut_final][match1[i]]
            sim = d2[match2[i]]
            angsep = angseps[i]
            f.write('%20s%10.2f%10.2f%10.2f%10.2f%12.2f%12.2f%10.2f\n'%(''.join(uga['NAME'].split()), uga['TS'], sim['SIG'], uga['ra'], uga['dec'], uga['modulus'], sim['modulus'], angsep))
        f.close()

        t = TexTable(8, justs='lccccccc', comments='', caption='')
        t.add_header_row(['Name', 'TS', 'SIG', r"$\alpha_{2000}$", r"$\delta_{2000}$", r"$m - M$", r"$m - M$", "Angular Separation"])
        t.add_header_row(['', '(ugali)', '(simple)', '(deg)', '(deg)', '(ugali)', '(simple)', '(deg)'])
        t.add_data(np.genfromtxt('textfiles/remains_{}_both.txt'.format(args.survey), skip_header=1, dtype=object).T)
        t.print_table('tables/remains_{}_both.tex'.format(args.survey))
        subprocess.call("pdflatex -output-directory tables tables/remains_{}_both.tex".format(args.survey).split())

    if args.no_textfile:
        raise SystemExit(0)


### Signal Detection

# Detections of known satellites
# As a .txt file
def print_detections(known_dwarf_catalog, write='w'):
    known_dwarfs = associate.catalogFactory(known_dwarf_catalog)
    match_candidate, match_known_dwarfs, angsep = known_dwarfs.match(d['RA'], d['DEC'], coord='cel', tol=0.2)

    f = open('textfiles/signal_{}_{}.txt'.format(args.survey, args.alg), write)
    f.write('%25s%15s%15s%15s%10s%10s\n'%('name', SIG, 'modulus', 'angsep', 'bit', 'cut'))

    for ii in range(0, len(known_dwarfs.data)):
        if ii in match_known_dwarfs:
            index = np.argmin(np.fabs(ii - match_known_dwarfs))
            sig = d[SIG][match_candidate[index]]
            a = angsep[index]
            modulus = d['MODULUS'][match_candidate[index]]
            bit = mask[pix[match_candidate[index]]]
            wascut = '' if cut_bulk[match_candidate[index]] else 'cut'
        else:
            sig = 0.
            a = 0.
            modulus= 0.
            bit = mask[ugali.utils.healpix.angToPix(4096, known_dwarfs.data['ra'][ii], known_dwarfs.data['dec'][ii], nest=True)]
            wascut = ''
        if not (bit & 0b10000): # Don't bother writing results for things outside the footprint
            f.write('%25s%15.2f%15.2f%15.3f%10i%10s\n'%(known_dwarfs.data['name'][ii], sig, modulus, a, bit, wascut))
    f.close()

print_detections("McConnachie15", 'w')
print_detections("ExtraDwarfs", 'a')

# As a LaTeX deluxetable
def signal_table():
    ncols = 5
    justs = 'lcccc'
    caption = '' #'\___caption' for paper
    comments = '' #'\___comments' for paper

    t = TexTable(ncols, justs=justs, caption=caption, comments=comments)
    t.add_header_row(['Name', SIG, r'$m - M$', 'Angular Separation', ''])
    t.add_header_row(['','','','(deg)', ''])

    def get_signal(signal, known_dwarf_catalog):
        known_dwarfs = associate.catalogFactory(known_dwarf_catalog)
        match_candidate, match_known_dwarfs, angsep = known_dwarfs.match(d['RA'], d['DEC'], coord='cel', tol=0.2)

        rows = []
        for ii in range(0, len(known_dwarfs.data)):
            if ii in match_known_dwarfs:
                index = np.argmin(np.fabs(ii - match_known_dwarfs))
                sig = d[SIG][match_candidate[index]]
                a = angsep[index]
                modulus = d['MODULUS'][match_candidate[index]]
                bit = mask[pix[match_candidate[index]]]
                wascut = '' if cut_bulk[match_candidate[index]] else 'cut'
            else:
                sig = 0.
                a = 0.
                modulus= 0.
                bit = mask[ugali.utils.healpix.angToPix(4096, known_dwarfs.data['ra'][ii], known_dwarfs.data['dec'][ii], nest=True)]
                wascut = ''
            if not (bit & 0b10000): # Don't bother writing results for things outside the footprint
                signal.append([known_dwarfs.data['name'][ii], sig, modulus, a, wascut])

    signal = []
    get_signal(signal, 'McConnachie15')
    get_signal(signal, 'ExtraDwarfs')

    galaxy_names = []
    cut_duplicates = []
    for galaxy in signal:
        cut_duplicates.append(galaxy[0] not in galaxy_names)
        galaxy_names.append(galaxy[0])
    
    t.add_data((np.array(signal, dtype=object)[np.array(cut_duplicates)]).T)
    t.print_table('tables/signal_{}_{}.tex'.format(args.survey, args.alg))
    subprocess.call("pdflatex -output-directory tables tables/signal_{}_{}.tex".format(args.survey, args.alg).split())
    
signal_table()

# If a known satellite matches something in another catalog (you can tell from the bit), use this to find out what it's matching
def what_matches(name, known_dwarf_catalog='McConnachie15'):
    known_dwarfs = associate.catalogFactory(known_dwarf_catalog)
    obj = known_dwarfs[known_dwarfs['name'] == name]
    ra, dec = obj['ra'], obj['dec']
    external_cat_list = ['Harris96', 'Corwen04', 'Nilson73', 'Webbink85', 'Kharchenko13', 'Bica08', 'WEBDA14', 'ExtraClusters']
    print '%15s%20s%10s%10s'%('catalog', 'name', 'radius', 'angsep')
    for cat in external_cat_list:
        catalog = associate.catalogFactory(cat)
        whocares, matches, angseps = catalog.match(ra, dec, coord='cel', tol=6.0)
        for i in range(len(matches)):
            index = matches[i]
            angsep = angseps[i]
            radius = catalog[index]['radius']
            print '%15s%20s%10.3f%10.3f'%(cat, catalog[index]['name'], radius, angsep) 


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
pylab.hist(d[SIG][cut_ebv & cut_footprint & cut_modulus & cut_associate & cut_bsc], bins=bins, color='black', histtype='step', cumulative=-1, label='above & no external association') # = cut_bulk
pylab.hist(d[SIG][cut_bulk & cut_dwarfs], bins=bins, color='purple', histtype='step', cumulative=-1, label='above & no dwarf assocation')
pylab.hist(d[SIG][cut_bulk & cut_dwarfs & cut_cross], bins=bins, color='darkturquoise', histtype='step', cumulative=-1, label='above & found by both algorithms') 
pylab.legend(loc='upper right')
pylab.xlabel(SIG)
pylab.ylabel('Cumulative Counts')
pylab.savefig('diagnostic_plots/significance_distribution_{}_{}.png'.format(args.survey, args.alg), bbox_inches='tight')

# Skymap of candidates
infile_dust = 'ebv_sfd98_fullres_nside_4096_nest_equatorial.fits.gz'
ebv_map = ugali.utils.healpix.read_map(infile_dust, nest=True)

pylab.figure()
healpy.mollview(ebv_map, xsize=1600, min=0., max=0.5, cmap='binary', title='{0} > {1} Hotspots'.format(SIG, min_sig), unit='E(B-V)', nest=True)
healpy.graticule()
healpy.projscatter(d['RA'][cut_sig], d['DEC'][cut_sig], lonlat=True, c='red', marker='o', edgecolor='none', s=2, vmax=0.5)
healpy.projscatter(d['RA'][cut_sig & cut_ebv & cut_footprint], d['DEC'][cut_sig & cut_ebv & cut_footprint], lonlat=True, c='blue', marker='o', edgecolor='none', s=2, vmax=0.5)
healpy.projscatter(d['RA'][cut_sig & cut_ebv & cut_footprint & cut_modulus], d['DEC'][cut_sig & cut_ebv & cut_footprint & cut_modulus], lonlat=True, c='green', marker='o', edgecolor='none', s=2, vmax=0.5)
healpy.projscatter(d['RA'][cut_final & ~cut_cross], d['DEC'][cut_final & ~cut_cross], lonlat=True, c='purple', marker='o', edgecolor='none', s=10, vmax=0.5)
healpy.projscatter(d['RA'][cut_final & cut_cross], d['DEC'][cut_final & cut_cross], lonlat=True, c='darkturquoise', marker='*', edgecolor='none', s=40, vmax=0.5)
pylab.savefig('diagnostic_plots/significance_map_{}_{}.png'.format(args.survey, args.alg), bbox_inches='tight')

# Modulus distribution of remainers
pylab.figure()
pylab.hist(d['MODULUS'][cut_final], bins=np.arange(14.25, 25.25, 0.5), cumulative=False)
pylab.xlabel('m-M')
pylab.ylabel('Counts')
pylab.title('Unassociated Hotspots')
pylab.savefig('diagnostic_plots/modulus_distribution_{}_{}.png'.format(args.survey, args.alg), bbox_inches='tight')

#print "Passed cuts:", sum(cut_bulk & cut_dwarfs)
#print "Passed sig:", sum(cut_sig)
#print "Passed final:", sum(cut_final)
#print "Passed cross:", sum(cut_cross)
