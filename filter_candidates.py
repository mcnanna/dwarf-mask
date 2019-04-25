#!/usr/bin/env python

import subprocess
import numpy as np
import healpy
import astropy.io.fits as pyfits
import pylab
pylab.ion()
from matplotlib.lines import Line2D
import argparse
from textable import TexTable

import ugali.utils.projector
import ugali.utils.healpix
import associate

import pdb

p = argparse.ArgumentParser()
p.add_argument('survey', help="'ps1' or 'des'")
p.add_argument('alg', help = "'u'/'ugali' or 's'/'simple'")
p.add_argument('--sig', type=float, help="Minimum significance/TS threshold.")
p.add_argument('--no_cross', action='store_true', help="Don't cross-check between the two algorithms. Use this is you only have the candidates for one algorithm")
p.add_argument('--no_fitsfile', action='store_true', help=argparse.SUPPRESS)
args = p.parse_args()
if 'u' in args.alg:
    args.alg = 'ugali'
elif 's' in args.alg:
    args.alg = 'simple'

subprocess.call('mkdir -p fits_files tables diagnostic_plots'.split())

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
if args.survey == 'ps1' and args.alg == 'ugali':
    cut_footprint &= np.where(mask[pix] & 0b100000, False, True) # Add failures to footprint cut 

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


# Write sorted remains to .fits and make a TeX table
pyfits.writeto('fits_files/remains_{}_{}.fits'.format(args.survey, args.alg), np.sort(d[cut_final], order=SIG)[::-1], overwrite=True)

justs = 'cccc'
header_row1 = [SIG, r"$\alpha_{2000}$", r"$\delta_{2000}$", r"$m - M$"]
header_row2 = ['', '(deg)', '(deg)', '']
data_headers = [SIG, 'ra', 'dec', 'modulus']
sigfigs = [3, 5, 4, 3]
if args.alg == 'ugali':
    justs = 'l' + justs
    header_row1 = ['Name'] + header_row1
    header_row2 = [''] + header_row2
    data_headers = ['name'] + data_headers
    sigfigs = [0] + sigfigs
t = TexTable(len(justs), justs=justs, comments='', caption='')
t.add_header_row(header_row1)
t.add_header_row(header_row2)
data = np.sort(d[cut_final], order=SIG)[::-1]
t.add_data([data[header] for header in data_headers], sigfigs=sigfigs)
t.print_table('tables/remains_{}_{}.tex'.format(args.survey, args.alg))
subprocess.call("pdflatex -output-directory tables tables/remains_{}_{}.tex".format(args.survey, args.alg).split())


### Signal Detection

# Detections of known satellites
signal = []
for known_dwarf_catalog in ['McConnachie15', 'ExtraDwarfs']:
    known_dwarfs = associate.catalogFactory(known_dwarf_catalog)
    match_candidate, match_known_dwarfs, angsep = known_dwarfs.match(d['RA'], d['DEC'], coord='cel', tol=0.2)

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
            name = known_dwarfs.data['name'][ii]
            if (len(signal) == 0) or (name not in np.array(signal)[:, 0]): # Try to avoid duplicates from the mutliple catalogs
                signal.append((known_dwarfs.data['name'][ii], sig, modulus, a, wascut, bit))

dtype = [('name','|S18'),(SIG, float),('modulus',float),('angsep',float),('cut','|S3'),('bit',int)]
signal = np.array(signal, dtype=dtype)
pyfits.writeto('fits_files/signal_{}_{}.fits'.format(args.survey, args.alg), signal, overwrite=True)

justs = 'lcccc'
t = TexTable(len(justs), justs=justs, caption='', comments='')
t.add_header_row(['Name', SIG, r'$m - M$', 'Angular Separation', ''])
t.add_header_row(['','','','(deg)', ''])
sigfigs=[0, 3, 3, 2, 0]
t.add_data([signal[header] for header in signal.dtype.names[:-1]], sigfigs=sigfigs)
t.print_table('tables/signal_{}_{}.tex'.format(args.survey, args.alg))
subprocess.call("pdflatex -output-directory tables tables/signal_{}_{}.tex".format(args.survey, args.alg).split())


### Combine results of two algorithsm

if not args.no_cross:
    # Load other alg's results from fits created above

    # THIS ASSUMES DEFAULT SIGNIFICANCES DEFINED TOWARDS THE TOP. IF YOU WANT TO USE DIFFERENT SIGS, RUN ONCE
    # USING args.no_cross THEN RE-RUN ONCE THE NECESSARY .fits FILES ARE CREATED
    if args.alg == 'simple':
        try:
            d2 = pyfits.open('fits_files/remains_{}_ugali.fits'.format(args.survey))[1].data
            signal2 = pyfits.open('fits_files/signal_{}_ugali.fits'.format(args.survey))[1].data
        except IOError:
            subprocess.call('python filter_candidates.py {} ugali --no_fitsfile'.format(args.survey).split())
            d2 = pyfits.open('fits_files/remains_{}_ugali.fits'.format(args.survey))[1].data
            signal2 = pyfits.open('fits_files/signal_{}_ugali.fits'.format(args.survey))[1].data
    if args.alg == 'ugali':
        try:
            d2 = pyfits.open('fits_files/remains_{}_simple.fits'.format(args.survey))[1].data
            signal2 = pyfits.open('fits_files/signal_{}_simple.fits'.format(args.survey))[1].data
        except IOError:
            subprocess.call('python filter_candidates.py {} simple --no_fitsfile'.format(args.survey).split())
            d2 = pyfits.open('fits_files/remains_{}_simple.fits'.format(args.survey))[1].data
            signal2 = pyfits.open('fits_files/signal_{}_simple.fits'.format(args.survey))[1].data

    match1, match2, angseps = ugali.utils.projector.match(d['ra'][cut_final], d['dec'][cut_final], d2['ra'], d2['dec'], tol=0.2)
    matches = d[cut_final][match1]
    cut_cross = np.array([d[i]['ra'] in matches['ra'] for i in range(len(d))])


    # Remains
    if args.alg == 'ugali': # This 'if' is just so that it only writes this file once, although writing it twice isn't a big deal
        uga = d[cut_final][match1]
        sim = d2[match2]

        dtype=[('name','|S12'),('TS',float),('sig',float),('ra',float),('dec',float),('mod_ugali',float),('mod_simple',float),('angsep',float)]
        both = np.array([(uga['name'][i], uga['TS'][i], sim['sig'][i], uga['ra'][i], uga['dec'][i], uga['modulus'][i], sim['modulus'][i], angseps[i]) for i in range(len(angseps))], dtype=dtype)
        both = np.sort(both, order='TS')[::-1]
        pyfits.writeto('fits_files/remains_{}_both.fits'.format(args.survey), both, overwrite=True)

        justs = 'lccccccc'
        t = TexTable(len(justs), justs=justs, comments='', caption='')
        t.add_header_row(['Name', 'TS', 'SIG', r"$\alpha_{2000}$", r"$\delta_{2000}$", r"$m - M$", r"$m - M$", "Angular Separation"])
        t.add_header_row(['', '(ugali)', '(simple)', '(deg)', '(deg)', '(ugali)', '(simple)', '(deg)'])
        sigfigs = [0, 3, 3, 5, 4, 3, 3, 2]
        t.add_data([both[header] for header in both.dtype.names], sigfigs=sigfigs)
        t.print_table('tables/remains_{}_both.tex'.format(args.survey))
        subprocess.call("pdflatex -output-directory tables tables/remains_{}_both.tex".format(args.survey).split())

        # Signal
        uga = signal
        sim = signal2
        if len(uga) >= len(sim):
            s1, s2 = uga, sim
            uga1 = True
            sortkey = 'TS'
        else:
            s1, s2 = sim, uga
            uga1 = False
            sortkey = 'sig'

        combined_signal = []
        # This for loop has to be done because sometimes the signal arrays aren't the same length
        for i in range(len(s1)):
            name = s1[i]['name']
            for j in range(len(s2)):
                if s2[j]['name'] == name:
                    if uga1:
                        combined_signal.append((name, uga['TS'][i], sim['sig'][j], uga['modulus'][i], sim['modulus'][j], uga['angsep'][i], sim['angsep'][j]))
                        break
                    else:
                        combined_signal.append((name, uga['TS'][j], sim['sig'][i], uga['modulus'][j], sim['modulus'][i], uga['angsep'][j], sim['angsep'][i]))
                        break
            else:
                if uga1:
                    combined_signal.append((name, s1['TS'][i], 0.0, s1['modulus'][i], 0.0, s1['angsep'][i], 0.0))
                else:
                    combined_signal.append((name, 0.0, s1['sig'][i], 0.0, s1['modulus'][i], 0.0, s1['angsep'][i]))
                
        dtype=[('name','|S18'),('TS',float),('sig',float),('mod_ugali',float),('mod_simple',float),('angsep_ugali',float),('angsep_simple',float)]
        combined_signal = np.sort(np.array(combined_signal, dtype=dtype), order=sortkey)[::-1]
        pyfits.writeto('fits_files/signal_{}_both.fits'.format(args.survey), combined_signal, overwrite=True)
        
        justs = 'lcccccc'
        t = TexTable(len(justs), justs=justs, comments='', caption='')
        t.add_header_row(['Name', 'TS', 'SIG', r"$m - M$", r"$m - M$", "Ang. Sep.", "Ang. Sep."])
        t.add_header_row(['', '(ugali)', '(simple)', '(ugali)', '(simple)', '(deg, ugali)', '(deg, simple)'])
        sigfigs = [0, 3, 3, 5, 4, 3, 3, 2, 2]
        t.add_data([combined_signal[header] for header in combined_signal.dtype.names], sigfigs=sigfigs)
        t.print_table('tables/signal_{}_both.tex'.format(args.survey))
        subprocess.call("pdflatex -output-directory tables tables/signal_{}_both.tex".format(args.survey).split())

    if args.no_fitsfile:
        raise SystemExit(0)


# If a known satellite matches something in another catalog (you can tell from the bit), use this to find out what it's matching
def what_matches(name, known_dwarf_catalog='McConnachie15'):
    known_dwarfs = associate.catalogFactory(known_dwarf_catalog)
    obj = known_dwarfs[known_dwarfs['name'] == name]
    ra, dec = obj['ra'], obj['dec']
    external_cat_list = ['Harris96', 'Corwen04', 'Nilson73', 'Webbink85', 'Kharchenko13', 'Bica08', 'WEBDA14', 'ExtraClusters', 'McConnachie15', 'ExtraDwarfs']
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
fig = pylab.figure()
ax = fig.add_subplot(111)
if args.alg == 'simple':
    bins = np.arange(5., 40.25, 0.5)
elif args.alg == 'ugali':
    bins = np.logspace(np.log10(25.), np.log10(200000.), num=70)
    pylab.xscale('log')
pylab.yscale('log')
ax.hist(d[SIG], bins=bins, color='red', histtype='step', cumulative=-1, label='All')
ax.hist(d[SIG][cut_ebv & cut_footprint], bins=bins, color='blue', histtype='step', cumulative=-1, label='E(B-V) < 0.2 mag & in footprint')
ax.hist(d[SIG][cut_ebv & cut_footprint & cut_modulus], bins=bins, color='green', histtype='step', cumulative=-1, label='above & (m - M) < {}'.format(21.75 if args.survey == 'ps1' else 23.5))
ax.hist(d[SIG][cut_ebv & cut_footprint & cut_modulus & cut_bsc], bins=bins, color='orange', histtype='step', cumulative=-1, label='above & no bright star assoc.') 
ax.hist(d[SIG][cut_ebv & cut_footprint & cut_modulus & cut_associate & cut_bsc], bins=bins, color='black', histtype='step', cumulative=-1, label='above & no catalog assoc.') # = cut_bulk
ax.hist(d[SIG][cut_bulk & cut_dwarfs], bins=bins, color='purple', histtype='step', cumulative=-1, label='above & no known dwarf assoc.')
if not args.no_cross:
    ax.hist(d[SIG][cut_bulk & cut_dwarfs & cut_cross], bins=bins, color='darkturquoise', histtype='step', cumulative=-1, label='above & found by both algorithms') 
handles, labels = ax.get_legend_handles_labels()
new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles]
pylab.legend(loc='upper right', handles=new_handles, labels=labels, prop={'size':8})
pylab.title('Candidate dwarf galaxies ({}, {})'.format('Pan-STARRS' if args.survey == 'ps1' else 'DES', args.alg))
pylab.xlabel(SIG)
pylab.ylabel('Cumulative Count')
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
if not args.no_cross:
    healpy.projscatter(d['RA'][cut_final & ~cut_cross], d['DEC'][cut_final & ~cut_cross], lonlat=True, c='purple', marker='o', edgecolor='none', s=10, vmax=0.5)
    healpy.projscatter(d['RA'][cut_final & cut_cross], d['DEC'][cut_final & cut_cross], lonlat=True, c='darkturquoise', marker='*', edgecolor='none', s=40, vmax=0.5)
else:
    healpy.projscatter(d['RA'][cut_final], d['DEC'][cut_final], lonlat=True, c='purple', marker='o', edgecolor='none', s=10, vmax=0.5)
pylab.savefig('diagnostic_plots/significance_map_{}_{}.png'.format(args.survey, args.alg), bbox_inches='tight')

# Modulus distribution of remainers
pylab.figure()
pylab.hist(d['MODULUS'][cut_final], bins=np.arange(14.25, 25.25, 0.5), cumulative=False)
pylab.xlabel('m-M')
pylab.ylabel('Count')
pylab.title('Unassociated Hotspots')
pylab.savefig('diagnostic_plots/modulus_distribution_{}_{}.png'.format(args.survey, args.alg), bbox_inches='tight')

print "Passed cuts:", sum(cut_bulk & cut_dwarfs)
print "Passed sig:", sum(cut_sig)
print "Passed final:", sum(cut_final)
if not args.no_cross:
    print "Passed cross:", sum(cut_cross)
