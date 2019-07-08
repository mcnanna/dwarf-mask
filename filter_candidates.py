#!/usr/bin/env python

import subprocess
import numpy as np
import healpy as hp
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
plt.ion()
from matplotlib.lines import Line2D
import argparse

import ugali.utils.projector
import ugali.utils.healpix
import associate
import make_nice_tables

version = 1.0
fiducialsigdict = {('ps1', 'simple'):6, ('ps1', 'ugali'):6**2, ('des', 'simple'):6, ('des', 'ugali'):6**2}

def custom_sort(signal, order=None):
    # Put highest to lowest
    sorted_ = np.sort(signal, order=order)[::-1] 

    # Put np.nan on the bottom
    if order is not None:
        if type(order) != list:
            order = [order]

        to_bottom = []
        for i in range(len(sorted_)):
            vals = []
            for key in order:
                vals.append(sorted_[i][key])

            if True in map(np.isnan, vals):
                # Might be high sig detection
                highsig = False
                for val in vals:
                    if not np.isnan(val) and val > 10:
                        highsig = True

                if not highsig:
                    to_bottom.append(i)
    
        sorted_ = sorted_[[j for j in range(len(sorted_)) if j not in to_bottom] + to_bottom]

    return sorted_


class Candidates:
    def __init__(self, survey, alg, infile=None, cross=True, threshold=None, sigdict=fiducialsigdict, quiet=False):
        if 'd' in survey.lower():
            self.survey = 'des'
        elif 'p' in survey.lower():
            self.survey = 'ps1'
        else:
            raise ValueError("survey must be 'des' or 'ps1'")

        if 's' in alg.lower():
            self.alg = 'simple'
        elif 'u' in alg.lower():
            self.alg = 'ugali'
        else:
            raise ValueError("alg must be 'ugali' or 'simple'")
        self.cross = cross
        if threshold is not None:
            self.threshold = threshold
        else:
            self.threshold = sigdict[(self.survey, self.alg)]


        if infile is None:
            infile = 'candidates/candidate_list_{0}_{1}.fits'.format(self.survey, self.alg)
        r = pyfits.open(infile)
        d_original = r[1].data
        r.close()

        # Minimum significance threshold (not used right now)
        self.data = d_original#[d_original['SIG'] > 5.] [d_original['TS'] > 25?]

        if self.alg == 'simple':
            self.SIG = 'SIG'
        elif self.alg == 'ugali':
            self.SIG = 'TS'

        ### Define and apply cuts

        # Consolidate nearby peaks iteratively
        while True:
            match_1, match_2, angsep = ugali.utils.projector.match(self.data['ra'], self.data['dec'], self.data['ra'], self.data['dec'], tol=0.2, nnearest=2)
            #index_exclude = np.where(self.data[self.SIG][match_1] > self.data[self.SIG][match_2], match_2, match_1)
            index_exclude = []
            for j in range(len(angsep)):
                i1 = match_1[j]
                i2 = match_2[j]
                sig1 = self.data[self.SIG][i1]
                sig2 = self.data[self.SIG][i2]
                if sig1 > sig2:
                    index_exclude.append(i2)
                elif sig2 > sig1:
                    index_exclude.append(i1)
                elif sig1 == sig2:
                    index_exclude.append(max(i1, i2)) # This is arbitrary
            index_exclude = np.array(index_exclude)
            if len(index_exclude) == 0:
                break
            cut_consolidate = np.tile(True, len(self.data))
            cut_consolidate[index_exclude] = False
            self.data = self.data[cut_consolidate]

        # Geometric cuts
        pix = ugali.utils.healpix.angToPix(4096, self.data['ra'], self.data['dec'], nest=True)
        mask = ugali.utils.healpix.read_map('healpix_mask_{}_v{}.fits.gz'.format(self.survey, version), nest=True)

        self.cut_ebv = np.where(mask[pix] & 0b00001, False, True)
        self.cut_associate = np.where(mask[pix] & 0b00010, False, True)
        self.cut_dwarfs = np.where(mask[pix] & 0b00100, False, True)
        self.cut_bsc = np.where(mask[pix] & 0b01000, False, True)
        self.cut_footprint = np.where(mask[pix] & 0b10000, False, True)
        if self.survey == 'ps1' and self.alg == 'ugali':
            self.cut_footprint &= np.where(mask[pix] & 0b100000, False, True) # Add failures to footprint cut 
        #self.cut_footprint &= np.where(mask[pix] & 0b1000000, False, True) # Add artifacts to footprint cut

        # Other cuts (modulus, size, shape)
        if self.survey == 'ps1':
            self.cut_modulus = (self.data['MODULUS'] < 21.75)
        elif self.survey == 'des':
            self.cut_modulus = (self.data['MODULUS'] < 23.5)

        # Significance cut
        self.cut_sig = self.data[self.SIG] > self.threshold

        # Combine cuts which should filter everything but dwarfs
        self.cut_bulk = self.cut_ebv & self.cut_footprint & self.cut_modulus & self.cut_associate & self.cut_bsc

        # Final cut also filters known satellites and low sig detections
        self.cut_final = self.cut_bulk & self.cut_dwarfs & self.cut_sig


        ### Signal Detection

        lvdb = pyfits.open('lvdb_v2.fits')[1].data
        pdet_dict = np.load('pdet_{}.npy'.format(self.survey.upper())).item()
        signal = []
        for known_dwarf_catalog in ['McConnachie15', 'ExtraDwarfs']:
            known_dwarfs = associate.catalogFactory(known_dwarf_catalog)
            for known_dwarf in known_dwarfs:
                tol = known_dwarf['radius']
                if np.isnan(tol) or tol < 0.2: # Minimum matching radius of 0.2
                    tol = 0.2
                match_candidate, match_known_dwarf, angseps = ugali.utils.projector.match(self.data['RA'], self.data['DEC'], [known_dwarf['ra']], [known_dwarf['dec']], tol=tol)

                name = known_dwarf['name']
                ra = known_dwarf['ra']
                dec = known_dwarf['dec']
                if len(match_known_dwarf) > 0:
                    # keep highest significance match (Usually only 1 match. Exceptions: LMC, Sagittarius dSph, Bootes III, Crater II (simple only))
                    mx = np.argmax(self.data[match_candidate][self.SIG])
                    idx = match_candidate[mx]
                    sig = self.data[self.SIG][idx]
                    a = angseps[0]*60.
                    modulus = self.data['MODULUS'][idx]
                    bit = mask[pix[idx]]
                    wascut = '' if self.cut_bulk[idx] else 'cut'
                else:
                    sig = np.nan
                    a = np.nan
                    modulus = np.nan
                    #modulus_actual = np.nan
                    #distance = np.nan
                    #rhalf = np.nan
                    #M_V = np.nan
                    bit = mask[ugali.utils.healpix.angToPix(4096, known_dwarf['ra'], known_dwarf['dec'], nest=True)]
                    wascut = ''

                # Look for object in Local Volume Database
                try:
                    lv = lvdb[lvdb['name'] == name][0]
                except IndexError:
                    modulus_actual = np.nan
                    distance = np.nan
                    rhalf_obs = np.nan
                    rhalf_pc = np.nan
                    M_V = np.nan
                    ref = r"\ldots"
                else:
                    modulus_actual = lv['distance_modulus']
                    distance = lv['distance_kpc']
                    rhalf_obs = lv['rhalf']
                    rhalf_pc = (distance*1000.) * np.radians(rhalf_obs/60.)
                    M_V = lv['m_v']
                    ref = lv['structure_ref']
                    if ref == 'None':
                        ref = lv['distance_ref']

                # Get p_det from dictionary (make sure to update if the sigs ever chnge)
                try:
                    pdet = pdet_dict[name][0] # pdet_dict[name] is an array with a single element, hence the [0]
                except KeyError:
                    pdet = np.nan

                if (not (bit & 0b10000)) or (not np.isnan(sig)): # Don't bother writing results for non-detections outside of the footprint
                    if (len(signal) == 0) or (name not in np.array(signal)[:, 0]): # Try to avoid duplicates from the multiple catalogs
                        signal.append((name, sig, pdet, ra, dec, modulus, modulus_actual, distance, rhalf_obs, rhalf_pc, M_V, a, wascut, bit, ref))

        dtype = [('name','|S18'),(self.SIG, float),('pdet', float),
                ('ra',float),('dec',float),
                ('modulus',float),('modulus_actual',float),('distance',float),('rhalf_obs',float),('rhalf_pc',float),('M_V',float),
                ('angsep',float),('cut','|S3'),('bit',int),
                ('ref','|S24')]
        self.signal = custom_sort(np.array(signal, dtype=dtype), order=self.SIG)


        
        ### Combine results of two algorithsm

        if self.cross:
            # Make sure you pass a sigdict if you want to compare significances other than the defaults in fiducialsigdict
            if self.alg == 'simple':
                cands2 = Candidates(self.survey, 'ugali', cross=False, sigdict=sigdict, quiet=True)
            if self.alg == 'ugali':
                cands2 = Candidates(self.survey, 'simple', cross=False, sigdict=sigdict, quiet=True)
            data2 = cands2.data[cands2.cut_final]
            self._signal2 = cands2.signal

            match1, match2, angseps = ugali.utils.projector.match(self.data['ra'][self.cut_final], self.data['dec'][self.cut_final], data2['ra'], data2['dec'], tol=0.2)
            self.matches = self.data[self.cut_final][match1]
            self._matches2 = data2[match2]
            self._angseps = angseps
            self.cut_cross = np.array([(self.data[i]['ra'], self.data[i]['dec']) in zip(self.matches['ra'], self.matches['dec']) for i in range(len(self.data))])

        if not quiet:    
            print
            print self.survey, self.alg
            print "Passed cuts:", sum(self.cut_bulk & self.cut_dwarfs)
            print "Passed sig:", sum(self.cut_sig)
            print "Passed final:", sum(self.cut_final)
            if self.cross:
                print "Passed cross:", sum(self.cut_cross)
            print


    # Write sorted remains to .fits and make a TeX table
    def write_remains(self, table=True):
        subprocess.call('mkdir -p fits_files'.split())
        pyfits.writeto('fits_files/remains_{}_{}.fits'.format(self.survey, self.alg), np.sort(self.data[self.cut_final], order=self.SIG)[::-1], overwrite=True)
        if table:
            subprocess.call('mkdir -p tables'.split())
            make_nice_tables.remains_table('tables/remains_{}_{}.tex'.format(self.survey, self.alg), np.sort(self.data[self.cut_final], order=self.SIG)[::-1], alg=self.alg)
            subprocess.call("pdflatex -interaction nonstopmode -output-directory tables tables/remains_{}_{}.tex".format(self.survey, self.alg).split())

    def write_cross_remains(self, table=True):
        subprocess.call('mkdir -p fits_files'.split())
        if not self.cross:
            raise Exception("No cross data, make sure you call with cross=True")

        if self.alg == 'ugali':
            uga = self.matches
            sim = self._matches2
        elif self.alg == 'simple':
            sim = self.matches
            uga = self._matches2
        angseps = self._angseps

        dtype=[('name','|S12'),('TS',float),('SIG',float),('ra',float),('dec',float),('mod_ugali',float),('mod_simple',float),('angsep',float)]
        both = np.array([(uga['name'][i], uga['TS'][i], sim['SIG'][i], uga['ra'][i], uga['dec'][i], uga['modulus'][i], sim['modulus'][i], angseps[i]*60.) for i in range(len(angseps))], dtype=dtype)
        both = np.sort(both, order=['TS','SIG'])[::-1]
        pyfits.writeto('fits_files/remains_{}_both.fits'.format(self.survey), both, overwrite=True)
        if table:
            subprocess.call('mkdir -p tables'.split())
            make_nice_tables.remains_table('tables/remains_{}_both.tex'.format(self.survey), both, alg='both')
            subprocess.call("pdflatex -interaction nonstopmode -output-directory tables tables/remains_{}_both.tex".format(self.survey).split())


    # Write signal to .fits file and make a Tex table
    def write_signal(self):
        pyfits.writeto('fits_files/signal_{}_{}.fits'.format(self.survey, self.alg), self.signal, overwrite=True)

    def write_cross_signal(self, table=True):
        if not self.cross:
            raise Exception("No cross data, make sure you call with cross=True")
        
        if self.alg == 'ugali':
            uga = self.signal
            sim = self._signal2
        elif self.alg == 'simple':
            sim = self.signal
            uga = self._signal2

        combined_signal = []
        for i in range(len(uga)):
            name = uga[i]['name']
            for j in range(len(sim)):
                if sim[j]['name'] == name:
                    combined_signal.append(
                            (name, uga['TS'][i], sim['SIG'][j], uga['pdet'][i],
                                uga['ra'][i], uga['dec'][i],
                                uga['modulus'][i], sim['modulus'][j],
                                uga['modulus_actual'][i], uga['rhalf_obs'][i], 
                                uga['distance'][i], uga['rhalf_pc'][i], uga['M_V'][i], 
                                uga['angsep'][i], sim['angsep'][j],
                                uga['ref'][i]))
                    sim = np.delete(sim, j)
                    break
            else:
                combined_signal.append(
                        (name, uga['TS'][i], np.nan, uga['pdet'][i],
                            uga['ra'][i], uga['dec'][i], 
                            uga['modulus'][i], np.nan,
                            uga['modulus_actual'][i], uga['rhalf_obs'][i], 
                            uga['distance'][i], uga['rhalf_pc'][i], uga['M_V'][i], 
                            uga['angsep'][i], np.nan,
                            uga['ref'][i]))

        for j in range(len(sim)):
                combined_signal.append(
                        (sim['name'][j], np.nan, sim['SIG'][j], sim['pdet'][j],
                            sim['ra'][j], sim['dec'][j],
                            np.nan, sim['modulus'][j],
                            sim['modulus_actual'][j], sim['rhalf_obs'][j], 
                            sim['distance'][j], sim['rhalf_pc'][j], sim['M_V'][j],
                            np.nan, sim['angsep'][j],
                            sim['ref'][j]))
        
        dtype=[('name','|S18'),('TS',float),('SIG',float),('pdet',float),
                ('ra',float),('dec',float),
                ('mod_ugali',float),('mod_simple',float),
                ('mod_actual',float),('rhalf_obs',float),
                ('distance',float),('rhalf_pc',float),('M_V',float),
                ('angsep_ugali',float),('angsep_simple',float),
                ('ref','|S24')]
        combined_signal = custom_sort(np.array(combined_signal, dtype=dtype), order=['TS', 'SIG'])
        pyfits.writeto('fits_files/signal_{}_both.fits'.format(self.survey), combined_signal, overwrite=True)
        if table:
            make_nice_tables.signal_table('tables/signal_{}.tex'.format(self.survey), combined_signal)
            subprocess.call("pdflatex -interaction nonstopmode -output-directory tables tables/signal_{}.tex".format(self.survey).split())


    ### Diagnostic plots
    def sighist(self, legend=True, title=True, text=True):
        fig = plt.figure(figsize=(6.0,5.0))
        ax = fig.add_subplot(111)
        if self.alg == 'simple':
            bins = np.arange(0., 40.0, 0.5) # n=81 bins
        elif self.alg == 'ugali':
            bins = np.logspace(np.log10(10.), np.log10(230500.), num=81)
            plt.xscale('log')
        plt.yscale('log')
        ax.hist(self.data[self.SIG], bins=bins, color='red', histtype='step', cumulative=-1, label='All')
        ax.hist(self.data[self.SIG][self.cut_ebv & self.cut_footprint], bins=bins, color='blue', histtype='step', cumulative=-1, label= r'$E(B-V) < 0.2$ mag & in footprint')
        ax.hist(self.data[self.SIG][self.cut_ebv & self.cut_footprint & self.cut_modulus], bins=bins, color='green', histtype='step', cumulative=-1, label=r'above & $m - M < {}$'.format(21.75 if self.survey == 'ps1' else 23.5))
        ax.hist(self.data[self.SIG][self.cut_ebv & self.cut_footprint & self.cut_modulus & self.cut_bsc], bins=bins, color='orange', histtype='step', cumulative=-1, label='above & no bright star assoc.') 
        ax.hist(self.data[self.SIG][self.cut_ebv & self.cut_footprint & self.cut_modulus & self.cut_associate & self.cut_bsc], bins=bins, color='black', histtype='step', cumulative=-1, label='above & no catalog assoc.') # = self.cut_bulk
        ax.hist(self.data[self.SIG][self.cut_bulk & self.cut_dwarfs], bins=bins, color='purple', histtype='step', cumulative=-1, label='above & no known dwarf assoc.')
        if self.cross:
            ax.hist(self.data[self.SIG][self.cut_bulk & self.cut_dwarfs & self.cut_cross], bins=bins, color='darkturquoise', histtype='step', cumulative=-1, label='above & found by both algorithms') 
        if legend:
            handles, labels = ax.get_legend_handles_labels()
            new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles]
            plt.legend(loc='upper right', handles=new_handles, labels=labels, prop={'size':8})
        if title:
            plt.title('SimpleBinner' if self.alg == 'simple' else 'Ugali', fontdict={'family': 'monospace'})
        else:
            plt.title(' ')
        #plt.title('{}'.format('Pan-STARRS' if self.survey == 'ps1' else 'DES', 'SimpleBinner' if self.alg == 'simple' else 'Ugali'))
        if text:
            plt.text(0.9, 0.55, self.survey.upper(), transform=ax.transAxes, horizontalalignment='left', fontsize=12, weight='bold')
        plt.xlabel(self.SIG)
        plt.ylabel('Cumulative Count')

        ax.axvline(self.threshold, color='gray', linestyle='-', linewidth=1.0)

        subprocess.call('mkdir -p diagnostic_plots'.split())
        plt.savefig('diagnostic_plots/significance_distribution_{}_{}.png'.format(self.survey, self.alg), bbox_inches='tight')

    
    def sigmap(self):
        infile_dust = 'ebv_sfd98_fullres_nside_4096_nest_equatorial.fits.gz'
        ebv_map = ugali.utils.healpix.read_map(infile_dust, nest=True)

        #plt.figure()
        hp.mollview(ebv_map, xsize=1600, min=0., max=0.5, cmap='binary', title='{0} > {1} Hotspots'.format(self.SIG, self.threshold), unit='E(B-V)', nest=True)
        hp.graticule()
        hp.projscatter(self.data['RA'][self.cut_sig], self.data['DEC'][self.cut_sig], lonlat=True, c='red', marker='o', edgecolor='none', s=2, vmax=0.5)
        hp.projscatter(self.data['RA'][self.cut_sig & self.cut_ebv & self.cut_footprint], self.data['DEC'][self.cut_sig & self.cut_ebv & self.cut_footprint], lonlat=True, c='blue', marker='o', edgecolor='none', s=2, vmax=0.5)
        hp.projscatter(self.data['RA'][self.cut_sig & self.cut_ebv & self.cut_footprint & self.cut_modulus], self.data['DEC'][self.cut_sig & self.cut_ebv & self.cut_footprint & self.cut_modulus], lonlat=True, c='green', marker='o', edgecolor='none', s=2, vmax=0.5)
        if self.cross:
            hp.projscatter(self.data['RA'][self.cut_final & ~self.cut_cross], self.data['DEC'][self.cut_final & ~self.cut_cross], lonlat=True, c='purple', marker='o', edgecolor='none', s=10, vmax=0.5)
            hp.projscatter(self.data['RA'][self.cut_final & self.cut_cross], self.data['DEC'][self.cut_final & self.cut_cross], lonlat=True, c='darkturquoise', marker='*', edgecolor='none', s=40, vmax=0.5)
        else:
            hp.projscatter(self.data['RA'][self.cut_final], self.data['DEC'][self.cut_final], lonlat=True, c='purple', marker='o', edgecolor='none', s=10, vmax=0.5)
        subprocess.call('mkdir -p diagnostic_plots'.split())
        plt.savefig('diagnostic_plots/significance_map_{}_{}.png'.format(self.survey, self.alg), bbox_inches='tight')

    def modhist(self):
        plt.figure()
        plt.hist(self.data['MODULUS'][self.cut_final], bins=np.arange(14.25, 25.25, 0.5), cumulative=False)
        plt.xlabel('m-M')
        plt.ylabel('Count')
        plt.title('Unassociated Hotspots')
        subprocess.call('mkdir -p diagnostic_plots'.split())
        plt.savefig('diagnostic_plots/modulus_distribution_{}_{}.png'.format(self.survey, self.alg), bbox_inches='tight')


    def doitall(self, table=True, legend=True, title=True, text=True):
        self.write_remains(table)
        self.write_signal()
        if self.cross:
            self.write_cross_remains(table)
            self.write_cross_signal(table)
        self.sighist(legend, title, text)
        self.sigmap()
        self.modhist()


# If a known satellite matches something in another catalog (you can tell from the bit), use this to find out what it's matching
def what_matches(name, known_dwarf_catalog='McConnachie15', tol=6.0):
    known_dwarfs = associate.catalogFactory(known_dwarf_catalog)
    obj = known_dwarfs[known_dwarfs['name'] == name]
    ra, dec = obj['ra'], obj['dec']
    external_cat_list = ['Harris96', 'Corwen04', 'Nilson73', 'Webbink85', 'Kharchenko13', 'Bica08', 'WEBDA14', 'ExtraClusters', 'McConnachie15', 'ExtraDwarfs', 'ExtraStructures']
    print '%15s%20s%10s%10s'%('catalog', 'name', 'radius', 'angsep')
    for cat in external_cat_list:
        catalog = associate.catalogFactory(cat)
        whocares, matches, angseps = catalog.match(ra, dec, coord='cel', tol=tol)
        for i in range(len(matches)):
            index = matches[i]
            angsep = angseps[i]*60.
            radius = catalog[index]['radius']
            print '%15s%20s%10.3f%10.3f'%(cat, catalog[index]['name'], radius, angsep) 


if __name__ == '__main__':
    TABLES = True
    if TABLES:
        subprocess.call('mkdir -p tables'.split())
    subprocess.call('mkdir -p fits_files diagnostic_plots'.split())

    p = argparse.ArgumentParser()
    p.add_argument('-s', '--survey', help="'ps1' or 'des'")
    p.add_argument('-a', '--alg', help = "'u'/'ugali' or 's'/'simple'")
    p.add_argument('--sig', type=float, help="Minimum significance/TS threshold.")
    p.add_argument('--no_cross', action='store_true', help="Don't cross-check between the two algorithms. Use this is you only have the candidates for one algorithm")
    args = p.parse_args()

    if (args.survey is not None) and (args.alg is not None):
        cands = Candidates(args.survey, args.alg, cross = (not args.no_cross), threshold = args.sig)
        cands.doitall(table=TABLES)
    else:
        print 'nada'

