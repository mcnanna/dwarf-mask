#!/usr/bin/env python
import subprocess
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
plt.ion()
from matplotlib.lines import Line2D

from filter_candidates import Candidates
import make_nice_tables

# Define significances
conservativesigdict = {
        ('ps1','simple'):7.,
        ('ps1','ugali'):7**2,
        ('des','simple'):7.,
        ('des','ugali'):8**2
        }

fiducialsigdict = {
        ('ps1','simple'):6.,
        ('ps1','ugali'):6**2,
        ('des','simple'):6,
        ('des','ugali'):6**2
        }

# Filter candidates for each survey/alg combination
sigdict = fiducialsigdict
des_ugali = Candidates('des', 'ugali', sigdict=sigdict)
des_simple = Candidates('des', 'simple', sigdict=sigdict)
ps1_ugali = Candidates('ps1', 'ugali', sigdict=sigdict)
ps1_simple = Candidates('ps1', 'simple', sigdict=sigdict)

# Make fits files, tables, and plots for each
des_ugali.doitall()
des_simple.doitall()
ps1_ugali.doitall()
ps1_simple.doitall()

# Combine signals into one big table
signal_des = pyfits.open("fits_files/signal_des_both.fits")[1].data
signal_ps1 = pyfits.open("fits_files/signal_ps1_both.fits")[1].data
make_nice_tables.signal_table("tables/signal.tex", signal_des, signal_ps1)
subprocess.call("pdflatex -interaction nonstopmode -output-directory tables tables/signal.tex".split())

# Combine remains into one big table
remains_des = pyfits.open("fits_files/remains_des_both.fits")[1].data
remains_ps1 = pyfits.open("fits_files/remains_ps1_both.fits")[1].data
make_nice_tables.remains_table("tables/remains.tex", remains_des, remains_ps1, alg='both')
subprocess.call("pdflatex -interaction nonstopmode -output-directory tables tables/remains.tex".split())

# Combine sighists into one plot
def sighist(cands, ax, legend=True, title=True, text=True, xlabel=True, ylabel=True):
    fiducial_threshold = fiducialsigdict[(cands.survey, cands.alg)]
    conservative_threshold = conservativesigdict[(cands.survey, cands.alg)]

    if cands.alg == 'simple':
        sigs = cands.data['SIG']
        bins = np.arange(0., 40., 0.5) # step=0.5, sigs max out at 37.5
        axis_label = 'SIG'
    elif cands.alg == 'ugali':
        sigs = np.sqrt(cands.data['TS'])
        bins = np.arange(3.0, 481, 0.5) # Lowest sqrt(TS) is sqrt(10) = 3.2
        axis_label = r'$\sqrt{\mathrm{TS}}$'
        ax.set_xlim(left=None, right=51.0)

        fiducial_threshold = np.sqrt(fiducial_threshold)
        conservative_threshold = np.sqrt(conservative_threshold)

    ax.set_ylim(0.5, 2.5*10**5)
    ax.set_yscale('log')
    ax.hist(sigs, bins=bins, color='red', histtype='step', cumulative=-1, label='All')
    ax.hist(sigs[cands.cut_ebv & cands.cut_footprint], bins=bins, color='blue', histtype='step', cumulative=-1, label= r'In footprint & $E(B-V) < 0.2$ mag')
    #ax.hist(sigs[cands.cut_ebv & cands.cut_footprint & cands.cut_modulus], bins=bins, color='darkturquoise', histtype='step', cumulative=-1, label=r'above & $m - M < {}$'.format(21.75 if cands.survey == 'ps1' else 23.5))
    #ax.hist(sigs[cands.cut_ebv & cands.cut_footprint & cands.cut_modulus & cands.cut_bsc], bins=bins, color='orange', histtype='step', cumulative=-1, label='above & no bright star assoc.') 
    ax.hist(sigs[cands.cut_ebv & cands.cut_footprint & cands.cut_modulus & cands.cut_associate & cands.cut_bsc], bins=bins, color='green', histtype='step', cumulative=-1, label='above & no catalog assoc.') # = cands.cut_bulk
    ax.hist(sigs[cands.cut_bulk & cands.cut_dwarfs], bins=bins, color='purple', histtype='step', cumulative=-1, label='above & no known dwarf assoc.')
    if cands.cross:
        ax.hist(sigs[cands.cut_bulk & cands.cut_dwarfs & cands.cut_cross], bins=bins, color='black', histtype='step', cumulative=-1, label='above & found by both algorithms') 
    if legend:
        handles, labels = ax.get_legend_handles_labels()
        new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles]
        ax.legend(loc='upper right', handles=new_handles, labels=labels, prop={'size':8})
    if title:
        ax.set_title('SimpleBinner' if cands.alg == 'simple' else 'ugali', fontdict={'family': 'monospace'})
    else:
        ax.set_title(' ')
    if text:
        ax.set_ylabel(' '*8+ cands.survey.upper(), rotation=0, horizontalalignment='right') # Without extra spaces, it overlaps the subplot outline
        ax.yaxis.set_label_position("right")
    if xlabel:
        ax.set_xlabel(axis_label)
    if ylabel:
        ax.set_ylabel('Cumulative Count')

    ax.axvline(fiducial_threshold, color='gray', linestyle='-', linewidth=1.0)
    ax.axvline(conservative_threshold, color='gray', linestyle='--', linewidth=1.0)
    
fig, axes = plt.subplots(nrows=2, ncols=2, sharex='col', sharey='all', figsize=(9, 7.5))
sighist(des_simple, axes[0,0], title=True, text=False, legend=False, xlabel=False, ylabel=True)
sighist(des_ugali, axes[0,1], title=True, text=True, legend=True, xlabel=False, ylabel=False)
sighist(ps1_simple, axes[1,0], title=False, text=False, legend=False, xlabel=True, ylabel=True)
sighist(ps1_ugali, axes[1,1], title=False, text=True, legend=False, xlabel=True, ylabel=False)
plt.subplots_adjust(wspace=0.1, hspace=0.1)
plt.savefig('diagnostic_plots/significance_distribution_all.png', bbox_inches='tight')

