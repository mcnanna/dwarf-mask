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
highsigdict = {
        ('ps1','simple'):6.,
        ('ps1','ugali'):80.,
        ('des','simple'):7.,
        ('des','ugali'):50.
        }

lowsigdict = {
        ('ps1','simple'):5.,
        ('ps1','ugali'):50.,
        ('des','simple'):7.,
        ('des','ugali'):30.
        }

# Filter candidates for each survey/alg combination
sigdict = highsigdict
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


# Combine sighists into one plot
# This is pretty specialized to get it into a good format for the paper
def sighist(cands, ax, legend=True, title=True, text=True, xlabel=True, ylabel=True):
    if cands.alg == 'simple':
        bins = np.arange(5., 40.25, 0.5)
    elif cands.alg == 'ugali':
        bins = np.logspace(np.log10(10.), np.log10(200000.), num=70)
        ax.set_xscale('log')
    ax.set_ylim(0.7, 4*10**4)
    ax.set_yscale('log')
    ax.hist(cands.data[cands.SIG], bins=bins, color='red', histtype='step', cumulative=-1, label='All')
    ax.hist(cands.data[cands.SIG][cands.cut_ebv & cands.cut_footprint], bins=bins, color='blue', histtype='step', cumulative=-1, label= r'In footprint & $E(B-V) < 0.2$ mag')
    #ax.hist(cands.data[cands.SIG][cands.cut_ebv & cands.cut_footprint & cands.cut_modulus], bins=bins, color='darkturquoise', histtype='step', cumulative=-1, label=r'above & $m - M < {}$'.format(21.75 if cands.survey == 'ps1' else 23.5))
    #ax.hist(cands.data[cands.SIG][cands.cut_ebv & cands.cut_footprint & cands.cut_modulus & cands.cut_bsc], bins=bins, color='orange', histtype='step', cumulative=-1, label='above & no bright star assoc.') 
    ax.hist(cands.data[cands.SIG][cands.cut_ebv & cands.cut_footprint & cands.cut_modulus & cands.cut_associate & cands.cut_bsc], bins=bins, color='green', histtype='step', cumulative=-1, label='above & no catalog assoc.') # = cands.cut_bulk
    ax.hist(cands.data[cands.SIG][cands.cut_bulk & cands.cut_dwarfs], bins=bins, color='purple', histtype='step', cumulative=-1, label='above & no known dwarf assoc.')
    if cands.cross:
        ax.hist(cands.data[cands.SIG][cands.cut_bulk & cands.cut_dwarfs & cands.cut_cross], bins=bins, color='black', histtype='step', cumulative=-1, label='above & found by both algorithms') 
    if legend:
        handles, labels = ax.get_legend_handles_labels()
        new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles]
        ax.legend(loc='upper right', handles=new_handles, labels=labels, prop={'size':8})
    if title:
        ax.set_title('SimpleBinner' if cands.alg == 'simple' else 'Ugali', fontdict={'family': 'monospace'})
    else:
        ax.set_title(' ')
    if text:
        ax.set_ylabel(' '*8+ cands.survey.upper(), rotation=0, horizontalalignment='right') # Without extra spaces, it overlaps the subplot outline
        ax.yaxis.set_label_position("right")
    if xlabel:
        ax.set_xlabel(cands.SIG)
    if ylabel:
        ax.set_ylabel('Cumulative Count')

fig, axes = plt.subplots(nrows=2, ncols=2, sharex='col', sharey='all', figsize=(9, 7.5))
sighist(des_simple, axes[0,0], title=True, text=False, legend=False, xlabel=False, ylabel=True)
sighist(des_ugali, axes[0,1], title=True, text=True, legend=True, xlabel=False, ylabel=False)
sighist(ps1_simple, axes[1,0], title=False, text=False, legend=False, xlabel=True, ylabel=True)
sighist(ps1_ugali, axes[1,1], title=False, text=True, legend=False, xlabel=True, ylabel=False)
plt.subplots_adjust(wspace=0.1, hspace=0.1)
plt.savefig('histograms.png', bbox_inches='tight')
