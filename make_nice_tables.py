#!/usr/bin/env python

import numpy as np
from textable import TexTable

def signal_table(outname, *signals):
    #if len(signals) > 1:
    #    nonoverlap = np.array([i for i in range(len(signals[1])) if signals[1][i]['name'] not in signals[0]['name']])
    #    signal[1] = signal[1][nonoverlap]

    def reduce_signal(signal):
        cut = []
        keepers = ['Indus 1']
        losers = ['LMC', 'Sagittarius dIrr', 'Sagittarius dSph', 'Leo A', 'Leo P', 'Leo T', 'Triangulum', 'WLM', 'Pegasus dIrr']
        # Remove all losers and Andromeda satellites and objects with an integer in their name, unless they're in keepers
        for name in signal['name']:
            if name in keepers:
                cut.append(True)
            elif name in losers:
                cut.append(False)
            elif 'And' in name:
                cut.append(False)
            else:
                noint = True
                for char in name:
                    try:
                        int(char)
                    except ValueError:
                        pass
                    else:
                        noint = False
                        break
                cut.append(noint)
        cut = np.array(cut)

        return np.sort(signal[cut], order = ['TS','SIG'])[::-1] # Sorting may be redundant with how the list is orignally sorted
                
    signals = map(reduce_signal, signals)

    justs = 'lcccccccccc'
    t = TexTable(len(justs), justs=justs, comments="\\knowncomments", caption="\\knowncaption", notes="\\knownnotes", fontsize="\\tiny", doc=True)
    t.add_header_row(['Name', 'TS', 'SIG', r"$m - M$", r"$m - M$", r"$m - M$", "Distance", r"$r_h$", r"$M_V$", "Ang. Sep.", "Ang. Sep."])
    t.add_header_row(['', '(ugali)', '(simple)', '(ugali)', '(simple)', '(acutal)', '(kpc)', r"($'$)", '(mag)', r"($'$, ugali)", r"($'$, simple)"])
    sigfigs = [0, 3, 3, 3, 3, 3, 3, 2, 3, 2, 2]
    t.add_data([signals[0][header] for header in signals[0].dtype.names], sigfigs=sigfigs)
    for signal in signals[1:]:
        t.add_data(['hline']*len(justs))
        t.add_data([signal[header] for header in signal.dtype.names], sigfigs=sigfigs)
    t.print_table(outname)
    # A hack to put in \hline
    with open(outname, 'r+') as f:
        lines = f.read().splitlines()
        newlines = [(line if 'hline' not in line else '\hline') for line in lines]
        f.seek(0)
        f.write('\n'.join(newlines))
        f.truncate()


def remains_table(outname, remains, alg='simple'):
    SIG = 'TS' if alg == 'ugali' else 'SIG'

    justs = 'cccc'
    header_row1 = [SIG, r"$\alpha_{2000}$", r"$\delta_{2000}$", r"$m - M$"]
    header_row2 = ['', '(deg)', '(deg)', '']
    data_headers = [SIG, 'ra', 'dec', 'modulus']
    sigfigs = [3, 5, 4, 3]

    if alg == 'ugali':
        justs = 'l' + justs
        header_row1 = ['Name'] + header_row1
        header_row2 = [''] + header_row2
        data_headers = ['name'] + data_headers
        sigfigs = [0] + sigfigs

    if alg == 'both':
        justs = 'lc' + justs + 'cc'
        header_row1 = ['Name', 'TS'] + header_row1 + [r"$m - M$", "Angular Separation"]
        header_row2 = ['', '(ugali)', '(simple)', '(deg)', '(deg)', '(ugali)', '(simple)', r"($'$)"]
        data_headers = remains.dtype.names
        sigfigs = [0, 3] + sigfigs + [3, 2]
    
    t = TexTable(len(justs), justs=justs, comments="\\candidatecomments", caption="\\candidatecaption", notes="\\knownnotes", fontsize="\\tiny", doc=True)
    t.add_header_row(header_row1)
    t.add_header_row(header_row2)
    t.add_data([remains[header] for header in data_headers], sigfigs=sigfigs)
    t.print_table(outname)
