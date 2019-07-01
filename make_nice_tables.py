#!/usr/bin/env python

import numpy as np
from textable import TexTable

def signal_table(outname, *signals):
    #if len(signals) > 1:
    #    nonoverlap = np.array([i for i in range(len(signals[1])) if signals[1][i]['name'] not in signals[0]['name']])
    #    signal[1] = signal[1][nonoverlap]

    def reduce_signal(signal):
        if ('TS' in signal.dtype.names) and ('SIG' in signal.dtype.names): # If it has both sig types, sqrt(TS) for comparison
            signal['TS'] = np.sqrt(signal['TS'])

        cut = []
        keepers = ['Segue 1', 'Segue 2', 'Willman 1']
        losers = ['LMC', 'Sagittarius dIrr', 'Sagittarius dSph', 'Leo A', 'Leo P', 'Leo T', 'Triangulum', 'WLM', 'Pegasus dIrr', 'Tucana', 'Cetus', 'Sextans B', 'Phoenix', 'Sextans A', 'Antlia', 'Aquarius']
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

        return signal[cut]
                
    signals = map(reduce_signal, signals)

    ref_dict = {}
    i = 1
    for signal in signals:
        for ref in signal['ref']:
            if ref in ['None', r"\ldots"]:
                continue

            elif ref not in ref_dict:
                ref_dict[ref] = i
                i += 1

    for signal in signals:
        for line in signal:
            ref = line['ref']
            if ref in ref_dict:
                line['ref'] = ref_dict[line['ref']]

    ref_dict_str = ("({}) {}, "*len(ref_dict)).format(*np.array([(ref_dict[key], key) for key in ref_dict.keys()]).flatten())

    #justs = 'lcccccccccl'
    justs = 'lccccccccl'
    t = TexTable(len(justs), justs=justs, comments="\\knowncomments" + " References: " + ref_dict_str, caption="\\knowncaption", notes="\\knownnotes", fontsize="\\scriptsize", doc=True)
    #t.add_header_row(['({})'.format(i+1) for i in range(len(justs))])
    #t.add_header_row(['Name', r'$\sqrt{\mathrm{TS}}$', 'SIG', r"$P_{\rm det}$", "RA", "Dec", "$m - M$", "Distance", r"$r_h$", r"$M_V$", "Ref."])
    #t.add_header_row(['', '', '', '', '(deg)', '(deg)', '', '(kpc)', r"($'$)", '(mag)', ''])
    #roundings = [0, 1, 1, 2, 2, 2, 1, 0, 1, 1, 0]
    t.add_header_row(['({})'.format(i+1) for i in range(len(justs))])
    t.add_header_row(['Name', r'$\sqrt{\mathrm{TS}}$', 'SIG', "RA", "Dec", "$m - M$", "Distance", r"$r_h$", r"$M_V$", "Ref."])
    t.add_header_row(['', '', '', '(deg)', '(deg)', '', '(kpc)', r"($'$)", '(mag)', ''])
    roundings = [0, 1, 1, 2, 2, 1, 0, 1, 1, 0]
    ignore = ['mod_ugali', 'mod_simple', 'angsep_ugali', 'angsep_simple']
    t.add_data([signals[0][header] for header in signals[0].dtype.names if header not in ignore], sigfigs=roundings)
    for signal in signals[1:]:
        t.add_data(['hline']*len(justs))
        t.add_data([signal[header] for header in signal.dtype.names if header not in ignore], sigfigs=roundings)
    t.print_table(outname)
    # A hack to put in \hline
    with open(outname, 'r+') as f:
        lines = f.read().splitlines()
        newlines = [(line if 'hline' not in line else '\hline') for line in lines]
        f.seek(0)
        f.write('\n'.join(newlines))
        f.truncate()


def remains_table(outname, *remains, **alg):
    alg = alg.pop('alg', 'simple') # Default alg is 'simple'

    SIG = 'TS' if alg == 'ugali' else 'SIG'

    justs = 'cccc'
    header_row1 = [SIG, r"$\alpha_{2000}$", r"$\delta_{2000}$", r"$m - M$"]
    header_row2 = ['', '(deg)', '(deg)', '']
    data_headers = [SIG, 'ra', 'dec', 'modulus']
    #sigfigs = [3, 5, 4, 3]
    roundings = [1, 2, 2, 1]

    if alg == 'ugali':
        justs = 'l' + justs
        header_row1 = ['Name'] + header_row1
        header_row2 = [''] + header_row2
        data_headers = ['Name'] + data_headers
        #sigfigs = [0] + sigfigs
        roundings = [0] + roundings

    if alg == 'both':
        for remain in remains:
            if 'TS' in remain.dtype.names:
                remain['TS'] = np.sqrt(remain['TS'])

        justs = 'lc' + justs + 'cc'
        header_row1 = ['Name', r'$\sqrt{\mathrm{TS}}$'] + header_row1 + [r"$m - M$", "Angular Separation"]
        header_row2 = ['', '(ugali)', '(simple)', '(deg)', '(deg)', '(ugali)', '(simple)', r"($'$)"]
        data_headers = remains[0].dtype.names
        #sigfigs = [0, 3] + sigfigs + [3, 2]
        roundings = [0, 1] + roundings + [1, 2]
    
    t = TexTable(len(justs), justs=justs, comments="\\candidatecomments", caption="\\candidatecaption", notes="\\knownnotes", fontsize="\\tiny", doc=True)
    t.add_header_row(header_row1)
    t.add_header_row(header_row2)
    t.add_data([remains[0][header] for header in data_headers], sigfigs=roundings)
    for remain in remains[1:]:
        t.add_data(['hline']*len(justs))
        t.add_data([remain[header] for header in data_headers], sigfigs=roundings)
    t.print_table(outname)
    # A hack to put in \hline
    with open(outname, 'r+') as f:
        lines = f.read().splitlines()
        newlines = [(line if 'hline' not in line else '\hline') for line in lines]
        f.seek(0)
        f.write('\n'.join(newlines))
        f.truncate()
