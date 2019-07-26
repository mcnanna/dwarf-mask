#!/usr/bin/env python

import numpy as np
from textable import TexTable

def signal_table(outname, *signals, **kwargs):
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

    justs = 'lcccccccccccc'
    
    comments = kwargs['comments'] if 'comments' in kwargs.keys() else "\\knowncomments"
    #comments += " References: " + ref_dict_str
    caption = kwargs['caption'] if 'caption' in kwargs.keys() else "\\knowncaption"
    notes = kwargs['notes'] if 'notes' in kwargs.keys() else "\\knownnotes"

    t = TexTable(len(justs), justs=justs, comments=comments, caption=caption, notes=notes, fontsize="\\scriptsize", doc=True)
    t.add_header_row(['({})'.format(i+1) for i in range(len(justs))])
    t.add_header_row(
            ['Name', r'$\sqrt{\mathrm{TS}}$', 'SIG', r"$P_{\rm det}$",
                "RA", "Dec",
                "$m - M$", r"$a_h$", r"$\epsilon$",
                "Distance", r"$M_V$", r"$r_{1/2}$", 
                "Ref."])
    t.add_header_row(
            ['', '', '', '',
                r'($\deg$)', r'($\deg$)', 
                '', r"($'$)", '',
                '(kpc)', '(mag)', '(pc)',
                ''])
    roundings = [0, 1, 1, 2,
            2, 2, 
            1, 1, 2,
            0, 1, 0,
            0]
    #ignore = ['mod_ugali', 'mod_simple', 'angsep_ugali', 'angsep_simple']
    #t.add_data([signals[0][header] for header in signals[0].dtype.names if header not in ignore], sigfigs=roundings)
    data_headers = ['name', 'TS', 'SIG', 'pdet', 'ra', 'dec', 'mod_actual', 'ah', 'ellipticity', 'distance', 'M_V', 'r12', 'ref']
    t.add_data([signals[0][header] for header in data_headers], sigfigs=roundings)
    for signal in signals[1:]:
        t.add_data(['hline']*len(justs))
        #t.add_data([signal[header] for header in signal.dtype.names if header not in ignore], sigfigs=roundings)
        t.add_data([signal[header] for header in data_headers], sigfigs=roundings)
    t.print_table(outname)
    # A hack to put in \hline
    with open(outname, 'r+') as f:
        lines = f.read().splitlines()
        newlines = [(line if 'hline' not in line else '\hline') for line in lines]
        f.seek(0)
        f.write('\n'.join(newlines))
        f.truncate()

    # A hack to add spaces after columns
    if 'spaced' in kwargs.keys():
        spaced = kwargs['spaced']
        gaps = [(spaced[i] - spaced[i-1]) if i>0 else spaced[i] for i in range(len(spaced))]
        
        with open(outname, 'r+') as f:
            lines = f.read().splitlines()
            for i in range(len(lines)):
                # Add spaces between justifications
                if justs in lines[i]:
                    idx = lines[i].find(justs)
                    newline = lines[i][:idx]
                    for gap in gaps:
                        newline += lines[i][idx:idx+gap]
                        newline += r' @{\hspace{0.3in}} '
                        idx += gap
                    newline += lines[i][idx:]
                    lines[i] = newline
                # Add spaces in columns heads
                elif 'colhead' in lines[i]:
                    heads = lines[i].split(' & ')
                    for col in spaced:
                        heads[col-1] += r'\hspace{0.25in}'
                    lines[i] = (' & ').join(heads)

            f.seek(0)
            f.write('\n'.join(lines))
            f.truncate()

    # A hack to add tablenotemarks. This is not currently configureable, you just have to edit this code
    marks_des = [('Tucana IV', 'a')]
    marks_ps1 = [('Leo I', 'a'), ('Leo II', 'a'), ('Leo IV', 'b'), ('Columba I', 'bc'), ('Bootes III', 'd'), ('Bootes IV', 'b')]
    marks_both = [('Tucana IV', 'a'), ('Leo I', 'b'), ('Leo II', 'b'), ('Leo IV', 'c'), ('Columba I', 'cd'), ('Bootes III', 'e'), ('Bootes IV', 'c')]

    with open(outname, 'r+') as f:
        if 'des' in outname:
            marks = marks_des
        elif 'ps1' in outname:
            marks = marks_ps1
        else:
            marks = marks_both

        contents = f.read()
        for sat, mark in marks:
            idx = contents.index(sat)
            with_mark = sat + '\\tablenotemark{{{0}}}'.format(mark)
            contents = contents[:idx] + with_mark + contents[idx+len(sat):]

        f.seek(0)
        f.write(contents)
        f.truncate()

def remains_table(outname, *remains, **alg):
    alg = alg.pop('alg', 'simple') # Default alg is 'simple'

    SIG = 'TS' if alg == 'ugali' else 'SIG'

    justs = 'cccc'
    header_row1 = [SIG, r"$\alpha_{2000}$", r"$\delta_{2000}$", r"$m - M$"]
    header_row2 = ['', '(deg)', '(deg)', '']
    data_headers = [SIG, 'ra', 'dec', 'modulus']
    roundings = [1, 2, 2, 1]

    if alg == 'ugali':
        justs = 'l' + justs
        header_row1 = ['Name'] + header_row1
        header_row2 = [''] + header_row2
        data_headers = ['name'] + data_headers
        roundings = [0] + roundings

    if alg == 'both':
        for remain in remains:
            if 'TS' in remain.dtype.names:
                remain['TS'] = np.sqrt(remain['TS'])

        justs = 'lc' + justs + 'cc'
        header_row1 = ['Name', r'$\sqrt{\mathrm{TS}}$'] + header_row1 + [r"$m - M$", "Angular Separation"]
        header_row2 = ['', '(ugali)', '(simple)', '(deg)', '(deg)', '(ugali)', '(simple)', r"($'$)"]
        data_headers = remains[0].dtype.names
        roundings = [0, 1] + roundings + [1, 2]
    if 'RA' in remains[0].dtype.names:
        data_headers = [header.upper() for header in data_headers]
    
    t = TexTable(len(justs), justs=justs, comments="\\candidatecomments", caption="\\candidatecaption", notes="\\knownnotes", fontsize="\\scriptsize", doc=True)
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
