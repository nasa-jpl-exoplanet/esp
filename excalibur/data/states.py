'''Data Database Products View'''

# Heritage code shame:
# pylint: disable=too-many-locals,too-many-statements

# -- IMPORTS -- ------------------------------------------------------

import dawgie

import excalibur
from excalibur.util.plotters import save_plot_toscreen

import numpy as np
import matplotlib.pyplot as plt
import scipy
import math


# ------------- ------------------------------------------------------
# -- CALIBRATE -- ----------------------------------------------------
class CalibrateSV(dawgie.StateVector):
    '''data.calibration view'''

    def __init__(self, name):
        '''__init__ ds'''
        self.__name = name
        self._version_ = dawgie.VERSION(1, 1, 2)
        self['data'] = excalibur.ValuesDict()
        self['STATUS'] = excalibur.ValuesList()
        self['STATUS'].append(False)
        return

    def name(self):
        '''name ds'''
        return self.__name

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''view ds'''
        if len(self['STATUS']) == 2:
            if 'Spitzer' in self.__name:
                pass
            elif "JWST" in self.__name:
                scores = [
                    self['data']['SCORES'][k] for k in self['data']['SCORES']
                ]
                labels = self['data']['SCORES'].keys()
                myfig = plt.figure(figsize=(12, 9))
                plt.plot(scores, 'o--')
                plt.xticks(
                    ticks=np.arange(len(labels)),
                    labels=labels,
                    fontsize=16,
                    rotation=45,
                )
                plt.yticks(fontsize=16)
                plt.ylabel('DN', fontsize=20)
                save_plot_toscreen(myfig, visitor)
                pass
            else:
                data = self['data']
                timing = np.array(
                    [d for d, i in zip(data['TIME'], data['IGNORED']) if not i]
                )
                dispersion = np.array(
                    [
                        d
                        for d, i in zip(data['DISPERSION'], data['IGNORED'])
                        if not i
                    ]
                )
                shift = np.array(
                    [d for d, i in zip(data['SHIFT'], data['IGNORED']) if not i]
                )
                spec = np.array(
                    [
                        d
                        for d, i in zip(data['SPECTRUM'], data['IGNORED'])
                        if not i
                    ]
                )
                photoc = np.array(
                    [
                        d
                        for d, i in zip(data['PHT2CNT'], data['IGNORED'])
                        if not i
                    ]
                )
                wave = np.array(
                    [d for d, i in zip(data['WAVE'], data['IGNORED']) if not i]
                )
                errspec = np.array(
                    [
                        d
                        for d, i in zip(data['SPECERR'], data['IGNORED'])
                        if not i
                    ]
                )
                allignore = np.array(data['IGNORED'])
                allindex = np.arange(len(data['LOC']))
                vrange = data['VRANGE']
                torder = np.argsort(timing)
                allerr = []
                for s, e, w in zip(spec, errspec, wave):
                    select = (w > vrange[0]) & (w < vrange[1])
                    allerr.extend(e[select] / np.sqrt(s[select]))
                    pass
                allerr = np.array(allerr)
                select = np.isfinite(allerr)
                allerr = allerr[select]
                allerr = allerr[allerr > 0.9]
                vldi = np.random.uniform()
                vldi = int(vldi * (np.sum(~allignore.astype(bool)) - 1))
                vldi = allindex[~allignore.astype(bool)][vldi]

                strignore = (
                    str(np.nansum(allignore)) + ' / ' + str(len(allignore))
                )
                visitor.add_declaration('IGNORED: ' + strignore)

                myfig = plt.figure()
                plt.imshow(data['MEXP'][vldi])
                plt.colorbar()
                plt.title('Frame Index: ' + str(vldi))
                save_plot_toscreen(myfig, visitor)

                if np.sum(allignore) > 0:
                    gbgi = np.random.uniform()
                    gbgi = int(gbgi * (np.sum(allignore) - 1))
                    gbgi = allindex[allignore][gbgi]
                    myfig = plt.figure()
                    plt.imshow(data['MEXP'][gbgi])
                    plt.colorbar()
                    plt.title(
                        'Frame Index: '
                        + str(gbgi)
                        + ' Ignored: '
                        + data['TRIAL'][gbgi]
                    )
                    save_plot_toscreen(myfig, visitor)
                    pass

                myfig = plt.figure()
                for spectrum in data['SPECTRUM']:
                    plt.plot(spectrum)
                plt.xlabel('Pixel Number')
                plt.ylabel('Stellar Spectra [Counts]')
                save_plot_toscreen(myfig, visitor)

                myfig = plt.figure()
                for w, p, s in zip(wave, photoc, spec):
                    select = (w > vrange[0]) & (w < vrange[1])
                    plt.plot(w[select], s[select] / p[select])
                    pass
                plt.xlabel('Wavelength [$\\mu$m]')
                plt.ylabel('Stellar Spectra [Photons]')
                save_plot_toscreen(myfig, visitor)

                myfig = plt.figure()
                num_bins = int(math.sqrt(len(allerr)) / 10)
                plt.hist(allerr, bins=num_bins)
                plt.xlabel('Error Distribution [Noise Model Units]')
                mean = np.nanmean(allerr)
                stdev = np.nanstd(allerr)
                skew = scipy.stats.skew(allerr)
                stats_summary = (
                    "Mean: "
                    + str(round(mean, 3))
                    + "\nStandard Deviation: "
                    + str(round(stdev, 3))
                    + "\nSkewness: "
                    + str(round(skew, 3))
                )
                ax = plt.gca()
                plt.text(
                    0.95,
                    0.95,
                    stats_summary,
                    horizontalalignment='right',
                    verticalalignment='top',
                    transform=ax.transAxes,
                )
                save_plot_toscreen(myfig, visitor)

                myfig = plt.figure()
                plt.plot(dispersion[torder], 'o')
                plt.xlabel('Time Ordered Frame Number')
                plt.ylabel('Dispersion [Angstroms/Pixel]')
                plt.ylim(data['DISPLIM'][0], data['DISPLIM'][1])
                save_plot_toscreen(myfig, visitor)

                myfig = plt.figure()
                plt.plot(shift[torder] - np.nanmin(shift), 'o')
                plt.xlabel('Time Ordered Frame Number')
                plt.ylabel('Relative Shift [Pixels]')
                save_plot_toscreen(myfig, visitor)
        return


# -------------- -----------------------------------------------------
# -- TIMING -- -------------------------------------------------------
class TimingSV(dawgie.StateVector):
    '''data.timing view'''

    def __init__(self, name):
        '''1.2.0: GMR: Creates top level keys for each instrument'''
        self._version_ = dawgie.VERSION(1, 2, 0)
        self.__name = name
        # TOP LEVEL KEYS
        self['STATUS'] = excalibur.ValuesList()
        self['STATUS'].append(False)
        self['EXT'] = excalibur.ValuesDict()

        self['data'] = excalibur.ValuesDict()
        self['transit'] = excalibur.ValuesList()
        self['eclipse'] = excalibur.ValuesList()
        self['phasecurve'] = excalibur.ValuesList()
        return

    def name(self):
        '''name ds'''
        return self.__name

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''view ds'''
        if self['STATUS'][-1]:
            for p in self['data'].keys():
                if 'Spitzer' in self.__name or 'JWST' in self.__name:
                    vlabels = ['TRANSIT', 'ECLIPSE', 'PHASE CURVE']
                    hlabels = ['PLANET: ' + p, 'VISIT NUMBER']
                    table = visitor.add_table(
                        clabels=hlabels, rows=len(vlabels)
                    )
                    table.get_cell(0, 0).add_primitive(vlabels[0])
                    table.get_cell(0, 1).add_primitive(
                        self['data'][p]['transit']
                    )
                    table.get_cell(1, 0).add_primitive(vlabels[1])
                    table.get_cell(1, 1).add_primitive(
                        self['data'][p]['eclipse']
                    )
                    table.get_cell(2, 0).add_primitive(vlabels[2])
                    table.get_cell(2, 1).add_primitive(
                        self['data'][p]['phasecurve']
                    )
                else:
                    vlabels = ['TRANSIT', 'ECLIPSE', 'PHASE CURVE']
                    hlabels = ['PLANET: ' + p, 'VISIT NUMBER']
                    table = visitor.add_table(
                        clabels=hlabels, rows=len(vlabels)
                    )
                    table.get_cell(0, 0).add_primitive(vlabels[0])
                    table.get_cell(0, 1).add_primitive(
                        self['data'][p]['transit']
                    )
                    table.get_cell(1, 0).add_primitive(vlabels[1])
                    table.get_cell(1, 1).add_primitive(
                        self['data'][p]['eclipse']
                    )
                    table.get_cell(2, 0).add_primitive(vlabels[2])
                    table.get_cell(2, 1).add_primitive(
                        self['data'][p]['phasecurve']
                    )

                    tmetod = self['data'][p]['tmetod']
                    thro = self['data'][p]['thro']
                    thrs = self['data'][p]['thrs']
                    whereo = self['data'][p]['whereo']
                    wherev = self['data'][p]['wherev']
                    phase = self['data'][p]['phase']
                    ignore = self['data'][p]['ignore']
                    dvis = self['data'][p]['dvisits']
                    phsto = phase.copy()[self['data'][p]['ordt']]
                    ignto = ignore.copy()[self['data'][p]['ordt']]
                    dvisto = dvis.copy()[self['data'][p]['ordt']]

                    myfig = plt.figure()
                    plt.plot(phsto, 'k.')
                    plt.plot(np.arange(phsto.size)[~ignto], phsto[~ignto], 'bo')
                    for i in wherev:
                        plt.axvline(i, ls='--', color='r')
                    for i in whereo:
                        plt.axvline(i, ls='-.', color='g')
                    plt.xlim(0, phsto.size - 1)
                    plt.ylim(-0.5, 0.5)
                    plt.xlabel('Time index')
                    plt.ylabel('Orbital Phase [2pi rad]')
                    save_plot_toscreen(myfig, visitor)

                    myfig = plt.figure()
                    plt.plot(tmetod, 'o')
                    plt.plot(tmetod * 0 + thro, 'r--')
                    plt.plot(tmetod * 0 + thrs, 'g-.')
                    for i in wherev:
                        plt.axvline(i, ls='--', color='r')
                    for i in whereo:
                        plt.axvline(i, ls='-.', color='g')
                    plt.xlim(0, tmetod.size - 1)
                    plt.xlabel('Time index')
                    plt.ylabel('Frame Separation [Days]')
                    plt.semilogy()
                    save_plot_toscreen(myfig, visitor)

                    if np.max(dvis) > np.max(self['data'][p]['visits']):
                        myfig = plt.figure()
                        plt.plot(dvisto, 'o')
                        plt.xlim(0, tmetod.size - 1)
                        plt.xlabel('Time index')
                        plt.ylabel('Double Scan Visit Number')
                        save_plot_toscreen(myfig, visitor)

        return


# ------------ -------------------------------------------------------
