'''transit algorithms ds'''

# Heritage code shame:
# pylint: disable=duplicate-code
# pylint: disable=too-many-arguments,too-many-branches,too-many-locals,too-many-nested-blocks,too-many-positional-arguments,too-many-statements

# -- IMPORTS -- ------------------------------------------------------
import dawgie
import dawgie.context

import numexpr

import logging

from collections import defaultdict
import numpy as np

import excalibur.transit.core as trncore
import excalibur.transit.starspots as trnspots
import excalibur.transit.states as trnstates

from excalibur.util.checksv import checksv

import excalibur.data as dat
import excalibur.data.algorithms as datalg
import excalibur.runtime as rtime
import excalibur.runtime.algorithms as rtalg
import excalibur.runtime.binding as rtbind
import excalibur.system as sys
import excalibur.system.algorithms as sysalg

from importlib import import_module as fetch  # avoid cicular dependencies

log = logging.getLogger(__name__)

numexpr.ncores = 1  # this is actually a performance enhancer!

# ------------- ------------------------------------------------------
# -- ALGO RUN OPTIONS -- ---------------------------------------------
# FILTERS
fltrs = [str(fn) for fn in rtbind.filter_names.values()]


# ---------------------- ---------------------------------------------
# -- ALGORITHMS -- ---------------------------------------------------
class Normalization(dawgie.Algorithm):
    '''
    G. ROUDIER: Light curve normalization by Out Of Transit data
    '''

    def __init__(self):
        self._version_ = trncore.normversion()
        self._type = 'transit'
        self.__cal = datalg.Calibration()
        self.__tme = datalg.Timing()
        self.__rt = rtalg.Autofill()
        self.__fin = sysalg.Finalize()
        self.__out = [trnstates.NormSV(fltr) for fltr in fltrs]
        return

    def name(self):
        '''Database name for subtask extension'''
        return 'normalization'

    def previous(self):
        '''Input State Vectors: data.calibration, data.timing, system.finalize'''
        return [
            dawgie.ALG_REF(dat.task, self.__cal),
            dawgie.ALG_REF(dat.task, self.__tme),
            dawgie.ALG_REF(sys.task, self.__fin),
        ] + self.__rt.refs_for_proceed()

    def state_vectors(self):
        '''Output State Vectors: transit.normalization'''
        return self.__out

    def run(self, ds, ps):
        '''Top level algorithm call'''

        svupdate = []
        vfin, sfin = checksv(self.__fin.sv_as_dict()['parameters'])
        for fltr in self.__rt.sv_as_dict()['status']['allowed_filter_names']:
            # stop here if it is not a runtime target
            self.__rt.proceed(fltr)

            update = False
            vcal, scal = checksv(self.__cal.sv_as_dict()[fltr])
            vtme, stme = checksv(self.__tme.sv_as_dict()[fltr])
            if vcal and vtme and vfin:
                log.warning(
                    '--< %s NORMALIZATION: %s >--', self._type.upper(), fltr
                )
                update = self._norm(
                    self.__cal.sv_as_dict()[fltr],
                    self.__tme.sv_as_dict()[fltr],
                    self.__fin.sv_as_dict()['parameters'],
                    fltrs.index(fltr),
                )
                pass
            else:
                errstr = [m for m in [scal, stme, sfin] if m is not None]
                self._failure(errstr[0])
                pass
            if update:
                svupdate.append(self.__out[fltrs.index(fltr)])
            pass
        self.__out = svupdate
        if self.__out:
            ds.update()
        else:
            raise dawgie.NoValidOutputDataError(
                f'No output created for {self._type.upper()}.{self.name()}'
            )
        return

    def _norm(self, cal, tme, fin, index):
        '''Core code call'''
        if 'Spitzer' in fltrs[index]:
            normed = trncore.norm_spitzer(
                cal, tme, fin, self.__out[index], self._type
            )
            pass
        elif 'NIRISS' in fltrs[index]:
            # Not compatible with new calibration. To be updated.
            # normed = trncore.norm_jwst_niriss(cal, tme, fin, self.__out[index], self._type)
            normed = False
            pass
        elif 'NIRSPEC' in fltrs[index]:
            normed = trncore.norm_jwst(
                cal,
                tme,
                fin,
                fltrs[index],
                self.__out[index],
                self._type,
                debug=False,
            )
            pass
        else:
            normed = trncore.norm(
                cal,
                tme,
                fin,
                fltrs[index],
                self.__out[index],
                self._type,
                verbose=False,
            )
            pass
        return normed

    def _failure(self, errstr):
        '''Failure log'''
        log.warning('--< %s NORMALIZATION: %s >--', self._type.upper(), errstr)
        return

    pass


class WhiteLight(dawgie.Algorithm):
    '''
    G. ROUDIER: See inheritance and CI5 thread with A NIESSNER
    for __init__() method and class attributes
    https://github-fn.jpl.nasa.gov/EXCALIBUR/esp/pull/86
    '''

    def __init__(self, nrm=Normalization()):
        '''__init__ ds'''
        self._version_ = trncore.wlversion()
        self._type = 'transit'
        self._nrm = nrm
        self.__fin = sysalg.Finalize()
        self.__rt = rtalg.Autofill()
        self.__out = [trnstates.WhiteLightSV(fltr) for fltr in fltrs]
        self.__out.append(trnstates.WhiteLightSV('HST'))
        return

    def name(self):
        '''Database name for subtask extension'''
        return 'whitelight'

    def previous(self):
        '''Input State Vectors: transit.normalization, system.finalize'''
        return [
            dawgie.ALG_REF(fetch('excalibur.transit').task, self._nrm),
            dawgie.ALG_REF(sys.task, self.__fin),
            dawgie.V_REF(
                rtime.task,
                self.__rt,
                self.__rt.sv_as_dict()['status'],
                'spectrum_steps',
            ),
        ] + self.__rt.refs_for_proceed()

    def state_vectors(self):
        '''Output State Vectors: transit.whitelight'''
        return self.__out

    def run(self, ds, ps):
        '''Top level algorithm call'''

        svupdate = []
        fin = self.__fin.sv_as_dict()['parameters']
        vfin, sfin = checksv(fin)
        # MERGE PROTOTYPE
        if self._type == "transit":
            allnormdata = []
            allfilters = []
            hstfltrs = [
                'HST-WFC3-IR-G141-SCAN',
                'HST-WFC3-IR-G102-SCAN',
                'HST-STIS-CCD-G750L-STARE',
                'HST-STIS-CCD-G430L-STARE',
            ]
            for fltr in self.__rt.sv_as_dict()['status'][
                'allowed_filter_names'
            ]:
                if fltr in hstfltrs:
                    # stop here if it is not a runtime target
                    self.__rt.proceed(fltr)

                    update = False
                    try:
                        nrm = self._nrm.sv_as_dict()[fltr]
                    except KeyError:
                        break
                    vnrm, snrm = checksv(nrm)
                    if vnrm and vfin:
                        log.warning(
                            '--< %s MERGING: %s >--', self._type.upper(), fltr
                        )
                        allnormdata.append(nrm)
                        allfilters.append(fltr)
                        update = True
                    else:
                        errstr = [m for m in [snrm, sfin] if m is not None]
                        self._failure(errstr[0])
                    pass
                pass
            if allnormdata:
                try:
                    update = self._hstwhitelight(
                        allnormdata,
                        fin,
                        self.__rt.sv_as_dict()['status'][
                            'spectrum_steps'
                        ].value(),
                        self.__out[-1],
                        allfilters,
                    )
                    if update:
                        svupdate.append(self.__out[-1])
                    pass
                except TypeError:
                    log.warning(
                        '>-- HST orbit solution failed %s', self._type.upper()
                    )
                    pass
                pass
        # FILTER LOOP
        for fltr in self.__rt.sv_as_dict()['status']['allowed_filter_names']:
            update = False
            nrm = self._nrm.sv_as_dict()[fltr]
            vnrm, snrm = checksv(nrm)
            if vnrm and vfin:
                log.warning(
                    '--< %s WHITE LIGHT: %s >--', self._type.upper(), fltr
                )
                update = self._whitelight(
                    nrm,
                    fin,
                    self.__rt.sv_as_dict()['status']['spectrum_steps'].value(),
                    self.__out[fltrs.index(fltr)],
                    fltr,
                )
            else:
                errstr = [m for m in [snrm, sfin] if m is not None]
                self._failure(errstr[0])
            if update:
                svupdate.append(self.__out[fltrs.index(fltr)])
        self.__out = svupdate
        if self.__out:
            ds.update()
        else:
            raise dawgie.NoValidOutputDataError(
                f'No output created for {self._type.upper()}.{self.name()}'
            )
        return

    def _hstwhitelight(self, nrm, fin, chain_length, out, fltr):
        '''Core code call for merged HST data'''
        # chain_length = 10
        # print('chain length in hstwhitelight',chain_length)

        wl = trncore.hstwhitelight(
            nrm,
            fin,
            out,
            fltr,
            self._type,
            chainlen=chain_length,
            verbose=False,
        )
        return wl

    def _whitelight(self, nrm, fin, chain_length, out, fltr):
        '''Core code call'''
        # chain_length = 10
        # print('chain length in whitelight',chain_length)

        if 'Spitzer' in fltr:
            wl = trncore.lightcurve_spitzer(
                nrm, fin, out, self._type, fltr, self.__out[-1]
            )
        elif 'JWST' in fltr:
            wl = trncore.lightcurve_jwst_niriss(
                nrm, fin, out, self._type, fltr, self.__out[-1], method='ns'
            )
        else:
            wl = trncore.whitelight(
                nrm,
                fin,
                out,
                fltr,
                self._type,
                self.__out[-1],
                chainlen=chain_length,
                verbose=False,
                parentprior=True,
            )
            pass
        return wl

    def _failure(self, errstr):
        '''Failure log'''
        log.warning('--< %s WHITE LIGHT: %s >--', self._type.upper(), errstr)
        return

    pass


class Spectrum(dawgie.Algorithm):
    '''
    G. ROUDIER: See inheritance and CI5 thread with A NIESSNER for __init__() method and class attributes https://github-fn.jpl.nasa.gov/EXCALIBUR/esp/pull/86
    '''

    def __init__(self, nrm=Normalization(), wht=WhiteLight()):
        '''__init__ ds'''
        self._version_ = trncore.spectrumversion()
        self._type = 'transit'
        self.__fin = sysalg.Finalize()
        self._nrm = nrm
        self.__rt = rtalg.Autofill()
        self._wht = wht
        self.__out = [trnstates.SpectrumSV(fltr) for fltr in fltrs]
        if sum('HST' in fltr for fltr in fltrs) > 1:
            self.__out.append(trnstates.SpectrumSV('STIS-WFC3'))
            pass
        return

    def name(self):
        '''Database name for subtask extension'''
        return 'spectrum'

    def previous(self):
        '''Input State Vectors: system.finalize, transit.normalization,
        transit.whitelight'''
        return [
            dawgie.ALG_REF(sys.task, self.__fin),
            dawgie.ALG_REF(fetch('excalibur.transit').task, self._nrm),
            dawgie.ALG_REF(fetch('excalibur.transit').task, self._wht),
            dawgie.V_REF(
                rtime.task,
                self.__rt,
                self.__rt.sv_as_dict()['status'],
                'spectrum_steps',
            ),
        ] + self.__rt.refs_for_proceed()

    def state_vectors(self):
        '''Output State Vectors: transit.spectrum'''
        return self.__out

    def run(self, ds, ps):
        '''Top level algorithm call'''

        svupdate = []
        vfin, sfin = checksv(self.__fin.sv_as_dict()['parameters'])

        for fltr in self.__rt.sv_as_dict()['status']['allowed_filter_names']:
            # stop here if it is not a runtime target
            self.__rt.proceed(fltr)

            update = False
            vnrm, snrm = checksv(self._nrm.sv_as_dict()[fltr])
            vwht, swht = checksv(self._wht.sv_as_dict()[fltr])
            if vfin and vnrm and vwht:
                log.warning('--< %s SPECTRUM: %s >--', self._type.upper(), fltr)
                update = self._spectrum(
                    self.__fin.sv_as_dict()['parameters'],
                    self._nrm.sv_as_dict()[fltr],
                    self._wht.sv_as_dict()[fltr],
                    self.__rt.sv_as_dict()['status']['spectrum_steps'].value(),
                    self.__out[fltrs.index(fltr)],
                    fltr,
                )
            else:
                errstr = [m for m in [sfin, snrm, swht] if m is not None]
                self._failure(errstr[0])
            if update:
                svupdate.append(self.__out[fltrs.index(fltr)])

        merg = trncore.hstspectrum(self.__out, fltrs)
        log.warning('--< %s SPECTRUM MERGED: %s >--', self._type.upper(), merg)
        if merg:
            svupdate.append(self.__out[-1])

        self.__out = (
            svupdate  # it will take all the elements that are not empty
        )
        if self.__out:
            ds.update()
        else:
            raise dawgie.NoValidOutputDataError(
                f'No output created for {self._type.upper()}.{self.name()}'
            )
        return

    def _spectrum(self, fin, nrm, wht, chain_length, out, fltr):
        '''Core code call'''
        # chain_length = 10
        # print('chain length in spectrum',chain_length)

        if "Spitzer" in fltr:
            s = trncore.spitzer_spectrum(wht, out, fltr)
        elif "JWST" in fltr:
            s = trncore.jwst_niriss_spectrum(
                nrm, fin, out, self._type, wht, method='lm'
            )
        else:
            s = trncore.spectrum(
                fin,
                nrm,
                wht,
                out,
                fltr,
                self._type,
                chainlen=chain_length,
                verbose=False,
            )
            pass
        return s

    def _failure(self, errstr):
        '''Failure log'''
        log.warning('--< %s SPECTRUM: %s >--', self._type.upper(), errstr)
        return

    pass


class StarSpots(dawgie.Algorithm):
    '''
    include Viktor's starspot model
    '''

    def __init__(self, spc=Spectrum()):
        '''__init__ ds'''
        self._version_ = trncore.spectrumversion()
        self._type = 'transit'
        self.__fin = sysalg.Finalize()
        self._spc = spc
        self.__rt = rtalg.Autofill()
        self.__out = [trnstates.StarspotSV(fltr) for fltr in fltrs]
        return

    def name(self):
        '''Database name for subtask extension'''
        return 'starspots'

    def previous(self):
        '''Input State Vectors: system.finalize, transit.normalization,
        transit.whitelight'''
        return [
            dawgie.ALG_REF(sys.task, self.__fin),
            dawgie.ALG_REF(fetch('excalibur.transit').task, self._spc),
        ] + self.__rt.refs_for_proceed()

    def state_vectors(self):
        '''Output State Vectors: transit.starspots'''
        return self.__out

    def run(self, ds, ps):
        '''Top level algorithm call'''

        svupdate = []
        vfin, sfin = checksv(self.__fin.sv_as_dict()['parameters'])

        # print('STARSPOTS: allowed filters',
        #    self.__rt.sv_as_dict()['status']['allowed_filter_names'])

        # for fltr in ['HST-WFC3-IR-G141-SCAN']:  # for debugging, just run G141
        for fltr in self.__rt.sv_as_dict()['status']['allowed_filter_names']:
            # stop here if it is not a runtime target
            self.__rt.proceed(fltr)

            update = False
            vspc, sspc = checksv(self._spc.sv_as_dict()[fltr])
            if vfin and vspc:
                # log.warning(
                #     '--< %s STARSPOTS: %s >--', self._type.upper(), fltr
                # )
                update = self._starspots(
                    self.__fin.sv_as_dict()['parameters'],
                    self._spc.sv_as_dict()[fltr],
                    self.__out[fltrs.index(fltr)],
                )
            else:
                errstr = [m for m in [sfin, sspc] if m is not None]
                self._failure(errstr[0])
            if update:
                svupdate.append(self.__out[fltrs.index(fltr)])

        self.__out = svupdate
        if self.__out:
            ds.update()
        else:
            raise dawgie.NoValidOutputDataError(
                f'No output created for {self._type.upper()}.{self.name()}'
            )
        return

    @staticmethod
    def _starspots(fin, spc, out):
        '''Core code call'''
        spotmodel = trnspots.starspots(fin, spc, out)
        return spotmodel

    def _failure(self, errstr):
        '''Failure log'''
        log.warning('--< %s STARSPOTS: %s >--', self._type.upper(), errstr)
        return

    pass


class Population(dawgie.Analyzer):
    '''population ds'''

    def __init__(self):
        '''__init__ ds'''
        self._version_ = dawgie.VERSION(1, 0, 2)
        self.__out = [trnstates.PopulationSV(fltr) for fltr in fltrs]
        return

    def feedback(self):
        '''feedback ds'''
        return []

    def name(self):
        '''Database name for subtask extension'''
        return 'population'

    def traits(self) -> [dawgie.SV_REF, dawgie.V_REF]:
        '''Input State Vectors'''
        return [
            *[
                dawgie.SV_REF(
                    fetch('excalibur.transit').task,
                    Spectrum(),
                    Spectrum().state_vectors()[i],
                )
                for i in range(len(Spectrum().state_vectors()))
            ],
            *[
                dawgie.SV_REF(
                    fetch('excalibur.transit').task,
                    WhiteLight(),
                    WhiteLight().state_vectors()[i],
                )
                for i in range(len(WhiteLight().state_vectors()))
            ],
        ]

    def state_vectors(self):
        '''Output State Vectors: transit.population'''
        return self.__out

    def run(self, aspects: dawgie.Aspect):
        '''Top level algorithm call'''
        # now handle IM parameter distribution
        banned_params = ['rprs']
        sv_prefix = 'transit.spectrum.'
        wl_prefix = 'transit.whitelight.'
        for fltr in fltrs:
            im_bins = defaultdict(lambda: defaultdict(list))
            wl_im_bins = defaultdict(lambda: defaultdict(list))
            for trgt in aspects:
                svname = sv_prefix + fltr
                wlname = wl_prefix + fltr
                if svname not in aspects[trgt] or wlname not in aspects[trgt]:
                    continue
                tr_data = aspects[trgt][svname]
                wl_data = aspects[trgt][wlname]
                for pl in tr_data['data']:
                    # verify SV succeeded for target
                    if (
                        tr_data['STATUS'][-1]
                        or 'MCTRACE' in tr_data['data'][pl]
                    ):
                        # logic for saving spectrum IM parameters
                        if 'MCTRACE' in tr_data['data'][pl]:
                            if pl in wl_data['data']:
                                num_visits = len(wl_data['data'][pl]['visits'])
                            else:  # always at least one visit
                                num_visits = 1
                            trace = tr_data['data'][pl]['MCTRACE']
                            for row in trace:
                                for key in row:
                                    keyparts = key.split('__')
                                    if (
                                        len(keyparts) > 1
                                        and int(keyparts[1]) >= num_visits
                                    ):
                                        continue
                                    param_name = keyparts[0]
                                    if param_name in banned_params:
                                        continue
                                    param_val = np.nanmedian(row[key])
                                    im_bins[param_name]['values'].append(
                                        param_val
                                    )
                                    bins = 60
                                    if 'G430' in fltr:
                                        bins = 30
                                    elif 'G750' in fltr:
                                        bins = 40
                                    elif 'G102' in fltr:
                                        bins = 30
                                    im_bins[param_name]['bins'] = bins
                    if pl not in wl_data['data']:
                        continue
                    if (
                        wl_data['STATUS'][-1]
                        or 'mctrace' in wl_data['data'][pl]
                    ):
                        # logic for saving whitelight IM parameters
                        if 'mctrace' in wl_data['data'][pl]:
                            if 'visits' not in wl_data['data'][pl]:
                                continue
                            num_visits = len(wl_data['data'][pl]['visits'])
                            trace = wl_data['data'][pl]['mctrace']
                            for key in trace:
                                keyparts = key.split('__')
                                if (
                                    len(keyparts) > 1
                                    and int(keyparts[1]) >= num_visits
                                ):
                                    continue
                                param_name = keyparts[0]
                                if param_name in banned_params:
                                    continue
                                param_val = np.nanmedian(trace[key])
                                wl_im_bins[param_name]['values'].append(
                                    param_val
                                )

            idx = fltrs.index(fltr)
            self.__out[idx]['data']['IMPARAMS'] = dict(im_bins)
            self.__out[idx]['data']['wl'] = {}
            self.__out[idx]['data']['wl']['imparams'] = dict(wl_im_bins)
            self.__out[idx]['STATUS'].append(True)

        aspects.ds().update()
        return

    pass
