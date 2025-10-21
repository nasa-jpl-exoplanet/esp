'''data algorithms dc'''

# Heritage code shame:
# pylint: disable=too-many-arguments,too-many-locals,too-many-positional-arguments

# -- IMPORTS -- ------------------------------------------------------
import logging

import dawgie
import dawgie.context

import excalibur

import excalibur.data.core as datcore
import excalibur.data.states as datstates

import excalibur.target as trg
import excalibur.target.states as trgstates
import excalibur.target.algorithms as trgalg

import excalibur.system as sys
import excalibur.system.algorithms as sysalg

import excalibur.runtime.algorithms as rtalg
import excalibur.runtime.binding as rtbind

from excalibur.util.checksv import checksv

from importlib import import_module as fetch  # avoid cicular dependencies

log = logging.getLogger(__name__)

# ------------- ------------------------------------------------------
# -- ALGO RUN OPTIONS -- ---------------------------------------------
fltrs = [str(fn) for fn in rtbind.filter_names.values()]


# ---------------------- ---------------------------------------------
# -- ALGORITHMS -- ---------------------------------------------------
class Collect(dawgie.Algorithm):
    '''
    G. ROUDIER: Data collection by filters
    '''

    def __init__(self):
        '''__init__ ds'''
        self._version_ = datcore.collectversion()
        self.__create = trgalg.Create()
        self.__rt = rtalg.Autofill()
        self.__scrape = trgalg.Scrape()
        self.__out = trgstates.FilterSV('frames')
        return

    def name(self):
        '''Database name for subtask extension'''
        return 'collect'

    def previous(self):
        '''Input State Vectors: target.create, target.scrape'''
        return [
            dawgie.ALG_REF(trg.analysis, self.__create),
            dawgie.ALG_REF(trg.task, self.__scrape),
        ] + self.__rt.refs_for_validity()

    def state_vectors(self):
        '''Output State Vectors: data.collect'''
        return [self.__out]

    def run(self, ds, ps):
        '''Top level algorithm call'''

        # stop here if it is not a runtime target
        if not self.__rt.is_valid():
            log.info('--< DATA.%s: not a valid target >--', self.name().upper())
            pass
        # 7/2/24 we want to always collect all filters, so don't use proceed(fltr) here
        else:
            update = False
            create = self.__create.sv_as_dict()
            scrape = self.__scrape.sv_as_dict()['databases']
            valid, errstring = checksv(scrape)
            if valid:
                for key in create['filters'].keys():
                    self.__out[key] = create['filters'][key]
                    pass
                # for fltr in self.__rt.sv_as_dict()['status']['allowed_filter_names']:
                for fltr in fltrs:
                    ok = self._collect(fltr, scrape, self.__out)
                    update = update or ok
                    pass
                if update:
                    _ = excalibur.lagger()
                    ds.update()
                    pass
                else:
                    self._raisenoout(self.name())
            else:
                self._failure(errstring)
            pass
        return

    @staticmethod
    def _raisenoout(myname):
        '''No output exception'''
        raise dawgie.NoValidOutputDataError(
            f'No output created for DATA.{myname}'
        )

    @staticmethod
    def _collect(fltr, scrape, out):
        '''Core code call'''
        collected = datcore.collect(fltr, scrape, out)
        return collected

    @staticmethod
    def _failure(errstr):
        '''Failure log'''
        log.warning('--< DATA COLLECT: %s >--', errstr)
        return

    pass


class Timing(dawgie.Algorithm):
    '''
    G. ROUDIER: Categorize data into 3 science purposes: TRANSIT, ECLIPSE, PHASECURVE
    '''

    def __init__(self):
        '''__init__ ds'''
        self._version_ = datcore.timingversion()
        self.__fin = sysalg.Finalize()
        self.__col = Collect()
        self.__rt = rtalg.Autofill()
        self.__out = [datstates.TimingSV(fltr) for fltr in fltrs]
        return

    def name(self):
        '''Database name for subtask extension'''
        return 'timing'

    def previous(self):
        '''Input State Vectors: system.finalize, data.collect'''
        return [
            dawgie.ALG_REF(sys.task, self.__fin),
            dawgie.ALG_REF(fetch('excalibur.data').task, self.__col),
        ] + self.__rt.refs_for_proceed()

    def state_vectors(self):
        '''Output State Vectors: data.timing'''
        return self.__out

    def run(self, ds, ps):
        '''Top level algorithm call'''
        update = False
        fin = self.__fin.sv_as_dict()['parameters']
        vfin, efin = checksv(fin)
        col = self.__col.sv_as_dict()['frames']
        vcol, ecol = checksv(col)
        svupdate = []
        if vfin and vcol:
            log.info('>-- DATA COLLECT: \n\t%s', list(col['activefilters']))
            letmethrough = self.__rt.sv_as_dict()['status'][
                'allowed_filter_names'
            ]
            log.info('>-- RUNTIME: \n\t%s', letmethrough)
            letmethrough = [
                k for k in letmethrough if k in col['activefilters']
            ]
            log.info('>-- DATA TIMING: \n\t%s', letmethrough)
            for fltr in letmethrough:
                # stop here if it is not a runtime target
                self.__rt.proceed(fltr)
                update = self._timing(
                    fin,
                    fltr,
                    col['activefilters'][fltr],
                    self.__out[fltrs.index(fltr)],
                )
                if update:
                    svupdate.append(self.__out[fltrs.index(fltr)])
                pass
            pass
        else:
            message = [m for m in [efin, ecol] if m is not None]
            self._failure(message[0])
            pass
        self.__out = svupdate
        if self.__out:
            _ = excalibur.lagger()
            ds.update()
            pass
        else:
            self._raisenoout(self.name())
        return

    @staticmethod
    def _raisenoout(myname):
        '''No output exception'''
        raise dawgie.NoValidOutputDataError(
            f'No output created for DATA.{myname}'
        )

    @staticmethod
    def _timing(fin, fltr, colin, out):
        '''Core code call'''
        log.info('--< DATA TIMING: %s >--', fltr)
        chunked = datcore.timing(fin, fltr, colin, out)
        return chunked

    @staticmethod
    def _failure(errstr):
        '''Failure log'''
        log.warning('--< DATA TIMING: %s >--', errstr)
        return

    pass


class Calibration(dawgie.Algorithm):
    '''
    G. ROUDIER: Data re-calibration and reduction
    '''

    def __init__(self):
        '''__init__ ds'''
        self._version_ = dawgie.VERSION(1, 4, 4)
        self.__fin = sysalg.Finalize()
        self.__col = Collect()
        self.__rt = rtalg.Autofill()
        self.__tim = Timing()
        self.__out = [datstates.CalibrateSV(fltr) for fltr in fltrs]
        return

    def name(self):
        '''Database name for subtask extension'''
        return 'calibration'

    def previous(self):
        '''Input State Vectors: system.finalize, data.collect, data.timing'''
        return [
            dawgie.ALG_REF(sys.task, self.__fin),
            dawgie.ALG_REF(fetch('excalibur.data').task, self.__col),
            dawgie.ALG_REF(fetch('excalibur.data').task, self.__tim),
        ] + self.__rt.refs_for_proceed()

    def state_vectors(self):
        '''Output State Vectors: data.calibration'''
        return self.__out

    def run(self, ds, ps):
        '''Top level algorithm call'''

        update = False
        cll = self.__col.sv_as_dict()['frames']
        vcll, ecll = checksv(cll)
        fin = self.__fin.sv_as_dict()['parameters']
        vfin, sfin = checksv(fin)
        svupdate = []

        for fltr in self.__rt.sv_as_dict()['status']['allowed_filter_names']:
            if fltr in cll['activefilters']:
                # stop here if it is not a runtime target
                self.__rt.proceed(fltr)
                tim = self.__tim.sv_as_dict()[fltr]
                vtim, etim = checksv(tim)
                if vfin and vcll and vtim:
                    # This code needs repair by moving out to config
                    # Dirty secret in repr(self)
                    update = self._calib(
                        fin,
                        cll['activefilters'][fltr],
                        tim,
                        repr(self).split('.')[1],
                        fltr,
                        self.__out[fltrs.index(fltr)],
                    )
                    if update:
                        svupdate.append(self.__out[fltrs.index(fltr)])
                        pass
                    pass
                else:
                    message = [m for m in [sfin, ecll, etim] if m is not None]
                    self._failure(message[0])
                    pass
                pass
            pass
        self.__out = svupdate
        if self.__out:
            _ = excalibur.lagger()
            ds.update()
            pass
        else:
            self._raisenoout(self.name())
        return

    @staticmethod
    def _raisenoout(myname):
        '''No output exception'''
        raise dawgie.NoValidOutputDataError(
            f'No output created for DATA.{myname}'
        )

    @staticmethod
    def _calib(fin, cll, tim, tid, flttype, out):
        '''Core code call'''
        log.info('--< DATA CALIBRATION: %s >--', flttype)
        caled = False
        if 'SCAN' in flttype:
            caled = datcore.scancal(cll, tim, tid, flttype, out)
            pass
        if 'G430' in flttype:
            caled = datcore.stiscal_G430L(fin, cll, tim, tid, flttype, out)
            pass
        if 'G750' in flttype:
            caled = datcore.stiscal_G750L(fin, cll, tim, tid, flttype, out)
            pass
        if 'Spitzer' in flttype:
            caled = datcore.spitzercal(cll, out)
            pass
        if 'JWST' in flttype:
            caled = datcore.jwstcal(fin, tim, flttype, out, verbose=False)
            pass
        return caled

    @staticmethod
    def _failure(errstr):
        '''Failure log'''
        log.warning('--< DATA CALIBRATION: %s >--', errstr)
        return

    pass


# ---------------- ---------------------------------------------------
