'''target algorithms ds'''

# -- IMPORTS -- ------------------------------------------------------
import os
import logging

import dawgie
import dawgie.context

import excalibur

import excalibur.runtime as rtime
import excalibur.runtime.algorithms as rtalg
import excalibur.runtime.binding as rtbind

from excalibur.util.checksv import checksv

import excalibur.target.core as trgcore
import excalibur.target.edit as trgedit
import excalibur.target.monitor as trgmonitor
import excalibur.target.states as trgstates

from importlib import import_module as fetch  # avoid cicular dependencies

log = logging.getLogger(__name__)

fltrs = [str(fn) for fn in rtbind.filter_names.values()]

# ------------- ------------------------------------------------------
# -- ALGO RUN OPTIONS -- ---------------------------------------------
# GENERATE DATABASE IDs
GEN_IDS = True
# NEXSCI QUERY
WEB = 'https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query='
# DATA ON DISK
diskloc = os.path.join(excalibur.context['data_dir'], 'sci')
# MAST MIRRORS
QUERY_FORM = 'https://archive.stsci.edu/hst/search.php?target='
MIRROR_1 = 'http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub/HSTCA/'
MIRROR_2 = 'http://archives.esac.esa.int/ehst-sl-server/servlet/data-action?ARTIFACT_ID='
# MAST API
DURL = 'https://mast.stsci.edu/api/v0.1/Download/file?'
# HSTURL = 'https://mast.stsci.edu/search/hst/api/v0.1/retrieve_product?'


# ---------------------- ---------------------------------------------
# -- ALGORITHMS -- ---------------------------------------------------
class Alert(dawgie.Analyzer):
    '''alert ds'''

    def __init__(self):
        '''version 1.2.0 has new 'web' url for Exoplanet Archive queries'''
        self._version_ = dawgie.VERSION(1, 2, 0)
        self.__out = trgstates.AlertSV()
        return

    def feedback(self):
        '''feedback ds'''
        return [
            dawgie.V_REF(
                fetch('excalibur.target').analysis, self, self.__out, 'known'
            ),
            dawgie.V_REF(
                fetch('excalibur.target').analysis, self, self.__out, 'table'
            ),
        ]

    def name(self):
        '''Database name for subtask extension'''
        return 'alert_from_variations_of'

    def traits(self) -> [dawgie.SV_REF, dawgie.V_REF]:
        '''traits ds'''
        return [
            dawgie.V_REF(
                fetch('excalibur.target').regress,
                Regress(),
                Regress().state_vectors()[0],
                'last',
            )
        ]

    def state_vectors(self):
        '''Output State Vectors: target.alert_from_variations_of'''
        return [self.__out]

    def run(self, aspects: dawgie.Aspect):
        '''Top level algorithm call'''
        c, k, t = trgmonitor.alert(
            aspects, self.__out['known'], self.__out['table']
        )
        self.__out['changes'].clear()
        self.__out['known'].clear()
        self.__out['table'].clear()
        self.__out['changes'].extend(c)
        self.__out['known'].extend(k)
        self.__out['table'].extend(t)
        excalibur.lagger()
        aspects.ds().update()
        return

    pass


class Create(dawgie.Analyzer):
    '''Creates a list of targets from edit.py'''

    def __init__(self):
        '''__init__ ds'''
        self._version_ = trgedit.createversion()
        self.__rt = rtalg.Create()
        self.__out = [
            trgstates.TargetSV('starIDs'),
            trgstates.FilterSV('filters'),
        ]
        return

    def name(self):
        '''Database name for subtask extension'''
        return 'create'

    def traits(self) -> [dawgie.SV_REF, dawgie.V_REF]:
        '''Aspect Input: None (it s the clutch)'''
        return []

    def state_vectors(self):
        '''Output State Vectors: target.create'''
        return self.__out

    def run(self, aspects: dawgie.Aspect):
        '''Top level algorithm call'''

        runtime = self.__rt.sv_as_dict()['controls']
        # print('runtime', runtime)
        runtime_params = trgcore.TargetCreateParams(
            num_reruns=runtime['testcerb_Nrepeats'].value(),
        )
        # print('runtime_params in target.create alg', runtime_params)

        trgcore.scrapeids(
            aspects.ds(), runtime_params, self.__out[0], WEB, gen_ids=GEN_IDS
        )
        update = trgcore.createfltrs(self.__out[1])
        if update:
            excalibur.lagger()
            aspects.ds().update()
        else:
            raise dawgie.NoValidOutputDataError(
                f'No output created for TARGET.{self.name()}'
            )
        return

    pass


class Autofill(dawgie.Algorithm):
    '''Fills mandatory info to get the ball rolling'''

    def __init__(self):
        '''__init__ ds'''
        self._version_ = trgcore.autofillversion()
        self.__create = Create()
        self.__rt = rtalg.Autofill()
        self.__out = trgstates.TargetSV('parameters')
        return

    def name(self):
        '''Database name for subtask extension'''
        return 'autofill'

    def previous(self):
        '''Input State Vectors: target.create'''
        return [
            dawgie.ALG_REF(fetch('excalibur.target').analysis, self.__create),
            dawgie.V_REF(
                rtime.task,
                self.__rt,
                self.__rt.sv_as_dict()['status'],
                'isValidTarget',
            ),
        ]

    def state_vectors(self):
        '''Output State Vectors: target.autofill'''
        return [self.__out]

    def run(self, ds, ps):
        '''Top level algorithm call'''
        update = False

        # FIXMEE: this code needs repaired by moving out to config
        target = repr(self).split('.')[1]

        # stop here if it is not a runtime target
        if not self.__rt.is_valid():
            log.warning(
                '--< TARGET.%s: %s not a valid target >--',
                target,
                self.name().upper(),
            )
            pass
        else:
            crt = self.__create.sv_as_dict()
            valid, errstring = checksv(crt['starIDs'])
            if valid and (target in crt['starIDs']['starID']):
                log.warning('--< TARGET AUTOFILL: %s >--', target)
                update = self._autofill(crt, target)
                pass
            else:
                self._failure(errstring)
                pass
            if update:
                _ = excalibur.lagger()
                ds.update()
                pass
            else:
                raise dawgie.NoValidOutputDataError(
                    f'!!! No SV created for TARGET.{self.name()} !!!'
                )
            pass
        return

    def _autofill(self, crt, thistarget):
        '''Core code call'''
        # currently we are running this on all filters, not just the available ones
        solved = trgcore.autofill(crt, thistarget, self.__out, fltrs)
        return solved

    @staticmethod
    def _failure(errstr):
        '''Failure log'''
        if errstr is None:
            errstr = 'TARGET NOT EXPECTED'
            pass
        log.warning('--< TARGET AUTOFILL: %s >--', errstr)
        return

    pass


class Scrape(dawgie.Algorithm):
    '''
    Download data or ingest data from disk
    '''

    def __init__(self):
        '''__init__ ds'''
        self._version_ = trgcore.scrapeversion()
        self.__autofill = Autofill()
        self.__rt = rtalg.Autofill()
        self.__out = trgstates.DatabaseSV('databases')
        return

    def name(self):
        '''Database name for subtask extension'''
        return 'scrape'

    def previous(self):
        '''Input State Vectors: target.autofill'''
        return [
            dawgie.ALG_REF(fetch('excalibur.target').task, self.__autofill),
            dawgie.V_REF(
                rtime.task,
                self.__rt,
                self.__rt.sv_as_dict()['status'],
                'isValidTarget',
            ),
        ]

    def state_vectors(self):
        '''Output State Vectors: target.scrape'''
        return [self.__out]

    def run(self, ds, ps):
        '''Top level algorithm call'''

        # stop here if it is not a runtime target
        if not self.__rt.is_valid():
            log.warning(
                '--< TARGET.%s: not a valid target >--', self.name().upper()
            )
        else:
            var_autofill = self.__autofill.sv_as_dict()['parameters']
            valid, errstring = checksv(var_autofill)
            if valid:
                log.warning(
                    '--< TARGET SCRAPE: %s >--', repr(self).split('.')[1]
                )
                self._scrape(var_autofill, self.__out)
            else:
                self._failure(errstring)
            # GMR: always update.
            # Sims / proceed() do not require data nor full set of system parameters.
            _ = excalibur.lagger()
            ds.update()
        return

    @staticmethod
    def _scrape(tfl, out):
        '''Core code call'''
        dbs = os.path.join(dawgie.context.data_dbs, 'mast')
        if not os.path.exists(dbs):
            os.makedirs(dbs)
            pass
        # Download from MAST
        umast = trgcore.mastapi(
            tfl,
            out,
            dbs,
            download_url=DURL,
            hst_url=DURL,
            verbose=False,
        )
        # Data on DISK
        # udisk gets prioritized over umast for duplicates
        # Ditching that. We have no idea what we have on disk right now
        # udisk = trgcore.disk(tfl, out, diskloc, dbs)
        udisk = False
        return udisk or umast

    @staticmethod
    def _failure(errstr):
        '''Failure log'''
        log.warning('--< TARGET SCRAPE: %s >--', errstr)
        return

    pass


class Regress(dawgie.Regression):
    '''regress ds'''

    def __init__(self):
        '''__init__ ds'''
        self._version_ = dawgie.VERSION(1, 0, 2)
        self.__out = trgstates.MonitorSV()
        return

    def feedback(self):
        '''feedback ds'''
        return [
            dawgie.V_REF(
                fetch('excalibur.target').regress, self, self.__out, 'planet'
            ),
            dawgie.V_REF(
                fetch('excalibur.target').regress, self, self.__out, 'runid'
            ),
        ]

    def name(self):
        '''Database name for subtask extension'''
        return 'variations_of'

    def run(self, ps: int, timeline: dawgie.Timeline):
        '''Top level algorithm call'''
        last, outlier = trgmonitor.regress(
            self.__out['planet'], self.__out['runid'], timeline
        )
        self.__out['last'].clear()
        self.__out['last'].update(last)
        self.__out['outlier'].clear()
        self.__out['outlier'].extend(outlier)
        timeline.ds().update()
        return

    def state_vectors(self):
        '''Output State Vectors: target.variations_of'''
        return [self.__out]

    def variables(self) -> [dawgie.SV_REF, dawgie.V_REF]:
        '''variables ds'''
        return [
            dawgie.SV_REF(
                fetch('excalibur.target').task,
                Autofill(),
                Autofill().state_vectors()[0],
            )
        ]

    pass


# -------------------------------------------------------------------
