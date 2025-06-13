'''cerberus algorithms ds'''

# Heritage code shame:
# pylint: disable=too-many-arguments,too-many-branches,too-many-locals,too-many-positional-arguments

# -- IMPORTS -- ------------------------------------------------------
import dawgie
import dawgie.context

import numexpr

import logging

import excalibur
import excalibur.system as sys
import excalibur.system.algorithms as sysalg
import excalibur.ancillary as anc
import excalibur.ancillary.algorithms as ancillaryalg
import excalibur.runtime as rtime
import excalibur.runtime.algorithms as rtalg
import excalibur.runtime.binding as rtbind
import excalibur.transit as trn
import excalibur.transit.algorithms as trnalg
from excalibur import ariel
import excalibur.ariel.algorithms as arielalg
import excalibur.cerberus.core as crbcore
import excalibur.cerberus.states as crbstates
from excalibur.util.checksv import checksv

from importlib import import_module as fetch  # avoid cicular dependencies

log = logging.getLogger(__name__)

numexpr.ncores = 1  # this is actually a performance enhancer!

fltrs = [str(fn) for fn in rtbind.filter_names.values()]


# ----------------------- --------------------------------------------
# -- ALGORITHMS -- ---------------------------------------------------
class XSLib(dawgie.Algorithm):
    '''Cross Section Library'''

    def __init__(self):
        '''__init__ ds'''
        self._version_ = crbcore.myxsecsversion()
        self.__spc = trnalg.Spectrum()
        self.__arielsim = arielalg.SimSpectrum()
        self.__rt = rtalg.Autofill()
        self.__out = [crbstates.XslibSv(fltr) for fltr in fltrs]
        return

    def name(self):
        '''Database name for subtask extension'''
        return 'xslib'

    def previous(self):
        '''Input State Vectors: transit.spectrum'''
        return [
            dawgie.ALG_REF(trn.task, self.__spc),
            dawgie.ALG_REF(ariel.task, self.__arielsim),
        ] + self.__rt.refs_for_proceed()

    def state_vectors(self):
        '''Output State Vectors: cerberus.xslib'''
        return self.__out

    def run(self, ds, ps):
        '''Top level algorithm call'''

        svupdate = []
        # just one filter, while debugging:
        # for fltr in ['HST-WFC3-IR-G141-SCAN']:
        for fltr in self.__rt.sv_as_dict()['status']['allowed_filter_names']:
            # stop here if it is not a runtime target
            self.__rt.proceed(fltr)
            update = False

            if fltr == 'Ariel-sim':
                sv = self.__arielsim.sv_as_dict()['parameters']
                vspc, sspc = checksv(sv)
                sspc = 'Ariel-sim spectrum not found'
            elif fltr in self.__spc.sv_as_dict().keys():
                sv = self.__spc.sv_as_dict()[fltr]
                vspc, sspc = checksv(sv)
            else:
                vspc = False
                sspc = 'This filter doesnt have a spectrum: ' + fltr

            if vspc:
                log.warning('--< CERBERUS XSLIB: %s >--', fltr)
                update = self._xslib(sv, fltrs.index(fltr))
            else:
                errstr = [m for m in [sspc] if m is not None]
                self._failure(errstr[0])

            if update:
                svupdate.append(self.__out[fltrs.index(fltr)])
        self.__out = svupdate
        if self.__out:
            _ = excalibur.lagger()
            ds.update()
            pass
        else:
            raise dawgie.NoValidOutputDataError(
                f'No output created for CERBERUS.{self.name()}'
            )
        return

    def _xslib(self, spc, index):
        '''Core code call'''
        cs = crbcore.myxsecs(spc, self.__out[index], verbose=False)
        return cs

    @staticmethod
    def _failure(errstr):
        '''Failure log'''
        log.warning('--< CERBERUS XSLIB: %s >--', errstr)
        return

    pass


class Atmos(dawgie.Algorithm):
    '''Atmospheric retrievial'''

    def __init__(self):
        '''__init__ ds'''
        self._version_ = crbcore.atmosversion()
        self.__spc = trnalg.Spectrum()
        self.__fin = sysalg.Finalize()
        self.__xsl = XSLib()
        self.__arielsim = arielalg.SimSpectrum()
        self.__rt = rtalg.Autofill()
        self.__out = [crbstates.AtmosSv(fltr) for fltr in fltrs]
        return

    def name(self):
        '''Database name for subtask extension'''
        return 'atmos'

    def previous(self):
        '''Input State Vectors: transit.spectrum, system.finalize, cerberus.xslib'''
        return [
            dawgie.ALG_REF(trn.task, self.__spc),
            dawgie.ALG_REF(sys.task, self.__fin),
            dawgie.ALG_REF(fetch('excalibur.cerberus').task, self.__xsl),
            dawgie.ALG_REF(ariel.task, self.__arielsim),
            dawgie.V_REF(
                rtime.task,
                self.__rt,
                self.__rt.sv_as_dict()['status'],
                'cerberus_steps',
            ),
            dawgie.V_REF(
                rtime.task,
                self.__rt,
                self.__rt.sv_as_dict()['status'],
                'cerberus_atmos_fitCloudParameters',
            ),
            dawgie.V_REF(
                rtime.task,
                self.__rt,
                self.__rt.sv_as_dict()['status'],
                'cerberus_atmos_fitT',
            ),
            dawgie.V_REF(
                rtime.task,
                self.__rt,
                self.__rt.sv_as_dict()['status'],
                'cerberus_atmos_fitCtoO',
            ),
            dawgie.V_REF(
                rtime.task,
                self.__rt,
                self.__rt.sv_as_dict()['status'],
                'cerberus_atmos_fitNtoO',
            ),
        ] + self.__rt.refs_for_proceed()

    def state_vectors(self):
        '''Output State Vectors: cerberus.atmos'''
        return self.__out

    def run(self, ds, ps):
        '''Top level algorithm call'''

        vfin, sfin = checksv(self.__fin.sv_as_dict()['parameters'])
        if sfin:
            sfin = 'Missing system params!'

        svupdate = []
        # just one filter, while debugging:
        # for fltr in ['HST-WFC3-IR-G141-SCAN']:
        # for fltr in ['Ariel-sim']:
        for fltr in self.__rt.sv_as_dict()['status']['allowed_filter_names']:
            # stop here if it is not a runtime target
            self.__rt.proceed(fltr)

            update = False
            if fltr in self.__xsl.sv_as_dict():
                vxsl, sxsl = checksv(self.__xsl.sv_as_dict()[fltr])
                if sxsl:
                    sxsl = fltr + ' missing XSL'
            else:
                vxsl, sxsl = (False, fltr + ' missing XSL')

            if fltr == 'Ariel-sim':
                sv = self.__arielsim.sv_as_dict()['parameters']
                vspc, sspc = checksv(sv)
                sspc = 'Ariel-sim spectrum not found'
            elif fltr in self.__spc.sv_as_dict().keys():
                sv = self.__spc.sv_as_dict()[fltr]
                vspc, sspc = checksv(sv)
            else:
                vspc = False
                sspc = 'This filter doesnt have a spectrum: ' + fltr

            if vfin and vxsl and vspc:
                log.warning('--< CERBERUS ATMOS: %s >--', fltr)
                runtime = self.__rt.sv_as_dict()['status']

                runtime_params = crbcore.CerbAtmosParams(
                    MCMC_chain_length=runtime['cerberus_steps'].value(),
                    MCMC_chains=runtime['cerberus_chains'].value(),
                    MCMC_sliceSampler=runtime['cerberus_atmos_sliceSampler'],
                    fitCloudParameters=runtime[
                        'cerberus_atmos_fitCloudParameters'
                    ],
                    fitT=runtime['cerberus_atmos_fitT'],
                    fitCtoO=runtime['cerberus_atmos_fitCtoO'],
                    fitNtoO=runtime['cerberus_atmos_fitNtoO'],
                    nlevels=runtime['cerberus_atmos_crbmodel_nlevels'].value(),
                    solrad=runtime['cerberus_atmos_crbmodel_solrad'].value(),
                    Hsmax=runtime['cerberus_atmos_crbmodel_Hsmax'].value(),
                    lbroadening=runtime['cerberus_atmos_crbmodel_lbroadening'],
                    lshifting=runtime['cerberus_atmos_crbmodel_lshifting'],
                    isothermal=runtime['cerberus_atmos_crbmodel_isothermal'],
                    boundTeq=runtime['cerberus.atmos.bounds.Teq'],
                    boundAbundances=runtime['cerberus.atmos.bounds.abundances'],
                    boundCTP=runtime['cerberus.atmos.bounds.CTP'],
                    boundHLoc=runtime['cerberus.atmos.bounds.HLoc'],
                    boundHScale=runtime['cerberus.atmos.bounds.HScale'],
                    boundHThick=runtime['cerberus.atmos.bounds.HThick'],
                )
                # print('runtime params',runtime_params)
                update = self._atmos(
                    self.__fin.sv_as_dict()['parameters'],
                    self.__xsl.sv_as_dict()[fltr],
                    sv,
                    runtime_params,
                    fltrs.index(fltr),
                    fltr,
                )
            else:
                errstr = [m for m in [sfin, sspc, sxsl] if m is not None]
                self._failure(errstr[0])
            if update:
                svupdate.append(self.__out[fltrs.index(fltr)])
        self.__out = svupdate
        if self.__out:
            _ = excalibur.lagger()
            ds.update()
            pass
        else:
            raise dawgie.NoValidOutputDataError(
                f'No output created for CERBERUS.{self.name()}'
            )
        return

    def _atmos(self, fin, xsl, spc, runtime_params, index, fltr):
        '''Core code call'''

        mcmc_chains = runtime_params.MCMC_chains
        mcmc_chain_length = runtime_params.MCMC_chain_length
        # print('MCMC_chain_length', mcmc_chain_length)
        # mcmc_chain_length = 1000
        # mcmc_chain_length = 10
        # print('MCMC_chain_length', mcmc_chain_length)
        log.info(
            ' calling atmos from cerb-alg-atmos  chain len=%d',
            mcmc_chain_length,
        )
        am = crbcore.atmos(
            fin,
            xsl,
            spc,
            runtime_params,
            self.__out[index],
            fltr,
            Nchains=mcmc_chains,
            chainlen=mcmc_chain_length,
            verbose=False,
        )  # singlemod='TEC' after chainlen
        return am

    @staticmethod
    def _failure(errstr):
        '''Failure log'''
        log.warning('--< CERBERUS ATMOS: %s >--', errstr)
        return

    pass


# ---------------- ---------------------------------------------------


class Results(dawgie.Algorithm):
    '''
    Plot the best-fit spectrum, to see how well it fits the data
    Plot the corner plot, to see how well each parameter is constrained
    '''

    def __init__(self):
        '''__init__ ds'''
        self._version_ = crbcore.resultsversion()
        self.__fin = sysalg.Finalize()
        self.__anc = ancillaryalg.Estimate()
        self.__xsl = XSLib()
        self.__atm = Atmos()
        self.__rt = rtalg.Autofill()
        self.__out = [crbstates.ResSv(fltr) for fltr in fltrs]
        return

    def name(self):
        '''Database name for subtask extension'''
        return 'results'

    def previous(self):
        '''Input State Vectors: cerberus.atmos'''
        return [
            dawgie.ALG_REF(sys.task, self.__fin),
            dawgie.ALG_REF(anc.task, self.__anc),
            dawgie.ALG_REF(fetch('excalibur.cerberus').task, self.__xsl),
            dawgie.ALG_REF(fetch('excalibur.cerberus').task, self.__atm),
        ] + self.__rt.refs_for_proceed()

    def state_vectors(self):
        '''Output State Vectors: cerberus.results'''
        return self.__out

    def run(self, ds, ps):
        '''Top level algorithm call'''

        svupdate = []
        vfin, sfin = checksv(self.__fin.sv_as_dict()['parameters'])
        vanc, sanc = checksv(self.__anc.sv_as_dict()['parameters'])

        update = False
        if vfin and vanc:
            # available_filters = self.__xsl.sv_as_dict().keys()
            # available_filters = self.__atm.sv_as_dict().keys()
            # print('available_filters',available_filters)
            # allowed_filters = self.__rt.sv_as_dict()['status']['allowed_filter_names']
            # print('allowed filters in cerb.results',allowed_filters)

            # just one filter, while debugging:
            # for fltr in ['HST-WFC3-IR-G141-SCAN']:
            # for fltr in ['Ariel-sim']:
            for fltr in self.__rt.sv_as_dict()['status'][
                'allowed_filter_names'
            ]:
                # stop here if it is not a runtime target
                self.__rt.proceed(fltr)

                vxsl, sxsl = checksv(self.__xsl.sv_as_dict()[fltr])
                vatm, satm = checksv(self.__atm.sv_as_dict()[fltr])
                if vxsl and vatm:
                    log.warning('--< CERBERUS RESULTS: %s >--', fltr)
                    # FIXMEE: this code needs repaired by moving out to config (Geoff added)

                    runtime = self.__rt.sv_as_dict()['status']

                    runtime_params = crbcore.CerbResultsParams(
                        nrandomwalkers=runtime['cerberus_results_nrandomwalkers'].value(),
                        randomseed=runtime['cerberus_results_randomseed'].value(),
                    )

                    update = self._results(
                        repr(self).split('.')[1],  # this is the target name
                        fltr,
                        runtime_params,
                        self.__fin.sv_as_dict()['parameters'],
                        self.__anc.sv_as_dict()['parameters'],
                        self.__xsl.sv_as_dict()[fltr]['data'],
                        self.__atm.sv_as_dict()[fltr]['data'],
                        fltrs.index(fltr),
                    )
                    if update:
                        svupdate.append(self.__out[fltrs.index(fltr)])
                else:
                    errstr = [m for m in [sxsl, satm] if m is not None]
                    self._failure(errstr[0])
        else:
            errstr = [m for m in [sanc, sfin] if m is not None]
            self._failure(errstr[0])

        self.__out = svupdate
        if self.__out:
            _ = excalibur.lagger()
            ds.update()
            pass
        else:
            raise dawgie.NoValidOutputDataError(
                f'No output created for CERBERUS.{self.name()}'
            )
        return

    def _results(self, trgt, fltr, fin, ancil, xsl, atm, index):
        '''Core code call'''
        resout = crbcore.results(
            trgt, fltr, fin, ancil, xsl, atm, self.__out[index], verbose=False
        )
        return resout

    @staticmethod
    def _failure(errstr):
        '''Failure log'''
        log.warning('--< CERBERUS RESULTS: %s >--', errstr)
        return

    pass


# ---------------- ---------------------------------------------------


class Analysis(dawgie.Analyzer):
    '''analysis ds'''

    def __init__(self):
        '''__init__ ds'''
        self._version_ = (
            crbcore.resultsversion()
        )  # same version number as results
        # self.__fin = sysalg.finalize()
        # self.__xsl = xslib()
        # self.__atm = atmos()
        # self.__out = crbstates.AnalysisSv('retrievalCheck')
        self.__out = [crbstates.AnalysisSv(fltr) for fltr in fltrs]
        return

    # def previous(self):
    #    '''Input State Vectors: cerberus.atmos'''
    #        return [dawgie.ALG_REF(sys.task, self.__fin)]

    def feedback(self):
        '''feedback ds'''
        return []

    def name(self):
        '''Database name for subtask extension'''
        return 'analysis'

    def traits(self) -> [dawgie.SV_REF, dawgie.V_REF]:
        '''traits ds'''
        return [
            dawgie.SV_REF(fetch('excalibur.cerberus').task, Atmos(), sv)
            for sv in Atmos().state_vectors()
        ]

    def state_vectors(self):
        '''Output State Vectors: cerberus.analysis'''
        return self.__out

    def run(self, aspects: dawgie.Aspect):
        '''Top level algorithm call'''

        svupdate = []
        if len(aspects) == 0:
            log.warning('--< CERBERUS ANALYSIS: contains no targets >--')
        else:
            # determine which filters have results from cerb.atmos (in aspects)
            #  (you have to loop through all targets, since filters vary by target)
            fwr = []
            for trgt in aspects:
                for fltr in fltrs:
                    if (fltr not in fwr) and (
                        'cerberus.atmos.' + fltr in aspects[trgt]
                    ):
                        # print('This filter exists in the cerb.atmos aspect:',fltr,trgt)
                        fwr.append(fltr)
            if not fwr:
                log.warning(
                    '--< CERBERUS ANALYSIS: NO FILTERS WITH ATMOS DATA!!!>--'
                )

            # filtersWithResults=['Ariel-sim']  # just one filter, while debugging
            # filtersWithResults=['HST-WFC3-IR-G141-SCAN']  # just one filter, while debugging

            # only consider filters that have cerb.atmos results loaded in as an aspect
            for fltr in fwr:
                # if 'cerberus.atmos.'+fltr not in aspects[trgt]:
                #    log.warning('--< CERBERUS ANALYSIS: %s not found IMPOSSIBLE!!!!>--', fltr)
                # else:
                log.warning('--< CERBERUS ANALYSIS: %s  >--', fltr)
                update = self._analysis(aspects, fltr, fltrs.index(fltr))
                if update:
                    svupdate.append(self.__out[fltrs.index(fltr)])
        self.__out = svupdate
        if self.__out:
            aspects.ds().update()
        else:
            raise dawgie.NoValidOutputDataError(
                f'No output created for CERBERUS.{self.name()}'
            )
        return

    def _analysis(self, aspects, fltr, index):
        '''Core code call'''
        analysisout = crbcore.analysis(
            aspects, fltr, self.__out[index], verbose=False
        )
        return analysisout

    @staticmethod
    def _failure(errstr):
        '''Failure log'''
        log.warning('--< CERBERUS ANALYSIS: %s >--', errstr)
        return


# -------------------------------------------------------------------
class Release(dawgie.Algorithm):
    '''Format release products Roudier et al. 2021'''

    def __init__(self):
        '''__init__ ds'''
        self._version_ = crbcore.rlsversion()
        self.__fin = sysalg.Finalize()
        # self.__atm = Atmos()
        self.__out = [crbstates.RlsSv(fltr) for fltr in fltrs]
        return

    def name(self):
        '''Database name for subtask extension'''
        return 'release'

    def previous(self):
        '''Input State Vectors: cerberus.atmos'''
        return []
        #    dawgie.ALG_REF(sys.task, self.__fin),
        #    dawgie.ALG_REF(fetch('excalibur.cerberus').task, self.__atm),
        # ]

    def state_vectors(self):
        '''Output State Vectors: cerberus.release'''
        return self.__out

    def run(self, ds, ps):
        '''Top level algorithm call'''
        svupdate = []
        vfin, sfin = checksv(self.__fin.sv_as_dict()['parameters'])
        fltr = 'HST-WFC3-IR-G141-SCAN'
        update = False
        if vfin:
            log.warning(
                '--< CERBERUS RELEASE: %s %s >--',
                fltr,
                repr(self).split('.')[1],
            )
            update = self._release(
                repr(self).split('.')[1],  # this is the target name
                self.__fin.sv_as_dict()['parameters'],
                fltrs.index(fltr),
            )
            pass
        else:
            errstr = [m for m in [sfin] if m is not None]
            self._failure(errstr[0])
            pass
        if update:
            svupdate.append(self.__out[fltrs.index(fltr)])
        self.__out = svupdate
        if self.__out:
            _ = excalibur.lagger()
            ds.update()
            pass
        else:
            raise dawgie.NoValidOutputDataError(
                f'No output created for CERBERUS.{self.name()}'
            )
        return

    def _release(self, trgt, fin, index):
        '''Core code call'''
        rlsout = crbcore.release(trgt, fin, self.__out[index], verbose=False)
        return rlsout

    @staticmethod
    def _failure(errstr):
        '''Failure log'''
        log.warning('--< CERBERUS RELEASE: %s >--', errstr)
        return

    pass


# ---------------- ---------------------------------------------------
