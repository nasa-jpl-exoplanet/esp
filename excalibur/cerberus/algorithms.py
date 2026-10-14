'''cerberus algorithms ds'''

# Heritage code shame:
# pylint: disable=too-many-arguments,too-many-branches,too-many-locals,too-many-positional-arguments,too-many-statements,too-many-nested-blocks

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
import excalibur.runtime.algorithms as rtalg
import excalibur.runtime.binding as rtbind
import excalibur.transit as trn
import excalibur.transit.algorithms as trnalg
from excalibur import ariel
import excalibur.ariel.algorithms as arielalg
import excalibur.cerberus.core as crbcore
import excalibur.cerberus.states as crbstates
from excalibur.util.checksv import checksv

from excalibur.target.targetlists import get_target_lists

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

        target = repr(self).split('.')[1]
        svupdate = []

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

            runtime = self.__rt.sv_as_dict()['status']

            # for Ariel targets, option to only do the actually Tier-2 targets
            targetlistcheck = True
            only_these_planets = []
            if (
                fltr == 'Ariel-sim'
                and runtime['cerberus_arielsample_tier'].value() == 2
            ):
                alltargetlists = get_target_lists()
                targetlist = alltargetlists['ariel_stars_tier2']
                if target not in targetlist:
                    targetlistcheck = False

                if targetlistcheck:
                    planetlist = alltargetlists['ariel_planets_tier2']
                    for planet in planetlist:
                        if planet.startswith(target + ' '):
                            only_these_planets.append(planet[-1])
                # print('only these planets', only_these_planets)

            if vspc and targetlistcheck:
                log.info('--< CERBERUS XSLIB: %s  %s >--', fltr, target)
                update = self._xslib(sv, runtime, only_these_planets, fltr)
            else:
                if targetlistcheck:
                    errstr = [m for m in [sspc] if m is not None]
                else:
                    errstr = ['not in the Ariel target list']
                self._failure(errstr[0], target)

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

    def _xslib(self, spc, runtime, only_these_planets, fltr):
        '''Core code call'''
        if 'JWST' in fltr:
            cs = crbcore.jwstwxs(
                spc,
                runtime,
                self.__out[fltrs.index(fltr)],
                verbose=False,
            )
            pass
        else:
            cs = crbcore.myxsecs(
                spc,
                runtime,
                self.__out[fltrs.index(fltr)],
                only_these_planets=only_these_planets,
                verbose=False,
            )
            pass
        return cs

    @staticmethod
    def _failure(errstr, target):
        '''Failure log'''
        log.warning('--< CERBERUS XSLIB: %s  %s >--', errstr, target)
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
        return (
            [
                dawgie.ALG_REF(trn.task, self.__spc),
                dawgie.ALG_REF(sys.task, self.__fin),
                dawgie.ALG_REF(fetch('excalibur.cerberus').task, self.__xsl),
                dawgie.ALG_REF(ariel.task, self.__arielsim),
            ]
            + self.__rt.trigger('cerberus')
            + self.__rt.refs_for_proceed()
        )

    def state_vectors(self):
        '''Output State Vectors: cerberus.atmos'''
        return self.__out

    def run(self, ds, ps):
        '''Top level algorithm call'''

        target = repr(self).split('.')[1]

        vfin, sfin = checksv(self.__fin.sv_as_dict()['parameters'])
        if sfin:
            sfin = 'Missing system params!'

        runtime = self.__rt.sv_as_dict()['status']

        svupdate = []
        # for fltr in ['Ariel-sim']:
        for fltr in self.__rt.sv_as_dict()['status']['allowed_filter_names']:
            # stop here if it is not a runtime target
            self.__rt.proceed(fltr)

            update = False
            # XSLIB CHECK
            if fltr in self.__xsl.sv_as_dict():
                vxsl, sxsl = checksv(self.__xsl.sv_as_dict()[fltr])
                if sxsl:
                    sxsl = fltr + ' missing XSL'
            else:
                vxsl, sxsl = (False, fltr + ' missing XSL')
                pass
            # INPUT SPECTRUM CHECK
            if fltr == 'Ariel-sim':
                sv = self.__arielsim.sv_as_dict()['parameters']
                vspc, sspc = checksv(sv)
                sspc = 'Ariel-sim spectrum not found'
                pass
            elif fltr in self.__spc.sv_as_dict().keys():
                sv = self.__spc.sv_as_dict()[fltr]
                if 'data' in sv and 'target' not in sv['data']:
                    sv['data']['target'] = 'HST spectrum needs target name'
                    pass
                vspc, sspc = checksv(sv)
                pass
            else:
                vspc = False
                sspc = 'This filter doesnt have a spectrum: ' + fltr
                pass

            # for Ariel targets, option to only do the Tier-2 targets
            targetlistcheck = True
            only_these_planets = []
            if (
                fltr == 'Ariel-sim'
                and runtime['cerberus_arielsample_tier'].value() == 2
            ):
                alltargetlists = get_target_lists()
                targetlist = alltargetlists['ariel_stars_tier2']
                if target not in targetlist:
                    targetlistcheck = False

                if targetlistcheck:
                    planetlist = alltargetlists['ariel_planets_tier2']
                    for planet in planetlist:
                        if planet.startswith(target + ' '):
                            only_these_planets.append(planet[-1])

            if vfin and vxsl and vspc and targetlistcheck:
                log.info('--< CERBERUS ATMOS: %s  %s >--', fltr, target)

                update = self._atmos(
                    self.__fin.sv_as_dict()['parameters'],
                    self.__xsl.sv_as_dict()[fltr],
                    sv,
                    runtime,
                    only_these_planets,
                    fltr,
                )
                pass
            else:
                if targetlistcheck:
                    errstr = [m for m in [sfin, sspc, sxsl] if m is not None]
                    pass
                else:
                    errstr = ['not in the Ariel target list']
                    pass
                self._failure(errstr[0], target)
                pass
            if update:
                svupdate.append(self.__out[fltrs.index(fltr)])
                pass
            pass
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

    def _atmos(self, fin, xsl, spc, runtime, only_these_planets, fltr):
        '''
        Core code call
        '''
        log.info(
            '--< CERBERUS ATMOS: Chain length %d >--',
            runtime['cerberus_steps'].value(),
        )
        if 'JWST' in fltr:
            am = crbcore.jwstatmos(
                self.__fin.sv_as_dict()['parameters'],
                xsl,
                spc,
                self.__rt.sv_as_dict()['status'],
                self.__out[fltrs.index(fltr)],
                verbose=False,
            )
            pass
        else:
            am = crbcore.atmos(
                fin,
                xsl,
                spc,
                runtime,
                self.__out[fltrs.index(fltr)],
                fltr,
                only_these_planets=only_these_planets,
                Nchains=runtime['cerberus_chains'].value(),
                chainlen=runtime['cerberus_steps'].value(),
                verbose=False,
            )
            pass
        return am

    @staticmethod
    def _failure(errstr, target):
        '''Failure log'''
        log.warning('--< CERBERUS ATMOS: %s  %s >--', errstr, target)
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

        target = repr(self).split('.')[1]

        svupdate = []
        vfin, sfin = checksv(self.__fin.sv_as_dict()['parameters'])
        vanc, sanc = checksv(self.__anc.sv_as_dict()['parameters'])

        update = False
        if vfin and vanc:
            runtime = self.__rt.sv_as_dict()['status']

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

                # for Ariel targets, option to only do the actually Tier-2 targets
                targetlistcheck = True
                only_these_planets = []
                if (
                    fltr == 'Ariel-sim'
                    and runtime['cerberus_arielsample_tier'].value() == 2
                ):
                    alltargetlists = get_target_lists()
                    targetlist = alltargetlists['ariel_stars_tier2']
                    if target not in targetlist:
                        targetlistcheck = False

                    if targetlistcheck:
                        planetlist = alltargetlists['ariel_planets_tier2']
                        for planet in planetlist:
                            if planet.startswith(target + ' '):
                                only_these_planets.append(planet[-1])
                    # print('only these planets', only_these_planets)

                if vxsl and vatm and targetlistcheck:
                    log.info('--< CERBERUS RESULTS: %s  %s >--', fltr, target)

                    update = self._results(
                        repr(self).split('.')[1],  # this is the target name
                        fltr,
                        runtime,
                        only_these_planets,
                        self.__fin.sv_as_dict()['parameters'],
                        self.__anc.sv_as_dict()['parameters'],
                        self.__xsl.sv_as_dict()[fltr]['data'],
                        self.__atm.sv_as_dict()[fltr]['data'],
                        fltrs.index(fltr),
                    )
                    if update:
                        svupdate.append(self.__out[fltrs.index(fltr)])
                else:
                    if targetlistcheck:
                        errstr = [m for m in [sxsl, satm] if m is not None]
                    else:
                        errstr = ['not in the Ariel target list']
                    self._failure(errstr[0], target)
        else:
            errstr = [m for m in [sanc, sfin] if m is not None]
            self._failure(errstr[0], target)

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

    def _results(
        self,
        trgt,
        fltr,
        runtime,
        only_these_planets,
        fin,
        ancil,
        xsl,
        atm,
        index,
    ):
        '''Core code call'''
        resout = crbcore.results(
            trgt,
            fltr,
            runtime,
            fin,
            ancil,
            xsl,
            atm,
            self.__out[index],
            only_these_planets=only_these_planets,
            verbose=False,
        )
        return resout

    @staticmethod
    def _failure(errstr, target):
        '''Failure log'''
        log.warning('--< CERBERUS RESULTS: %s  %s >--', errstr, target)
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
        self.__rt = rtalg.Autofill()
        # self.__rtc = rtalg.Create()
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

            # fwr = ['Ariel-sim']  # just one filter, while debugging
            # fwr =['HST-WFC3-IR-G141-SCAN']  # just one filter, while debugging

            # only consider filters that have cerb.atmos results loaded in as an aspect
            for fltr in fwr:
                # if 'cerberus.atmos.'+fltr not in aspects[trgt]:
                #    log.warning('--< CERBERUS ANALYSIS: %s not found IMPOSSIBLE!!!!>--', fltr)

                # (asdf: this is still not working)
                runtime = self.__rt.sv_as_dict()['status']
                # print('runtime old way',runtime)
                # runtime2 = self.__rtc.sv_as_dict()['status']
                # print('runtime old2 way',runtime2)

                # this might actually work now that it's not a tuple
                # if runtime['ariel_simspectrum_tier'].value() == None:
                #    runtime['ariel_simspectrum_tier'].value() = 2
                # print('runtime', runtime)

                log.info('--< CERBERUS ANALYSIS: %s  >--', fltr)
                update = self._analysis(
                    aspects, fltr, runtime, fltrs.index(fltr)
                )
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

    def _analysis(self, aspects, fltr, runtime, index):
        '''Core code call'''
        analysisout = crbcore.analysis(
            aspects, fltr, runtime, self.__out[index], verbose=False
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
            log.info(
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
