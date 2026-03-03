'''gemli algorithms ds'''

# after updates are made, check if there's still duplicate code!
# Heritage code shame:
# pylint: disable=duplicate-code
# pylint: disable=too-many-arguments,too-many-branches,too-many-locals,too-many-positional-arguments,too-many-statements,too-many-nested-blocks,too-many-instance-attributes

# -- IMPORTS -- ------------------------------------------------------
import dawgie
import dawgie.context

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
import excalibur.cerberus as crb
import excalibur.cerberus.algorithms as crbalg
import excalibur.cerberus.core as crbcore
import excalibur.gemli.core as gemlicore
import excalibur.gemli.states as gemlistates
from excalibur.util.checksv import checksv

from excalibur.target.targetlists import get_target_lists

from importlib import import_module as fetch  # avoid cicular dependencies

log = logging.getLogger(__name__)

fltrs = [str(fn) for fn in rtbind.filter_names.values()]


# ---------------- ---------------------------------------------------


class MLfit(dawgie.Algorithm):
    '''
    Determine the best fit spectrum via Machine Learning
    Plot the best-fit spectrum, to see how well it fits the data
    '''

    def __init__(self):
        '''__init__ ds'''
        self._version_ = gemlicore.mlfitversion()
        self.__fin = sysalg.Finalize()
        self.__anc = ancillaryalg.Estimate()
        self.__xsl = crbalg.XSLib()
        self.__atm = crbalg.Atmos()
        self.__spc = trnalg.Spectrum()
        self.__arielsim = arielalg.SimSpectrum()
        self.__rt = rtalg.Autofill()
        self.__out = [gemlistates.MLfitSv(fltr) for fltr in fltrs]
        return

    def name(self):
        '''Database name for subtask extension'''
        return 'mlfit'

    def previous(self):
        '''Input State Vectors: gemli.atmos'''
        return [
            dawgie.ALG_REF(sys.task, self.__fin),
            dawgie.ALG_REF(anc.task, self.__anc),
            dawgie.ALG_REF(crb.task, self.__xsl),
            dawgie.ALG_REF(crb.task, self.__atm),
            dawgie.ALG_REF(trn.task, self.__spc),
            dawgie.ALG_REF(ariel.task, self.__arielsim),
        ] + self.__rt.refs_for_proceed()

    def state_vectors(self):
        '''Output State Vectors: gemli.mlfit'''
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
            runtime_params = crbcore.CerbResultsParams(
                nrandomwalkers=runtime[
                    'cerberus_results_nrandomwalkers'
                ].value(),
                randomseed=runtime['cerberus_results_randomseed'].value(),
                knownspecies=runtime[
                    'cerberus_crbmodel_HITEMPmolecules'
                ].molecules,
                cialist=runtime['cerberus_crbmodel_HITRANmolecules'].molecules,
                xmollist=runtime['cerberus_crbmodel_EXOMOLmolecules'].molecules,
                nlevels=runtime['cerberus_crbmodel_nlevels'].value(),
                Hsmax=runtime['cerberus_crbmodel_Hsmax'].value(),
                solrad=runtime['cerberus_crbmodel_solrad'].value(),
                cornerBins=runtime['cerberus_plotters_cornerBins'].value(),
                lbroadening=runtime['cerberus_crbmodel_lbroadening'],
                lshifting=runtime['cerberus_crbmodel_lshifting'],
                isothermal=runtime['cerberus_crbmodel_isothermal'],
            )

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

                if fltr == 'Ariel-sim':
                    spectrum = self.__arielsim.sv_as_dict()['parameters']
                    vspc, sspc = checksv(spectrum)
                    sspc = 'Ariel-sim spectrum not found'
                elif fltr in self.__spc.sv_as_dict().keys():
                    spectrum = self.__spc.sv_as_dict()[fltr]
                    if 'data' in spectrum and 'target' not in spectrum['data']:
                        spectrum['data'][
                            'target'
                        ] = 'HST spectrum needs target name'
                    vspc, sspc = checksv(spectrum)
                else:
                    vspc = False
                    sspc = 'This filter doesnt have a spectrum: ' + fltr

                # for Ariel targets, option to only do the actually Tier-2 targets
                targetlistcheck = True
                only_these_planets = []
                if (
                    fltr == 'Ariel-sim'
                    and runtime['cerberus_arielsample_tier'].value() == 2
                ):
                    alltargetlists = get_target_lists()
                    targetlist = alltargetlists['ariel_Nov2024_2years']
                    if target not in targetlist:
                        targetlistcheck = False

                    if targetlistcheck:
                        planetlist = alltargetlists[
                            'ariel_Nov2024_2years_withPlanetletters'
                        ]
                        for planet in planetlist:
                            if planet.startswith(target + ' '):
                                only_these_planets.append(planet[-1])
                    # print('only these planets', only_these_planets)

                if vxsl and vatm and vspc and targetlistcheck:
                    log.info('--< GEMLI MLFIT: %s  %s >--', fltr, target)

                    update = self._mlfit(
                        repr(self).split('.')[1],  # this is the target name
                        fltr,
                        runtime_params,
                        only_these_planets,
                        self.__fin.sv_as_dict()['parameters'],
                        self.__anc.sv_as_dict()['parameters'],
                        self.__xsl.sv_as_dict()[fltr]['data'],
                        self.__atm.sv_as_dict()[fltr]['data'],
                        spectrum,
                        fltrs.index(fltr),
                    )
                    if update:
                        svupdate.append(self.__out[fltrs.index(fltr)])
                else:
                    if targetlistcheck:
                        errstr = [
                            m for m in [sxsl, satm, sspc] if m is not None
                        ]
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
                f'No output created for GEMLI.{self.name()}'
            )
        return

    def _mlfit(
        self,
        trgt,
        fltr,
        runtime_params,
        fin,
        ancil,
        xsl,
        crbatmos,
        spectrum,
        index,
        only_these_planets,
    ):
        '''Core code call'''
        mlfitout = gemlicore.mlfit(
            trgt,
            fltr,
            runtime_params,
            fin,
            ancil,
            xsl,
            crbatmos,
            spectrum,
            self.__out[index],
            only_these_planets=only_these_planets,
            verbose=False,
        )
        return mlfitout

    @staticmethod
    def _failure(errstr, target):
        '''Failure log'''
        log.warning('--< GEMLI MLFIT: %s  %s >--', errstr, target)
        return

    pass


# ---------------- ---------------------------------------------------


class Analysis(dawgie.Analyzer):
    '''analysis ds'''

    def __init__(self):
        '''__init__ ds'''
        # use same version number as results
        self._version_ = gemlicore.mlfitversion()
        # self.__fin = sysalg.finalize()
        # self.__xsl = xslib()
        # self.__atm = atmos()
        # self.__out = gemlistates.AnalysisSv('retrievalCheck')
        self.__rt = rtalg.Autofill()
        # self.__rtc = rtalg.Create()
        self.__out = [gemlistates.AnalysisSv(fltr) for fltr in fltrs]
        return

    # def previous(self):
    #    '''Input State Vectors: gemli.atmos'''
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
            dawgie.SV_REF(fetch('excalibur.gemli').task, MLfit(), sv)
            for sv in MLfit().state_vectors()
        ]

    def state_vectors(self):
        '''Output State Vectors: gemli.analysis'''
        return self.__out

    def run(self, aspects: dawgie.Aspect):
        '''Top level algorithm call'''

        svupdate = []
        if len(aspects) == 0:
            log.warning('--< GEMLI ANALYSIS: contains no targets >--')
        else:
            # determine which filters have results from cerb.atmos (in aspects)
            #  (you have to loop through all targets, since filters vary by target)
            fwr = []
            for trgt in aspects:
                for fltr in fltrs:
                    if (fltr not in fwr) and (
                        'gemli.atmos.' + fltr in aspects[trgt]
                    ):
                        # print('This filter exists in the cerb.atmos aspect:',fltr,trgt)
                        fwr.append(fltr)
            if not fwr:
                log.warning(
                    '--< GEMLI ANALYSIS: NO FILTERS WITH ATMOS DATA!!!>--'
                )

            # fwr = ['Ariel-sim']  # just one filter, while debugging
            # fwr =['HST-WFC3-IR-G141-SCAN']  # just one filter, while debugging

            # only consider filters that have cerb.atmos results loaded in as an aspect
            for fltr in fwr:
                # if 'gemli.atmos.'+fltr not in aspects[trgt]:
                #    log.warning('--< GEMLI ANALYSIS: %s not found IMPOSSIBLE!!!!>--', fltr)

                # (asdf: this is still not working)
                runtime = self.__rt.sv_as_dict()['status']
                # print('runtime old way',runtime)
                # runtime2 = self.__rtc.sv_as_dict()['status']
                # print('runtime old2 way',runtime2)

                # RUNTIME DOESNT WORK YET FOR ASPECTS!!
                runtime_params = crbcore.CerbAnalysisParams(
                    # tier=runtime['ariel_simspectrum_tier'].value(),
                    tier=2,
                    boundTeq=runtime['cerberus_atmos_bounds_Teq'],
                    boundAbundances=runtime['cerberus_atmos_bounds_abundances'],
                    boundCTP=runtime['cerberus_atmos_bounds_CTP'],
                    boundHLoc=runtime['cerberus_atmos_bounds_HLoc'],
                    boundHScale=runtime['cerberus_atmos_bounds_HScale'],
                    boundHThick=runtime['cerberus_atmos_bounds_HThick'],
                )
                # if runtime_params.tier == None:
                #    runtime_params.tier = 2  # no dice. it's too tupley
                # print('runtime', runtime_params)

                log.info('--< GEMLI ANALYSIS: %s  >--', fltr)
                update = self._analysis(
                    aspects, fltr, runtime_params, fltrs.index(fltr)
                )
                if update:
                    svupdate.append(self.__out[fltrs.index(fltr)])
        self.__out = svupdate
        if self.__out:
            aspects.ds().update()
        else:
            raise dawgie.NoValidOutputDataError(
                f'No output created for GEMLI.{self.name()}'
            )
        return

    def _analysis(self, aspects, fltr, runtime_params, index):
        '''Core code call'''
        analysisout = gemlicore.analysis(
            aspects, fltr, runtime_params, self.__out[index], verbose=False
        )
        return analysisout

    @staticmethod
    def _failure(errstr):
        '''Failure log'''
        log.warning('--< GEMLI ANALYSIS: %s >--', errstr)
        return


# -------------------------------------------------------------------
