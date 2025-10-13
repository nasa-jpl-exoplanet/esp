'''selftest algorithms ds'''

# Heritage code shame:
# pylint: disable=invalid-name,duplicate-code,
# pylint: disable=too-many-arguments,too-many-locals,too-many-positional-arguments,

# -- IMPORTS -- ------------------------------------------------------
import numexpr

import dawgie

# import dawgie.context

import excalibur
import excalibur.system as sys
import excalibur.system.algorithms as sysalg
import excalibur.ancillary as anc
import excalibur.ancillary.algorithms as ancillaryalg

import excalibur.ariel.core as arielcore
import excalibur.ariel.states as arielstates
import excalibur.cerberus.core as crbcore
import excalibur.cerberus.states as crbstates

import excalibur.selftest.core as selftestcore

import excalibur.runtime as rtime
import excalibur.runtime.algorithms as rtalg
import excalibur.runtime.binding as rtbind
from excalibur.util.checksv import checksv

from importlib import import_module as fetch  # avoid cicular dependencies

import logging


log = logging.getLogger(__name__)

numexpr.ncores = 1  # this is actually a performance enhancer!

fltrs = [str(fn) for fn in rtbind.filter_names.values()]


# ------------- ------------------------------------------------------
# -- ALGORITHMS -- ---------------------------------------------------
class SimSpectrum(dawgie.Algorithm):
    '''
    Create a simulated ariel spectrum for selftest
    '''

    def __init__(self):
        '''__init__ ds'''
        self._version_ = dawgie.VERSION(1, 0, 0)
        self.__rt = rtalg.Autofill()
        self.__system_finalize = sysalg.Finalize()
        self.__ancillary = ancillaryalg.Estimate()
        self.__out = arielstates.SimSpectrumSV('parameters')
        return

    def name(self):
        '''Database name for subtask extension'''
        return 'sim_spectrum'

    def previous(self):
        '''Input State Vectors: system.finalize'''
        return [
            dawgie.ALG_REF(sys.task, self.__system_finalize),
            dawgie.ALG_REF(anc.task, self.__ancillary),
        ] + self.__rt.refs_for_validity()

    def state_vectors(self):
        '''Output State Vectors: selftest.sim_spectrum'''
        return [self.__out]

    def run(self, ds, ps):
        '''Top level algorithm call'''

        target = repr(self).split('.')[1]

        # stop here if it is not a runtime target
        if not self.__rt.is_valid() or not target.startswith('test'):
            log.info(
                '--< SELFTEST.%s: not a valid target >--', self.name().upper()
            )

        else:
            update = False

            system_dict = self.__system_finalize.sv_as_dict()['parameters']
            sysvalid, errstring = checksv(system_dict)
            ancil_dict = self.__ancillary.sv_as_dict()['parameters']
            ancvalid, errstring = checksv(ancil_dict)
            if sysvalid and ancvalid:
                runtime = self.__rt.sv_as_dict()['status']
                runtime_params = arielcore.ArielParams(
                    tier=runtime['ariel_simspectrum_tier'].value(),
                    SNRfactor=runtime[
                        'ariel_simspectrum_SNRadjustment'
                    ].value(),
                    randomSeed=runtime['ariel_simspectrum_randomseed'].value(),
                    randomCloudProperties=runtime[
                        'ariel_simspectrum_randomCloudProperties'
                    ],
                    thorngrenMassMetals=runtime[
                        'ariel_simspectrum_thorngrenMassMetals'
                    ],
                    chachanMassMetals=runtime[
                        'ariel_simspectrum_chachanMassMetals'
                    ],
                    includeMetallicityDispersion=runtime[
                        'ariel_simspectrum_includeMetallicityDispersion'
                    ],
                    metallicityDispersion=runtime[
                        'ariel_simspectrum_metallicityDispersion'
                    ].value(),
                    CtoOaverage=runtime[
                        'ariel_simspectrum_CtoOaverage'
                    ].value(),
                    CtoOdispersion=runtime[
                        'ariel_simspectrum_CtoOdispersion'
                    ].value(),
                    knownspecies=runtime[
                        'cerberus_crbmodel_HITEMPmolecules'
                    ].molecules,
                    cialist=runtime[
                        'cerberus_crbmodel_HITRANmolecules'
                    ].molecules,
                    xmollist=runtime[
                        'cerberus_crbmodel_EXOMOLmolecules'
                    ].molecules,
                    nlevels=runtime['cerberus_crbmodel_nlevels'].value(),
                    solrad=runtime['cerberus_crbmodel_solrad'].value(),
                    Hsmax=runtime['cerberus_crbmodel_Hsmax'].value(),
                    lbroadening=runtime['cerberus_crbmodel_lbroadening'],
                    lshifting=runtime['cerberus_crbmodel_lshifting'],
                    isothermal=runtime['cerberus_crbmodel_isothermal'],
                )
                update = self._sim_spectrum(
                    target,
                    system_dict,
                    ancil_dict,
                    runtime_params,
                    self.__out,
                )
            else:
                self._failure(errstring)
            if update:
                _ = excalibur.lagger()
                ds.update()
                pass
            elif sysvalid and ancvalid:
                raise dawgie.NoValidOutputDataError(
                    f'No output created for SELFTEST.{self.name()}'
                )
        return

    @staticmethod
    def _sim_spectrum(target, system_dict, ancil_dict, runtime_params, out):
        '''Core code call'''
        filled = arielcore.simulate_spectra(
            target, system_dict, ancil_dict, runtime_params, out
        )
        return filled

    @staticmethod
    def _failure(errstr):
        '''Failure log'''
        log.warning('--< SELFTEST SIM_SPECTRUM: %s >--', errstr)
        return


# ----------------------- --------------------------------------------
class XSLib(dawgie.Algorithm):
    '''Cross Section Library'''

    def __init__(self):
        '''__init__ ds'''
        self._version_ = crbcore.myxsecsversion()
        self.__arielsim = SimSpectrum()
        self.__rt = rtalg.Autofill()
        self.__out = [crbstates.XslibSv(fltr) for fltr in fltrs]
        return

    def name(self):
        '''Database name for subtask extension'''
        return 'xslib'

    def previous(self):
        '''Input State Vectors: transit.spectrum'''
        return [
            dawgie.ALG_REF(fetch('excalibur.selftest').task, self.__arielsim),
        ] + self.__rt.refs_for_proceed()

    def state_vectors(self):
        '''Output State Vectors: cerberus.xslib'''
        return self.__out

    def run(self, ds, ps):
        '''Top level algorithm call'''

        svupdate = []
        # for fltr in self.__rt.sv_as_dict()['status']['allowed_filter_names']:
        for fltr in ['Ariel-sim']:

            # stop here if it is not a runtime target
            self.__rt.proceed(fltr)
            update = False

            if fltr == 'Ariel-sim':
                sv = self.__arielsim.sv_as_dict()['parameters']
                vspc, sspc = checksv(sv)
                sspc = 'Ariel-sim spectrum not found'
            else:
                vspc = False
                sspc = 'This filter doesnt have a spectrum: ' + fltr

            if vspc:
                log.info('--< SELFTEST XSLIB: %s >--', fltr)

                runtime = self.__rt.sv_as_dict()['status']
                runtime_params = crbcore.CerbXSlibParams(
                    knownspecies=runtime[
                        'cerberus_crbmodel_HITEMPmolecules'
                    ].molecules,
                    cialist=runtime[
                        'cerberus_crbmodel_HITRANmolecules'
                    ].molecules,
                    xmollist=runtime[
                        'cerberus_crbmodel_EXOMOLmolecules'
                    ].molecules,
                    nlevels=runtime['cerberus_crbmodel_nlevels'].value(),
                    solrad=runtime['cerberus_crbmodel_solrad'].value(),
                    Hsmax=runtime['cerberus_crbmodel_Hsmax'].value(),
                    lbroadening=runtime['cerberus_crbmodel_lbroadening'],
                    lshifting=runtime['cerberus_crbmodel_lshifting'],
                )

                update = self._xslib(sv, runtime_params, fltrs.index(fltr))
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
                f'No output created for SELFTEST.{self.name()}'
            )
        return

    def _xslib(self, spc, runtime_params, index):
        '''Core code call'''
        cs = crbcore.myxsecs(
            spc, runtime_params, self.__out[index], verbose=False
        )
        return cs

    @staticmethod
    def _failure(errstr):
        '''Failure log'''
        log.warning('--< SELFTEST XSLIB: %s >--', errstr)
        return

    pass


class Atmos(dawgie.Algorithm):
    '''Atmospheric retrievial'''

    def __init__(self):
        '''__init__ ds'''
        self._version_ = crbcore.atmosversion()
        self.__fin = sysalg.Finalize()
        self.__xsl = XSLib()
        self.__arielsim = SimSpectrum()
        self.__rt = rtalg.Autofill()
        self.__out = [crbstates.AtmosSv(fltr) for fltr in fltrs]
        return

    def name(self):
        '''Database name for subtask extension'''
        return 'atmos'

    def previous(self):
        '''Input State Vectors: transit.spectrum, system.finalize, cerberus.xslib'''
        return [
            dawgie.ALG_REF(sys.task, self.__fin),
            dawgie.ALG_REF(fetch('excalibur.selftest').task, self.__xsl),
            dawgie.ALG_REF(fetch('excalibur.selftest').task, self.__arielsim),
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
        # for fltr in self.__rt.sv_as_dict()['status']['allowed_filter_names']:
        for fltr in ['Ariel-sim']:

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
            else:
                vspc = False
                sspc = 'This filter doesnt have a spectrum: ' + fltr

            if vfin and vxsl and vspc:
                log.info('--< SELFTEST ATMOS: %s >--', fltr)

                runtime = self.__rt.sv_as_dict()['status']
                runtime_params = crbcore.CerbAtmosParams(
                    MCMC_chain_length=runtime['cerberus_steps'].value(),
                    # MCMC_chain_length=33,
                    MCMC_chains=runtime['cerberus_chains'].value(),
                    MCMC_sliceSampler=runtime['cerberus_atmos_sliceSampler'],
                    fitCloudParameters=runtime[
                        'cerberus_atmos_fitCloudParameters'
                    ],
                    cornerBins=runtime['cerberus_plotters_cornerBins'].value(),
                    fitT=runtime['cerberus_atmos_fitT'],
                    fitCtoO=runtime['cerberus_atmos_fitCtoO'],
                    fitNtoO=runtime['cerberus_atmos_fitNtoO'],
                    fitmolecules=runtime[
                        'cerberus_crbmodel_fitmolecules'
                    ].molecules,
                    knownspecies=runtime[
                        'cerberus_crbmodel_HITEMPmolecules'
                    ].molecules,
                    cialist=runtime[
                        'cerberus_crbmodel_HITRANmolecules'
                    ].molecules,
                    xmollist=runtime[
                        'cerberus_crbmodel_EXOMOLmolecules'
                    ].molecules,
                    nlevels=runtime['cerberus_crbmodel_nlevels'].value(),
                    solrad=runtime['cerberus_crbmodel_solrad'].value(),
                    Hsmax=runtime['cerberus_crbmodel_Hsmax'].value(),
                    lbroadening=runtime['cerberus_crbmodel_lbroadening'],
                    lshifting=runtime['cerberus_crbmodel_lshifting'],
                    isothermal=runtime['cerberus_crbmodel_isothermal'],
                    boundTeq=runtime['cerberus_atmos_bounds_Teq'],
                    boundAbundances=runtime['cerberus_atmos_bounds_abundances'],
                    boundCTP=runtime['cerberus_atmos_bounds_CTP'],
                    boundHLoc=runtime['cerberus_atmos_bounds_HLoc'],
                    boundHScale=runtime['cerberus_atmos_bounds_HScale'],
                    boundHThick=runtime['cerberus_atmos_bounds_HThick'],
                )
                # print('runtime params in selftest.alg', runtime_params)

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
                f'No output created for SELFTEST.{self.name()}'
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
        chemistrymodel = 'TEC'
        am = crbcore.atmos(
            fin,
            xsl,
            spc,
            runtime_params,
            self.__out[index],
            fltr,
            Nchains=mcmc_chains,
            chainlen=mcmc_chain_length,
            singlemod=chemistrymodel,
            verbose=False,
        )
        return am

    @staticmethod
    def _failure(errstr):
        '''Failure log'''
        log.warning('--< SELFTEST ATMOS: %s >--', errstr)
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
            dawgie.ALG_REF(fetch('excalibur.selftest').task, self.__xsl),
            dawgie.ALG_REF(fetch('excalibur.selftest').task, self.__atm),
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
            # for fltr in self.__rt.sv_as_dict()['status'][
            #    'allowed_filter_names'
            # ]:
            for fltr in ['Ariel-sim']:
                # stop here if it is not a runtime target
                self.__rt.proceed(fltr)

                vxsl, sxsl = checksv(self.__xsl.sv_as_dict()[fltr])
                vatm, satm = checksv(self.__atm.sv_as_dict()[fltr])

                if vxsl and vatm:
                    log.info('--< SELFTEST RESULTS: %s >--', fltr)

                    runtime = self.__rt.sv_as_dict()['status']
                    runtime_params = crbcore.CerbResultsParams(
                        nrandomwalkers=runtime[
                            'cerberus_results_nrandomwalkers'
                        ].value(),
                        randomseed=runtime[
                            'cerberus_results_randomseed'
                        ].value(),
                        knownspecies=runtime[
                            'cerberus_crbmodel_HITEMPmolecules'
                        ].molecules,
                        cialist=runtime[
                            'cerberus_crbmodel_HITRANmolecules'
                        ].molecules,
                        xmollist=runtime[
                            'cerberus_crbmodel_EXOMOLmolecules'
                        ].molecules,
                        nlevels=runtime['cerberus_crbmodel_nlevels'].value(),
                        Hsmax=runtime['cerberus_crbmodel_Hsmax'].value(),
                        solrad=runtime['cerberus_crbmodel_solrad'].value(),
                        cornerBins=runtime[
                            'cerberus_plotters_cornerBins'
                        ].value(),
                        lbroadening=runtime['cerberus_crbmodel_lbroadening'],
                        lshifting=runtime['cerberus_crbmodel_lshifting'],
                        isothermal=runtime['cerberus_crbmodel_isothermal'],
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
                f'No output created for SELFTEST.{self.name()}'
            )
        return

    def _results(self, trgt, fltr, runtime_params, fin, ancil, xsl, atm, index):
        '''Core code call'''
        resout = crbcore.results(
            trgt,
            fltr,
            runtime_params,
            fin,
            ancil,
            xsl,
            atm,
            self.__out[index],
            verbose=False,
        )
        return resout

    @staticmethod
    def _failure(errstr):
        '''Failure log'''
        log.warning('--< SELFTEST RESULTS: %s >--', errstr)
        return

    pass


# ---------------- ---------------------------------------------------


class Analysis(dawgie.Analyzer):
    '''analysis ds'''

    def __init__(self):
        '''__init__ ds'''
        self._version_ = crbcore.resultsversion()
        # self.__rt = rtalg.Create()
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
            dawgie.SV_REF(fetch('excalibur.selftest').task, Atmos(), sv)
            for sv in Atmos().state_vectors()
        ]

    def state_vectors(self):
        '''Output State Vectors: selftest.analysis'''
        return self.__out

    def run(self, aspects: dawgie.Aspect):
        '''Top level algorithm call'''

        svupdate = []
        if len(aspects) == 0:
            log.warning('--< SELFTEST ANALYSIS: contains no targets >--')
        else:
            filtersWithResults = ['Ariel-sim']  # just consider Ariel-sim

            for fltr in filtersWithResults:
                log.info('--< SELFTEST ANALYSIS: %s  >--', fltr)

                # runtime is not actually needed in analysis.
                #  prior range is loaded in from previous state vector

                # runtime = self.__rt.sv_as_dict()['composite']['controls']
                # these two lines should be the same thing I think
                # runtime = self.__rt.sv_as_dict()['controls']
                # print('runtime new way',runtime)
                # import pdb; pdb.set_trace()

                # runtime_params = selftestcore.SelftestAnalysisParams(
                #    tier=runtime['ariel_simspectrum_tier'].value(),
                #    boundTeq=runtime['cerberus_atmos_bounds_Teq'],
                #    boundAbundances=runtime['cerberus_atmos_bounds_abundances'],
                #    boundCTP=runtime['cerberus_atmos_bounds_CTP'],
                #    boundHLoc=runtime['cerberus_atmos_bounds_HLoc'],
                #    boundHScale=runtime['cerberus_atmos_bounds_HScale'],
                #    boundHThick=runtime['cerberus_atmos_bounds_HThick'],
                # )
                # print()
                # print('runtimeparams in selftest.alg',runtime_params)
                # print()
                # update = self._analysis(aspects, fltr, runtime_params, fltrs.index(fltr))
                chemistrymodel = 'TEC'
                update = self._analysis(
                    aspects, fltr, chemistrymodel, fltrs.index(fltr)
                )
                if update:
                    svupdate.append(self.__out[fltrs.index(fltr)])
        self.__out = svupdate
        if self.__out:
            aspects.ds().update()
        else:
            raise dawgie.NoValidOutputDataError(
                f'No output created for SELFTEST.{self.name()}'
            )
        return

    # def _analysis(self, aspects, fltr, runtime_params, index):
    def _analysis(self, aspects, fltr, chemistrymodel, index):
        '''Core code call'''
        # aspects, fltr, runtime_params, self.__out[index], verbose=False
        analysisout = selftestcore.analysis(
            aspects, fltr, chemistrymodel, self.__out[index], verbose=False
        )
        return analysisout

    @staticmethod
    def _failure(errstr):
        '''Failure log'''
        log.warning('--< SELFTEST ANALYSIS: %s >--', errstr)
        return


# -------------------------------------------------------------------
