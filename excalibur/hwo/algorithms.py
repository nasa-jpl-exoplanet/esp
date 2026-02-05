'''hwo algorithms ds'''

# -- IMPORTS -- ------------------------------------------------------
import logging

import dawgie

import excalibur
import excalibur.hwo.core as hwocore
import excalibur.hwo.states as hwostates
import excalibur.runtime.algorithms as rtalg
import excalibur.system as sys
import excalibur.system.algorithms as sysalg
import excalibur.ancillary as anc
import excalibur.ancillary.algorithms as ancalg

from excalibur.util.checksv import checksv

log = logging.getLogger(__name__)


# ------------- ------------------------------------------------------
# -- ALGORITHMS -- ---------------------------------------------------
class SimSpectrum(dawgie.Algorithm):
    '''
    Create a simulated HWO spectrum
    '''

    def __init__(self):
        '''__init__ ds'''
        self._version_ = dawgie.VERSION(1, 0, 0)
        self.__rt = rtalg.Autofill()
        self.__system_finalize = sysalg.Finalize()
        self.__ancillary_estimate = ancalg.Estimate()
        self.__out = hwostates.SimSpectrumSV('parameters')
        return

    def name(self):
        '''Database name for subtask extension'''
        return 'sim_spectrum'

    def previous(self):
        '''Input State Vectors: system.finalize'''
        return [
            dawgie.ALG_REF(sys.task, self.__system_finalize),
            dawgie.ALG_REF(anc.task, self.__ancillary_estimate),
        ] + self.__rt.refs_for_validity()

    def state_vectors(self):
        '''Output State Vectors: hwo.sim_spectrum'''
        return [self.__out]

    def run(self, ds, ps):
        '''Top level algorithm call'''

        # stop here if it is not a runtime target
        if not self.__rt.is_valid():
            log.info('--< HWO.%s: not a valid target >--', self.name().upper())

        else:
            update = False

            system_dict = self.__system_finalize.sv_as_dict()['parameters']
            sysvalid, errstring = checksv(system_dict)

            ancil_dict = self.__ancillary_estimate.sv_as_dict()['parameters']
            ancvalid, errstring = checksv(ancil_dict)

            target = repr(self).split('.')[1]
            if target.startswith('test'):
                sysvalid = False
                errstring = 'Do not run hwo-sim for test targets'

            if sysvalid and ancvalid:
                runtime = self.__rt.sv_as_dict()['status']
                runtime_params = hwocore.HWOparams(
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
                    CtoOdaSilva=runtime['ariel_simspectrum_CtoOdaSilva'],
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
                    f'No output created for HWO.{self.name()}'
                )
        return

    @staticmethod
    def _sim_spectrum(target, system_dict, ancil_dict, runtime_params, out):
        '''Core code call'''
        filled = hwocore.simulate_spectra(
            target, system_dict, ancil_dict, runtime_params, out
        )
        return filled

    @staticmethod
    def _failure(errstr):
        '''Failure log'''
        log.warning('--< HWO SIM_SPECTRUM: %s >--', errstr)
        return
