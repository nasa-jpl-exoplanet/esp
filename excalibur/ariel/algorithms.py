'''ariel algorithms ds'''

# -- IMPORTS -- ------------------------------------------------------
import logging

import dawgie

import excalibur
import excalibur.ariel.core as arielcore
import excalibur.ariel.states as arielstates
import excalibur.runtime as rtime
import excalibur.runtime.algorithms as rtalg
import excalibur.system as sys
import excalibur.system.algorithms as sysalg

from excalibur.util.checksv import checksv

log = logging.getLogger(__name__)


# ------------- ------------------------------------------------------
# -- ALGORITHMS -- ---------------------------------------------------
class SimSpectrum(dawgie.Algorithm):
    '''
    Create a simulated Ariel spectrum
    '''

    def __init__(self):
        '''__init__ ds'''
        self._version_ = dawgie.VERSION(1, 0, 0)
        self.__rt = rtalg.Autofill()
        self.__system_finalize = sysalg.Finalize()
        self.__out = arielstates.SimSpectrumSV('parameters')
        return

    def name(self):
        '''Database name for subtask extension'''
        return 'sim_spectrum'

    def previous(self):
        '''Input State Vectors: system.finalize'''
        return [
            dawgie.ALG_REF(sys.task, self.__system_finalize),
        ] + self.__rt.refs_for_validity()

    def state_vectors(self):
        '''Output State Vectors: ariel.sim_spectrum'''
        return [self.__out]

    def run(self, ds, ps):
        '''Top level algorithm call'''

        # stop here if it is not a runtime target
        if not self.__rt.is_valid():
            log.warning(
                '--< ARIEL.%s: not a valid target >--', self.name().upper()
            )

        else:
            update = False

            system_dict = self.__system_finalize.sv_as_dict()['parameters']
            valid, errstring = checksv(system_dict)

            target = repr(self).split('.')[1]
            if target.startswith('test'):
                valid = False
                errstring = 'Do not run ariel-sim for test targets'

            if valid:
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
                    runtime_params,
                    self.__out,
                )
            else:
                self._failure(errstring)
            if update:
                _ = excalibur.lagger()
                ds.update()
                pass
            elif valid:
                raise dawgie.NoValidOutputDataError(
                    f'No output created for ARIEL.{self.name()}'
                )
        return

    @staticmethod
    def _sim_spectrum(target, system_dict, runtime_params, out):
        '''Core code call'''
        filled = arielcore.simulate_spectra(
            target, system_dict, runtime_params, out
        )
        return filled

    @staticmethod
    def _failure(errstr):
        '''Failure log'''
        log.warning('--< ARIEL SIM_SPECTRUM: %s >--', errstr)
        return
