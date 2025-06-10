'''runtime configuration algorithms'''

import logging

import dawgie
import dawgie.context
import dawgie.db
import os

import excalibur

from . import core
from . import states
from importlib import import_module as fetch  # avoid cyclic-import

log = logging.getLogger(__name__)


class Autofill(dawgie.Algorithm):
    '''Breaks the levers and knobs global table to target specific subtables'''

    def __init__(self):
        '''init autofill'''
        self._version_ = dawgie.VERSION(1, 0, 1)
        self.__parent = Create()
        self.__status = states.StatusSV()

    def is_valid(self):
        '''convenience function'''
        return self.__status['isValidTarget']

    def name(self):
        '''database name'''
        return 'autofill'

    def previous(self):
        '''state vectors needed for graph

        The pipeline never really runs this task. Instead, it is called from
        runtime.create where the table and target name are given
        '''
        return [
            dawgie.SV_REF(
                fetch('excalibur.runtime').analysis,
                self.__parent,
                self.__parent.sv_as_dict()['composite'],
            )
        ]

    def proceed(self, ext: str = None):
        '''convenience function'''
        return self.__status.proceed(ext)

    def refs_for_proceed(self) -> [dawgie.V_REF]:
        '''return list for StatusSV.proceed() to work'''
        return [
            dawgie.V_REF(
                fetch('excalibur.runtime').task,
                self,
                self.__status,
                'allowed_filter_names',
            ),
            dawgie.V_REF(
                fetch('excalibur.runtime').task,
                self,
                self.__status,
                'isValidTarget',
            ),
            dawgie.V_REF(
                fetch('excalibur.runtime').task,
                self,
                self.__status,
                'runTarget',
            ),
        ]

    def refs_for_validity(self) -> [dawgie.V_REF]:
        '''return minimum list for StatusSV to determine if target is valid'''
        return [
            dawgie.V_REF(
                fetch('excalibur.runtime').task,
                self,
                self.__status,
                'isValidTarget',
            )
        ]

    def run(self, ds, ps):
        '''isolate target specific information from the global table'''
        table_index = ds._tn()  # pylint: disable=protected-access
        core.isolate(
            self.__status, self.__parent.sv_as_dict()['composite'], table_index
        )
        _ = excalibur.lagger()
        ds.update()
        return

    def state_vectors(self):
        '''state vectors generated from this algorithm'''
        return [self.__status]

    pass


class Create(dawgie.Analyzer):
    '''Read the configuration file then turn it into state vectors'''

    def __init__(self):
        '''init the create process'''
        self._version_ = dawgie.VERSION(1, 0, 1)
        self.__table = [
            states.FilterSV(),
            states.PymcSV('cerberuschains'),
            states.PymcSV('cerberuschainlen'),
            states.PymcSV('spectrumchains'),
            states.PymcSV('spectrumchainlen'),
            states.ControlsSV(),
            states.TargetsSV('run_only'),
            states.TargetsSV('sequester'),
        ]
        self.__table.append(states.CompositeSV(self.__table))

    @staticmethod
    def _do(arg):
        fetch('excalibur.runtime.bot').TaskTeam(*arg[0], **arg[1]).do()

    def name(self) -> str:
        '''database name'''
        return 'create'

    def run(self, aspects: dawgie.Aspect) -> None:
        '''load the configuration file then process it'''
        try:
            log.info('starting load of config')
            core.load(self.sv_as_dict(), dawgie.db.targets())
            log.info('updating state vector')
            aspects.ds().update()
        except FileNotFoundError as e:
            log.exception(e)
            raise dawgie.AbortAEError(
                f'The environment variable {core.ENV_NAME} points to the non-existent file: {os.environ[core.ENV_NAME]}'
            ) from e
        except KeyError as e:
            log.exception(e)
            raise dawgie.AbortAEError(
                f'The environment variable {core.ENV_NAME} must be defined.'
            ) from e
        except Exception as e:
            log.exception(e)
            raise dawgie.AbortAEError(
                f'The contents pointed to by {core.ENV_NAME} cannot be parsed. Try validating {os.environ[core.ENV_NAME]} with xmllint.'
            ) from e
        return

    def state_vectors(self) -> [dawgie.StateVector]:
        '''produce the configuration file as simple set of state vectors'''
        return self.__table

    def traits(self):
        '''no traits are required'''
        return []

    def where(self):
        '''make sure it always runs on the local cluster for disk access'''
        return dawgie.Distribution.cluster

    pass
