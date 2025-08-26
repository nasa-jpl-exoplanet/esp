'''algorithms to process metric data'''

import logging

import dawgie
import dawgie.context
import dawgie.db
import os

from . import states

log = logging.getLogger(__name__)


class Performance(dawgie.Analyzer):
    '''Convert the internal metric data to highlight performance elements'''

    def __init__(self):
        '''init the performance process'''
        self._version_ = dawgie.VERSION(1, 0, 0)
        self._fn = os.path.join(dawgie.context.data_per, 'known_metrics.pkl')
        self._ruse = states.CpuAndMem()
        self._what = states.Accomplished()

    def name(self) -> str:
        '''database name'''
        return 'performance'

    def run(self, aspects: dawgie.Aspect) -> None:
        '''load dawgie.db.metrics() then process it to states.CpuAndMem'''
        metrics = dawgie.db.metrics()
        self._ruse.fill(metrics)
        self._what.fill(metrics)
        aspects.ds().update()

    def state_vectors(self) -> [dawgie.StateVector]:
        '''products of this aspect that is fully dynamically created'''
        return [self._ruse, self._what]

    def traits(self):
        '''no traits are required'''
        return []

    def where(self):
        '''this task should always run on the local cluster'''
        return dawgie.Distribution.cluster
