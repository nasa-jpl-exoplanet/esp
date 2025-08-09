'''algorithms to process metric data'''

import logging

import dawgie
import dawgie.db

from . import states

log = logging.getLogger(__name__)


class Performance(dawgie.Analyzer):
    '''Convert the internal metric data to highlight performance elements'''

    def __init__(self):
        '''init the performance process'''
        self._version_ = dawgie.VERSION(1, 0, 0)
        self._svs = [states.CpuAndMem('test', [])]

    def name(self) -> str:
        '''database name'''
        return 'performance'

    def run(self, aspects: dawgie.Aspect) -> None:
        '''load dawgie.db.metrics() then process it to states.CpuAndMem'''
        metrics = dawgie.db.metrics()
        table = {}
        for md in metrics:
            name = '.'.join([md.task, md.alg_name])
            known = table.get(name, [])
            known.append(md)
            table[name] = known
        self._svs = [states.CpuAndMem(*i) for i in table.items()]
        aspects.ds().update()

    def state_vectors(self) -> [dawgie.StateVector]:
        '''products of this aspect that is fully dynamically created'''
        return self._svs

    def traits(self):
        '''no traits are required'''
        return []

    def where(self):
        '''this task should always run on the local cluster'''
        return dawgie.Distribution.cluster
