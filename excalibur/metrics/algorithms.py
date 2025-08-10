'''algorithms to process metric data'''

import logging

import dawgie
import dawgie.context
import dawgie.db
import os
import pickle

from . import states

log = logging.getLogger(__name__)


class Performance(dawgie.Analyzer):
    '''Convert the internal metric data to highlight performance elements'''

    def __init__(self):
        '''init the performance process'''
        self._version_ = dawgie.VERSION(1, 0, 0)
        self._fn = os.path.join(dawgie.context.data_per, 'known_metrics.pkl')
        self._svs = self._load()

    def _load(self) -> []:
        result = [states.CpuAndMem('test', [])]
        if os.path.isfile(self._fn):
            with open(self._fn, 'br') as file:
                known = pickle.load(file)
            result = [states.CpuAndMem(name, []) for name in known]
        else:
            log.warning('Could not read the file: %s', self._fn)
        return result

    def _save(self, known: []):
        if os.path.isdir(os.path.dirname(self._fn)):
            with open(self._fn, 'bw') as file:
                pickle.dump(list(known), file)
        else:
            log.warning('Could not write the file: %s', self._fn)

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
        self._save(table.keys())
        self._svs = [states.CpuAndMem(*i) for i in table.items()]
        tsk = aspects.ds()._bot()  # pylint: disable=protected-access
        for sv in self._svs:
            for vn, v in sv.items():
                dawgie.db.update(tsk, self, sv, vn, v)
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
