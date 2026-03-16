'''gemli bot ds'''

# -- IMPORTS -- ------------------------------------------------------
import dawgie

import excalibur.gemli.algorithms as gemlialg


# ------------- ------------------------------------------------------
# -- A&A -- ----------------------------------------------------------
class Actor(dawgie.Task):
    '''Actor ds'''

    def list(self) -> [dawgie.Task]:
        '''Subtasks top level ordered call'''
        return [
            gemlialg.MLfit(),
        ]

    pass


# -------------------------------------------------------------------
class Agent(dawgie.Analysis):
    '''Agent ds'''

    def list(self) -> [dawgie.Analyzer]:
        '''list ds'''
        return [gemlialg.Analysis()]

    pass


# -------------------------------------------------------------------
