'''The actors/agents that do something'''

import dawgie
from . import algorithms


class AnalysisTeam(dawgie.Analysis):
    '''Analytical team'''

    def list(self) -> [dawgie.Analyzer]:
        '''list of analysis to be done'''
        return [algorithms.Create()]

    pass


class TaskTeam(dawgie.Task):
    '''Task team'''

    def __init__(self, *args, **kwds):
        '''override task without knowing anything about it'''
        dawgie.Task.__init__(self, *args, **kwds)

    def list(self) -> [dawgie.Task]:
        '''list of tasks to perform'''
        return [algorithms.Autofill()]

    pass
