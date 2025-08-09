'''The actors/agents that do something'''

import dawgie
from . import algorithms


class AnalysisTeam(dawgie.Analysis):
    '''Analytical team'''

    def list(self) -> [dawgie.Analyzer]:
        '''list of analysis to be done'''
        return [algorithms.Performance()]

    pass
