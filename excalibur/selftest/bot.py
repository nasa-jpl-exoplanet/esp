'''selftest bot ds'''

#  Careful: Commenting out tasks below will mess up the pipeline during startup
#          (tasks listed here are compiled into the algorithm/task tree)

# -- IMPORTS ---------------------------------------------------------
import dawgie

import excalibur.selftest.algorithms as selftestalg


# --------------------------------------------------------------------
# -- A&A -------------------------------------------------------------
class Actor(dawgie.Task):
    '''Actor ds'''

    def list(self):
        '''Subtasks top level ordered call'''
        return [
            selftestalg.SimSpectrum(),
            selftestalg.XSLib(),
            selftestalg.Atmos(),
            selftestalg.Results(),
        ]

    pass


# -------------------------------------------------------------------
class Agent(dawgie.Analysis):
    '''Agent ds'''

    def list(self) -> [dawgie.Analyzer]:
        '''list ds'''
        return [selftestalg.Analysis()]

    pass


# -------------------------------------------------------------------
