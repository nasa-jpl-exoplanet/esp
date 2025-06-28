'''testcerb bot ds'''

#  Careful: Commenting out tasks below will mess up the pipeline during startup
#          (tasks listed here are compiled into the algorithm/task tree)

# -- IMPORTS ---------------------------------------------------------
import dawgie

import excalibur.testcerb.algorithms as testcerbalg


# --------------------------------------------------------------------
# -- A&A -------------------------------------------------------------
class Actor(dawgie.Task):
    '''Actor ds'''

    def list(self):
        '''Subtasks top level ordered call'''
        return [
            testcerbalg.SimSpectrum(),
            testcerbalg.XSLib(),
            testcerbalg.Atmos(),
            testcerbalg.Results(),
        ]

    pass


# -------------------------------------------------------------------
class Agent(dawgie.Analysis):
    '''Agent ds'''

    def list(self) -> [dawgie.Analyzer]:
        '''list ds'''
        return [testcerbalg.Analysis()]

    pass


# -------------------------------------------------------------------
