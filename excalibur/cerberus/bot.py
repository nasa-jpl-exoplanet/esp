'''cerberus bot ds'''

# -- IMPORTS -- ------------------------------------------------------
import dawgie

import excalibur.cerberus.algorithms as crbalg


# ------------- ------------------------------------------------------
# -- A&A -- ----------------------------------------------------------
class Actor(dawgie.Task):
    '''Actor ds'''

    def list(self) -> [dawgie.Task]:
        '''Subtasks top level ordered call'''
        return [crbalg.XSLib(), crbalg.Atmos(), crbalg.Results(), crb.Release()]

    pass


# -------------------------------------------------------------------
class Agent(dawgie.Analysis):
    '''Agent ds'''

    def list(self) -> [dawgie.Analyzer]:
        '''list ds'''
        return [crbalg.Analysis()]

    pass


# -------------------------------------------------------------------
