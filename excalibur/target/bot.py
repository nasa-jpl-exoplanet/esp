'''target bot ds'''

# -- IMPORTS -- ------------------------------------------------------
import dawgie

import excalibur.target.algorithms as trgalg


# ------------- ------------------------------------------------------
# -- A&A -- ----------------------------------------------------------
class Agent(dawgie.Analysis):
    '''Agent ds'''

    def list(self) -> [dawgie.Analyzer]:
        '''list ds'''
        # return [trgalg.Alert(), trgalg.Create()]
        return []

    pass


class Actor(dawgie.Task):
    '''Actor ds'''

    def list(self) -> [dawgie.Task]:
        '''list ds'''
        # return [trgalg.Autofill(), trgalg.Scrape()]
        return []

    pass


class Regress(dawgie.Regress):
    '''Regress ds'''

    def list(self) -> [dawgie.Regression]:
        '''list ds'''
        # return [trgalg.Regress(), trgalg.TargetScrapeRegression()]
        return [trgalg.TargetScrapeRegression()]

    pass


# --------- ----------------------------------------------------------
# pp_exec.sh target 'WASP-69'
