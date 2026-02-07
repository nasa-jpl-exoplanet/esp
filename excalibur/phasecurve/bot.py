'''phasecurve bot ds'''

# -- IMPORTS -- ------------------------------------------------------
import dawgie

import excalibur.phasecurve.algorithms as phcalg


# ------------- ------------------------------------------------------
# -- A&A -- ----------------------------------------------------------
class Actor(dawgie.Task):
    '''Actor ds'''

    def list(self) -> [dawgie.Task]:
        '''Subtasks top level ordered call'''
        return [
            phcalg.pcnormalization(),
            phcalg.pcwhitelight(),
            phcalg.pcflaredetection(),
        ]

    pass


# --------- ----------------------------------------------------------
