'''hwo bot ds'''

#  Careful: Commenting out tasks below will mess up the pipeline during startup
#          (tasks listed here are compiled into the algorithm/task tree)

# -- IMPORTS ---------------------------------------------------------
import dawgie

import excalibur.hwo.algorithms as hwoalg


# --------------------------------------------------------------------
# -- A&A -------------------------------------------------------------
class Actor(dawgie.Task):
    '''Actor ds'''

    def list(self):
        '''Subtasks top level ordered call'''
        return [
            hwoalg.SimSpectrum(),
        ]


# --------- ----------------------------------------------------------
