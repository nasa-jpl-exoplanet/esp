'''system __main__ ds'''

# -- IMPORTS -- ------------------------------------------------------
import dawgie
import dawgie.db
import dawgie.security

from excalibur.util.main import main_start

import excalibur.system.bot

# ------------- ------------------------------------------------------
if __name__ == "__main__":
    rid, tn = main_start()

    if tn in ['', '__all__']:
        NAME = 'population'
        subtasks = excalibur.system.bot.Agent('system', 4, rid)
        pass
    else:
        NAME = ['finalize', 'validate', None][-1]  # -1 to run them all
        subtasks = excalibur.system.bot.Actor('system', 4, rid, tn)
        pass

    subtasks.do(NAME)
    dawgie.db.close()
    dawgie.security.finalize()
    pass
