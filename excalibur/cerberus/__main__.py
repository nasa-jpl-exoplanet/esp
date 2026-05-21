'''cerberus __main__ ds'''

# -- IMPORTS -- ------------------------------------------------------
import dawgie
import dawgie.db
import dawgie.security

from excalibur.util.main import main_start

import excalibur.cerberus.bot

# ------------- ------------------------------------------------------

if __name__ == "__main__":
    rid, tn = main_start()

    if tn in ['', '__all__']:
        NAME = 'analysis'
        subtasks = excalibur.cerberus.bot.Agent('cerberus', 4, rid)
    else:
        NAME = ['atmos', 'results', 'xslib', None][-1]  # -1 to run them all
        subtasks = excalibur.cerberus.bot.Actor('cerberus', 4, rid, tn)
        pass

    subtasks.do(NAME)
    dawgie.db.close()
    dawgie.security.finalize()
    pass
