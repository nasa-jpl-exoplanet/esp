'''gemli __main__ ds'''

# -- IMPORTS -- ------------------------------------------------------
import dawgie
import dawgie.db
import dawgie.security

from excalibur.util.main import main_start

import excalibur.gemli

# ------------- ------------------------------------------------------
if __name__ == "__main__":
    rid, tn = main_start()

    if tn in ['', '__all__']:
        NAME = 'analysis'
        subtasks = excalibur.gemli.analysis('gemli', 4, rid)
    else:
        NAME = ['mlfit', None][-1]  # -1 to run them all
        subtasks = excalibur.gemli.task('gemli', 4, rid, tn)
        pass

    subtasks.do(NAME)
    dawgie.db.close()
    dawgie.security.finalize()
    pass
