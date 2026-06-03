'''ancillary __main__ ds'''

# -- IMPORTS -- ------------------------------------------------------
import dawgie
import dawgie.db
import dawgie.security

from excalibur.util.main import main_start

import excalibur.ancillary

# ------------- ------------------------------------------------------

if __name__ == "__main__":
    rid, tn = main_start()
    if tn in ['', '__all__']:
        NAME = 'population'
        subtasks = excalibur.ancillary.analysis('ancillary', 4, rid)
        pass
    else:
        NAME = 'estimate'
        subtasks = excalibur.ancillary.task('ancillary', 4, rid, tn)
        pass
    subtasks.do(NAME)
    dawgie.db.close()
    dawgie.security.finalize()
    pass
