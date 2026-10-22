'''selftest __main__ ds'''

# -- IMPORTS -- ------------------------------------------------------
import dawgie
import dawgie.db
import dawgie.security

from excalibur.util.main import main_start

import excalibur.selftest

# ------------- ------------------------------------------------------
if __name__ == "__main__":

    rid, tn = main_start()

    if tn in ['', '__all__']:
        NAME = 'analysis'
        subtasks = excalibur.selftest.analysis('selftest', 4, rid)
    else:
        NAME = ['simspectrum', 'xslib', 'atmos', 'results', None][
            -1
        ]  # -1 to run them all
        subtasks = excalibur.selftest.task('selftest', 4, rid, tn)
        pass

    subtasks.do(NAME)
    dawgie.db.close()
    dawgie.security.finalize()
    pass
