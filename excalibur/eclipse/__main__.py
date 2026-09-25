'''eclipse __main__ ds'''

# -- IMPORTS -- ------------------------------------------------------
import dawgie
import dawgie.db
import dawgie.security

from excalibur.util.main import main_start

import excalibur.eclipse

# ------------- ------------------------------------------------------
if __name__ == "__main__":
    rid, tn = main_start()

    if tn in ['', '__all__']:
        pass
    else:
        NAME = ['normalize', 'whitelight', 'spectrum', None][
            -1
        ]  # -1 to run them all
        subtasks = excalibur.eclipse.task('eclipse', 4, rid, tn)
        subtasks.do(NAME)
        dawgie.db.close()
        dawgie.security.finalize()
        pass
    pass
