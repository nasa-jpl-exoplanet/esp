'''phasecurve __main__ ds'''

# -- IMPORTS -- ------------------------------------------------------
import os

import dawgie
import dawgie.db
import dawgie.security

from excalibur.util.main import main_start

import excalibur.phasecurve.bot

# ------------- ------------------------------------------------------
if __name__ == "__main__":
    rid, tn = main_start()

    if tn in ["", "__all__"]:
        pass
    else:
        name = os.environ.get("PHASECURVE_SUBTASK")
        if name in (None, "", "all", "None"):
            name = None
        subtasks = excalibur.phasecurve.bot.Actor("phasecurve", 4, rid, tn)
        subtasks.do(name)
        pass

    dawgie.db.close()
    dawgie.security.finalize()
    pass
