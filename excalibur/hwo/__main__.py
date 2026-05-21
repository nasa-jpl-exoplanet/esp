'''hwo __main__ ds'''

# -- IMPORTS -- ------------------------------------------------------
import dawgie
import dawgie.db
import dawgie.security

from excalibur.util.main import main_start

import excalibur.hwo.bot

# ------------- ------------------------------------------------------
if __name__ == "__main__":
    rid, tn = main_start()

    excalibur.hwo.bot.Actor('hwo', 4, rid, tn).do()

    dawgie.db.close()
    dawgie.security.finalize()
    pass
