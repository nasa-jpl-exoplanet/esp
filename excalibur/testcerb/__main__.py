'''testcerb __main__ ds'''

# -- IMPORTS -- ------------------------------------------------------
import dawgie
import dawgie.db
import dawgie.security

from excalibur.util.main import main_start

import excalibur.testcerb.bot

# ------------- ------------------------------------------------------

rid, tn = main_start()

excalibur.testcerb.bot.Actor('testcerb', 4, rid, tn).do()

dawgie.db.close()
dawgie.security.finalize()
