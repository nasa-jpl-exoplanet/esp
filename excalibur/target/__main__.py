'''target __main__ ds'''

# -- IMPORTS -- ------------------------------------------------------

import dawgie
import dawgie.db
import dawgie.security

from excalibur.util.main import main_start

import excalibur.target.bot


# ------------- ------------------------------------------------------

rid, tn = main_start()

if tn in ['', '__all__']:
    name = ['create', None][-1]  # -1 to run them all
    subtasks = excalibur.target.bot.Agent('target', 4, rid)
    pass
elif rid == 0:
    name = ['scrape_regression', None][-1]  # -1 to run them all
    subtasks = excalibur.target.bot.Regress('target', 4, tn)
else:
    name = ['autofill', 'scrape', None][-1]  # -1 to run them all
    subtasks = excalibur.target.bot.Actor('target', 4, rid, tn)
    pass

subtasks.do(name)
dawgie.db.close()
dawgie.security.finalize()
