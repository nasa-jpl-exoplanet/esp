'''ancillary __main__ ds'''

# -- IMPORTS -- ------------------------------------------------------
import dawgie
import dawgie.db
import dawgie.security

from excalibur.util.main import main_start

import excalibur.classifier.bot

# ------------- ------------------------------------------------------
if __name__ == "__main__":
    rid, tn = main_start()

    if tn in ['', '__all__']:
        NAME = 'summarize_flags'
        subtasks = excalibur.classifier.bot.Agent('classifier', 4, rid)
    else:
        NAME = ['inference', 'flags', None][
            -1
        ]  # 0 is kicked off list, -1 to run them all
        subtasks = excalibur.classifier.bot.Actor('classifier', 4, rid, tn)
        pass

    subtasks.do(NAME)
    dawgie.db.close()
    dawgie.security.finalize()
    pass
