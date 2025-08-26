'''routinely process excalibur metric data to make it more accessible'''

# -- IMPORTS -- ------------------------------------------------------
import dawgie
import excalibur.console.bot

from importlib import import_module as fetch  # avoid cyclic-import

# ------------- ------------------------------------------------------

DAWGIE_IGNORE = False


def analysis(prefix: str, ps_hint: int = 0, runid: int = -1):
    '''console is global or an aspect'''
    return excalibur.console.bot.AnalysisTeam(prefix, ps_hint, runid)


def events():
    return [
        dawgie.schedule(
            analysis,
            fetch('excalibur.console.algorithms').Performance(),
            boot=True,
        )
    ]
