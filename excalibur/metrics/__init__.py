'''routinely process excalibur metric data to make it more accessible'''

# -- IMPORTS -- ------------------------------------------------------
import excalibur.metrics.bot

# ------------- ------------------------------------------------------
DAWGIE_IGNORE = False

import dawgie

from importlib import import_module as fetch  # avoid cyclic-import

def analysis(prefix: str, ps_hint: int = 0, runid: int = -1):
    '''metrics are global or an aspect'''
    return excalibur.metrics.bot.AnalysisTeam(prefix, ps_hint, runid)

def events():
    return [dawgie.schedule(analysis, fetch('excalibur.metrics.algorithms').Performance(), boot=True)
