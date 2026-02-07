'''
HWO simulates HWO observations
'''

# -- IMPORTS -- ------------------------------------------------------
import excalibur.hwo.bot as hwobot

# ------------- ------------------------------------------------------

DAWGIE_IGNORE = False


def task(
    prefix: str, ps_hint: int = 0, runid: int = -1, target: str = '__none__'
):
    '''Factory'''
    return hwobot.Actor(prefix, ps_hint, runid, target)
