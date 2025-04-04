'''ancillary __init__ ds'''

# -- IMPORTS -- ------------------------------------------------------
import excalibur.ancillary.bot as ancbot

# ------------- ------------------------------------------------------

DAWGIE_IGNORE = False


def analysis(prefix: str, ps_hint: int = 0, runid: int = -1):
    '''analysis ds'''
    return ancbot.Agent(prefix, ps_hint, runid)


def task(
    prefix: str, ps_hint: int = 0, runid: int = -1, target: str = '__none__'
):
    '''task ds'''
    return ancbot.Actor(prefix, ps_hint, runid, target)
