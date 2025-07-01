'''
SELFTEST runs a series of retrievals, to test bias, precision, and consistency
'''

# -- IMPORTS -- ------------------------------------------------------
import excalibur.selftest.bot as selftestbot

# ------------- ------------------------------------------------------

DAWGIE_IGNORE = False


def analysis(prefix: str, ps_hint: int = 0, runid: int = -1):
    '''analysis ds'''
    return selftestbot.Agent(prefix, ps_hint, runid)


def task(
    prefix: str, ps_hint: int = 0, runid: int = -1, target: str = '__none__'
):
    '''Factory'''
    return selftestbot.Actor(prefix, ps_hint, runid, target)
