'''
TESTCERB runs a series of retrievals, to test bias, precision, and consistency
'''

# -- IMPORTS -- ------------------------------------------------------
import excalibur.testcerb.bot as testcerbbot

# ------------- ------------------------------------------------------

DAWGIE_IGNORE = False


def analysis(prefix: str, ps_hint: int = 0, runid: int = -1):
    '''analysis ds'''
    return testcerbbot.Agent(prefix, ps_hint, runid)


def task(
    prefix: str, ps_hint: int = 0, runid: int = -1, target: str = '__none__'
):
    '''Factory'''
    return testcerbbot.Actor(prefix, ps_hint, runid, target)
