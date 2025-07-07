'''First clutch of the pipeline, generates target list and data scraping'''

# -- IMPORTS -- ------------------------------------------------------
import excalibur.target.bot as trgbot

# ------------- ------------------------------------------------------
DAWGIE_IGNORE = False


def analysis(prefix: str, ps_hint: int = 0, runid: int = -1):
    '''analysis (aspect) ds'''
    return trgbot.Agent(prefix, ps_hint, runid)


def regress(prefix: str, ps_hint: int = 0, target: str = '__none__'):
    '''regression ds'''
    return trgbot.Regress(prefix, ps_hint, target)


def task(
    prefix: str, ps_hint: int = 0, runid: int = -1, target: str = '__none__'
):
    '''task ds'''
    return trgbot.Actor(prefix, ps_hint, runid, target)