'''
SYSTEM manages the astrophysical parameters of the target observed
- VALIDATE parameters from target.autofill, report missing parameters
- FINALIZE parameters using overwriter.py function ppar()
'''

# -- IMPORTS -- ------------------------------------------------------
import excalibur.system.bot as sysbot

# ------------- ------------------------------------------------------
DAWGIE_IGNORE = False


def task(
    prefix: str, ps_hint: int = 0, runid: int = -1, target: str = '__none__'
):
    '''Factory'''
    return sysbot.Actor(prefix, ps_hint, runid, target)


def analysis(prefix: str, ps_hint: int = 0, runid: int = -1):
    '''analysis ds'''
    return sysbot.Agent(prefix, ps_hint, runid)
