'''cerberus __init__ ds'''

import dawgie

DAWGIE_IGNORE = False

# pylint: disable=duplicate-code


def analysis(
    prefix: str, ps_hint: int = 0, runid: int = -1
) -> dawgie.FactoryPlaceholder[dawgie.base.Analysis]:
    raise NotImplementedError('placeholder until dawgie monkey patches me')


def events() -> dawgie.FactoryPlaceholder[list[dawgie.EVENT]]:
    raise NotImplementedError('placeholder until dawgie monkey patches me')


def regress(
    prefix: str, ps_hint: int = 0, target: str = '__none__'
) -> dawgie.FactoryPlaceholder[dawgie.base.Regress]:
    raise NotImplementedError('placeholder until dawgie monkey patches me')


def task(
    prefix: str, ps_hint: int = 0, runid: int = -1, target: str = '__none__'
) -> dawgie.FactoryPlaceholder[dawgie.base.Task]:
    raise NotImplementedError('placeholder until dawgie monkey patches me')


# pylint: enable=duplicate-code
