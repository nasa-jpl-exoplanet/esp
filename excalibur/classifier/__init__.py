'''classifier __init__ ds'''

import dawgie

# Do not enable this module until you find the blame commit that creates this
# comment and put back all of the torch elements that were removed
DAWGIE_IGNORE = True

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
