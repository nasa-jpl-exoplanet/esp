'''phasecurve flare_det_context ds'''

# Heritage code shame:
# pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals

import excalibur

from collections import namedtuple

# -- GLOBAL CONTEXT FOR PYMC DETERMINISTICS ---------------------------------------------

CONTEXT = namedtuple(
    'CONTEXT',
    [
        'data',
        'sigma',
    ],
)


def ctxtinit():
    ctxt = CONTEXT(
        data=None,
        sigma=None,
    )
    return ctxt


def ctxtupdt(
    data=None,
    sigma=None,
):
    '''
    Update global context for pymc deterministics
    '''
    excalibur.phasecurve.flare_det_context.ctxt = CONTEXT(
        data=data,
        sigma=sigma,
    )

    return
