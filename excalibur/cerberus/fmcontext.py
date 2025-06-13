'''cerberus fmcontext ds'''

# Heritage code shame:
# pylint: disable=too-many-arguments,too-many-positional-arguments

import excalibur

from collections import namedtuple


# -- GLOBAL CONTEXT FOR PYMC DETERMINISTICS ---------------------------------------------

CONTEXT = namedtuple(
    'CONTEXT',
    [
        'cleanup',
        'model',
        'planet',
        'solidr',
        'orbp',
        'tspectrum',
        'xsl',
        'spc',
        'modparlbl',
        'hzlib',
        'fixedParams',
        'mcmcdat',
        'mcmcsig',
        'nodeshape',
        'forwardmodel',
        'nlevels',
        'solrad',
        'Hsmax',
        'lbroadening',
        'lshifting',
        'isothermal',
    ],
)


def ctxtinit():
    ctxt = CONTEXT(
        cleanup=None,
        model=None,
        planet=None,
        solidr=None,
        orbp=None,
        tspectrum=None,
        xsl=None,
        spc=None,
        modparlbl=None,
        hzlib=None,
        fixedParams=None,
        mcmcdat=None,
        mcmcsig=None,
        nodeshape=None,
        forwardmodel=None,
        nlevels=None,
        solrad=None,
        Hsmax=None,
        lbroadening=None,
        lshifting=None,
        isothermal=None,
    )
    return ctxt


def ctxtupdt(
    runtime=None,
    cleanup=None,
    model=None,
    planet=None,
    solidr=None,
    orbp=None,
    tspectrum=None,
    xsl=None,
    spc=None,
    modparlbl=None,
    hzlib=None,
    fixed_params=None,
    mcmcdat=None,
    mcmcsig=None,
    nodeshape=None,
    forwardmodel=None,
):
    '''
    G. ROUDIER: Update global context for pymc deterministics
    '''
    # sys.modules[__name__].ctxt = CONTEXT(
    excalibur.cerberus.forward_model.ctxt = CONTEXT(
        cleanup=cleanup,
        model=model,
        planet=planet,
        solidr=solidr,
        orbp=orbp,
        tspectrum=tspectrum,
        xsl=xsl,
        spc=spc,
        modparlbl=modparlbl,
        hzlib=hzlib,
        fixedParams=fixed_params,
        mcmcdat=mcmcdat,
        mcmcsig=mcmcsig,
        nodeshape=nodeshape,
        forwardmodel=forwardmodel,
        nlevels=runtime.nlevels,
        solrad=runtime.solrad,
        Hsmax=runtime.Hsmax,
        lbroadening=runtime.lbroadening,
        lshifting=runtime.lshifting,
        isothermal=runtime.isothermal,
    )

    excalibur.util.tensor.ctxt = excalibur.cerberus.forward_model.ctxt

    return
