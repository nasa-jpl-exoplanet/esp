'''cerberus fmcontext ds'''

# Heritage code shame:
# pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals

import excalibur

from collections import namedtuple


# -- GLOBAL CONTEXT FOR PYMC DETERMINISTICS ---------------------------------------------

CONTEXT = namedtuple(
    'CONTEXT',
    [
        'cleanup',
        'model',
        'planet',
        'rp0',
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
        'knownspecies',
        'cialist',
        'xmollist',
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
        rp0=None,
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
        knownspecies=None,
        cialist=None,
        xmollist=None,
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
    rp0=None,
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
        rp0=rp0,
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
        knownspecies=runtime.knownspecies,
        cialist=runtime.cialist,
        xmollist=runtime.xmollist,
        nlevels=runtime.nlevels,
        solrad=runtime.solrad,
        Hsmax=runtime.Hsmax,
        lbroadening=runtime.lbroadening,
        lshifting=runtime.lshifting,
        isothermal=runtime.isothermal,
    )

    excalibur.util.tensor.ctxt = excalibur.cerberus.forward_model.ctxt

    return
