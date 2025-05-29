'''cerberus fmcontext ds'''

import excalibur

from collections import namedtuple


# -- GLOBAL CONTEXT FOR PYMC DETERMINISTICS ---------------------------------------------

CONTEXT = namedtuple(
    'CONTEXT',
    [
        'cleanup',
        'model',
        'p',
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
    ],
)

# ctxt = CONTEXT(
#    cleanup=None,
#    model=None,
#    p=None,
#    solidr=None,
#    orbp=None,
#    tspectrum=None,
#    xsl=None,
#    spc=None,
#    modparlbl=None,
#    hzlib=None,
#    fixedParams=None,
#    mcmcdat=None,
#    mcmcsig=None,
#    nodeshape=None,
#    forwardmodel=None,
# )


def ctxtupdt(
    cleanup=None,
    model=None,
    p=None,
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
        p=p,
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
    )
    return
