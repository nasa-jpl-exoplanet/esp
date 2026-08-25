'''cerberus fmcontext ds'''

# Heritage code shame:
# pylint: disable=invalid-name
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
        'chemistry',
        'fixedParams',
        'mcmcdat',
        'mcmcsig',
        'nodeshape',
        'forwardmodel',
        'knownspecies',
        'cialist',
        'xmollist',
        'atomlist',
        'nlevels',
        'solrad',
        'Hsmax',
        'isothermal',
        'atom_xsec',
        'interp_tea',
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
        chemistry=None,
        fixedParams=None,
        mcmcdat=None,
        mcmcsig=None,
        nodeshape=None,
        forwardmodel=None,
        knownspecies=None,
        cialist=None,
        xmollist=None,
        atomlist=None,
        nlevels=None,
        solrad=None,
        Hsmax=None,
        isothermal=None,
        atom_xsec=None,
        interp_tea=None,
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
    chemistry=None,
    fixed_params=None,
    mcmcdat=None,
    mcmcsig=None,
    nodeshape=None,
    forwardmodel=None,
    atom_xsec=None,
    interp_tea=None,
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
        chemistry=chemistry,
        fixedParams=fixed_params,
        mcmcdat=mcmcdat,
        mcmcsig=mcmcsig,
        nodeshape=nodeshape,
        forwardmodel=forwardmodel,
        knownspecies=runtime.knownspecies,
        cialist=runtime.cialist,
        xmollist=runtime.xmollist,
        atomlist=runtime.atomlist,
        nlevels=runtime.nlevels,
        solrad=runtime.solrad,
        Hsmax=runtime.Hsmax,
        isothermal=runtime.isothermal,
        atom_xsec=atom_xsec,
        interp_tea=interp_tea,
    )

    excalibur.util.tensor.ctxt = excalibur.cerberus.forward_model.ctxt

    return
