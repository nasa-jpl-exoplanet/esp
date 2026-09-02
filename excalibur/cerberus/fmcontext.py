'''cerberus fmcontext ds'''

# Heritage code shame:
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals

import excalibur

from collections import namedtuple

# -- GLOBAL CONTEXT FOR PYMC DETERMINISTICS --
ndctx = {
    'atom_xsec'=None,
    'atomlist'=None,
    'chemistry'=None,
    'cialist'=None,
    'cleanup'=None,
    'fixedParams'=None,
    'forwardmodel'=None,
    'hitemplist'=None,
    'Hsmax'=None,
    'hzlib'=None,
    'interp_tea'=None,
    'isothermal'=None,
    'mcmcdat'=None,
    'mcmcsig'=None,
    'mcmcwav'=None,
    'model'=None,
    'modparlbl'=None,
    'nlevels'=None,
    'nodeshape'=None,
    'offsetthr'=None,
    'orbp'=None,
    'planet'=None,
    'rp0'=None,
    'runtime'=None,
    'solrad'=None,
    'spc'=None,
    'tspectrum'=None,
    'xmollist'=None,
    'xsl'=None,
}

def dctxupdt(dct={}):
    '''
    GMR: Rewriting this as it was intended to start with
    '''
    if not dct:  # INIT
        dctx = ndctx
        pass
    else:  # UPDATE
        dctx = excalibur.cerberus.forward_model.dctx
        for k in dct:
            dctx[k] = dct[k]
            pass
        pass
    excalibur.cerberus.forward_model.dctx = dctx
    excalibur.util.tensor.dctx = dctx
    return

# -- -------------------------------------- --

CONTEXT = namedtuple(
    'CONTEXT',
    [
        'atom_xsec',
        'atomlist',
        'chemistry',
        'cialist',
        'cleanup',
        'fixedParams',
        'forwardmodel',
        'hitemplist',
        'Hsmax',
        'hzlib',
        'interp_tea',
        'isothermal',
        'mcmcdat',
        'mcmcsig',
        'mcmcwav',
        'model',
        'modparlbl',
        'nlevels',
        'nodeshape',
        'offsetthr',
        'orbp',
        'planet',
        'rp0',
        'runtime',
        'solrad',
        'spc',
        'tspectrum',
        'xmollist',
        'xsl',
    ],
)

def ctxtinit():
    '''
    GMR: Init context variables
    Old method ditch that someday
    '''
    ctxt = CONTEXT(
        atom_xsec=None,
        atomlist=None,
        chemistry=None,
        cialist=None,
        cleanup=None,
        fixedParams=None,
        forwardmodel=None,
        hitemplist=None,
        Hsmax=None,
        hzlib=None,
        interp_tea=None,
        isothermal=None,
        mcmcdat=None,
        mcmcsig=None,
        mcmcwav=None,
        model=None,
        modparlbl=None,
        nlevels=None,
        nodeshape=None,
        offsetthr=None,
        orbp=None,
        planet=None,
        rp0=None,
        runtime=None,
        solrad=None,
        spc=None,
        tspectrum=None,
        xmollist=None,
        xsl=None,
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
        hitemplist=runtime.hitemplist,
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
