'''cerberus fmcontext ds'''

# Heritage code shame:
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals

import excalibur

from collections import namedtuple

# GLOBAL CONTEXT FOR PYMC DETERMINISTICS
# GMR: Labels are defined here and nowhere else, overcomplicated for retrocomp
ndctx = {
    'atom_xsec': None,
    'atomlist': None,
    'chemistry': None,
    'cialist': None,
    'cleanup': None,
    'fixedParams': None,
    'forwardmodel': None,
    'hitemplist': None,
    'Hsmax': None,
    'hzlib': None,
    'interp_tea': None,
    'isothermal': None,
    'mcmcdat': None,
    'mcmcsig': None,
    'mcmcwav': None,
    'model': None,
    'modparlbl': None,
    'nlevels': None,
    'nodeshape': None,
    'offsetthr': None,
    'orbp': None,
    'planet': None,
    'priors': None,
    'rp0': None,
    'runtime': None,
    'solrad': None,
    'spc': None,
    'tspectrum': None,
    'xmollist': None,
    'xsl': None,
}

CONTEXT = namedtuple('CONTEXT', ndctx.keys())


def dctxupdt(dct=None, freeze=False):
    '''
    GMR: Rewriting this as it was intended to start with
    '''
    if dct is None:  # INIT
        dctx = ndctx
        pass
    else:  # UPDATE
        dctx = excalibur.cerberus.forward_model.dctx
        for k in dct:
            dctx[k] = dct[k]
            pass
        pass
    excalibur.cerberus.forward_model.dctx = dctx
    # GMR: Immutables are no good at creation for us, use dicts.
    # Should freeze context before sampling and use namedtuples in forward model.
    if freeze:
        ctxt = CONTEXT(**dctx)
        excalibur.cerberus.forward_model.ctxt = ctxt
        # GMR: We do not want to duplicate heavy context interpolators
        # Clean that up someday
        excalibur.util.tensor.ctxt = ctxt
        pass
    return dctx


# --------------------------------
# -- DITCH WHAT S BELOW SOMEDAY --


def ctxtinit():
    '''
    GMR: Init context variables
    '''
    ctxt = CONTEXT(**ndctx)
    return ctxt
