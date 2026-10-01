'''transit fmcontext ds'''

import excalibur

from collections import namedtuple

ndctx = {
    'ecc': None,  # eccentricity (to be cleaned up we have orbp)
    'fixedpars': None,  # parameters not fit in mcmc
    'ginc': None,  # inclination (to be cleaned up we have orbp)
    'lclds': None,  # limb darkening coefficients
    'LETHE': None,  # LETHE interpolators
    'mcmcdat': None,  # data for MCMC
    'mcmcsig': None,  # errors on data for MCMC
    'nodeshape': None,  # list of dimensions of each added MCMC node
    'observatory': None,  # HST / JWST
    'orbp': None,  # orbital parameters
    'period': None,  # period (to be cleaned up we have orbp)
    'ref_IM': None,  # zero of x axis for IM formulation
    'selectfit': None,  # valid data point selection
    'modelwrapper': None,  # whitelight / spectrum
    'smaors': None,  # semi major (to be cleaned up we have orbp)
    'time': None,  # time of samples
    'tmjd': None,  # mid transit time (to be cleaned up we have orbp)
    'visits': None,  # list of available visits
    'orbp': None,  # HST
    'g1': None,  # HST
    'g2': None,  # HST
    'g3': None,  # HST
    'g4': None,  # HST
    'orbits': None,  # HST
    'ttv': None,  # HST
}

CONTEXT = namedtuple('CONTEXT', ndctx.keys())


def dctxupdt(dct=None, freeze=False):
    '''
    GMR: Globals for pymc
        dctxt = dctxupdt()  # INIT
        ...
        if whatever: dctxt['key'] = value  # EXAMPLE CONDITIONAL UPDATE
        ...
        dctxt = dctxupdt(dct=dctx, freeze=True)  # STORE AS IMMUTABLE IN ctxt
    '''
    if dct is None:  # INIT
        dctxt = ndctx
        pass
    else:  # UPDATE
        dctxt = excalibur.transit.core.dctxt
        for k in dct:
            dctxt[k] = dct[k]
            pass
        pass
    excalibur.transit.core.dctxt = dctxt
    if freeze:  # IMMUTABLES
        excalibur.transit.core.ctxt = CONTEXT(**dctxt)
        pass
    return dctxt
