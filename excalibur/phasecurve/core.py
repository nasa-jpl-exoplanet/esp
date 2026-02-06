'''phasecurve core ds'''

# Heritage code shame:
# pylint: disable=duplicate-code
# pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals,too-many-statements

# -- IMPORTS -- ------------------------------------------------------
import copy
import numpy as np

import logging

import excalibur.system.core as syscore
import excalibur.transit.core as trncore
import excalibur.util.monkey_patch  # side effects # noqa: F401 # pylint: disable=unused-import
from excalibur.util.plotters import (
    save_plot_myfit,
    save_plot_tosv,
    plot_residual_fft,
)
from excalibur.phasecurve.plotters import plot_phasecurve

from collections import namedtuple

log = logging.getLogger(__name__)

CONTEXT = namedtuple(
    'CONTEXT',
    [
        'alt',
        'ald',
        'allz',
        'orbp',
        'commonoim',
        'ecc',
        'g1',
        'g2',
        'g3',
        'g4',
        'ootoindex',
        'ootorbits',
        'orbits',
        'period',
        'selectfit',
        'smaors',
        'time',
        'tmjd',
        'ttv',
        'valid',
        'visits',
        'aos',
        'avi',
    ],
)
ctxt = CONTEXT(
    alt=None,
    ald=None,
    allz=None,
    orbp=None,
    commonoim=None,
    ecc=None,
    g1=None,
    g2=None,
    g3=None,
    g4=None,
    ootoindex=None,
    ootorbits=None,
    orbits=None,
    period=None,
    selectfit=None,
    smaors=None,
    time=None,
    tmjd=None,
    ttv=None,
    valid=None,
    visits=None,
    aos=None,
    avi=None,
)


def ctxtupdt(
    alt=None,
    ald=None,
    allz=None,
    orbp=None,
    commonoim=None,
    ecc=None,
    g1=None,
    g2=None,
    g3=None,
    g4=None,
    ootoindex=None,
    ootorbits=None,
    orbits=None,
    period=None,
    selectfit=None,
    smaors=None,
    time=None,
    tmjd=None,
    ttv=None,
    valid=None,
    visits=None,
    aos=None,
    avi=None,
):
    '''
    G. ROUDIER: Update global context for pymc deterministics
    '''
    trncore.ctxt = CONTEXT(
        alt=alt,
        ald=ald,
        allz=allz,
        orbp=orbp,
        commonoim=commonoim,
        ecc=ecc,
        g1=g1,
        g2=g2,
        g3=g3,
        g4=g4,
        ootoindex=ootoindex,
        ootorbits=ootorbits,
        orbits=orbits,
        period=period,
        selectfit=selectfit,
        smaors=smaors,
        time=time,
        tmjd=tmjd,
        ttv=ttv,
        valid=valid,
        visits=visits,
        aos=aos,
        avi=avi,
    )
    return


# ----------------- --------------------------------------------------
# -- NORMALIZATION -- ------------------------------------------------
def norm_spitzer(cal, tme, fin, _ext, out, selftype):
    '''
    K. PEARSON: aperture selection, remove nans, remove zeros, 3 sigma clip time series
    '''
    normed = False
    priors = fin['priors'].copy()

    planetloop = [
        pnet
        for pnet in tme['data'].keys()
        if (pnet in priors.keys()) and tme['data'][pnet][selftype]
    ]

    for p in planetloop:
        out['data'][p] = {}

        keys = [
            'TIME',
            'WX',
            'WY',
            'FRAME',
        ]
        for k in keys:
            out['data'][p][k] = np.array(cal['data'][k])

        # is set later during aperture selection
        out['data'][p]['PHOT'] = np.zeros(len(cal['data']['TIME']))
        out['data'][p]['NOISEPIXEL'] = np.zeros(len(cal['data']['TIME']))

        # remove nans
        nanmask = np.isnan(out['data'][p]['PHOT'])
        for k in out['data'][p].keys():
            nanmask = nanmask | np.isnan(out['data'][p][k])

        for k in out['data'][p].keys():
            out['data'][p][k] = out['data'][p][k][~nanmask]

        # time order things
        ordt = np.argsort(out['data'][p]['TIME'])
        for k in out['data'][p].keys():
            out['data'][p][k] = out['data'][p][k][ordt]
        cflux = np.array(cal['data']['PHOT'])[~nanmask][ordt]
        cnpp = np.array(cal['data']['NOISEPIXEL'])[~nanmask][ordt]

        # 3 sigma clip flux time series
        phase = (out['data'][p]['TIME'] - fin['priors'][p]['t0']) / fin[
            'priors'
        ][p]['period']
        badmask = np.zeros(out['data'][p]['TIME'].shape).astype(bool)
        for i in np.unique(tme['data'][p][selftype]):
            # mask out phase curve data (previous eclipse | transit | eclipse)
            omask = (
                (phase > i - 0.75)
                & (phase < i + 0.75)
                & ~np.isnan(out['data'][p]['TIME'])
            )

            if omask.sum() == 0:
                continue

            dt = np.nanmean(np.diff(out['data'][p]['TIME'][omask])) * 24 * 60
            ndt = int(7 / dt) * 2 + 1
            # aperture selection
            stds = []
            for j in range(cflux.shape[1]):
                stds.append(np.nanstd(trncore.sigma_clip(cflux[omask, j], ndt)))

            bi = np.argmin(stds)
            out['data'][p]['PHOT'][omask] = cflux[omask, bi]
            out['data'][p]['NOISEPIXEL'][omask] = cnpp[omask, bi]

            # sigma clip and remove nans
            photmask = np.isnan(
                trncore.sigma_clip(out['data'][p]['PHOT'][omask], ndt)
            )
            xmask = np.isnan(
                trncore.sigma_clip(out['data'][p]['WX'][omask], ndt)
            )
            ymask = np.isnan(
                trncore.sigma_clip(out['data'][p]['WY'][omask], ndt)
            )
            nmask = np.isnan(
                trncore.sigma_clip(out['data'][p]['NOISEPIXEL'][omask], ndt)
            )

            badmask[omask] = photmask | xmask | ymask | nmask

        # remove outliers
        for k in out['data'][p].keys():
            out['data'][p][k] = out['data'][p][k][~badmask]

        # pass information along
        out['data'][p]['transit'] = tme['data'][p]['transit']
        out['data'][p]['eclipse'] = tme['data'][p]['eclipse']
        out['data'][p]['phasecurve'] = tme['data'][p]['phasecurve']

        if out['data'][p][selftype]:
            normed = True
            out['STATUS'].append(True)

    return normed


def phasecurve_spitzer(nrm, fin, out, selftype, fltr):
    '''
    K. PEARSON: modeling of phase curves
    '''
    wl = False
    priors = fin['priors'].copy()
    ssc = syscore.ssconstants()
    planetloop = list(nrm['data'].keys())

    for p in planetloop:

        out['data'][p] = []

        # extract data based on phase
        if selftype == 'transit':
            phase = (nrm['data'][p]['TIME'] - fin['priors'][p]['t0']) / fin[
                'priors'
            ][p]['period']
        elif selftype == 'eclipse':
            priors = fin['priors']
            w = priors[p].get('omega', 0)
            tme = priors[p]['t0'] + priors[p]['period'] * 0.5 * (
                1 + priors[p]['ecc'] * (4.0 / np.pi) * np.cos(np.deg2rad(w))
            )
            phase = (nrm['data'][p]['TIME'] - tme) / fin['priors'][p]['period']
        elif selftype == 'phasecurve':
            phase = (nrm['data'][p]['TIME'] - fin['priors'][p]['t0']) / fin[
                'priors'
            ][p]['period']
        else:
            log.warning(
                'PHASECURVE phasecurve_spitzer: UNKNOWN DATA TYPE (%s)',
                selftype,
            )
            phase = []

        # loop through epochs
        ec = 0  # event counter
        for event in nrm['data'][p][selftype]:
            print('processing event:', event)

            # compute phase + priors
            smaors = priors[p]['sma'] / priors['R*'] / ssc['Rsun/AU']
            # smaors_up = (priors[p]['sma']+priors[p]['sma_uperr'])/(priors['R*']-abs(priors['R*_lowerr']))/ssc['Rsun/AU']
            # smaors_lo = (priors[p]['sma']-abs(priors[p]['sma_lowerr']))/(priors['R*']+priors['R*_uperr'])/ssc['Rsun/AU']
            priors[p]['ars'] = smaors

            tmid = priors[p]['t0'] + event * priors[p]['period']
            rprs = (priors[p]['rp'] * 7.1492e7) / (priors['R*'] * 6.955e8)
            # inc_lim = 90 - np.rad2deg(np.arctan((priors[p]['rp'] * ssc['Rjup/Rsun'] + priors['R*']) / (priors[p]['sma']/ssc['Rsun/AU'])))
            w = priors[p].get('omega', 0)

            # mask out data by event type
            pmask = (phase > event - 0.65) & (
                phase < event + 0.65
            )  # should check for eccentric orbits

            if pmask.sum() == 0:
                continue

            # extract aperture photometry data
            subt = nrm['data'][p]['TIME'][pmask]
            aper = nrm['data'][p]['PHOT'][pmask]
            aper_err = np.sqrt(aper)

            fpfs = trncore.eclipse_ratio(priors, p, fltr)
            edepth = fpfs * rprs**2

            tpars = {
                # Star
                'T*': priors['T*'],
                # transit
                'rprs': rprs,
                'ars': smaors,
                'tmid': tmid,
                'per': priors[p]['period'],
                'inc': priors[p]['inc'],
                'omega': priors['b'].get('omega', 0),
                'ecc': priors['b']['ecc'],
                # limb darkening (nonlinear - exotethys - pylightcurve)
                'u0': priors[p].get('u0', 0),
                'u1': priors[p].get('u1', 0),
                'u2': priors[p].get('u2', 0),
                'u3': priors[p].get('u3', 0),
                # phase curve amplitudes
                'c0': 0,
                'c1': edepth * 0.25,
                'c2': 0,
                'c3': 0,
                'c4': 0,
                'fpfs': fpfs,
            }

            # gather detrending parameters
            wxa = nrm['data'][p]['WX'][pmask]
            wya = nrm['data'][p]['WY'][pmask]
            npp = nrm['data'][p]['NOISEPIXEL'][pmask]
            syspars = np.array([wxa, wya, npp]).T

            # only report params with bounds, all others will be fixed to initial value
            mybounds = {
                'rprs': [0.5 * rprs, 1.5 * rprs],
                'tmid': [tmid - 0.01, tmid + 0.01],
                'inc': [tpars['inc'] - 3, max(90, tpars['inc'] + 3)],
                'fpfs': [0, fpfs * 3],
                # 'omega':[priors[p]['omega']-25,priors[p]['omega']+25],
                # 'ecc': [0,priors[p]['ecc']+0.1],
                # 'c0':[-0.025,0.025], gets set in code
                'c1': [0, 0.5 * edepth],
                'c2': [-0.15 * edepth, 0.15 * edepth],
            }

            # 10 minute time scale
            nneighbors = int(10.0 / 24.0 / 60.0 / np.mean(np.diff(subt)))
            nneighbors = min(300, nneighbors)
            print(" N neighbors:", nneighbors)
            print(" N datapoints:", len(subt))

            myfit = trncore.pc_fitter(
                subt,
                aper,
                aper_err,
                tpars,
                mybounds,
                syspars,
                neighbors=nneighbors,
                verbose=False,
            )

            # copy best fit parameters and uncertainties
            for k in myfit.bounds.keys():
                print(
                    f" {k} = {myfit.parameters[k]:.6f} +/- {myfit.errors[k]:.6f}"
                )
                pass
            out['data'][p].append({})
            out['data'][p][ec]['time'] = subt
            out['data'][p][ec]['flux'] = aper
            out['data'][p][ec]['err'] = aper_err
            out['data'][p][ec]['xcent'] = wxa
            out['data'][p][ec]['ycent'] = wya
            out['data'][p][ec]['npp'] = npp
            out['data'][p][ec]['wf'] = myfit.wf
            out['data'][p][ec]['model'] = myfit.model
            out['data'][p][ec]['transit'] = myfit.transit
            out['data'][p][ec]['residuals'] = myfit.residuals
            out['data'][p][ec]['detrended'] = myfit.detrended
            out['data'][p][ec]['filter'] = fltr
            out['data'][p][ec]['final_pars'] = copy.deepcopy(myfit.parameters)
            out['data'][p][ec]['final_errs'] = copy.deepcopy(myfit.errors)

            # 11/17/24 also save the MCMC results (for corner plot of the posteriors)
            out['data'][p][ec]['results'] = myfit.results

            out['data'][p][ec]['plot_bestfit'] = save_plot_myfit(
                myfit.plot_bestfit
            )
            out['data'][p][ec]['plot_posterior'] = save_plot_myfit(
                myfit.plot_posterior
            )
            out['data'][p][ec]['plot_pixelmap'] = save_plot_myfit(
                myfit.plot_pixelmap
            )

            out['data'][p][ec]['plot_residual_fft'] = plot_residual_fft(
                selftype,
                fltr,
                p,
                aper,
                subt,
                myfit,
            )

            ec += 1
            out['STATUS'].append(True)
            wl = True

            pass
        pass
    return wl


def flaredetection(whitelight, fin, out, selftype, fltr, verbose=False):
    '''
    find flares in a whitelight phasecurve
    '''
    processed = False

    planetloop = list(whitelight['data'].keys())
    for p in planetloop:
        out['data'][p] = []

        # print('# of events', len(whitelight['data'][p]))
        # loop through epochs
        ec = 0  # event counter
        for event in whitelight['data'][p]:
            # print('processing event:', event)
            # print('event keys', ec, whitelight['data'][p][ec].keys())

            whitelightdata = whitelight['data'][p][ec]
            # no need to pass in systemparam. fit values are in whitelightdata
            systemparam = fin['priors'][p]
            # (this isn't used; it's just to get around pylint checks)
            tmidAdjust = (
                systemparam['t0'] - whitelightdata['final_pars']['tmid']
            )

            randomTimeIndices = np.array(
                len(whitelightdata['time']) * np.random.rand(10), dtype=int
            )
            flaretimes = whitelightdata['time'][randomTimeIndices]
            flarephases = (
                flaretimes - whitelightdata['final_pars']['tmid']
            ) / whitelightdata['final_pars']['per']
            flarephases = flarephases - int(np.max(flarephases))

            # phase isn't defined yet. have to calculate it here
            # phasedtime = (
            #    whitelightdata['time'] - systemparam['t0']
            # ) / systemparam['period']
            # use the fit values, not the original system params
            phasedtime = (
                whitelightdata['time'] - whitelightdata['final_pars']['tmid']
            ) / whitelightdata['final_pars']['per']
            # this makes phase between -0.5 and 0.5
            # whitelightdata['phase'] = (phasedtime + 0.5) % 1 - 0.5
            # actually keep all the data in order; don't wrap around
            whitelightdata['phase'] = phasedtime - int(np.max(phasedtime))

            plot_lightcurve, _ = plot_phasecurve(
                whitelightdata, flarephases, verbose=verbose
            )

            out['data'][p].append({})
            out['data'][p][ec]['plot_lightcurve'] = save_plot_tosv(
                plot_lightcurve
            )

            ec += 1
            out['STATUS'].append(True)
            processed = True
            pass
        pass
    return processed
