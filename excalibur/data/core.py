'''data core ds'''

# Heritage code shame:
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments,too-many-branches,too-many-lines,too-many-locals,too-many-nested-blocks,too-many-positional-arguments,too-many-statements

# -- IMPORTS -- ------------------------------------------------------
import os
import glob
import logging

import dawgie
import dawgie.context

import excalibur
import excalibur.system.core as syscore
import excalibur.util.time
from excalibur.util import nerdclub

import lmfit as lm
import matplotlib.pyplot as plt
import time as raissatime
from ldtk import LDPSetCreator, BoxcarFilter
import datetime

import scipy.interpolate as itp
import scipy.signal
import scipy.optimize as opt
from scipy.ndimage.measurements import label
from scipy.ndimage.morphology import (
    binary_dilation,
    binary_closing,
    binary_erosion,
    binary_fill_holes,
)

from photutils.aperture import aperture_photometry, CircularAperture
import pyvo as vo

import numpy as np
import numpy.polynomial.polynomial as poly

import astropy.io.fits as pyfits
import astropy.units
from astropy.modeling.models import BlackBody as astrobb
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.wcs import WCS

log = logging.getLogger(__name__)


# ----------------- --------------------------------------------------
# -- COLLECT DATA -- -------------------------------------------------
def collectversion():
    '''
    1.1.3: GMR: bugfix for JWST NIRSPEC filters
    '''
    return dawgie.VERSION(1, 1, 3)


def collect(name, scrape, out):
    '''
    G. ROUDIER: Filters data from target.scrape.databases according to active filters
    '''
    collected = False
    # For simulated instrument, there is no data to collect
    if not name.split('-')[1] == 'sim':
        obs, ins, det, fil, mod = name.split('-')
        for rootname in scrape['name'].keys():
            ok = scrape['name'][rootname]['observatory'] in [obs.strip()]
            ok = ok and (
                scrape['name'][rootname]['instrument'] in [ins.strip()]
            )
            ok = ok and (det.strip() in scrape['name'][rootname]['detector'])
            ok = ok and (scrape['name'][rootname]['filter'] in [fil.strip()])
            ok = ok and (scrape['name'][rootname]['mode'] in [mod.strip()])
            if ok:
                out['activefilters'][name]['ROOTNAME'].append(rootname)
                loc = (
                    scrape['name'][rootname]['md5']
                    + '_'
                    + scrape['name'][rootname]['sha']
                )
                out['activefilters'][name]['LOC'].append(loc)
                out['activefilters'][name]['TOTAL'].append(True)
                collected = True
                pass
            pass
        pass
    if collected:
        log.info(
            '--< DATA COLLECT: %s in %s >--',
            str(int(np.sum(out['activefilters'][name]['TOTAL']))),
            name,
        )
        out['STATUS'].append(True)
        pass
    else:
        log.info('--< DATA COLLECT: NO DATA in %s >--', name)
        out['activefilters'].pop(name, None)
        pass
    return collected


# ------------------ -------------------------------------------------
# -- TIMING -- -------------------------------------------------------
def timingversion():
    '''
    1.1.0: GMR: HST
    1.2.0: K. Pearson: Spitzer
    1.3.0: GMR: JWST
    '''
    return dawgie.VERSION(1, 3, 0)


def timing(force, ext, clc, out):
    '''
    Uses system orbital parameters to guide the dataset towards
    transit, eclipse or phasecurve tasks
    '''
    ssc = syscore.ssconstants()
    chunked = False
    priors = force['priors'].copy()
    dbs = os.path.join(dawgie.context.data_dbs, 'mast')
    data = {
        'LOC': [],
        'SCANANGLE': [],
        'TIME': [],
        'EXPLEN': [],
        'ISORTEXP': [],
    }
    # LOAD DATA ------------------------------------------------------
    if 'JWST' in ext:
        # Segmented datasets, length clc['LOC'] is the number of segments
        alldinm = []  # Flatten integration number
        alldtms = []  # Flatten times
        alldidr = []  # Flatten integration durations
        for loc in sorted(clc['LOC']):
            fullloc = os.path.join(dbs, loc)
            with pyfits.open(fullloc) as hdulist:
                for hdu in hdulist:
                    if 'PRIMARY' in hdu.name:
                        pass
                    elif "INT_TIMES" in hdu.name:
                        alldinm.extend(hdu.data['integration_number'])
                        alldtms.extend(hdu.data['int_mid_MJD_UTC'])
                        alldidr.extend(
                            hdu.data['int_end_MJD_UTC']
                            - hdu.data['int_start_MJD_UTC']
                        )
                        pass
                    pass
                pass
            pass
        isort = np.argsort(alldtms)
        data['ISORTEXP'].extend(isort)
        data['TIME'].extend(np.array(alldtms)[isort])  # [MJD-UTC]
        data['EXPLEN'].extend(np.array(alldidr)[isort])  # [days]
        data['LOC'].extend(sorted(clc['LOC']))
        pass
    elif 'Spitzer' in ext:
        for loc in sorted(clc['LOC']):
            fullloc = os.path.join(dbs, loc)
            with pyfits.open(fullloc) as hdulist:
                header0 = hdulist[0].header
                ftime = []
                exptime = []
                for hdu in hdulist:
                    start = hdu.header.get('MJD_OBS')
                    if hdu.data.ndim == 3:  # data cube
                        idur = hdu.header.get('ATIMEEND') - hdu.header.get(
                            'AINTBEG'
                        )
                        nimgs = hdu.data.shape[0]
                        dt = idur / nimgs / (24 * 60 * 60)
                        for i in range(nimgs):
                            ftime.append(start + dt * i)
                            exptime.append(dt)
                            pass
                        pass
                    else:
                        ftime.append(start)
                        exptime.append(hdu.header['EXPTIME'])
                        pass
                    pass
                data['LOC'].append(loc)
                data['TIME'].extend(ftime)
                data['EXPLEN'].extend(exptime)
                pass
            pass
        pass
    elif ('WFC3' in ext) and ('SCAN' in ext):
        for loc in sorted(clc['LOC']):
            fullloc = os.path.join(dbs, loc)
            with pyfits.open(fullloc) as hdulist:
                header0 = hdulist[0].header
                if 'SCAN_ANG' in header0:
                    data['SCANANGLE'].append(header0['SCAN_ANG'])
                elif 'PA_V3' in header0:
                    data['SCANANGLE'].append(header0['PA_V3'])
                else:
                    data['SCANANGLE'].append(666)
                ftime = []
                for fits in hdulist:
                    if (fits.size != 0) and ('DELTATIM' in fits.header.keys()):
                        ftime.append(float(fits.header['ROUTTIME']))
                        pass
                    pass
                data['LOC'].append(loc)
                data['TIME'].append(np.nanmean(ftime))
                data['EXPLEN'].append(header0['EXPTIME'])
                pass
            pass
        pass
    elif 'STARE' in ext:
        for loc in sorted(clc['LOC']):
            fullloc = os.path.join(dbs, loc)
            with pyfits.open(fullloc) as hdulist:
                header0 = hdulist[0].header
                if 'SCAN_ANG' in header0:
                    scanangle = header0['SCAN_ANG']
                elif 'PA_V3' in header0:
                    scanangle = header0['PA_V3']
                else:
                    scanangle = 666
                ftime = []
                allexplen = []
                allloc = []
                for fits in hdulist:
                    if (fits.size != 0) and (fits.header['EXTNAME'] in ['SCI']):
                        ftime.append(float(fits.header['EXPEND']))
                        allexplen.append(float(fits.header['EXPTIME']))
                        allloc.append(fits.header['EXPNAME'])
                        pass
                    pass
                allscanangle = [scanangle] * len(allexplen)
                data['SCANANGLE'].extend(allscanangle)
                data['LOC'].extend(allloc)
                data['TIME'].extend(ftime)
                data['EXPLEN'].extend(allexplen)
                pass
            pass
        pass
    data['IGNORED'] = [False] * len(data['TIME'])
    time = np.array(data['TIME'].copy())
    ignore = np.array(data['IGNORED'].copy())
    exposlen = np.array(data['EXPLEN'].copy())
    scanangle = np.array(data['SCANANGLE'].copy())
    ordt = np.argsort(time)
    exlto = exposlen.copy()[ordt]
    tmeto = time.copy()[ordt]
    ignto = ignore.copy()[ordt]
    if 'HST' in ext:
        scato = scanangle.copy()[ordt]
    if tmeto.size > 1:
        timingplist = [
            p for p in priors['planets'] if p not in force['pignore']
        ]
        for p in timingplist:
            out['data'][p] = {}
            if 'JWST' in ext:  # JWST --------------------------------
                smaors = priors[p]['sma'] / priors['R*'] / ssc['Rsun/AU']
                tmjd = priors[p]['t0']
                if tmjd > 2400000.5:
                    tmjd -= 2400000.5
                z, phase = time2z(
                    time,
                    priors[p]['inc'],
                    tmjd,
                    smaors,
                    priors[p]['period'],
                    priors[p]['ecc'],
                )
                zto = z.copy()[ordt]
                phsto = phase.copy()[ordt]
                # VISIT NUMBERING
                visto = np.ones(time.size)
                vis = np.ones(time.size)
                # think of something else for phasecurves
                wherev = np.where(np.diff(phsto) < 0)[0]
                for index in wherev:
                    visto[index + 1 :] += 1
                # TRANSIT VISIT PHASECURVE
                out['data'][p]['transit'] = []
                out['data'][p]['eclipse'] = []
                out['data'][p]['phasecurve'] = []
                for v in set(visto):
                    selv = visto == v
                    trlim = 1e0
                    posphsto = phsto.copy()
                    posphsto[posphsto < 0] = posphsto[posphsto < 0] + 1e0
                    tecrit = abs(np.arcsin(trlim / smaors)) / (2e0 * np.pi)
                    select = abs(zto[selv]) < trlim
                    pcconde = False
                    if np.any(select) and (
                        (np.min(abs(posphsto[selv][select] - 0.5)) < tecrit)
                    ):
                        out['eclipse'].append(int(v))
                        out['data'][p]['eclipse'].append(int(v))
                        pcconde = True
                        pass
                    pccondt = False
                    if np.any(select) and (
                        np.min(abs(phsto[selv][select])) < tecrit
                    ):
                        out['transit'].append(int(v))
                        out['data'][p]['transit'].append(int(v))
                        pccondt = True
                        pass
                    if pcconde and pccondt:
                        out['phasecurve'].append(int(v))
                        out['data'][p]['phasecurve'].append(int(v))
                        pass
                    pass
                vis[ordt] = visto.astype(int)
                ignore[ordt] = ignto
                out['data'][p]['wherev'] = wherev
                out['data'][p]['visits'] = vis
                out['data'][p]['z'] = z
                out['data'][p]['phase'] = phase
                out['data'][p]['ordt'] = ordt
                out['data'][p]['ignore'] = ignore
                out['STATUS'].append(True)
                pass
            elif 'Spitzer' in ext:  # SPITZER ------------------------
                out['data'][p]['transit'] = []
                out['data'][p]['eclipse'] = []
                out['data'][p]['phasecurve'] = []
                vis = np.ones(time.size)

                smaors = priors[p]['sma'] / priors['R*'] / ssc['Rsun/AU']
                tdur = priors[p]['period'] / (np.pi) / smaors  # rough estimate
                pdur = tdur / priors[p]['period']

                # https://arxiv.org/pdf/1001.2010.pdf eq 33
                w = priors[p].get('omega', 0)
                tme = priors[p]['t0'] + priors[p]['period'] * 0.5 * (
                    1 + priors[p]['ecc'] * (4.0 / np.pi) * np.cos(np.deg2rad(w))
                )
                tm = priors[p]['t0']
                if tm > 2400000.5 and 'Spitzer' in ext:
                    tme -= 2400000.5
                    tm = priors[p]['t0'] - 2400000.5

                rsun = 6.955e8  # m
                au = 1.496e11  # m
                period = priors[p]['period']
                offset = 0.25
                tphase = (time - tm + offset * period) / period
                ephase = (time - tme + offset * period) / period
                pdur = (
                    2
                    * np.arctan(priors['R*'] * rsun / (priors[p]['sma'] * au))
                    / (2 * np.pi)
                )
                events = np.unique(np.floor(tphase))
                min_images = 200
                for e in events:  # loop through different years
                    tmask = (tphase > e + (offset - 2 * pdur)) & (
                        tphase < e + (offset + 2 * pdur)
                    )
                    mmask = (tphase > e + (offset - 2 * pdur + 0.25)) & (
                        tphase < e + (offset + 2 * pdur + 0.25)
                    )
                    emask = (ephase > e + (offset - 2 * pdur)) & (
                        ephase < e + (offset + 2 * pdur)
                    )
                    # emask2 = (ephase > (e-1) + (offset-2*pdur)) &
                    # (ephase < (e-1) + (offset+2*pdur))  # eclipse at prior orbit
                    # tmask2 = (tphase > (e+1) + (offset-2*pdur)) &
                    # (tphase < (e+1) + (offset+2*pdur) )  # transit at next orbit
                    if tmask.sum() > min_images:
                        out['data'][p]['transit'].append(e)
                    if emask.sum() > min_images:
                        out['data'][p]['eclipse'].append(e)
                    if (
                        (tmask.sum() > min_images)
                        and (emask.sum() > min_images)
                        and (mmask.sum() > min_images)
                    ):
                        out['data'][p]['phasecurve'].append(e)
                        pass
                    pass
                visto = np.floor(tphase)
                out['STATUS'].append(True)
                pass
            elif 'HST' in ext:  # HST --------------------------------
                smaors = priors[p]['sma'] / priors['R*'] / ssc['Rsun/AU']
                tmjd = priors[p]['t0']
                if tmjd > 2400000.5:
                    tmjd -= 2400000.5
                z, phase = time2z(
                    time,
                    priors[p]['inc'],
                    tmjd,
                    smaors,
                    priors[p]['period'],
                    priors[p]['ecc'],
                )
                zto = z.copy()[ordt]
                phsto = phase.copy()[ordt]
                tmetod = [np.diff(tmeto)[0]]
                tmetod.extend(list(np.diff(tmeto)))
                tmetod = np.array(tmetod)
                thrs = np.percentile(tmetod, 75)
                cftfail = tmetod > 3 * thrs
                if True in cftfail:
                    thro = np.percentile(tmetod[cftfail], 75)
                else:
                    thro = 0
                # THRESHOLDS
                rbtthr = 25e-1 * thrs  # HAT-P-11
                vstthr = 3e0 * thro
                # VISIT NUMBERING
                whereo = np.where(tmetod > rbtthr)[0]
                wherev = np.where(tmetod > vstthr)[0]
                visto = np.ones(tmetod.size)
                dvis = np.ones(tmetod.size)
                vis = np.ones(tmetod.size)
                for index in wherev:
                    visto[index:] += 1
                # DOUBLE SCAN VISIT RE NUMBERING
                dvisto = visto.copy()
                for v in set(visto):
                    selv = visto == v
                    vordsa = scato[selv].copy()
                    if len(set(vordsa)) > 1:
                        dvisto[visto > v] = dvisto[visto > v] + 1
                        dbthr = np.mean(list(set(vordsa)))
                        vdbvisto = dvisto[selv].copy()
                        vdbvisto[vordsa > dbthr] = vdbvisto[vordsa > dbthr] + 1
                        dvisto[selv] = vdbvisto
                        pass
                    pass
                # ORBIT NUMBERING
                orbto = np.ones(tmetod.size)
                orb = np.ones(tmetod.size)
                for v in set(visto):
                    selv = visto == v
                    if len(~ignto[selv]) < 4:
                        ignto[selv] = True
                    else:
                        select = np.where(tmetod[selv] > rbtthr)[0]
                        incorb = orbto[selv]
                        for indice in select:
                            incorb[indice:] = incorb[indice:] + 1
                        orbto[selv] = incorb
                        for o in set(orbto[selv]):
                            selo = orbto[selv] == o
                            if len(~ignto[selv][selo]) < 4:
                                visignto = ignto[selv]
                                visignto[selo] = True
                                ignto[selv] = visignto
                                pass
                            ref = np.median(exlto[selv][selo])
                            if len(set(exlto[selv][selo])) > 1:
                                rej = exlto[selv][selo] != ref
                                ovignto = ignto[selv][selo]
                                ovignto[rej] = True
                                ignto[selv][selo] = ovignto
                                pass
                            pass
                        pass
                    pass
                # TRANSIT VISIT PHASECURVE
                out['data'][p]['svntransit'] = []
                out['data'][p]['svneclipse'] = []
                out['data'][p]['svnphasecurve'] = []
                for v in set(visto):
                    selv = visto == v
                    trlim = 1e0
                    posphsto = phsto.copy()
                    posphsto[posphsto < 0] = posphsto[posphsto < 0] + 1e0
                    tecrit = abs(np.arcsin(trlim / smaors)) / (2e0 * np.pi)
                    select = abs(zto[selv]) < trlim
                    pcconde = False
                    if np.any(select) and (
                        np.min(abs(posphsto[selv][select] - 0.5)) < tecrit
                    ):
                        out['eclipse'].append(int(v))
                        out['data'][p]['svneclipse'].append(int(v))
                        pcconde = True
                        pass
                    pccondt = False
                    if np.any(select) and (
                        np.min(abs(phsto[selv][select])) < tecrit
                    ):
                        out['transit'].append(int(v))
                        out['data'][p]['svntransit'].append(int(v))
                        pccondt = True
                        pass
                    if pcconde and pccondt:
                        out['phasecurve'].append(int(v))
                        out['data'][p]['svnphasecurve'].append(int(v))
                        pass
                    pass
                out['data'][p]['transit'] = []
                out['data'][p]['eclipse'] = []
                out['data'][p]['phasecurve'] = []
                for v in set(dvisto):
                    selv = dvisto == v
                    trlim = 1e0
                    posphsto = phsto.copy()
                    posphsto[posphsto < 0] = posphsto[posphsto < 0] + 1e0
                    tecrit = abs(np.arcsin(trlim / smaors)) / (2e0 * np.pi)
                    select = abs(zto[selv]) < trlim
                    pcconde = False
                    if np.any(select) and (
                        np.min(abs(posphsto[selv][select] - 0.5)) < tecrit
                    ):
                        out['data'][p]['eclipse'].append(int(v))
                        pcconde = True
                        pass
                    pccondt = False
                    if np.any(select) and (
                        np.min(abs(phsto[selv][select])) < tecrit
                    ):
                        out['data'][p]['transit'].append(int(v))
                        pccondt = True
                        pass
                    if pcconde and pccondt:
                        out['data'][p]['phasecurve'].append(int(v))
                    pass
                vis[ordt] = visto.astype(int)
                orb[ordt] = orbto.astype(int)
                dvis[ordt] = dvisto.astype(int)
                ignore[ordt] = ignto
                out['data'][p]['tmetod'] = tmetod
                out['data'][p]['whereo'] = whereo
                out['data'][p]['wherev'] = wherev
                out['data'][p]['thrs'] = rbtthr
                out['data'][p]['thro'] = vstthr
                out['data'][p]['visits'] = vis
                out['data'][p]['orbits'] = orb
                out['data'][p]['dvisits'] = dvis
                out['data'][p]['z'] = z
                out['data'][p]['phase'] = phase
                out['data'][p]['ordt'] = ordt
                out['data'][p]['ignore'] = ignore
                out['STATUS'].append(True)
                pass
            log.info('>-- Planet: %s', p)
            log.info('--< Transit: %s', str(out['data'][p]['transit']))
            log.info('--< Eclipse: %s', str(out['data'][p]['eclipse']))
            log.info('--< Phase Curve: %s', str(out['data'][p]['phasecurve']))
            if (
                out['data'][p]['transit']
                or out['data'][p]['eclipse']
                or out['data'][p]['phasecurve']
            ):
                chunked = True
            pass
        pass
    return chunked


# ------------ -------------------------------------------------------
# -- JWST CALIBRATION -- ---------------------------------------------
def rampfits(raws, sb=False, nl=False, alldq=None, verbose=False):
    '''
    G. ROUDIER: Ramp fits
    '''
    out = raws.copy()
    alldexp = raws['alldexp'].copy()
    allsb = []
    alldqnl = []
    if sb:
        argsdict = {
            'progbar': verbose,
            'progsizemax': 35,
            'lbllen': 15,
            'proginprompt': True,
        }
        if nl:
            progbar = nerdclub.Progressbar(argsdict, '>-- SB + NL', alldexp)
            pass
        else:
            progbar = nerdclub.Progressbar(argsdict, '>-- SB', alldexp)
            pass
        for i in range(len(raws['alldexp'])):
            sbmap = jwstsbfiles(raws['alldet'][i], raws['allgrating'][i])
            alldexp[i] = raws['alldexp'][i] - sbmap
            allsb.append(sbmap)
            if nl:
                nlmap = jwstnlfiles(raws['alldet'][i])
                domain = raws['allstartpix'][i]
                coeffs = [
                    n[
                        domain[1] - 1 : domain[1] - 1 + sbmap.shape[0],
                        domain[0] - 1 : domain[0] - 1 + sbmap.shape[1],
                    ]
                    for n in nlmap[0]
                ]
                alldqnl.append(
                    ~np.bool(
                        nlmap[1][
                            domain[1] - 1 : domain[1] - 1 + sbmap.shape[0],
                            domain[0] - 1 : domain[0] - 1 + sbmap.shape[1],
                        ]
                    )
                )
                powers = np.arange(len(coeffs))
                alldexp[i] = np.sum(
                    [c * (alldexp[i] ** p) for c, p in zip(coeffs, powers)],
                    axis=0,
                )
                pass
            progbar.update()
            pass
        progbar.close()
        pass
    # TIMING FORMAT
    # TTYPE1  = 'integration_number'
    # TFORM1  = 'J       '
    # TTYPE2  = 'int_start_MJD_UTC'
    # TFORM2  = 'D       '
    # TTYPE3  = 'int_mid_MJD_UTC'
    # TFORM3  = 'D       '
    # TTYPE4  = 'int_end_MJD_UTC'
    # TFORM4  = 'D       '
    # TTYPE5  = 'int_start_BJD_TDB'
    # TFORM5  = 'D       '
    # TTYPE6  = 'int_mid_BJD_TDB'
    # TFORM6  = 'D       '
    # TTYPE7  = 'int_end_BJD_TDB'
    # TFORM7  = 'D       '
    exposures = []
    err = []
    dq = []  # valid data flag

    argsdict = {
        'progbar': verbose,
        'progsizemax': 35,
        'lbllen': 15,
        'proginprompt': True,
    }
    progbar = nerdclub.Progressbar(argsdict, '>-- RAMPFITS', alldexp)
    for i, _ in enumerate(alldexp):
        # VECTORIAL FIT
        dramps = [r.flatten() for r in alldexp[i]]
        # TIMESTAMPS
        deltat = (
            (raws['alltiming'][i][-1] - raws['alltiming'][i][-3]) * 24e0 * 36e2
        ) / (
            1e0 * len(alldexp[i])
        )  # [s]
        xfit = deltat * np.arange(len(dramps))
        # LINFIT
        fitresult = np.polyfit(xfit, np.array(dramps), 1, cov=True, full=False)
        exposures.append(fitresult[0][0].reshape(alldexp[i][0].shape))
        cov = fitresult[1][0, 0].reshape(alldexp[i][0].shape)
        err.append(np.sqrt(abs(cov)))
        dqmap = np.bool(err[-1])  # zero error means trouble
        dqmap = dqmap.astype(float)  # pylint: disable=E1101
        if sb:
            thrplus = np.nanpercentile(allsb[i], 50 + 68 / 2)
            sigma = thrplus - np.nanmedian(allsb[i])
            dqmap = (
                dqmap
                * (abs(allsb[i] - np.nanmedian(allsb[i])) < 3 * sigma).astype(
                    float
                )
                * alldq[i]
            )
            if nl:
                dqmap = dqmap * alldqnl[i]
                pass
            pass
        dqmap[dqmap < 1] = np.nan
        dq.append(dqmap)
        progbar.update()
        pass
    progbar.close()
    out['alldexp'] = exposures
    out['allerr'] = err
    out['alldq'] = dq
    return out


def getscore(expts):
    '''
    G. ROUDIER: Time variations of a patch of the frame
    Assuming no full frame, it is the upper band, 2040 wide
    '''
    patch = np.nanstd([expo[0:4, 4:-4] for expo in expts], axis=0)
    return np.nanmedian(patch)


def readfitsdata(
    loclist, dbs, raws=False, sb=False, nl=False, alldq=None, verbose=False
):
    '''
    G. ROUDIER: Creates a dictionnary of time series of
    data of interest for JWST datasets from fits files
    '''
    out = {
        'alldet': [],
        'alldexp': [],
        'allunits': [],
        'allerr': [],
        'alldq': [],
        'allwaves': [],
        'alltiming': [],
        'allfilter': [],
        'allgrating': [],
        'allslit': [],
        'allread': [],
        'allstartpix': [],
    }
    for loc in loclist:
        fullloc = os.path.join(dbs, loc)
        with pyfits.open(fullloc) as hdulist:
            nints = None
            for hdu in hdulist:
                if 'PRIMARY' in hdu.name:
                    nints = 1 + hdu.header['INTEND'] - hdu.header['INTSTART']
                    out['alldet'].extend(nints * [hdu.header['DETECTOR']])
                    out['allfilter'].extend(nints * [hdu.header['FILTER']])
                    out['allgrating'].extend(nints * [hdu.header['GRATING']])
                    out['allslit'].extend(nints * [hdu.header['FXD_SLIT']])
                    out['allread'].extend(nints * [hdu.header['READPATT']])
                    # (X,Y)
                    out['allstartpix'].extend(
                        nints
                        * [(hdu.header['SUBSTRT1'], hdu.header['SUBSTRT2'])]
                    )
                    pass
                elif 'SCI' in hdu.name:
                    fitsdata = np.empty(hdu.data.shape)
                    fitsdata[:] = hdu.data[:]
                    out['alldexp'].extend(fitsdata)
                    out['allunits'].extend(
                        [hdu.header['BUNIT']] * len(fitsdata)
                    )
                    del hdu.data
                    pass
                # <-- L2b data only
                elif 'ERR' in hdu.name:
                    fitsdata = np.empty(hdu.data.shape)
                    fitsdata[:] = hdu.data[:]
                    out['allerr'].extend(fitsdata)
                    del hdu.data
                    pass
                elif 'DQ' in hdu.name:
                    fitsdata = np.empty(hdu.data.shape)
                    fitsdata[:] = hdu.data[:]
                    out['alldq'].extend(fitsdata)
                    del hdu.data
                    pass
                elif ('WAVELENGTH' in hdu.name) and nints:
                    fitsdata = np.empty(hdu.data.shape)
                    fitsdata[:] = hdu.data[:]
                    out['allwaves'].extend(nints * [fitsdata])
                    del hdu.data
                    pass
                # -->
                elif 'INT_TIMES' in hdu.name:
                    fitsdata = hdu.data[:].copy()
                    out['alltiming'].extend(fitsdata)
                    del hdu.data
                    pass
                pass
            pass
        pass
    if raws:
        out = rampfits(out, sb=sb, nl=nl, alldq=alldq, verbose=verbose)
        pass
    return out


def jwstcal(fin, clc, tim, ext, out, verbose=False, fastdev=False):
    '''
    G. ROUDIER: Extracts and Wavelength calibrates JWST datasets
    '''
    dbs = os.path.join(dawgie.context.data_dbs, 'mast')
    # TEMP CI
    _ = tim
    # RAW VERSUS CALIBRATED LISTS
    rawloc = [
        loc
        for loc, n in zip(clc['LOC'], clc['ROOTNAME'])
        if n.endswith('uncal')
    ]
    calloc = [
        loc
        for loc, n in zip(clc['LOC'], clc['ROOTNAME'])
        if n.endswith('calints')
    ]
    out['data']['LOC'] = rawloc
    if fastdev:
        rawloc = rawloc[0:2]
        calloc = calloc[0:2]
        pass
    # DATASET
    rawdata = readfitsdata(rawloc, dbs, raws=True, verbose=verbose)
    caldata = readfitsdata(calloc, dbs, raws=False, verbose=verbose)
    # NRS1 STSCI clips ref.pixels and detector plate (28, 1271)
    # NRS2 STSCI (32, 2048)
    # RAWS (32, 2048)
    datatiming = np.array([t[5] for t in rawdata['alltiming']])  # [Days]
    allrexp = np.array(rawdata['alldexp'])
    alldet = np.array(rawdata['alldet'])
    allwaves = caldata['allwaves'].copy()

    # CALIBRATION STEPS
    allscores = {'0 RAW': getscore(allrexp)}
    # 1 - Data quality initialization
    alldq = np.array(rawdata['alldq'])
    for i, _ in enumerate(alldq):
        if alldq[i].shape != caldata['alldq'][i].shape:
            # Play with allstartpix
            pass
        else:
            alldq[i] = ~np.bool(caldata['alldq'][i])
            pass
        pass
    alldq = alldq.astype(float)
    alldq[alldq < 1] = np.nan
    allscores['1 DQ'] = getscore(allrexp * alldq)
    # 2 - Saturation check
    # Useless. Either the NL files are reliable or they re not.
    # There s no saturation for a detector. Only non linearities.
    # 2.1 - Superbias
    rawdata = readfitsdata(
        rawloc, dbs, raws=True, sb=True, alldq=alldq, verbose=verbose
    )
    allrexp = np.array(rawdata['alldexp'])
    alldq = np.array(rawdata['alldq'])
    allscores['2 SB'] = getscore(allrexp * alldq)
    # 2.2 - Linearity correction
    rawdata = readfitsdata(
        rawloc, dbs, raws=True, sb=True, nl=True, alldq=alldq, verbose=verbose
    )
    allrexp = np.array(rawdata['alldexp'])
    alldq = np.array(rawdata['alldq'])
    allscores['3 NL'] = getscore(allrexp * alldq)
    # 2.3 - Persistence correction [IM]
    # 2.4 - Dark subtraction [IM]
    # 3 - Reference pixel correction [Exposure Low Freq Noise]
    allrexp = lfnoise(allrexp, alldq, verbose=verbose)
    allscores['4 LFN'] = getscore(allrexp * alldq)
    out['data']['SCORES'] = allscores
    # 4 - Jump detection [that is a joke when we have 2, 3 or 4 groups...]

    if verbose:
        scores = [out['data']['SCORES'][k] for k in out['data']['SCORES']]
        labels = out['data']['SCORES'].keys()
        plt.figure(figsize=(12, 9))
        plt.plot(scores, 'o--')
        plt.xticks(
            ticks=np.arange(len(labels)),
            labels=labels,
            fontsize=16,
            rotation=45,
        )
        plt.yticks(fontsize=16)
        plt.ylabel('DN', fontsize=20)
        plt.show()
        pass

    reffile = jwstreffiles(ext)

    if 'NIRISS' in ext:
        Tstar = fin['priors']['T*']
        bbfunc = astrobb(Tstar * astropy.units.K)
        # NIRISS
        # reffile[0]: images of the 3 orders [296, 2088] +20 on each side
        # reffile[1]: calibrated wavelength, pixel XY, throughput for the 3 orders
        refwave = []  # Reference wavelength for each order
        refX = []  # X reference
        refY = []  # Y reference
        refT = []  # Throughput curve
        refB = []  # Blackbody model curve
        refS = []  # Blackbody*Throughput curve
        for order in np.arange(3):
            refwave.append(reffile[1][order]['WAVELENGTH'])
            refX.append(reffile[1][order]['X'] - 20.0)
            refY.append(reffile[1][order]['Y'] - 20.0)
            refT.append(reffile[1][order]['THROUGHPUT'])
            bbstar = bbfunc(refwave[-1] * astropy.units.micron)
            refB.append(bbstar)
            refS.append(refT[-1] * refB[-1] / np.nansum(refB[-1]))
            pass
        YY, XX = np.mgrid[0 : allrexp[0].shape[0], 0 : allrexp[0].shape[1]]
        # Mixing Matrix
        MM = list(reffile[0])
        # Transforms
        TMX = []
        TMY = []
        for thisrefX, thisrefS in zip(refX, refS):
            orderme = np.argsort(thisrefX)
            TMX.append(itp.CubicSpline(thisrefX[orderme], thisrefS[orderme]))
        # STOP TEST JWST
        _ = YY
        _ = XX
        _ = MM
        _ = TMY
        pass
    elif 'NIRSPEC' in ext:
        all1d = []
        all1dwave = []
        all1dvalid = []
        argsdict = {
            'progbar': verbose,
            'progsizemax': 35,
            'lbllen': 15,
            'proginprompt': True,
        }
        progbar = nerdclub.Progressbar(argsdict, '>-- WAVECAL', allrexp)
        for it, thisexp in enumerate(allrexp):
            if thisexp.shape != allwaves[it].shape:
                # NRS1: DUDES...
                # (32, 2048) versus (28, 1271)
                bs = np.zeros(thisexp.shape) * np.nan
                bs[:-4, -allwaves[it].shape[1] - 4 : -4] = allwaves[it]
                refwave = bs
                pass
            else:
                refwave = allwaves[it]
                pass
            wct = nirspeccal(thisexp, refwave)
            all1dvalid.append(wct[0])
            all1d.append(wct[1])
            all1dwave.append(wct[2])
            progbar.update()
            pass
        progbar.close()
        clean1D = cleanspec(all1d, all1dvalid, alldet)
        if verbose:
            plt.figure(figsize=(12, 9))
            for w, s, m in zip(all1dwave, all1d, all1dvalid):
                plt.plot(w[m], s[m])
                pass
            plt.figure(figsize=(12, 9))
            for w, s, m in zip(all1dwave, clean1D, all1dvalid):
                plt.plot(w[m], s[m])
                pass
            plt.show()
            pass
        out['STATUS'].append(True)
        out['data']['TIME'] = datatiming
        out['data']['DET'] = alldet
        out['data']['EXCLNUM'] = all1dvalid
        out['data']['IGNORED'] = [
            np.sum(x) < int(len(all1d[0]) / 2) for x in all1dvalid
        ]
        out['data']['SPECTRUM'] = all1d
        out['data']['WAVE'] = all1dwave
        pass
    return True


def cleanspec(all1d, all1dvalid, alldet):
    out = np.array(all1d.copy())
    valid = np.array(all1dvalid.copy())
    out[~valid] = np.nan
    for x in np.unique(alldet):
        select = alldet == x
        ref = np.nanmedian(out[select], axis=0)
        refup = np.nanpercentile(out[select], 50 - 34, axis=0)
        refdown = np.nanpercentile(out[select], 50 + 34, axis=0)
        detset = out[select].copy()
        detset[(out[select] > refup) | (out[select] < refdown)] = np.nan
        stdsample = np.nanstd(detset, axis=0)
        valid[select] = abs(out[select] - ref) < (3 * stdsample)
        pass
    out[~valid] = np.nan
    return out


def nirspeccal(thisexp, thosewaves):
    '''
    GMR: NIRSPEC wavelength calibration
    Assumes lfer removed data
    '''
    this1d = np.nansum(thisexp, axis=0)
    this1dwave = np.nanmedian(thosewaves, axis=0)
    select = np.isfinite(this1d)
    return (select, this1d, this1dwave)


def jwstdqfiles(thisdet):
    '''
    G. ROUDIER: Returns a mask for DQ init
    Source: https://jwst-crds.stsci.edu
    Local: /proj/sdp/data/cal/
    DQ definition: https://jwst-pipeline.readthedocs.io/en/stable/jwst/references_general/references_general.html#data-quality-flags
    Find a way to fetch the latest version
    '''
    local = ''
    fpath = ''
    if thisdet.startswith('NRS'):
        local = os.path.join(excalibur.context['data_cal'], 'NIRSPEC')
        pass
    if thisdet in ['NRS1']:
        fpath = os.path.join(local, 'jwst_nirspec_mask_0085.fits')
        pass
    if thisdet in ['NRS2']:
        fpath = os.path.join(local, 'jwst_nirspec_mask_0087.fits')
        pass
    with pyfits.open(fpath) as prfhdul:
        mask = []
        for hdu in prfhdul:
            if hdu.data is not None:
                fitsdata = np.empty(hdu.data.shape)
                fitsdata[:] = hdu.data[:]
                mask.append(fitsdata)
                del hdu.data
                pass
            pass
        pass
    # mask[1] contains bit definition
    return mask[0]  # Full array...


def jwstsbfiles(thisdet, thisgrating):
    '''
    G. ROUDIER: Returns a list of reference files for SuperBias sub
    Source: https://jwst-crds.stsci.edu
    Local: /proj/sdp/data/cal/
    Find a way to fetch the latest version
    '''
    local = ''
    fpath = ''
    if thisdet.startswith('NRS'):
        local = os.path.join(excalibur.context['data_cal'], 'NIRSPEC')
        pass
    if thisdet in ['NRS1'] and thisgrating in ['G395H']:
        fpath = os.path.join(local, 'jwst_nirspec_superbias_0479.fits')
        pass
    if thisdet in ['NRS2'] and thisgrating in ['G395H']:
        fpath = os.path.join(local, 'jwst_nirspec_superbias_0474.fits')
        pass
    with pyfits.open(fpath) as prfhdul:
        sb = []
        for hdu in prfhdul:
            if hdu.data is not None:
                fitsdata = hdu.data[:].copy()
                sb.append(fitsdata)
                del hdu.data
                pass
            pass
        pass
    return sb[0]


def jwstnlfiles(thisdet):
    '''
    G. ROUDIER: Returns a list of reference files for NL correction
    Source: https://jwst-crds.stsci.edu
    Local: /proj/sdp/data/cal/
    Find a way to fetch the latest version
    '''
    local = ''
    fpath = ''
    if thisdet.startswith('NRS'):
        local = os.path.join(excalibur.context['data_cal'], 'NIRSPEC')
        pass
    if thisdet in ['NRS1']:
        fpath = os.path.join(local, 'jwst_nirspec_linearity_0022.fits')
        pass
    if thisdet in ['NRS2']:
        fpath = os.path.join(local, 'jwst_nirspec_linearity_0025.fits')
        pass
    with pyfits.open(fpath) as prfhdul:
        nl = []
        for hdu in prfhdul:
            if hdu.data is not None:
                fitsdata = hdu.data[:].copy()
                nl.append(fitsdata)
                del hdu.data
                pass
            pass
        pass
    return nl  # [Coeff, DQ, DQdef, ?]


def jwstreffiles(thisext):
    '''
    G. ROUDIER: Returns a list of reference files for wavelegnth calibration
    Source: https://jwst-crds.stsci.edu
    Local: /proj/sdp/data/cal/
    '''
    thisdir = None
    thistrace = None
    if 'NIRISS' in thisext:
        # Note: do something smarter to get latest files
        # --< DATA CALIBRATION: JWST-NIRISS-NIS-CLEAR-GR700XD >--
        thisdir = 'NIRISS'
        thisprofile = 'jwst_niriss_specprofile_0022.fits'
        thistrace = 'jwst_niriss_spectrace_0023.fits'
        pass
    if 'NIRSPEC' in thisext:
        # --< DATA CALIBRATION: JWST-NIRSPEC-NRS-F290LP-G395H >--
        thisdir = 'NIRSPEC'
        thisprofile = 'jwst_nirspec_g395h_disp_20160902193401.fits'
        thistrace = thisprofile
        thisdir = None
        pass
    if thisdir is None:
        reffiles = None
    else:
        fpath = os.path.join(
            excalibur.context['data_cal'], thisdir, thisprofile
        )
        # Trace template for each order (2D ARRAY)
        with pyfits.open(fpath) as prfhdul:
            orders = []
            for hdu in prfhdul:
                if hdu.data is not None:
                    fitsdata = np.empty(hdu.data.shape)
                    fitsdata[:] = hdu.data[:]
                    orders.append(fitsdata)
                    del hdu.data
                    pass
                pass
            pass
        # Trace pixel (waves[0]['X'], waves[0]['Y']) and
        # associated calibrated wavelength (waves[0]['WAVELENGTH']) (1D ARRAYS)
        fpath = os.path.join(excalibur.context['data_cal'], thisdir, thistrace)
        with pyfits.open(fpath) as prfhdul:
            waves = []
            for hdu in prfhdul:
                if hdu.data is not None:
                    fitsdata = np.empty(hdu.data.shape)
                    fitsdata[:] = hdu.data[:]
                    waves.append(fitsdata)
                    del hdu.data
                    pass
                pass
            pass
        reffiles = [orders, waves]
        pass
    return reffiles


def lfnoise(xps, flg, verbose=False):
    '''
    G. ROUDIER: Low Frequency Noise Excess Removal
    '''
    argsdict = {
        'progbar': verbose,
        'progsizemax': 35,
        'lbllen': 15,
        'proginprompt': True,
    }
    progbar = nerdclub.Progressbar(argsdict, '>-- LFNER', xps)
    out = []
    for i, _ in enumerate(xps):
        lfcorr = xps[i].copy()
        lfcorr[4:-4] = np.nan
        lfr = np.nanmean(lfcorr, axis=0)
        select = ~np.isfinite(lfr)
        if np.sum(select):
            temp = flg[i]
            temp[:, select] = np.nan
            flg[i] = temp
            pass
        out.append(xps[i] - lfr)
        progbar.update()
        pass
    progbar.close()
    return np.array(out)


# ---------------------- ---------------------------------------------
# -- CALIBRATE SCAN DATA -- ------------------------------------------
def scancal(
    clc,
    tim,
    tid,
    flttype,
    out,
    emptythr=1e3,
    frame2png=False,
    debug=False,
):
    '''
    G. ROUDIER: Extracts and Wavelength calibrates WFC3 SCAN mode spectra
    '''
    # VISIT ------------------------------------------------------------------------------
    for pkey in tim['data'].keys():
        visits = np.array(tim['data'][pkey]['visits'])
    for pkey in tim['data'].keys():
        dvisits = np.array(tim['data'][pkey]['dvisits'])
    # DATA TYPE --------------------------------------------------------------------------
    arcsec2pix = dps(flttype)
    vrange = validrange(flttype)
    wvrng, disper, ldisp, udisp = fng(flttype)
    spectrace = np.round((np.max(wvrng) - np.min(wvrng)) / disper)
    # LOAD DATA --------------------------------------------------------------------------
    dbs = os.path.join(dawgie.context.data_dbs, 'mast')
    data = {
        'LOC': [],
        'EPS': [],
        'DISPLIM': [ldisp, udisp],
        'SCANRATE': [],
        'SCANLENGTH': [],
        'SCANANGLE': [],
        'EXP': [],
        'EXPERR': [],
        'EXPFLAG': [],
        'VRANGE': vrange,
        'TIME': [],
        'EXPLEN': [],
        'MIN': [],
        'MAX': [],
        'TRIAL': [],
    }
    for loc in sorted(clc['LOC']):
        fullloc = os.path.join(dbs, loc)
        with pyfits.open(fullloc) as hdulist:
            header0 = hdulist[0].header
            test = header0['UNITCORR']
            eps = False
            if test in ['COMPLETE', 'PERFORM']:
                eps = True
            data['EPS'].append(eps)
            if 'SCAN_RAT' in header0:
                data['SCANRATE'].append(header0['SCAN_RAT'])
            else:
                data['SCANRATE'].append(np.nan)
            if 'SCAN_LEN' in header0:
                data['SCANLENGTH'].append(header0['SCAN_LEN'])
            else:
                data['SCANLENGTH'].append(np.nan)
            if 'SCAN_ANG' in header0:
                data['SCANANGLE'].append(header0['SCAN_ANG'])
            elif 'PA_V3' in header0:
                data['SCANANGLE'].append(header0['PA_V3'])
            else:
                data['SCANANGLE'].append(666)
            frame = []
            errframe = []
            dqframe = []
            ftime = []
            fmin = []
            fmax = []
            for fits in hdulist:
                if (fits.size != 0) and ('DELTATIM' in fits.header.keys()):
                    fitsdata = np.empty(fits.data.shape)
                    fitsdata[:] = fits.data[:]
                    frame.append(fitsdata)
                    ftime.append(float(fits.header['ROUTTIME']))
                    fmin.append(float(fits.header['GOODMIN']))
                    fmax.append(float(fits.header['GOODMAX']))
                    del fits.data
                    pass
                if 'EXTNAME' in fits.header:
                    if fits.header['EXTNAME'] in ['ERR', 'DQ']:
                        fitsdata = np.empty(fits.data.shape)
                        fitsdata[:] = fits.data[:]
                        if fits.header['EXTNAME'] == 'ERR':
                            errframe.append(fitsdata)
                        if fits.header['EXTNAME'] == 'DQ':
                            dqframe.append(fitsdata)
                        del fits.data
                        pass
                    if eps and (fits.header['EXTNAME'] == 'TIME'):
                        xpsrl = np.array(float(fits.header['PIXVALUE']))
                        frame[-1] = frame[-1] * xpsrl
                        errframe[-1] = errframe[-1] * xpsrl
                        pass
                    pass
                pass
            data['LOC'].append(loc)
            data['EXP'].append(frame)
            data['EXPERR'].append(errframe)
            data['EXPFLAG'].append(dqframe)
            data['TIME'].append(ftime)
            data['EXPLEN'].append(header0['EXPTIME'])
            data['MIN'].append(fmin)
            data['MAX'].append(fmax)
            pass
        pass
    # MASK DATA --------------------------------------------------------------------------
    data['MEXP'] = data['EXP'].copy()
    data['MASK'] = data['EXPFLAG'].copy()
    data['IGNORED'] = [False] * len(data['LOC'])
    data['FLOODLVL'] = [np.nan] * len(data['LOC'])
    data['UP'] = [np.nan] * len(data['LOC'])
    data['DOWN'] = [np.nan] * len(data['LOC'])
    data['TRIAL'] = [''] * len(data['LOC'])
    for index, nm in enumerate(data['LOC']):
        maskedexp = []
        masks = []
        ignore = False
        for dd, ff in zip(data['EXP'][index], data['EXPFLAG'][index]):
            select = ff > 0
            if np.sum(select) > 0:
                dd[select] = np.nan
                if np.all(~np.isfinite(dd)):
                    data['TRIAL'][index] = 'Empty Subexposure'
                    ignore = True
                    pass
                else:
                    maskedexp.append(dd)
                pass
            else:
                maskedexp.append(dd)
            mm = np.isfinite(dd)
            masks.append(mm)
            pass
        if ignore:
            maskedexp = data['EXP'][index].copy()
        data['MEXP'][index] = maskedexp
        data['MASK'][index] = masks
        data['IGNORED'][index] = ignore
        pass
    # ALL FLOOD LEVELS -------------------------------------------------------------------
    for index, nm in enumerate(data['LOC']):
        ignore = data['IGNORED'][index]
        # MINKOWSKI ----------------------------------------------------------------------
        psdiff = np.diff(data['MEXP'][index][::-1].copy(), axis=0)
        floatsw = data['SCANLENGTH'][index] / arcsec2pix
        scanwdw = np.round(floatsw)
        if (scanwdw > psdiff[0].shape[0]) or (len(psdiff) < 2):
            scanwpi = np.round(floatsw / (len(psdiff)))
            pass
        else:
            scanwpi = np.round(floatsw / (len(psdiff) - 1))
        if scanwpi < 1:
            data['TRIAL'][index] = 'Subexposure Scan Length < 1 Pixel'
            ignore = True
            pass
        if not ignore:
            targetn = 0
            if tid in ['XO-2', 'HAT-P-1']:
                targetn = -1
            minlocs = []
            maxlocs = []
            floodlist = []
            for de, md in zip(psdiff[::-1], data['MIN'][index][::-1]):
                valid = np.isfinite(de)
                if np.nansum(~valid) > 0:
                    de[~valid] = 0
                select = de[valid] < md
                if np.nansum(select) > 0:
                    de[valid][select] = 0
                perfldlist = np.nanpercentile(de, np.arange(1001) / 1e1)
                perfldlist = np.diff(perfldlist)
                perfldlist[:100] = 0
                perfldlist[-1] = 0
                indperfld = list(perfldlist).index(np.max(perfldlist)) * 1e-1
                floodlist.append(np.nanpercentile(de, indperfld))
                pass
            fldthr = np.nanmax(floodlist)
            # CONTAMINATION FROM ANOTHER SOURCE IN THE UPPER FRAME -----------------------
            if tid in ['HAT-P-41']:
                fldthr /= 1.5
            if 'G102' in flttype:
                fldthr /= 3e0
            data['FLOODLVL'][index] = fldthr
            pass
        pass
    allfloodlvl = np.array(data['FLOODLVL'])
    for dv in set(dvisits):
        allfloodlvl[dvisits == dv] = np.nanmedian(allfloodlvl[dvisits == dv])
        pass
    data['FLOODLVL'] = allfloodlvl
    # ALL LIMITS  ------------------------------------------------------------------------
    for index, nm in enumerate(data['LOC']):
        ignore = data['IGNORED'][index]
        # MINKOWSKI FLOOD LEVEL ----------------------------------------------------------
        psdiff = np.diff(data['MEXP'][index][::-1].copy(), axis=0)
        floatsw = data['SCANLENGTH'][index] / arcsec2pix
        scanwdw = np.round(floatsw)
        if (scanwdw > psdiff[0].shape[0]) or (len(psdiff) < 2):
            scanwpi = np.round(floatsw / (len(psdiff)))
            pass
        else:
            scanwpi = np.round(floatsw / (len(psdiff) - 1))
        if scanwpi < 1:
            data['TRIAL'][index] = 'Subexposure Scan Length < 1 Pixel'
            ignore = True
            pass
        if not ignore:
            targetn = 0
            if tid in ['XO-2', 'HAT-P-1']:
                targetn = -1
            minlocs = []
            maxlocs = []
            for de, md in zip(psdiff.copy()[::-1], data['MIN'][index][::-1]):
                lmn, lmx = isolate(
                    de, md, spectrace, scanwpi, targetn, data['FLOODLVL'][index]
                )
                minlocs.append(lmn)
                maxlocs.append(lmx)
                pass
            # HEAVILY FLAGGED SCAN -------------------------------------------------------
            nanlocs = np.all(~np.isfinite(minlocs)) or np.all(
                ~np.isfinite(maxlocs)
            )
            almstare = scanwpi < 5
            if nanlocs or almstare:
                for de, md in zip(
                    psdiff.copy()[::-1], data['MIN'][index][::-1]
                ):
                    if (scanwpi / 3) < 2:
                        redscanwpi = scanwpi / 2
                    else:
                        redscanwpi = scanwpi / 3
                    lmn, lmx = isolate(
                        de,
                        md,
                        spectrace,
                        redscanwpi,
                        targetn,
                        data['FLOODLVL'][index],
                    )
                    minlocs.append(lmn)
                    maxlocs.append(lmx)
                    pass
                pass
            ignore = ignore or not (
                (np.any(np.isfinite(minlocs)))
                and (np.any(np.isfinite(maxlocs)))
            )
            if not ignore:
                minl = np.nanmin(minlocs)
                maxl = np.nanmax(maxlocs)
                # CONTAMINATION FROM ANOTHER SOURCE IN THE UPPER FRAME -------------------
                if (tid in ['HAT-P-41']) and ((maxl - minl) > 15):
                    minl = maxl - 15
                minl = max(minl, 10)
                maxl = min(psdiff[0].shape[0] - 10, maxl)
                pass
            else:
                minl = np.nan
                maxl = np.nan
                pass
            data['UP'][index] = minl
            data['DOWN'][index] = maxl
            pass
        pass
    allminl = np.array(data['UP'])
    allmaxl = np.array(data['DOWN'])
    for dv in set(dvisits):
        allminl[dv == dvisits] = np.nanpercentile(allminl[dv == dvisits], 2.5)
        allmaxl[dv == dvisits] = np.nanpercentile(allmaxl[dv == dvisits], 97.5)
        pass
    data['UP'] = allminl
    data['DOWN'] = allmaxl
    allscanlen = np.array(data['SCANLENGTH'])
    allignore = np.array(data['IGNORED'])
    alltrials = np.array(
        ['Exposure Scan Length Rejection'] * len(data['TRIAL'])
    )
    for v in set(visits):
        allfloodlvl[visits == v] = np.nanmedian(allfloodlvl[visits == v])
        visitign = allignore[visits == v]
        visittrials = alltrials[visits == v]
        select = allscanlen[visits == v] != np.nanmedian(
            allscanlen[visits == v]
        )
        visitign[select] = True
        visittrials[~select] = ''
        allignore[visits == v] = visitign
        alltrials[visits == v] = visittrials
        pass
    data['FLOODLVL'] = list(allfloodlvl)
    data['IGNORED'] = list(allignore)
    data['TRIAL'] = list(alltrials)
    ovszspc = False
    # BACKGROUND SUB AND ISOLATE ---------------------------------------------------------
    for index, nm in enumerate(data['LOC']):
        ignore = data['IGNORED'][index]
        psdiff = np.diff(data['MEXP'][index][::-1].copy(), axis=0)
        psminsel = np.array(data['MIN'][index]) < 0
        if True in psminsel:
            psmin = np.nansum(np.array(data['MIN'][index])[psminsel])
        else:
            psmin = np.nanmin(data['MIN'][index])
        minl = data['UP'][index]
        maxl = data['DOWN'][index]
        if not ignore:
            # BACKGROUND SUBTRACTION -----------------------------------------------------
            for eachdiff in psdiff:
                background = []
                for eachcol in eachdiff.T:
                    if True in np.isfinite(eachcol):
                        selfinite = np.isfinite(eachcol)
                        fineachcol = eachcol[selfinite]
                        test = fineachcol < psmin
                        if True in test:
                            fineachcol[test] = np.nan
                        eachcol[selfinite] = fineachcol
                        pass
                    nancounts = np.sum(~np.isfinite(eachcol))
                    thr = 1e2 * (1e0 - (scanwpi + nancounts) / eachcol.size)
                    if thr <= 0:
                        bcke = np.nan
                    else:
                        test = eachcol < np.nanpercentile(eachcol, thr)
                        if True in test:
                            bcke = np.nanmedian(eachcol[test])
                        else:
                            bcke = np.nan
                        pass
                    background.append(bcke)
                    pass
                background = np.array(
                    [np.array(background)] * eachdiff.shape[0]
                )
                eachdiff -= background
                eachdiff[: int(minl), :] = 0
                eachdiff[int(maxl) :, :] = 0
                pass
            # DIFF ACCUM -----------------------------------------------------------------
            thispstamp = np.nansum(psdiff, axis=0)
            thispstamp[thispstamp <= psmin] = np.nan
            thispstamp[thispstamp == 0] = np.nan
            if abs(spectrace - thispstamp.shape[1]) < 36:
                ovszspc = True
            # ISOLATE SCAN X -------------------------------------------------------------
            mltord = thispstamp.copy()
            targetn = 0
            if tid in ['XO-2']:
                targetn = -1
            minx, maxx = isolate(
                mltord,
                psmin,
                spectrace,
                scanwpi,
                targetn,
                data['FLOODLVL'][index],
                axis=0,
            )
            if np.isfinite(minx * maxx):
                minx -= 1.5 * 12
                maxx += 1.5 * 12
                if minx < 0:
                    minx = 5
                thispstamp[:, : int(minx)] = np.nan
                if maxx > (thispstamp.shape[1] - 1):
                    maxx = thispstamp.shape[1] - 5
                thispstamp[:, int(maxx) :] = np.nan
                if ((maxx - minx) < spectrace) and not ovszspc:
                    data['TRIAL'][index] = 'Could Not Find Full Spectrum'
                    ignore = True
                    pass
                pstamperr = np.array(data['EXPERR'][index].copy())
                select = ~np.isfinite(pstamperr)
                if np.nansum(select) > 0:
                    pstamperr[select] = 0
                pstamperr = np.sqrt(np.nansum(pstamperr**2, axis=0))
                select = ~np.isfinite(thispstamp)
                if np.nansum(select) > 0:
                    pstamperr[select] = np.nan
                pass
            else:
                garbage = data['EXP'][index][::-1].copy()
                thispstamp = np.sum(np.diff(garbage, axis=0), axis=0)
                pstamperr = thispstamp * np.nan
                data['TRIAL'][index] = 'Could Not Find X Edges'
                ignore = True
                pass
            pass
        else:
            garbage = data['EXP'][index][::-1].copy()
            thispstamp = np.sum(np.diff(garbage, axis=0), axis=0)
            pstamperr = thispstamp * np.nan
            if len(data['TRIAL'][index]) < 1:
                data['TRIAL'][index] = 'Could Not Find Y Edges'
                pass
            ignore = True
            pass
        data['MEXP'][index] = thispstamp
        data['TIME'][index] = np.nanmean(data['TIME'][index].copy())
        data['IGNORED'][index] = ignore
        data['EXPERR'][index] = pstamperr
        if debug:
            log.info('>-- %s / %s', str(index), str(len(data['LOC']) - 1))
        # PLOTS --------------------------------------------------------------------------
        if frame2png:
            if not os.path.exists('TEST'):
                os.mkdir('TEST')
            if not os.path.exists('TEST/' + tid):
                os.mkdir('TEST/' + tid)
            plt.figure()
            plt.title('Index: ' + str(index) + ' Ignored=' + str(ignore))
            plt.imshow(thispstamp)
            plt.colorbar()
            plt.savefig('TEST/' + tid + '/' + nm + '.png')
            plt.close()
            pass
        pass
    maxwasize = []
    for mexp in data['MEXP']:
        maxwasize.append(mexp.shape[1])
    maxwasize = np.nanmax(maxwasize)
    # SPECTRUM EXTRACTION ----------------------------------------------------------------
    data['SPECTRUM'] = [np.array([np.nan] * maxwasize)] * len(data['LOC'])
    data['SPECERR'] = [np.array([np.nan] * maxwasize)] * len(data['LOC'])
    data['NSPEC'] = [np.nan] * len(data['LOC'])
    for index, loc in enumerate(data['LOC']):
        floodlevel = data['FLOODLVL'][index]
        if floodlevel < emptythr:
            data['IGNORED'][index] = True
            data['TRIAL'][index] = 'Empty Frame'
            pass
        ignore = data['IGNORED'][index]
        if not ignore:
            frame = data['MEXP'][index].copy()
            frame = [line for line in frame if not np.all(~np.isfinite(line))]
            # OVERSIZED MASK -------------------------------------------------------------
            for line in frame:
                if np.nanmax(line) < floodlevel:
                    line *= np.nan
                if abs(spectrace - line.size) < 36:
                    ovszspc = True
                elif np.sum(np.isfinite(line)) < spectrace:
                    line *= np.nan
                pass
            frame = [line for line in frame if not np.all(~np.isfinite(line))]
            # SCAN RATE CORRECTION -------------------------------------------------------
            template = []
            for col in np.array(frame).T:
                if not np.all(~np.isfinite(col)):
                    template.append(np.nanmedian(col))
                else:
                    template.append(np.nan)
                pass
            template = np.array(template)
            for line in frame:
                errref = np.sqrt(abs(line)) / abs(template)
                line /= template
                refline = np.nanmedian(line)
                select = np.isfinite(line)
                minok = abs(line[select] - refline) < 3e0 * np.nanmin(
                    errref[select]
                )
                if np.nansum(minok) > 0:
                    alpha = np.nansum(minok) / np.nansum(line[select][minok])
                    line *= alpha
                    line *= template
                    pass
                else:
                    line *= np.nan
                pass
            frame = [line for line in frame if not np.all(~np.isfinite(line))]

            spectrum = []
            specerr = []
            nspectrum = []
            vtemplate = []
            for row in np.array(frame):
                if not np.all(~np.isfinite(row)):
                    vtemplate.append(np.nanmedian(row))
                else:
                    vtemplate.append(np.nan)
                pass
            vtemplate = np.array(vtemplate)
            for col in np.array(frame).T:
                ignorecol = False
                if not np.all(~np.isfinite(col)):
                    errref = np.sqrt(abs(np.nanmedian(col))) / abs(vtemplate)
                    ratio = col / vtemplate
                    refline = np.nanmedian(ratio)
                    select = np.isfinite(col)
                    ok = abs(ratio[select] - refline) < 3e0 * errref[select]
                    if np.nansum(ok) > 0:
                        alpha = np.nansum(ok) / np.nansum(ratio[select][ok])
                        valid = abs(
                            col[select] * alpha - vtemplate[select]
                        ) < 3 * np.sqrt(abs(vtemplate[select]))
                        pass
                    else:
                        valid = [False]
                    if np.nansum(valid) > 0:
                        spectrum.append(np.nanmedian(col[select][valid]))
                        specerr.append(np.nanstd(col[select][valid]))
                        nspectrum.append(np.nansum(valid))
                        pass
                    else:
                        ignorecol = True
                else:
                    ignorecol = True
                if ignorecol:
                    spectrum.append(np.nan)
                    specerr.append(np.nan)
                    nspectrum.append(0)
                    pass
                pass
            spectrum = np.array(spectrum)
            spectrum -= np.nanmin(spectrum)
            # EXCLUDE RESIDUAL GLITCHES
            seloutlrs = np.isfinite(template) & np.isfinite(spectrum)
            if True in seloutlrs:
                nanme = (
                    abs(spectrum[seloutlrs] - template[seloutlrs])
                    / template[seloutlrs]
                ) > 1e0
                if True in nanme:
                    spectrum[seloutlrs][nanme] = np.nan
                pass
            else:
                data['IGNORED'][index] = True
                data['TRIAL'][index] = 'Invalid Spectrum/Template'
                pass
            # TRUNCATED SPECTRUM
            # testspec = spectrum[np.isfinite(spectrum)]
            # if (np.all(testspec[-18:] > emptythr)) and not ovszspc:
            #    data['IGNORED'][index] = True
            #    data['TRIAL'][index] = 'Truncated Spectrum'
            #    pass
            data['SPECTRUM'][index] = np.array(spectrum)
            data['SPECERR'][index] = np.array(specerr)
            data['NSPEC'][index] = np.array(nspectrum)
            pass
        pass
    # WAVELENGTH CALIBRATION -------------------------------------------------------------
    wavett, tt = ag2ttf(flttype)
    if ovszspc:
        select = (wavett * 1e-4 < 1.68) & (wavett * 1e-4 > 1.09)
        wavett = wavett[select]
        tt = tt[select]
        pass
    scaleco = np.nanmax(tt) / np.nanmin(tt[tt > 0])
    data['PHT2CNT'] = [np.nan] * len(data['LOC'])
    data['WAVE'] = [np.array([np.nan] * maxwasize)] * len(data['LOC'])
    data['DISPERSION'] = [np.nan] * len(data['LOC'])
    data['SHIFT'] = [np.nan] * len(data['LOC'])
    data['BACKGROUND'] = [np.nan] * len(data['LOC'])
    spectralindex = []
    for index, loc in enumerate(data['LOC']):
        ignore = data['IGNORED'][index]
        if not ignore:
            spectrum = data['SPECTRUM'][index].copy()
            cutoff = np.nanmax(spectrum) / scaleco
            finitespec = spectrum[np.isfinite(spectrum)]
            test = finitespec < cutoff
            if True in test:
                finitespec[test] = np.nan
                spectrum[np.isfinite(spectrum)] = finitespec
                pass
            wave, disp, shift, si, bck = wavesol(
                abs(spectrum),
                tt,
                wavett,
                disper,
                ovszspc=ovszspc,
                bck=None,
            )
            if ldisp < disp < udisp:
                spectralindex.append(si)
            pass
        pass
    siv = np.nanmedian(spectralindex)
    for index, loc in enumerate(data['LOC']):
        ignore = data['IGNORED'][index]
        if not ignore:
            spectrum = data['SPECTRUM'][index].copy()
            cutoff = np.nanmax(spectrum) / scaleco
            finitespec = spectrum[np.isfinite(spectrum)]
            test = finitespec < cutoff
            if True in test:
                finitespec[test] = np.nan
                spectrum[np.isfinite(spectrum)] = finitespec
                pass
            wave, disp, shift, si, bck = wavesol(
                abs(spectrum),
                tt,
                wavett,
                disper,
                siv=siv,
                ovszspc=ovszspc,
                bck=None,
            )
            if (disp < ldisp) or (disp > udisp):
                data['TRIAL'][index] = 'Dispersion Out Of Bounds'
                ignore = True
                pass
            if (abs(disp - disper) < 1e-7) and not ovszspc:
                data['TRIAL'][index] = 'Dispersion Fit Failure'
                ignore = True
                pass
            pass
        if not ignore:
            liref = itp.interp1d(
                wavett * 1e-4, tt, bounds_error=False, fill_value=np.nan
            )
            phot2counts = liref(wave)
            data['PHT2CNT'][index] = phot2counts
            data['WAVE'][index] = wave  # MICRONS
            data['DISPERSION'][index] = disp  # ANGSTROMS/PIXEL
            data['SHIFT'][index] = shift * 1e4 / disp  # PIXELS
            data['BACKGROUND'][index] = bck
            data['SPECTRUM'][index] = data['SPECTRUM'][index] - bck
            pass
        else:
            data['WAVE'][index] = (data['SPECTRUM'][index]) * np.nan
        data['IGNORED'][index] = ignore
        pass
    allignore = data['IGNORED']
    allculprits = data['TRIAL']
    log.info(
        '>-- IGNORED: %s / %s', str(np.nansum(allignore)), str(len(allignore))
    )
    for index, ignore in enumerate(allignore):
        if ignore:
            log.info('>-- %s: %s', str(index), str(allculprits[index]))
        pass
    data.pop('EXP', None)
    data.pop('EXPFLAG', None)
    for k, v in data.items():
        out['data'][k] = v
    caled = not np.all(data['IGNORED'])
    if caled:
        out['STATUS'].append(True)
    return caled


# ------------------------- ------------------------------------------
# -- DETECTOR PLATE SCALE -- -----------------------------------------
def dps(flttype):
    '''
    G. ROUDIER: Detector plate scale

    http://www.stsci.edu/hst/wfc3/ins_performance/detectors
    http://www.stsci.edu/hst/stis/design/detectors
    http://www.stsci.edu/hst/stis/design/gratings
    '''
    detector = flttype.split('-')[2]
    fltr = flttype.split('-')[3]
    arcsec2pix = None
    if detector in ['IR']:
        arcsec2pix = 0.13
    if detector in ['UVIS']:
        arcsec2pix = 0.04
    if detector in ['CCD']:
        arcsec2pix = 0.05071
    if detector in ['FUV.MAMA']:
        if fltr in ['G140M']:
            arcsec2pix = np.sqrt(0.029 * 0.036)
        pass
    return arcsec2pix


# --------------------------------------------------------------------
# -- SCIENCE WAVELENGTH BAND -- --------------------------------------
def validrange(flttype):
    '''
    G. ROUDIER: Science wavelength band

    http://www.stsci.edu/hst/wfc3/documents/ISRs/WFC3-2009-18.pdf
    http://www.stsci.edu/hst/wfc3/documents/ISRs/WFC3-2009-17.pdf
    http://www.stsci.edu/hst/stis/design/gratings/documents/handbooks/currentIHB/c13_specref07.html
    http://www.stsci.edu/hst/stis/design/gratings/documents/handbooks/currentIHB/c13_specref16.html
    http://www.stsci.edu/hst/stis/design/gratings/documents/handbooks/currentIHB/c13_specref05.html
    '''
    fltr = flttype.split('-')[3]
    vrange = None
    if fltr in ['G141']:
        vrange = [1.12, 1.65]  # MICRONS
    if fltr in ['G102']:
        vrange = [0.80, 1.14]
    if fltr in ['G430L']:
        vrange = [0.30, 0.57]
    if fltr in ['G140M']:
        vrange = [0.12, 0.17]
    if fltr in ['G750L']:
        vrange = [0.53, 0.95]
    if fltr in ['3.6']:
        vrange = [3.1, 3.92]
    if fltr in ['4.5']:
        vrange = [3.95, 4.95]
    return vrange


# ----------------------------- --------------------------------------
# -- FILTERS AND GRISMS -- -------------------------------------------
def fng(flttype):
    '''
    G. ROUDIER: Filters and grisms

    http://www.stsci.edu/hst/wfc3/documents/handbooks/currentDHB/wfc3_dhb.pdf
    G102 http://www.stsci.edu/hst/wfc3/documents/ISRs/WFC3-2009-18.pdf
    G141 http://www.stsci.edu/hst/wfc3/documents/ISRs/WFC3-2009-17.pdf
    http://www.stsci.edu/hst/stis/design/gratings/documents/handbooks/currentIHB/c13_specref07.html
    http://www.stsci.edu/hst/stis/design/gratings/documents/handbooks/currentIHB/c13_specref16.html
    http://www.stsci.edu/hst/stis/design/gratings/documents/handbooks/currentIHB/c13_specref05.html
    '''
    fltr = flttype.split('-')[3]
    wvrng = None
    disp = None
    llim = None
    ulim = None
    if fltr == 'G141':
        wvrng = [1085e1, 17e3]  # Angstroms
        disp = 46.5  # Angstroms/Pixel
        llim = 45
        ulim = 47.7  # HD 189733 TRANSIT VS 47.5 STSCI
        pass
    if fltr == 'G102':
        wvrng = [8e3, 115e2]  # Angstroms
        disp = 24.5  # Angstroms/Pixel
        llim = 23.5
        ulim = 25
        pass
    if fltr == 'G430L':
        wvrng = [29e2, 57e2]  # Angstroms
        disp = 2.73  # Angstroms/Pixel
        llim = 2.6
        ulim = 2.9
        pass
    if fltr == 'G750L':
        wvrng = [524e1, 1027e1]  # Angstroms
        disp = 4.92  # Angstroms/Pixel
        llim = 4.91
        ulim = 4.93
        pass
    # if fltr == '3.6':
    #     wvrng = [3100e1, 3950e1 ]  # Angstroms
    #     disp = 0  # Angstroms/Pixel
    #     llim = 0
    #     ulim = 0
    #     pass
    # if fltr == '4.5':
    #     wvrng = [3900e1, 4950e1 ]  # Angstroms
    #     disp = 0  # Angstroms/Pixel
    #     llim = 0
    #     ulim = 0
    #     pass
    return wvrng, disp, llim, ulim


# ------------------------ -------------------------------------------
# -- ISOLATE -- ------------------------------------------------------
def isolate(
    thisdiff,
    psmin,
    spectrace,
    scanwdw,
    targetn,
    floodlevel,
    axis=1,
    stare=False,
):
    '''
    G. ROUDIER: Based on Minkowski functionnals decomposition algorithm
    '''
    valid = np.isfinite(thisdiff)
    if np.nansum(~valid) > 0:
        thisdiff[~valid] = 0
    select = thisdiff[valid] < psmin
    if np.nansum(select) > 0:
        thisdiff[valid][select] = 0
    flooded = thisdiff.copy()
    flooded[thisdiff < floodlevel] = 0
    eachprofile = np.nansum(flooded, axis=axis)
    eachpixels = np.arange(eachprofile.size)
    loc = eachpixels[eachprofile > 0]
    # ines mertz : changed the "if" statement from "if loc: " to "if loc.size:" since loc is now numpy.ndarray
    # if loc:
    if loc.size:
        diffloc = [0]
        diffloc.extend(np.diff(loc))
        minlocs = [np.nanmin(loc)]
        maxlocs = []
        cvcount = 0
        if axis > 0:
            thr = int(scanwdw - 6)
            if thr <= 0:
                thr = scanwdw
            thrw = int(scanwdw - 6)
            if thrw <= 0:
                thrw = scanwdw
            pass
        else:
            thr = 6
            thrw = int(spectrace - 12)
            pass
        for dl, c in zip(diffloc, np.arange(len(loc))):
            if (dl > thr) and (cvcount <= thrw):
                minlocs.pop(-1)
                minlocs.append(loc[c])
                pass
            if (dl > thr) and (cvcount > thrw):
                maxlocs.append(loc[c - 1])
                minlocs.append(loc[c])
                cvcount = 0
                pass
            if dl < thr:
                cvcount += 1
            pass
        maxlocs.append(loc[-1])
        mn = minlocs[targetn]
        mx = maxlocs[targetn]
        if (mx - mn) < thrw:
            mn = np.nan
            mx = np.nan
            pass
        if stare:
            mn = minlocs[targetn] - 1
            mx = maxlocs[targetn] + 1
            pass
        pass
    else:
        mn = np.nan
        mx = np.nan
        pass
    return mn, mx


# ------------- ------------------------------------------------------
# -- APERTURE AND FILTER TO TOTAL TRANSMISSION FILTER -- -------------
def ag2ttf(flttype):
    '''
    G ROUDIER: Aperture and filter to total transmission filter
    '''
    detector = flttype.split('-')[2]
    grism = flttype.split('-')[3]
    lightpath = ag2lp(detector, grism)
    mu, ttp = bttf(lightpath)
    if grism == 'G141':
        select = (mu > 1.055e4) & (mu < 1.75e4)
        mu = mu[select]
        ttp = ttp[select]
        pass
    if grism == 'G102':
        select = mu < 1.16e4
        mu = mu[select]
        ttp = ttp[select]
        pass
    if grism == 'G430L':
        select = (mu > 2.9e3) & (mu < 5.7394e3)
        mu = mu[select]
        ttp = ttp[select]
        pass
    if grism == 'G750L':
        select = (mu > 5.24e3) & (mu < 1.27e4)
        mu = mu[select]
        ttp = ttp[select]
        pass
    return mu, ttp


# ------------------------------------------------------ -------------
# -- APERTURE AND GRISM TO .FITS FILES -- ----------------------------
def ag2lp(detector, grism):
    '''
    G. ROUDIER: The first element of the returned list defines the default
    interpolation grid, filter/grism file is suited.

    ['Grism source', 'Refractive correction plate', 'Cold mask', 'Mirror 1', 'Mirror 2',
    'Fold mirror', 'Channel select mechanism', 'Pick off mirror', 'OTA']

    http://www.stsci.edu/hst/wfc3/documents/ISRs/WFC3-2011-05.pdf
    ftp://ftp.stsci.edu/cdbs/comp/ota/
    ftp://ftp.stsci.edu/cdbs/comp/wfc3/
    '''
    lightpath = []
    if grism == 'G141':
        lightpath.append('WFC3/wfc3_ir_g141_src_004_syn.fits')
    if grism == 'G102':
        lightpath.append('WFC3/wfc3_ir_g102_src_003_syn.fits')
    if grism == 'G430L':
        lightpath.append('STIS/stis_g430l_009_syn.fits')
    if grism == 'G750L':
        lightpath.append('STIS/stis_g750l_009_syn.fits')
    if detector == 'IR':
        lightpath.extend(
            [
                'WFC3/wfc3_ir_rcp_001_syn.fits',
                'WFC3/wfc3_ir_mask_001_syn.fits',
                'WFC3/wfc3_ir_mir1_001_syn.fits',
                'WFC3/wfc3_ir_mir2_001_syn.fits',
                'WFC3/wfc3_ir_fold_001_syn.fits',
                'WFC3/wfc3_ir_csm_001_syn.fits',
                'WFC3/wfc3_pom_001_syn.fits',
                'WFC3/hst_ota_007_syn.fits',
            ]
        )
        pass
    return lightpath


# --------------------------------------- ----------------------------
# -- BUILD TOTAL TRANSMISSION FILTER -- ------------------------------
def bttf(lightpath):
    '''
    G. ROUDIER: Builds total transmission filter
    '''
    ttp = 1e0
    muref = [np.nan]
    for name in lightpath:
        muref, t = loadcalf(name, muref)
        ttp *= t
    return muref, ttp


# ------------------------------------- ------------------------------
# -- WFC3 CAL FITS -- ------------------------------------------------
def loadcalf(name, muref, calloc=excalibur.context['data_cal']):
    '''
    G. ROUDIER: Loads optical element .fits calibration file
    '''
    # NOTEBOOK RUN
    # calloc = '/proj/sdp/data/cal'
    # ------------
    fitsfile = os.path.join(calloc, name)
    data = pyfits.getdata(fitsfile)
    muin = np.array(data.WAVELENGTH)
    tin = np.array(data.THROUGHPUT)
    if False in np.isfinite(muref):
        muref = muin
    f = itp.interp1d(muin, tin, bounds_error=False, fill_value=0)
    t = f(muref)
    return muref, t


# ------------------- ------------------------------------------------
# -- WAVELENGTH SOLUTION -- ------------------------------------------
def wavesol(
    spectrum,
    tt,
    wavett,
    disper,
    siv=None,
    fd=False,
    bck=None,
    fs=False,
    ovszspc=False,
):
    '''
    G. ROUDIER: Wavelength calibration on log10 spectrum to emphasize the
    edges, approximating the log(stellar spectrum) with a linear model
    '''
    mutt = wavett.copy()
    mutt /= 1e4
    xdata = np.arange(spectrum.size)
    test = tt == 0
    if True in test:
        tt[test] = np.nan
    logtt = np.log10(tt)
    test = spectrum == 0
    if True in test:
        spectrum[test] = np.nan
    logspec = np.log10(spectrum)
    wave = xdata * disper / 1e4
    select = np.isfinite(spectrum)
    minwave = np.nanmin(wave[select])
    select = np.isfinite(tt)
    reftt = np.nanmin(mutt[select])
    shift = reftt - minwave
    scale = np.nanmedian(logspec) - np.nanmedian(logtt)
    params = lm.Parameters()
    params.add('scale', value=scale)
    if bck is None:
        params.add('background', value=0e0, vary=False)
    else:
        params.add('background', value=bck)
    if siv is None:
        params.add('slope', value=1e-2)
    else:
        params.add('slope', value=siv, vary=False)
    if fd or ovszspc:
        params.add('disper', value=disper, vary=False)
    else:
        params.add('disper', value=disper)
    if fs:
        params.add('shift', value=shift, vary=False)
    else:
        params.add('shift', value=shift)
    out = lm.minimize(wcme, params, args=(logspec, mutt, logtt, False))
    disper = out.params['disper'].value
    shift = out.params['shift'].value
    slope = out.params['slope'].value
    scale = out.params['scale'].value
    background = out.params['background'].value
    wave = wcme(out.params, logspec, refmu=mutt, reftt=logtt)
    return wave, disper, shift + minwave, slope, background


# ------------------------- ------------------------------------------
# -- WAVELENGTH FIT FUNCTION -- --------------------------------------
def wcme(params, data, refmu=None, reftt=None, forward=True):
    '''
    G. ROUDIER: Wavelength calibration function for LMfit
    '''
    slope = params['slope'].value
    scale = params['scale'].value
    disper = params['disper'].value
    shift = params['shift'].value
    background = params['background'].value
    liref = itp.interp1d(refmu, reftt, bounds_error=False, fill_value=np.nan)
    wave = np.arange(data.size) * disper * 1e-4 + shift
    model = liref(wave) + scale + slope * wave
    select = (np.isfinite(model)) & (np.isfinite(data))
    d = np.log10(10 ** (data.copy()) - background)
    weights = np.ones(d.size)
    if np.sum(~select) > 0:
        model[~select] = 1e0
        d[~select] = 0e0
        weights[~select] = 1e2
        pass
    if forward:
        out = wave
    else:
        out = (d - model) / weights
    return out


# ----------------------------- --------------------------------------
# -- TIME TO Z -- ----------------------------------------------------
def time2z(
    time,
    ipct,
    tknot,
    sma,
    orbperiod,
    ecc,
    tperi=None,
    epsilon=1e-10,
):
    '''
    G. ROUDIER: Time samples in [Days] to separation in [R*]
    '''
    return excalibur.util.time.time2z(
        time,
        ipct,
        tknot,
        sma,
        orbperiod,
        ecc,
        tperi,
        epsilon,
        True,
    )


# --------------- ----------------------------------------------------
# -- CALIBRATE STARE DATA -- -----------------------------------------
def starecal(
    _fin,
    clc,
    tim,
    tid,
    flttype,
    out,
    emptythr=1e3,
    frame2png=False,
):
    '''
    G. ROUDIER: WFC3 STARE Calibration
    '''
    calibrated = False
    # VISIT ----------------------------------------------------------
    for pkey in tim['data'].keys():
        visits = np.array(tim['data'][pkey]['visits'])
    # DATA TYPE ------------------------------------------------------
    vrange = validrange(flttype)
    wvrng, disper, ldisp, udisp = fng(flttype)
    spectrace = np.round((np.max(wvrng) - np.min(wvrng)) / disper)
    # LOAD DATA ------------------------------------------------------
    dbs = os.path.join(dawgie.context.data_dbs, 'mast')
    data = {
        'LOC': [],
        'EPS': [],
        'DISPLIM': [ldisp, udisp],
        'SCANRATE': [],
        'SCANLENGTH': [],
        'SCANANGLE': [],
        'EXP': [],
        'EXPERR': [],
        'EXPFLAG': [],
        'VRANGE': vrange,
        'TIME': [],
        'EXPLEN': [],
        'MIN': [],
        'MAX': [],
        'TRIAL': [],
    }
    for loc in sorted(clc['LOC']):
        fullloc = os.path.join(dbs, loc)
        with pyfits.open(fullloc) as hdulist:
            header0 = hdulist[0].header
            eps = False
            test = header0['UNITCORR']
            if test in ['COMPLETE', 'PERFORM']:
                eps = True
            data['EPS'].append(eps)
            if 'SCAN_RAT' in header0:
                data['SCANRATE'].append(header0['SCAN_RAT'])
            else:
                data['SCANRATE'].append(np.nan)
            if 'SCAN_LEN' in header0:
                data['SCANLENGTH'].append(header0['SCAN_LEN'])
            else:
                data['SCANLENGTH'].append(np.nan)
            if 'SCAN_ANG' in header0:
                data['SCANANGLE'].append(header0['SCAN_ANG'])
            elif 'PA_V3' in header0:
                data['SCANANGLE'].append(header0['PA_V3'])
            else:
                data['SCANANGLE'].append(666)
            frame = []
            errframe = []
            dqframe = []
            ftime = []
            fmin = []
            fmax = []
            datahdu = [h for h in hdulist if h.size != 0]
            for fits in datahdu:
                if fits.header['EXTNAME'] in ['SCI']:
                    fitsdata = np.empty(fits.data.shape)
                    fitsdata[:] = fits.data[:]
                    if eps and ('DELTATIM' in fits.header):
                        frame.append(fitsdata * float(fits.header['DELTATIM']))
                        fmin.append(
                            float(fits.header['GOODMIN'])
                            * float(fits.header['DELTATIM'])
                        )
                        fmax.append(
                            float(fits.header['GOODMAX'])
                            * float(fits.header['DELTATIM'])
                        )
                        pass
                    elif eps and 'DELTATIM' not in fits.header:
                        frame.append(fitsdata * np.nan)
                        fmin.append(np.nan)
                        fmax.append(np.nan)
                        pass
                    else:
                        frame.append(fitsdata)
                        fmin.append(float(fits.header['GOODMIN']))
                        fmax.append(float(fits.header['GOODMAX']))
                        pass
                    ftime.append(float(fits.header['ROUTTIME']))
                    del fits.data
                    pass
                if fits.header['EXTNAME'] in ['ERR', 'DQ']:
                    fitsdata = np.empty(fits.data.shape)
                    fitsdata[:] = fits.data[:]
                    if fits.header['EXTNAME'] == 'ERR':
                        errframe.append(fitsdata)
                    if fits.header['EXTNAME'] == 'DQ':
                        dqframe.append(fitsdata)
                    del fits.data
                    pass
                pass
            data['LOC'].append(loc)
            data['EXP'].append(frame)
            data['EXPERR'].append(errframe)
            data['EXPFLAG'].append(dqframe)
            data['TIME'].append(ftime)
            data['EXPLEN'].append(header0['EXPTIME'])
            data['MIN'].append(fmin)
            data['MAX'].append(fmax)
            pass
        pass
    # MASK DATA ------------------------------------------------------
    data['MEXP'] = data['EXP'].copy()
    data['MASK'] = data['EXPFLAG'].copy()
    data['IGNORED'] = [False] * len(data['LOC'])
    data['TRUNCSPEC'] = [False] * len(data['LOC'])
    data['FLOODLVL'] = [np.nan] * len(data['LOC'])
    data['TRIAL'] = [''] * len(data['LOC'])
    # FLOOD LEVEL STABILIZATION --------------------------------------
    for index in enumerate(data['LOC']):
        ignore = data['IGNORED'][index[0]]
        # ISOLATE SCAN Y ---------------------------------------------
        sampramp = np.array(data['MEXP'][index[0]]).copy()
        for sutr in sampramp:
            select = ~np.isfinite(sutr)
            if True in select:
                sutr[select] = 0
            pass
        psdiff = sampramp[0] - sampramp[-1]
        psmin = np.nanmin(data['MIN'][index[0]])
        fldthr = np.nanpercentile(
            psdiff,
            1e2 * (1e0 - spectrace / (psdiff.shape[0] * psdiff.shape[1])),
        )
        data['FLOODLVL'][index[0]] = fldthr
        pass
    allfloodlvl = np.array(data['FLOODLVL'])
    for v in set(visits):
        select = visits == v
        allfloodlvl[select] = np.nanmedian(allfloodlvl[select])
        pass
    data['FLOODLVL'] = list(allfloodlvl)
    # DATA CUBE ------------------------------------------------------
    for index, nm in enumerate(data['LOC']):
        ignore = data['IGNORED'][index]
        # ISOLATE SCAN Y ---------------------------------------------
        sampramp = np.array(data['MEXP'][index]).copy()
        for sutr in sampramp:
            select = ~np.isfinite(sutr)
            if True in select:
                sutr[select] = 0
            pass
        psdiff = sampramp[0] - sampramp[-1]
        psmin = np.nanmin(data['MIN'][index])
        targetn = 0
        scanwpi = 1
        minlocs = []
        maxlocs = []
        fldthr = data['FLOODLVL'][index]
        for de in [psdiff]:
            lmn, lmx = isolate(
                de,
                psmin,
                spectrace,
                scanwpi,
                targetn,
                fldthr,
                stare=True,
            )
            minlocs.append(lmn)
            maxlocs.append(lmx)
            pass
        ignore = ignore or not (
            (np.any(np.isfinite(minlocs))) and (np.any(np.isfinite(maxlocs)))
        )
        if not ignore:
            minl = np.nanmin(minlocs)
            maxl = np.nanmax(maxlocs)
            minl -= 12
            maxl += 12
            if minl < 0:
                minl = 0
            else:
                psdiff[: int(minl), :] = np.nan
            if maxl > (psdiff[0].shape[0] - 1):
                maxl = psdiff[0].shape[0] - 1
            else:
                psdiff[int(maxl) :, :] = np.nan
            thispstamp = psdiff.copy()
            thispstamp[thispstamp <= psmin] = np.nan
            # ISOLATE SCAN X -----------------------------------------
            targetn = 0
            minx, maxx = isolate(
                thispstamp.copy(),
                psmin,
                spectrace,
                scanwpi,
                targetn,
                fldthr,
                axis=0,
                stare=True,
            )
            if np.isfinite(minx * maxx):
                minx -= 1.5 * 12
                maxx += 1.5 * 12
                if minx < 0:
                    minx = 5
                    data['TRUNCSPEC'][index] = True
                    pass
                if maxx > (thispstamp.shape[1] - 1):
                    maxx = thispstamp.shape[1] - 5
                    data['TRUNCSPEC'][index] = True
                    pass
                thispstamp[:, : int(minx)] = np.nan
                thispstamp[:, int(maxx) :] = np.nan
                pstamperr = np.sqrt(abs(thispstamp))
                pass
            else:
                garbage = data['EXP'][index][::-1].copy()
                thispstamp = np.sum(np.diff(garbage, axis=0), axis=0)
                pstamperr = thispstamp * np.nan
                data['TRIAL'][index] = 'Could Not Find X Edges'
                ignore = True
                pass
            pass
        else:
            garbage = data['EXP'][index][::-1].copy()
            thispstamp = np.sum(np.diff(garbage, axis=0), axis=0)
            pstamperr = thispstamp * np.nan
            data['TRIAL'][index] = 'Could Not Find Y Edges'
            ignore = True
            pass
        data['MEXP'][index] = thispstamp
        data['TIME'][index] = np.nanmax(data['TIME'][index].copy())
        data['IGNORED'][index] = ignore
        data['EXPERR'][index] = pstamperr
        log.info('>-- %s / %s', str(index), str(len(data['LOC']) - 1))
        # PLOTS ------------------------------------------------------
        if frame2png:
            if not os.path.exists('TEST'):
                os.mkdir('TEST')
            if not os.path.exists('TEST/' + tid):
                os.mkdir('TEST/' + tid)
            plt.figure()
            plt.title('Ignored = ' + str(ignore))
            plt.imshow(thispstamp)
            plt.colorbar()
            plt.savefig('TEST/' + tid + '/' + nm + '.png')
            plt.close()
            pass
        pass
    # SPECTRUM EXTRACTION --------------------------------------------
    data['SPECTRUM'] = [np.nan] * len(data['LOC'])
    data['SPECERR'] = [np.nan] * len(data['LOC'])
    data['NSPEC'] = [np.nan] * len(data['LOC'])
    for index, loc in enumerate(data['LOC']):
        floodlevel = data['FLOODLVL'][index]
        if floodlevel < emptythr:
            data['IGNORED'][index] = True
            data['TRIAL'][index] = 'Empty Frame'
            pass
        ignore = data['IGNORED'][index]
        if not ignore:
            frame = data['MEXP'][index].copy()
            spectrum = []
            specerr = []
            nspectrum = []
            for col in np.array(frame).T:
                valid = np.isfinite(col)
                if True in valid:
                    spectrum.append(np.nansum(col[valid]))
                    specerr.append(np.sqrt(abs(np.nansum(col[valid]))))
                    nspectrum.append(np.nansum(valid))
                    pass
                else:
                    spectrum.append(np.nan)
                    specerr.append(np.nan)
                    nspectrum.append(0)
                    pass
                pass
            spectrum = np.array(spectrum)
            if np.all(~np.isfinite(spectrum)):
                data['IGNORED'][index] = True
                data['TRIAL'][index] = 'NaN Spectrum'
                pass
            data['SPECTRUM'][index] = np.array(spectrum)
            data['SPECERR'][index] = np.array(specerr)
            data['NSPEC'][index] = np.array(nspectrum)
            pass
        pass
    # WAVELENGTH CALIBRATION -----------------------------------------
    wavett, tt = ag2ttf(flttype)
    scaleco = np.nanmax(tt) / np.nanmin(tt[tt > 0])
    data['PHT2CNT'] = [np.nan] * len(data['LOC'])
    data['WAVE'] = [np.nan] * len(data['LOC'])
    data['DISPERSION'] = [np.nan] * len(data['LOC'])
    data['SHIFT'] = [np.nan] * len(data['LOC'])
    spectralindex = []
    for index, loc in enumerate(data['LOC']):
        ignore = data['IGNORED'][index]
        ovszspc = False
        if data['TRUNCSPEC'][index]:
            ovszspc = True
            if 'G141' in flttype:
                wthr = 1.67 * 1e4
            select = wavett < wthr
            wavett = wavett[select]
            tt = tt[select]
            pass
        if not ignore:
            spectrum = data['SPECTRUM'][index].copy()
            cutoff = np.nanmax(spectrum) / scaleco
            spectrum[spectrum < cutoff] = np.nan
            spectrum = abs(spectrum)
            _w, d, _s, si, _bck = wavesol(
                spectrum, tt, wavett, disper, ovszspc=ovszspc
            )
            if ldisp < d < udisp:
                spectralindex.append(si)
            pass
        pass
    siv = np.nanmedian(spectralindex)
    for index, loc in enumerate(data['LOC']):
        ignore = data['IGNORED'][index]
        ovszspc = False
        if data['TRUNCSPEC'][index]:
            ovszspc = True
            if 'G141' in flttype:
                wthr = 1.67 * 1e4
            select = wavett < wthr
            wavett = wavett[select]
            tt = tt[select]
            pass
        if not ignore:
            spectrum = data['SPECTRUM'][index].copy()
            cutoff = np.nanmax(spectrum) / scaleco
            spectrum[spectrum < cutoff] = np.nan
            spectrum = abs(spectrum)
            wave, disp, shift, si, _bck = wavesol(
                spectrum, tt, wavett, disper, siv=siv, ovszspc=ovszspc
            )
            liref = itp.interp1d(
                wavett * 1e-4, tt, bounds_error=False, fill_value=np.nan
            )
            phot2counts = liref(wave)
            data['PHT2CNT'][index] = phot2counts
            data['WAVE'][index] = wave  # MICRONS
            data['DISPERSION'][index] = disp  # ANGSTROMS/PIXEL
            data['SHIFT'][index] = shift * 1e4 / disp  # PIXELS
            pass
        else:
            data['WAVE'][index] = (data['SPECTRUM'][index]) * np.nan
        data['IGNORED'][index] = ignore
        pass
    allignore = data['IGNORED']
    allculprits = data['TRIAL']
    allindex = np.arange(len(data['LOC']))
    log.info(
        '>-- IGNORED: %s / %s', str(np.nansum(allignore)), str(len(allignore))
    )
    for index in allindex:
        if allculprits[index]:
            log.info('>-- %s: %s', str(index), str(allculprits[index]))
            pass
        pass
    data.pop('EXP', None)
    data.pop('EXPFLAG', None)
    for k, v in data.items():
        out['data'][k] = v
    calibrated = not np.all(data['IGNORED'])
    if calibrated:
        out['STATUS'].append(True)
    return calibrated


# -------------------------- -----------------------------------------
# -- STIS CALIBRATION -- ---------------------------------------------
def stiscal_G750L(_fin, clc, tim, tid, flttype, out):
    '''
    R. ESTRELA: STIS .flt data extraction and wavelength calibration
    '''
    calibrated = False
    # VISIT NUMBERING --------------------------------------------------------------------
    for pkey in tim['data'].keys():
        visits = np.array(tim['data'][pkey]['dvisits'])
    # PHASE ------------------------------------------------------------------------------
    # for pkey in tim['data'].keys(): phase = np.array(tim['data'][pkey]['phase'])
    # OPTICS AND FILTER ------------------------------------------------------------------
    vrange = validrange(flttype)
    _wvrng, disp, ldisp, udisp = fng(flttype)
    # DATA FORMAT ------------------------------------------------------------------------
    dbs = os.path.join(dawgie.context.data_dbs, 'mast')
    data = {
        'LOC': [],
        'EPS': [],
        'DISPLIM': [ldisp, udisp],
        'SCANRATE': [],
        'SCANLENGTH': [],
        'SCANANGLE': [],
        'EXP': [],
        'EXPERR': [],
        'EXPFLAG': [],
        'VRANGE': vrange,
        'TIME': [],
        'EXPLEN': [],
        'MIN': [],
        'MAX': [],
        'TRIAL': [],
        'TIMEOBS': [],
        'DATEOBS': [],
    }
    # LOAD DATA --------------------------------------------------------------------------
    for loc in sorted(clc['LOC']):
        fullloc = os.path.join(dbs, loc)
        with pyfits.open(fullloc) as hdulist:
            header0 = hdulist[0].header
            if 'SCAN_ANG' in header0:
                scanangle = header0['SCAN_ANG']
            elif 'PA_V3' in header0:
                scanangle = header0['PA_V3']
            else:
                scanangle = 666
            allloc = []
            alltime = []
            allexplen = []
            alleps = []
            allexp = []
            allexperr = []
            allmask = []
            allmin = []
            allmax = []
            alldate = []
            alltimeobs = []
            for fits in hdulist:
                if (fits.size != 0) and (fits.header['EXTNAME'] == 'SCI'):
                    allloc.append(fits.header['EXPNAME'])
                    alltime.append(float(fits.header['EXPEND']))
                    allexplen.append(float(fits.header['EXPTIME']))
                    alldate.append(header0['TDATEOBS'])
                    alltimeobs.append(header0['TTIMEOBS'])
                    fitsdata = np.empty(fits.data.shape)
                    fitsdata[:] = fits.data[:]
                    test = fits.header['BUNIT']
                    eps = False
                    if test != 'COUNTS':
                        eps = True
                    alleps.append(eps)
                    allmin.append(float(fits.header['GOODMIN']))
                    allmax.append(float(fits.header['GOODMAX']))
                    # BINARIES
                    # GMR: Let's put that in the mask someday
                    nam = 'MIDPOINT'
                    if tid in ['HAT-P-1'] and nam in header0['TARGNAME']:
                        allexp.append(fitsdata[200:380, :])
                    else:
                        allexp.append(fitsdata)
                    del fits.data
                    pass
                if 'EXTNAME' in fits.header:
                    if fits.header['EXTNAME'] in ['ERR', 'DQ']:
                        fitsdata = np.empty(fits.data.shape)
                        fitsdata[:] = fits.data[:]
                        if fits.header['EXTNAME'] == 'ERR':
                            allexperr.append(fitsdata)
                        if fits.header['EXTNAME'] == 'DQ':
                            allmask.append(fitsdata)
                        del fits.data
                        pass
                    if eps:
                        eps2count = allexplen[-1] * float(header0['CCDGAIN'])
                        allexp[-1] = allexp[-1] * eps2count
                        allexperr[-1] = allexperr[-1] * eps2count
                        pass
                    pass
                pass
            allscanangle = [scanangle] * len(allloc)
            allscanlength = [1e0] * len(allloc)
            allscanrate = [0e0] * len(allloc)
            data['LOC'].extend(allloc)
            data['EPS'].extend(alleps)
            data['SCANRATE'].extend(allscanrate)
            data['SCANLENGTH'].extend(allscanlength)
            data['SCANANGLE'].extend(allscanangle)
            data['EXP'].extend(allexp)
            data['EXPERR'].extend(allexperr)
            data['TIME'].extend(alltime)
            data['EXPLEN'].extend(allexplen)
            data['MIN'].extend(allmin)
            data['MAX'].extend(allmax)
            data['TIMEOBS'].extend(alltimeobs)
            data['DATEOBS'].extend(alldate)
            pass
        pass
    data['MEXP'] = data['EXP'].copy()
    data['MASK'] = data['EXPFLAG'].copy()
    data['ALLDATEOBS'] = data['DATEOBS'].copy()
    data['ALLTIMEOBS'] = data['TIMEOBS'].copy()
    data['IGNORED'] = np.array([False] * len(data['LOC']))
    data['FLOODLVL'] = [np.nan] * len(data['LOC'])
    data['TRIAL'] = [''] * len(data['LOC'])
    data['SPECTRUM'] = [np.array([np.nan])] * len(data['LOC'])
    data['SPECERR'] = [np.array([np.nan])] * len(data['LOC'])
    data['TEMPLATE'] = [np.array([np.nan])] * len(data['LOC'])
    data['NSPEC'] = [1e0] * len(data['LOC'])
    # REJECT OUTLIERS IN EXPOSURE LENGTH -------------------------------------------------
    for v in set(visits):
        select = visits == v
        visitexplength = np.array(data['EXPLEN'])[select]
        visitignore = data['IGNORED'][select]
        ref = np.nanmedian(visitexplength)
        visitignore[visitexplength != ref] = True
        data['IGNORED'][select] = visitignore
        pass
    for index, ignore in enumerate(data['IGNORED']):
        # SELECT DATE AND TIME OF THE EXPOSURE FOR FLAT FRINGE SELECTION
        frame = data['MEXP'][index].copy()
        dateobs_exp = data['ALLDATEOBS'][index]
        timeobs_exp = data['ALLTIMEOBS'][index]
        tog_exp = dateobs_exp + ' ' + timeobs_exp
        time_exp = raissatime.mktime(
            datetime.datetime.strptime(tog_exp, "%Y-%m-%d %H:%M:%S").timetuple()
        )
        # LOAD FRINGE FLAT ---------------------------------------------------------------
        obs_name = clc['ROOTNAME'][0]
        name_sel = obs_name[:-5]
        lightpath_fringe = 'STIS/CCDFLAT/'
        calloc = excalibur.context['data_cal']
        filefringe = os.path.join(calloc, lightpath_fringe)
        if tid in ['HD 209458']:
            lightpath_fringe = 'STIS/CCDFLAT/h230851ao_pfl.fits'
            calloc = excalibur.context['data_cal']
            filefringe = os.path.join(calloc, lightpath_fringe)
            hdu = pyfits.open(filefringe)
            data_fringe = hdu[1].data
            err_fringe = hdu[2].data
            pass
        else:
            diff_list = []
            all_infile = []
            for infile in glob.glob(f'{filefringe}/{name_sel}*_flt.fits'):
                # print(' infile',infile)
                hdu = pyfits.open(infile)
                all_infile.append(infile)
                header_flat = hdu[0].header
                date_time = header_flat['TDATEOBS']
                hour_time = header_flat['TTIMEOBS']
                tog = date_time + ' ' + hour_time
                time_flat_s = raissatime.mktime(
                    datetime.datetime.strptime(
                        tog, "%Y-%m-%d %H:%M:%S"
                    ).timetuple()
                )
                diff = abs(time_exp - time_flat_s)
                diff_list.append(diff)
                pass
            # 6/11/24 GB
            #  There is a bug here where AU Mic and WASP-127 crash.
            #  The diff_list is empty, so the np.min(diff_list) crashes.
            cond_win = np.where(diff_list == np.min(diff_list))
            all_infile = np.array(all_infile)
            sel_flatfile = all_infile[cond_win][0]
            hdulist = pyfits.open(sel_flatfile)
            data_fringe = hdulist[4].data
            err_fringe = hdulist[5].data
            pass
        smooth_fringe = scipy.signal.medfilt(data_fringe, 7)
        sigma_fringe = np.median(err_fringe)
        bad_fringe = (np.abs(data_fringe - smooth_fringe) / sigma_fringe) > 2
        img_fringe = data_fringe.copy()
        img_fringe[bad_fringe] = smooth_fringe[bad_fringe]

        cont_data = img_fringe.copy()
        div_list = []
        for i in range(508, 515):
            pixels = np.arange(0, 1024, 1)
            coefs = poly.polyfit(pixels, cont_data[i, :], 11)
            ffit = poly.polyval(pixels, coefs)
            div = cont_data[i, :] / ffit
            div_list.append(div)
        #############
        # COSMIC RAY REJECTION IN THE 2D IMAGE
        img_cr = frame.copy()
        allframe_list = []
        for i in range(0, len(frame)):
            img_sm = scipy.signal.medfilt(img_cr[i, :], 9)
            # std = np.std(img_cr[i,:] - img_sm)
            std = np.std(img_sm)
            bad = np.abs(img_cr[i, :] - img_sm) > 3 * std
            line = img_cr[i, :]
            line[bad] = img_sm[bad]
            allframe_list.append(line)
            pass
        allframe = np.array(allframe_list)
        # APPLY FLAT FRINGE
        if not ignore:
            find_spec = np.where(allframe == np.max(allframe))
            spec_idx = find_spec[0][0]
            spec_idx_up = spec_idx + 4
            spec_idx_dwn = spec_idx - 3
            spec_idx_all = np.arange(spec_idx_dwn, spec_idx_up, 1)
            frame2 = allframe.copy()
            for i, flatnorm in zip(spec_idx_all, div_list):
                # frame_sel = allframe[i, :]
                # coefs_f = poly.polyfit(pixels, frame_sel, 12)
                # ffit_f = poly.polyval(pixels, coefs_f)
                frame2[i, 400:1023] = frame2[i, 400:1023] / flatnorm[400:1023]

                data['SPECTRUM'][index] = np.nansum(frame2, axis=0)
                data['SPECERR'][index] = np.sqrt(np.nansum(frame2, axis=0))
                pass
        else:
            data['SPECTRUM'][index] = np.nansum(frame, axis=0) * np.nan
            data['SPECERR'][index] = np.nansum(frame, axis=0) * np.nan
            data['TRIAL'][index] = 'Exposure Length Outlier'
            pass
        pass
    ####
    # MASK BAD PIXELS IN SPECTRUM --------------------------------------------------------
    for v in set(visits):
        select = (visits == v) & ~(data['IGNORED'])
        # for index, valid in enumerate(select):
        #    spec_ind = data['SPECTRUM'][index]
        #    if valid and 'G750' in flttype:
        #            b, a = signal.butter(3., 0.05)
        #            zi = signal.lfilter_zi(b, a)
        #            spec0 = spec_ind[500:1024]
        #            spec_ind[500:1024], _ = signal.lfilter(b, a, spec_ind[500:1024], zi=zi*spec0[0])
        specarray = np.array(
            [s for s, ok in zip(data['SPECTRUM'], select) if ok]
        )
        trans = np.transpose(specarray)
        template = np.nanmedian(trans, axis=1)
        # TEMPLATE MEDIAN 5 POINTS LOW PASS FILTER ---------------------------------------
        smootht = []
        smootht.extend([template[0]] * 2)
        smootht.extend(template)
        smootht.extend([template[-1]] * 2)
        for index in np.arange(len(template)):
            medianvalue = np.nanmedian(template[index : index + 5])
            smootht[2 + index] = medianvalue
            pass
        smootht = smootht[2:-2]

        template = np.array(smootht)
        for vindex, valid in enumerate(select):
            if valid:
                data['TEMPLATE'][vindex] = template
            pass
        pass
    wavett, tt = ag2ttf(flttype)
    # COSMIC RAYS REJECTION --------------------------------------------------------------
    data['PHT2CNT'] = [np.nan] * len(data['LOC'])
    data['WAVE'] = [np.array([np.nan])] * len(data['LOC'])
    data['DISPERSION'] = [np.nan] * len(data['LOC'])
    data['SHIFT'] = [np.nan] * len(data['LOC'])
    allsi = []
    for index, rejected in enumerate(data['IGNORED']):
        if not rejected:
            template = data['TEMPLATE'][index]
            spec = data['SPECTRUM'][index]
            temp_spec = spec / template
            ht25 = np.nanpercentile(temp_spec, 25)
            lt75 = np.nanpercentile(temp_spec, 75)
            std = np.std(temp_spec[(temp_spec > ht25) & (temp_spec < lt75)])
            # BAD PIXEL THRESHOLD --------------------------------------------------------
            bpthr = temp_spec > np.nanmedian(temp_spec) + 3e0 * std
            if True in bpthr:
                spec[bpthr] = np.nan
            # if 'G430' in flttype:
            #    cond_max = np.where(spec == np.nanmax(spec))
            #    spec[cond_max] = np.nan
            #    data['SPECTRUM'][index] = spec
            # else:
            #    data['SPECTRUM'][index] = spec

            # FIRST WAVESOL --------------------------------------------------------------
            # scaleco = np.nanmax(tt) / np.nanmin(tt[tt > 0])
            scaleco = 1e1
            if (
                np.sum(np.isfinite(spec)) > (spec.size / 2)
                and 'G750' in flttype
            ):
                wavecalspec = spec.copy()
                finitespec = np.isfinite(spec)
                nanme = spec[finitespec] < (np.nanmax(wavecalspec) / scaleco)
                if True in nanme:
                    finwavecalspec = wavecalspec[finitespec]
                    finwavecalspec[nanme] = np.nan
                    wavecalspec[finitespec] = finwavecalspec
                    pass
                selref = tt > (np.nanmax(tt) / scaleco)
                # dispersion is fixed, shift is not
                w, d, s, si, _bck = wavesol(
                    abs(wavecalspec),
                    tt[selref],
                    wavett[selref],
                    disp,
                    fd=True,
                    fs=False,
                )
                data['WAVE'][index] = w
                allsi.append(si)
                pass
            else:
                data['IGNORED'][index] = True
                data['TRIAL'][
                    index
                ] = 'Not Enough Valid Points In Extracted Spectrum'
                pass
            if 'G430' in flttype:
                pixel = np.arange(len(spec))
                w = (pixel + 1079.96) / 0.37
                data['WAVE'][index] = w * 1e-4
                liref = itp.interp1d(
                    wavett * 1e-4, tt, bounds_error=False, fill_value=np.nan
                )
                phot2counts = liref(w * 1e-4)
                data['PHT2CNT'][index] = phot2counts
                data['DISPERSION'][index] = 2.70
                data['SHIFT'][index] = -1079.96
                pass
            pass
        pass
    # SECOND Wavesol
    # SECOND SIGMA CLIPPING
    for index, rejected in enumerate(data['IGNORED']):
        if not rejected:
            # template = data['TEMPLATE'][index]
            # spec = data['SPECTRUM'][index]
            # temp_spec = spec/template
            # ht25 = np.nanpercentile(temp_spec,25)
            # lt75 = np.nanpercentile(temp_spec,75)
            # selfin = np.isfinite(temp_spec)
            # std1 = np.nanstd(temp_spec[selfin][(temp_spec[selfin] > ht25) & (temp_spec[selfin] < lt75)])
            # BAD PIXEL THRESHOLD --------------------------------------------------------
            # bpthr = temp_spec[selfin] > np.nanmedian(temp_spec) + 2e0*std1
            # temp_cut=template.copy()
            # spec_cut = spec.copy()
            # data['SPECTRUM'][index] = spec_cut

            # SECOND WAVESOL --------------------------------------------------------------
            # scaleco = np.nanmax(tt) / np.nanmin(tt[tt > 0])
            scaleco = 1e1

            if (
                np.sum(np.isfinite(spec)) > (spec.size / 2)
                and 'G750' in flttype
            ):
                wavecalspec = spec.copy()
                selfinspec = np.isfinite(spec)
                nanme = spec[selfinspec] < (np.nanmax(wavecalspec) / scaleco)
                if True in nanme:
                    wavefin = wavecalspec[selfinspec]
                    wavefin[nanme] = np.nan
                    wavecalspec[selfinspec] = wavefin
                    pass
                selref = tt > (np.nanmax(tt) / scaleco)
                # dispersion is fixed, shift is not
                w, d, s, si, _bck = wavesol(
                    abs(wavecalspec),
                    tt[selref],
                    wavett[selref],
                    disp,
                    siv=np.nanmedian(allsi),
                    fd=True,
                    fs=False,
                )
                liref = itp.interp1d(
                    wavett * 1e-4, tt, bounds_error=False, fill_value=np.nan
                )
                phot2counts = liref(w)
                data['PHT2CNT'][index] = phot2counts
                data['WAVE'][index] = w
                data['DISPERSION'][index] = d
                data['SHIFT'][index] = s
                pass
                # else:
                # data['IGNORED'][index] = True
                # data['TRIAL'][index] = 'Not Enough Valid Points In Extracted Spectrum'
                pass

            # WAVELENGTH CALIBRATION - G430L
        #             if 'G430' in flttype:
        #                 pixel = np.arange(len(spec))
        #                 w = (pixel + 1079.96)/0.37
        #                 data['WAVE'][index] = w*1e-4
        #                 liref = itp.interp1d(wavett*1e-4, tt, bounds_error=False, fill_value=np.nan)
        #                 phot2counts = liref(w*1e-4)
        #                 data['PHT2CNT'][index] = phot2counts
        #                 data['DISPERSION'][index] = 2.70
        #                 data['SHIFT'][index] = -1079.96
        #             pass
        pass
    allignore = data['IGNORED']
    allculprits = data['TRIAL']
    log.info(
        '>-- IGNORED: %s / %s', str(np.nansum(allignore)), str(len(allignore))
    )
    for index, ignore in enumerate(allignore):
        if ignore:
            log.info('>-- %s: %s', str(index), str(allculprits[index]))
        pass
    data.pop('EXP', None)
    for k, v in data.items():
        out['data'][k] = v
    calibrated = not np.all(data['IGNORED'])
    if calibrated:
        out['STATUS'].append(True)
    return calibrated


# ---------------------- ---------------------------------------------
# -------------------------- -----------------------------------------
# -- STIS CALIBRATION -- ---------------------------------------------
def stiscal_G430L(fin, clc, tim, tid, flttype, out):
    '''
    R. ESTRELA: STIS .flt data extraction and wavelength calibration
    '''
    calibrated = False
    # VISIT NUMBERING --------------------------------------------------------------------
    for pkey in tim['data'].keys():
        visits = np.array(tim['data'][pkey]['dvisits'])
    # PHASE ------------------------------------------------------------------------------
    # for pkey in tim['data'].keys():
    #    phase = np.array(tim['data'][pkey]['phase'])
    # OPTICS AND FILTER ------------------------------------------------------------------
    vrange = validrange(flttype)
    _wvrng, _disp, ldisp, udisp = fng(flttype)
    # DATA FORMAT ------------------------------------------------------------------------
    dbs = os.path.join(dawgie.context.data_dbs, 'mast')
    data = {
        'LOC': [],
        'EPS': [],
        'DISPLIM': [ldisp, udisp],
        'SCANRATE': [],
        'SCANLENGTH': [],
        'SCANANGLE': [],
        'EXP': [],
        'EXPERR': [],
        'EXPFLAG': [],
        'VRANGE': vrange,
        'TIME': [],
        'EXPLEN': [],
        'MIN': [],
        'MAX': [],
        'TRIAL': [],
    }
    # LOAD DATA --------------------------------------------------------------------------
    for loc in sorted(clc['LOC']):
        fullloc = os.path.join(dbs, loc)
        with pyfits.open(fullloc) as hdulist:
            header0 = hdulist[0].header
            if 'SCAN_ANG' in header0:
                scanangle = header0['SCAN_ANG']
            elif 'PA_V3' in header0:
                scanangle = header0['PA_V3']
            else:
                scanangle = 666
            allloc = []
            alltime = []
            allexplen = []
            alleps = []
            allexp = []
            allexperr = []
            allmask = []
            allmin = []
            allmax = []
            for fits in hdulist:
                if (fits.size != 0) and (fits.header['EXTNAME'] == 'SCI'):
                    allloc.append(fits.header['EXPNAME'])
                    alltime.append(float(fits.header['EXPEND']))
                    allexplen.append(float(fits.header['EXPTIME']))
                    fitsdata = np.empty(fits.data.shape)
                    fitsdata[:] = fits.data[:]
                    test = fits.header['BUNIT']
                    eps = False
                    if test != 'COUNTS':
                        eps = True
                    alleps.append(eps)
                    allmin.append(float(fits.header['GOODMIN']))
                    allmax.append(float(fits.header['GOODMAX']))
                    # BINARIES
                    # GMR: Let's put that in the mask someday
                    if tid in ['HAT-P-1']:
                        allexp.append(fitsdata[0:120, :])
                    else:
                        allexp.append(fitsdata)
                    del fits.data
                    pass
                if 'EXTNAME' in fits.header:
                    if fits.header['EXTNAME'] in ['ERR', 'DQ']:
                        fitsdata = np.empty(fits.data.shape)
                        fitsdata[:] = fits.data[:]
                        if fits.header['EXTNAME'] == 'ERR':
                            allexperr.append(fitsdata)
                        if fits.header['EXTNAME'] == 'DQ':
                            allmask.append(fitsdata)
                        del fits.data
                        pass
                    if eps:
                        eps2count = allexplen[-1] * float(header0['CCDGAIN'])
                        allexp[-1] = allexp[-1] * eps2count
                        allexperr[-1] = allexperr[-1] * eps2count
                        pass
                    pass
                pass
            allscanangle = [scanangle] * len(allloc)
            allscanlength = [1e0] * len(allloc)
            allscanrate = [0e0] * len(allloc)
            data['LOC'].extend(allloc)
            data['EPS'].extend(alleps)
            data['SCANRATE'].extend(allscanrate)
            data['SCANLENGTH'].extend(allscanlength)
            data['SCANANGLE'].extend(allscanangle)
            data['EXP'].extend(allexp)
            data['EXPERR'].extend(allexperr)
            data['TIME'].extend(alltime)
            data['EXPLEN'].extend(allexplen)
            data['MIN'].extend(allmin)
            data['MAX'].extend(allmax)
            pass
        pass
    data['MEXP'] = data['EXP'].copy()
    data['MASK'] = data['EXPFLAG'].copy()
    data['IGNORED'] = np.array([False] * len(data['LOC']))
    data['FLOODLVL'] = [np.nan] * len(data['LOC'])
    data['TRIAL'] = [''] * len(data['LOC'])
    data['SPECTRUM0'] = [np.array([np.nan])] * len(data['LOC'])
    data['SPECTRUM'] = [np.array([np.nan])] * len(data['LOC'])
    data['SPECTRUM_CLEAN'] = [np.array([np.nan])] * len(data['LOC'])
    data['SPECERR'] = [np.array([np.nan])] * len(data['LOC'])
    data['TEMPLATE'] = [np.array([np.nan])] * len(data['LOC'])
    data['PHT2CNT'] = [np.array([np.nan])] * len(data['LOC'])
    data['NSPEC'] = [1e0] * len(data['LOC'])
    # REJECT OUTLIERS IN EXPOSURE LENGTH -------------------------------------------------
    for v in set(visits):
        select = visits == v
        visitexplength = np.array(data['EXPLEN'])[select]
        visitignore = data['IGNORED'][select]
        ref = np.nanmedian(visitexplength)
        visitignore[visitexplength != ref] = True
        data['IGNORED'][select] = visitignore
        pass
    # COSMIC RAYS REJECTION - MEDIAN FILTER + SIGMA CLIPPING
    for index, ignore in enumerate(data['IGNORED']):
        # COSMIC RAY REJECTION IN THE 2D IMAGE
        frame = data['MEXP'][index].copy()
        img_cr = frame.copy()
        allframe_list = []
        for i in range(0, len(frame)):
            img_sm = scipy.signal.medfilt(img_cr[i, :], 9)
            std = np.std(img_cr[i, :] - img_sm)
            # std = np.std(img_sm)
            bad = np.abs(img_cr[i, :] - img_sm) > 2 * std
            line = img_cr[i, :]
            line[bad] = img_sm[bad]
            img_sm2 = scipy.signal.medfilt(line, 9)
            std2 = np.std(line - img_sm2)
            bad2 = np.abs(line - img_sm2) > 2 * std2
            line2 = line.copy()
            line2[bad2] = img_sm2[bad2]
            allframe_list.append(line2)
            pass
        allframe = np.array(allframe_list)
        if not ignore:
            data['SPECTRUM'][index] = np.nansum(allframe, axis=0)
            data['SPECERR'][index] = np.sqrt(np.nansum(allframe, axis=0))
            data['PHT2CNT'][index] = [np.nan] * len(frame[0])
            pass
        else:
            data['SPECTRUM'][index] = np.nansum(allframe, axis=0) * np.nan
            data['SPECERR'][index] = np.nansum(allframe, axis=0) * np.nan
            data['TRIAL'][index] = 'Exposure Length Outlier'
            data['PHT2CNT'][index] = [np.nan] * len(frame[0])
            pass
        pass
    pass
    wavett, tt = ag2ttf(flttype)
    if 'G430' in flttype:
        select = wavett > 0.29e4
        wavett = wavett[select]
        tt = tt[select]
        pass
    # MASK BAD PIXELS IN SPECTRUM --------------------------------------------------------
    for v in set(visits):
        select = (visits == v) & ~(data['IGNORED'])
        specarray = np.array(
            [s for s, ok in zip(data['SPECTRUM'], select) if ok]
        )
        trans = np.transpose(specarray)
        template = np.nanmedian(trans, axis=1)
        # TEMPLATE MEDIAN 5 POINTS LOW PASS FILTER ---------------------------------------
        smootht = []
        smootht.extend([template[0]] * 2)
        smootht.extend(template)
        smootht.extend([template[-1]] * 2)
        for index in np.arange(len(template)):
            medianvalue = np.nanmedian(template[index : index + 5])
            smootht[2 + index] = medianvalue
            pass
        smootht = smootht[2:-2]
        template = np.array(smootht)
        for vindex, valid in enumerate(select):
            if valid:
                data['TEMPLATE'][vindex] = template
            pass
        pass
    # COSMIC RAYS REJECTION
    # data['PHT2CNT'] = [np.nan]*len(data['LOC'])
    data['WAVE'] = [np.array([np.nan])] * len(data['LOC'])
    data['DISPERSION'] = [np.nan] * len(data['LOC'])
    data['SHIFT'] = [np.nan] * len(data['LOC'])

    set_wav = np.array([290, 570])

    def phoenix(set_wav):
        # PHOENIX MODELS
        filters = [BoxcarFilter('a', 300, 570)]  # Define your passbands
        feherr = np.sqrt(
            abs(fin['priors']['FEH*_uperr'] * fin['priors']['FEH*_lowerr'])
        )
        loggerr = np.sqrt(
            abs(fin['priors']['LOGG*_uperr'] * fin['priors']['LOGG*_lowerr'])
        )
        terr = np.sqrt(
            abs(fin['priors']['T*_uperr'] * fin['priors']['T*_lowerr'])
        )

        # 6/11/24 GB
        # bug for LP 791-18 when running LDPSetCreator called from stiscal_G430L()
        #  the parameters passed into it are fine:
        #   teff 2960.0 55.0
        #   logg 3.0512 0.09
        #   FEH* -0.09 0.19
        # print('calling LDPSetCreator')
        # print('   teff',fin['priors']['T*'], terr)
        # print('   logg',fin['priors']['b']['logg'], loggerr)
        # print('   FEH*',fin['priors']['FEH*'], feherr)
        # print('  filters',filters)
        loggerr = np.max(loggerr, 0.1)

        sc = LDPSetCreator(
            teff=(fin['priors']['T*'], terr),
            logg=(fin['priors']['b']['logg'], loggerr),
            z=(fin['priors']['FEH*'], feherr),
            filters=filters,
        )
        list_diff = []
        for thisfile in sc.files:
            hdul = pyfits.open(thisfile)
            teff = hdul[0].header['PHXTEFF']
            zz = hdul[0].header['PHXM_H']
            logg_str = hdul[0].header['PHXLOGG']
            diff1 = abs(fin['priors']['T*'] - teff)
            diff2 = abs(fin['priors']['FEH*'] - zz)
            diff3 = abs(fin['priors']['b']['logg'] - logg_str)
            diff_total = diff1 + diff2 + diff3
            list_diff.append(diff_total)
            pass
        cond_win = np.where(list_diff == np.min(list_diff))
        hdul2 = pyfits.open(
            sc.files[cond_win[0][0]]
        )  # 1 for HAT-p-26 and Hat-P-11, 3 for Hat-p-18, 1 for WASP-52, 4 for WASP-80
        data_all = hdul2[0].data
        wl0 = (
            hdul[0].header['crval1'] * 1e-1
        )  # defines the wavelength at pixel CRPIX1
        dwl = hdul[0].header['cdelt1'] * 1e-1  # Delta wavelength     [nm]
        nwl = hdul[0].header['naxis1']  # Number of wl samples
        wl = wl0 + np.arange(nwl) * dwl
        model = data_all[77]  # take only the last spectra of each fits
        # Average the spectra to get 1 spectrum model
        # new_spec=[]
        # trans_listdata = np.transpose(list_models)
        # for i in range(0,len(trans_listdata)):
        #    med_wav = np.mean(trans_listdata[i])
        #    new_spec.append(med_wav)
        # f_spec = itp.interp1d(wl, new_spec, bounds_error=False)
        # spec_sel = f_spec(wavett*0.1)
        # spec_sel_norm = spec_sel/np.max(spec_sel)
        cond_wav = np.where((wl > set_wav[0]) & (wl < set_wav[1]))
        wl_sel = wl[cond_wav]
        new_spec = np.array(model)
        spec_sm = scipy.signal.medfilt(new_spec, 9)
        new_spec_sel = spec_sm[cond_wav]
        mid, low, high, binsz = binnagem(wl_sel, 1024)
        # func_spec = scipy.interpolate.interp1d(wl_sel,new_spec_sel)
        # BINNING PHOENIX MODEL
        bin_spec = []
        for w_low, w_hi in zip(low, high):
            select = np.where((wl_sel > w_low) & (wl_sel < w_hi))
            # inte = scipy.integrate.quad(lambda x: func_spec(x), w_low, w_hi)
            inte = np.sum(
                new_spec_sel[select]
                * (wl_sel[select[0] + 1] - wl_sel[select[0]])
            )
            databin = inte / binsz[0]  # inte[0] if scipy.integrate
            bin_spec.append(databin)
        bin_spec = np.array(bin_spec)
        bin_spec = scipy.signal.medfilt(bin_spec, 5)
        return mid, bin_spec

    def chisqfunc(args, g_wav, bin_spec_norm, cond_mid, f, mid_ang):
        avar, bvar = args
        # 2/13/25 Geoff: the missing parameters are now passed in as an arg() for scipy.optimize
        #  but I'm not sure if the above syntax is correct.  try running this later..
        chisq = np.sum(
            (
                (g_wav * bin_spec_norm[cond_mid])
                - f(bvar + (avar) * mid_ang[cond_mid])
            )
            ** 2
        )
        return chisq

    dispersion_list = []
    for index, rejected in enumerate(data['IGNORED']):
        if not rejected:
            spec = data['SPECTRUM'][index]
            template = data['TEMPLATE'][index]
            temp_spec = spec / template
            ht25 = np.nanpercentile(temp_spec, 25)
            lt75 = np.nanpercentile(temp_spec, 75)
            std = np.std(temp_spec[(temp_spec > ht25) & (temp_spec < lt75)])
            # BAD PIXEL THRESHOLD --------------------------------------------------------
            bpthr = temp_spec > np.nanmedian(temp_spec) + 2e0 * std
            if True in bpthr:
                spec[bpthr] = np.nan
            # second
            temp_spec2 = spec / template
            ht25 = np.nanpercentile(temp_spec2, 25)
            lt75 = np.nanpercentile(temp_spec2, 75)
            selfin = np.isfinite(temp_spec2)
            std1 = np.nanstd(
                temp_spec2[selfin][
                    (temp_spec2[selfin] > ht25) & (temp_spec2[selfin] < lt75)
                ]
            )
            # BAD PIXEL THRESHOLD --------------------------------------------------------
            bpthr = temp_spec2[selfin] > np.nanmedian(temp_spec2) + 2e0 * std1
            spec_cut = spec.copy()
            data['SPECTRUM'][index] = spec_cut
            data['SPECTRUM_CLEAN'][index] = spec_cut
            # WAVELENGTH CALIBRATION -----------------------------------------------------
            disp_all = []
            if np.sum(np.isfinite(spec)) > (spec.size / 2):
                # wavecalspec = spec.copy()
                wavecalspec = spec[:-1]
                finitespec = np.isfinite(wavecalspec)
                spec_norm = wavecalspec[finitespec] / np.max(
                    wavecalspec[finitespec]
                )
                phoenix_model = phoenix(set_wav)
                bin_spec = phoenix_model[1]
                mid = phoenix_model[0]
                bin_spec_norm = bin_spec / np.max(bin_spec)
                # select=spec_norm > 1e-1
                x = np.arange(len(wavecalspec))
                x_finite = x[finitespec]
                th_norm = tt / np.max(tt)
                f = itp.interp1d(
                    x_finite, spec_norm, bounds_error=False, fill_value=0
                )
                # f_x = f(x)
                mid_ang = mid * 10
                g = itp.interp1d(
                    wavett, th_norm, bounds_error=False, fill_value=0
                )
                cond_mid = np.where(
                    (mid_ang >= wavett[0]) & (mid_ang <= wavett[-1])
                )
                g_wav = g(mid_ang[cond_mid])
                # model = scipy.signal.medfilt(g_wav*bin_spec_norm[cond_mid], 5)
                # wave = np.arange(spec.size)*disper*1e-4 + shift
                x0 = (1.0 / 2.72, -1000)
                result = opt.minimize(
                    chisqfunc,
                    x0,
                    args=(g_wav, bin_spec_norm, cond_mid, f, mid_ang),
                    method='Nelder-Mead',
                )
                d_frc = result.x[0]
                d = 1.0 / result.x[0]
                dispersion_list.append(d)
                s = result.x[1]
                calib_spec = f(s + (d_frc) * mid_ang[cond_mid])
                data['SPECTRUM'][index] = calib_spec * np.max(
                    wavecalspec[finitespec]
                )
                liref = itp.interp1d(
                    wavett, tt, bounds_error=False, fill_value=np.nan
                )
                phot2counts = liref(mid_ang[cond_mid])
                data['PHT2CNT'][index] = phot2counts
                data['WAVE'][index] = mid[cond_mid] * 0.001
                data['DISPERSION'][index] = d
                data['SHIFT'][index] = s
                err = data['SPECERR'][index]
                data['SPECERR'][index] = err[cond_mid]
                disp_all.append(np.median(dispersion_list))
                pass
            pass

    for v in set(visits):
        select = (visits == v) & ~(data['IGNORED'])
        for index, valid in enumerate(select):
            spec_all = np.array(data['SPECTRUM'][index])
            wave_all = np.array(data['WAVE'][index])
            phot_all = np.array(data['PHT2CNT'][index])
            err_all = np.array(data['SPECERR'][index])
            cond_wavcut = np.where(wave_all > 0.45)
            data['SPECTRUM'][index] = spec_all[cond_wavcut]
            data['WAVE'][index] = wave_all[cond_wavcut]
            data['PHT2CNT'][index] = phot_all[cond_wavcut]
            data['SPECERR'][index] = err_all[cond_wavcut]
            pass
        pass
    allignore = data['IGNORED']
    allculprits = data['TRIAL']
    log.info(
        '>-- IGNORED: %s / %s', str(np.nansum(allignore)), str(len(allignore))
    )
    for index, ignore in enumerate(allignore):
        if ignore:
            log.info('>-- %s: %s', str(index), str(allculprits[index]))
        pass
    data.pop('EXP', None)
    for k, v in data.items():
        out['data'][k] = v
    calibrated = not np.all(data['IGNORED'])
    if calibrated:
        out['STATUS'].append(True)
    return calibrated


# ---------------------- ---------------------------------------------
# -------------------------- -----------------------------------------


def stiscal_unified(fin, clc, tim, tid, flttype, out):
    '''
    R. ESTRELA: STIS .flt data extraction and wavelength calibration FILTERS G430L and G750L
    '''
    calibrated = False
    # VISIT NUMBERING --------------------------------------------------------------------
    for pkey in tim['data'].keys():
        visits = np.array(tim['data'][pkey]['dvisits'])
    # PHASE ------------------------------------------------------------------------------
    # for pkey in tim['data'].keys():
    #    phase = np.array(tim['data'][pkey]['phase'])
    # OPTICS AND FILTER ------------------------------------------------------------------
    vrange = validrange(flttype)
    _wvrng, _disp, ldisp, udisp = fng(flttype)
    # DATA FORMAT ------------------------------------------------------------------------
    dbs = os.path.join(dawgie.context.data_dbs, 'mast')
    data = {
        'LOC': [],
        'EPS': [],
        'DISPLIM': [ldisp, udisp],
        'SCANRATE': [],
        'SCANLENGTH': [],
        'SCANANGLE': [],
        'EXP': [],
        'EXPERR': [],
        'EXPFLAG': [],
        'VRANGE': vrange,
        'TIME': [],
        'EXPLEN': [],
        'MIN': [],
        'MAX': [],
        'TRIAL': [],
        'TIMEOBS': [],
        'DATEOBS': [],
    }
    # LOAD DATA --------------------------------------------------------------------------
    for loc in sorted(clc['LOC']):
        fullloc = os.path.join(dbs, loc)
        with pyfits.open(fullloc) as hdulist:
            header0 = hdulist[0].header
            if 'SCAN_ANG' in header0:
                scanangle = header0['SCAN_ANG']
            elif 'PA_V3' in header0:
                scanangle = header0['PA_V3']
            else:
                scanangle = 666
            allloc = []
            alltime = []
            allexplen = []
            alleps = []
            allexp = []
            allexperr = []
            allmask = []
            allmin = []
            allmax = []
            alldate = []
            alltimeobs = []
            for fits in hdulist:
                if (fits.size != 0) and (fits.header['EXTNAME'] == 'SCI'):
                    allloc.append(fits.header['EXPNAME'])
                    alltime.append(float(fits.header['EXPEND']))
                    allexplen.append(float(fits.header['EXPTIME']))
                    alldate.append(header0['TDATEOBS'])
                    alltimeobs.append(header0['TTIMEOBS'])
                    fitsdata = np.empty(fits.data.shape)
                    fitsdata[:] = fits.data[:]
                    test = fits.header['BUNIT']
                    eps = False
                    if test != 'COUNTS':
                        eps = True
                    alleps.append(eps)
                    allmin.append(float(fits.header['GOODMIN']))
                    allmax.append(float(fits.header['GOODMAX']))
                    # BINARIES
                    # GMR: Let's put that in the mask someday
                    if tid in ['HAT-P-1']:
                        allexp.append(fitsdata[0:120, :])
                    else:
                        allexp.append(fitsdata)
                    del fits.data
                    pass
                if 'EXTNAME' in fits.header:
                    if fits.header['EXTNAME'] in ['ERR', 'DQ']:
                        fitsdata = np.empty(fits.data.shape)
                        fitsdata[:] = fits.data[:]
                        if fits.header['EXTNAME'] == 'ERR':
                            allexperr.append(fitsdata)
                        if fits.header['EXTNAME'] == 'DQ':
                            allmask.append(fitsdata)
                        del fits.data
                        pass
                    if eps:
                        eps2count = allexplen[-1] * float(header0['CCDGAIN'])
                        allexp[-1] = allexp[-1] * eps2count
                        allexperr[-1] = allexperr[-1] * eps2count
                        pass
                    pass
                pass
            allscanangle = [scanangle] * len(allloc)
            allscanlength = [1e0] * len(allloc)
            allscanrate = [0e0] * len(allloc)
            data['LOC'].extend(allloc)
            data['EPS'].extend(alleps)
            data['SCANRATE'].extend(allscanrate)
            data['SCANLENGTH'].extend(allscanlength)
            data['SCANANGLE'].extend(allscanangle)
            data['EXP'].extend(allexp)
            data['EXPERR'].extend(allexperr)
            data['TIME'].extend(alltime)
            data['EXPLEN'].extend(allexplen)
            data['MIN'].extend(allmin)
            data['MAX'].extend(allmax)
            data['TIMEOBS'].extend(alltimeobs)
            data['DATEOBS'].extend(alldate)
            pass
        pass
    data['MEXP'] = data['EXP'].copy()
    data['MASK'] = data['EXPFLAG'].copy()
    data['ALLDATEOBS'] = data['DATEOBS'].copy()
    data['ALLTIMEOBS'] = data['TIMEOBS'].copy()
    data['IGNORED'] = np.array([False] * len(data['LOC']))
    data['FLOODLVL'] = [np.nan] * len(data['LOC'])
    data['TRIAL'] = [''] * len(data['LOC'])
    data['SPECTRUM0'] = [np.array([np.nan])] * len(data['LOC'])
    data['SPECTRUM'] = [np.array([np.nan])] * len(data['LOC'])
    data['SPECTRUM_CLEAN'] = [np.array([np.nan])] * len(data['LOC'])
    data['SPECERR'] = [np.array([np.nan])] * len(data['LOC'])
    data['TEMPLATE'] = [np.array([np.nan])] * len(data['LOC'])
    data['PHT2CNT'] = [np.array([np.nan])] * len(data['LOC'])
    data['NSPEC'] = [1e0] * len(data['LOC'])
    # REJECT OUTLIERS IN EXPOSURE LENGTH -------------------------------------------------
    for v in set(visits):
        select = visits == v
        visitexplength = np.array(data['EXPLEN'])[select]
        visitignore = data['IGNORED'][select]
        ref = np.nanmedian(visitexplength)
        visitignore[visitexplength != ref] = True
        data['IGNORED'][select] = visitignore
        pass
    # COSMIC RAYS REJECTION - MEDIAN FILTER + SIGMA CLIPPING
    for index, ignore in enumerate(data['IGNORED']):
        if 'G430L' in flttype:
            frame = data['MEXP'][index].copy()
            img_cr = frame.copy()
            allframe_list = []
            for i in range(0, len(frame)):
                img_sm = scipy.signal.medfilt(img_cr[i, :], 9)
                std = np.std(img_cr[i, :] - img_sm)
                # std = np.std(img_sm)
                bad = np.abs(img_cr[i, :] - img_sm) > 2 * std
                line = img_cr[i, :]
                line[bad] = img_sm[bad]
                img_sm2 = scipy.signal.medfilt(line, 9)
                std2 = np.std(line - img_sm2)
                bad2 = np.abs(line - img_sm2) > 2 * std2
                line2 = line.copy()
                line2[bad2] = img_sm2[bad2]
                allframe_list.append(line2)
                pass
            allframe = np.array(allframe_list)

        # FLAT FRINGE G750L
    for index, ignore in enumerate(data['IGNORED']):
        if 'G750L' in flttype:
            # SELECT DATE AND TIME OF THE EXPOSURE FOR FLAT FRINGE SELECTION
            frame = data['MEXP'][index].copy()
            dateobs_exp = data['ALLDATEOBS'][index]
            timeobs_exp = data['ALLTIMEOBS'][index]
            tog_exp = dateobs_exp + ' ' + timeobs_exp
            time_exp = raissatime.mktime(
                datetime.datetime.strptime(
                    tog_exp, "%Y-%m-%d %H:%M:%S"
                ).timetuple()
            )
            # LOAD FRINGE FLAT -------------------------------------------------------------------
            obs_name = clc['ROOTNAME'][0]
            name_sel = obs_name[:-5]
            lightpath_fringe = 'STIS/CCDFLAT/'
            calloc = excalibur.context['data_cal']
            filefringe = os.path.join(calloc, lightpath_fringe)
            if tid in ['HD 209458']:
                lightpath_fringe = 'STIS/CCDFLAT/h230851ao_pfl.fits'
                calloc = excalibur.context['data_cal']
                filefringe = os.path.join(calloc, lightpath_fringe)
                hdu = pyfits.open(filefringe)
                data_fringe = hdu[1].data
                err_fringe = hdu[2].data
                pass
            else:
                diff_list = []
                all_infile = []
                for infile in glob.glob(f'{filefringe}/{name_sel}*_flt.fits'):
                    hdu = pyfits.open(infile)
                    all_infile.append(infile)
                    header_flat = hdu[0].header
                    date_time = header_flat['TDATEOBS']
                    hour_time = header_flat['TTIMEOBS']
                    tog = date_time + ' ' + hour_time
                    time_flat_s = raissatime.mktime(
                        datetime.datetime.strptime(
                            tog, "%Y-%m-%d %H:%M:%S"
                        ).timetuple()
                    )
                    diff = abs(time_exp - time_flat_s)
                    diff_list.append(diff)
                    pass
                cond_win = np.where(diff_list == np.min(diff_list))
                all_infile = np.array(all_infile)
                sel_flatfile = all_infile[cond_win][0]
                hdulist = pyfits.open(sel_flatfile)
                data_fringe = hdulist[4].data
                err_fringe = hdulist[5].data
                pass
            smooth_fringe = scipy.signal.medfilt(data_fringe, 7)
            sigma_fringe = np.median(err_fringe)
            bad_fringe = (
                np.abs(data_fringe - smooth_fringe) / sigma_fringe
            ) > 2
            img_fringe = data_fringe.copy()
            img_fringe[bad_fringe] = smooth_fringe[bad_fringe]
            cont_data = img_fringe.copy()
            div_list = []
            for ll in range(508, 515):
                pixels = np.arange(0, 1024, 1)
                coefs = poly.polyfit(pixels, cont_data[ll, :], 11)
                ffit = poly.polyval(pixels, coefs)
                div = cont_data[ll, :] / ffit
                div_list.append(div)
                # COSMIC RAY REJECTION IN THE 2D IMAGE
            img_cr = frame.copy()
            allframe_list = []
            for i in range(0, len(frame)):
                img_sm = scipy.signal.medfilt(img_cr[i, :], 9)
                # std = np.std(img_cr[i,:] - img_sm)
                std = np.std(img_sm)
                bad = np.abs(img_cr[i, :] - img_sm) > 3 * std
                line = img_cr[i, :]
                line[bad] = img_sm[bad]
                allframe_list.append(line)
                pass
            allframe = np.array(allframe_list)
            # APPLY FLAT FRINGE
            if not ignore:
                find_spec = np.where(allframe == np.max(allframe))
                spec_idx = find_spec[0][0]
                spec_idx_up = spec_idx + 4
                spec_idx_dwn = spec_idx - 3
                spec_idx_all = np.arange(spec_idx_dwn, spec_idx_up, 1)
                frame2 = allframe.copy()
                for i, flatnorm in zip(spec_idx_all, div_list):
                    # frame_sel = allframe[i, :]
                    # coefs_f = poly.polyfit(pixels, frame_sel, 12)
                    # ffit_f = poly.polyval(pixels, coefs_f)
                    frame2[i, 400:1023] = (
                        frame2[i, 400:1023] / flatnorm[400:1023]
                    )
                    data['SPECTRUM'][index] = np.nansum(frame2, axis=0)
                    data['SPECERR'][index] = np.sqrt(np.nansum(frame2, axis=0))
                    pass
                pass
            else:
                data['SPECTRUM'][index] = np.nansum(frame, axis=0) * np.nan
                data['SPECERR'][index] = np.nansum(frame, axis=0) * np.nan
                data['TRIAL'][index] = 'Exposure Length Outlier'
                pass
            pass
        if 'G430' in flttype:
            if not ignore:
                data['SPECTRUM'][index] = np.nansum(allframe, axis=0)
                data['SPECERR'][index] = np.sqrt(np.nansum(allframe, axis=0))
                data['PHT2CNT'][index] = [np.nan] * len(frame[0])
                pass
            else:
                data['SPECTRUM'][index] = np.nansum(allframe, axis=0) * np.nan
                data['SPECERR'][index] = np.nansum(allframe, axis=0) * np.nan
                data['TRIAL'][index] = 'Exposure Length Outlier'
                data['PHT2CNT'][index] = [np.nan] * len(frame[0])
                pass
            pass
    pass
    wavett, tt = ag2ttf(flttype)
    if 'G430' in flttype:
        select = wavett > 0.29e4
        wavett = wavett[select]
        tt = tt[select]
        pass
    # MASK BAD PIXELS IN SPECTRUM ----------------------------------------------------
    for v in set(visits):
        select = (visits == v) & ~(data['IGNORED'])
        specarray = np.array(
            [s for s, ok in zip(data['SPECTRUM'], select) if ok]
        )
        trans = np.transpose(specarray)
        template = np.nanmedian(trans, axis=1)
        # TEMPLATE MEDIAN 5 POINTS LOW PASS FILTER ---------------------------------------
        smootht = []
        smootht.extend([template[0]] * 2)
        smootht.extend(template)
        smootht.extend([template[-1]] * 2)
        for index in np.arange(len(template)):
            medianvalue = np.nanmedian(template[index : index + 5])
            smootht[2 + index] = medianvalue
            pass
        smootht = smootht[2:-2]
        template = np.array(smootht)
        for vindex, valid in enumerate(select):
            if valid:
                data['TEMPLATE'][vindex] = template
            pass
        pass
    # COSMIC RAYS REJECTION
    # data['PHT2CNT'] = [np.nan]*len(data['LOC'])
    data['WAVE'] = [np.array([np.nan])] * len(data['LOC'])
    data['DISPERSION'] = [np.nan] * len(data['LOC'])
    data['SHIFT'] = [np.nan] * len(data['LOC'])

    if 'G430' in flttype:
        set_wav = np.array([290, 570])
    elif 'G750' in flttype:
        set_wav = np.array([524, 1027])
    else:
        log.warning('DATA stiscal_unified: UNKNOWN HST FILTER %s', flttype)
        set_wav = np.array([666, -666])

    def phoenix(set_wav):
        # PHOENIX MODELS
        filters = [BoxcarFilter('a', 550, 950)]  # Define your passbands
        feherr = np.sqrt(
            abs(fin['priors']['FEH*_uperr'] * fin['priors']['FEH*_lowerr'])
        )
        loggerr = np.sqrt(
            abs(fin['priors']['LOGG*_uperr'] * fin['priors']['LOGG*_lowerr'])
        )
        terr = np.sqrt(
            abs(fin['priors']['T*_uperr'] * fin['priors']['T*_lowerr'])
        )
        sc = LDPSetCreator(
            teff=(fin['priors']['T*'], terr),
            logg=(fin['priors']['b']['logg'], loggerr),
            z=(fin['priors']['FEH*'], feherr),
            filters=filters,
        )
        list_diff = []
        for thisfile in sc.files:
            hdul = pyfits.open(thisfile)
            teff = hdul[0].header['PHXTEFF']
            zz = hdul[0].header['PHXM_H']
            logg_str = hdul[0].header['PHXLOGG']
            diff1 = abs(fin['priors']['T*'] - teff)
            diff2 = abs(fin['priors']['FEH*'] - zz)
            diff3 = abs(fin['priors']['b']['logg'] - logg_str)
            diff_total = diff1 + diff2 + diff3
            list_diff.append(diff_total)
            pass
        cond_win = np.where(list_diff == np.min(list_diff))
        hdul2 = pyfits.open(
            sc.files[cond_win[0][0]]
        )  # 1 for HAT-p-26 and Hat-P-11, 3 for Hat-p-18, 1 for WASP-52, 4 for WASP-80
        data_all = hdul2[0].data
        wl0 = (
            hdul[0].header['crval1'] * 1e-1
        )  # defines the wavelength at pixel CRPIX1
        dwl = hdul[0].header['cdelt1'] * 1e-1  # Delta wavelength     [nm]
        nwl = hdul[0].header['naxis1']  # Number of wl samples
        wl = wl0 + np.arange(nwl) * dwl
        model = data_all[77]  # take only the last spectra of each fits
        # Average the spectra to get 1 spectrum model
        # new_spec=[]
        # trans_listdata = np.transpose(list_models)
        # for i in range(0,len(trans_listdata)):
        #    med_wav = np.mean(trans_listdata[i])
        #    new_spec.append(med_wav)
        # f_spec = itp.interp1d(wl, new_spec, bounds_error=False)
        # spec_sel = f_spec(wavett*0.1)
        # spec_sel_norm = spec_sel/np.max(spec_sel)
        cond_wav = np.where((wl > set_wav[0]) & (wl < set_wav[1]))
        wl_sel = wl[cond_wav]
        new_spec = np.array(model)
        if 'G430' in flttype:
            spec_sm = scipy.signal.medfilt(new_spec, 9)
        else:
            spec_sm = new_spec
        new_spec_sel = spec_sm[cond_wav]
        mid, low, high, binsz = binnagem(wl_sel, 1024)
        # func_spec = scipy.interpolate.interp1d(wl_sel,new_spec_sel)
        # BINNING PHOENIX MODEL
        bin_spec = []
        for w_low, w_hi in zip(low, high):
            select = np.where((wl_sel > w_low) & (wl_sel < w_hi))
            # inte = scipy.integrate.quad(lambda x: func_spec(x), w_low, w_hi)
            inte = np.sum(
                new_spec_sel[select]
                * (wl_sel[select[0] + 1] - wl_sel[select[0]])
            )
            databin = inte / binsz[0]  # inte[0] if scipy.integrate
            bin_spec.append(databin)
        bin_spec = np.array(bin_spec)
        if 'G430' in flttype:
            window = 5
        else:
            window = 5
        bin_spec = scipy.signal.medfilt(bin_spec, window)
        return mid, bin_spec

    def chisqfunc(args, g_wav, bin_spec_norm, cond_mid, f, mid_ang):
        avar, bvar, scvar = args
        # 2/13/25 Geoff: the missing parameters are now passed in as an arg() for scipy.optimize
        #  but I'm not sure if the above syntax is correct.  try running this later..
        chisq = np.sum(
            (
                (g_wav * bin_spec_norm[cond_mid]) * scvar
                - f(bvar + (avar) * mid_ang[cond_mid])
            )
            ** 2
        )
        return chisq

    # CALIBRATION
    dispersion_list = []
    for index, rejected in enumerate(data['IGNORED']):
        if not rejected:
            spec = data['SPECTRUM'][index]
            template = data['TEMPLATE'][index]
            temp_spec = spec / template
            ht25 = np.nanpercentile(temp_spec, 25)
            lt75 = np.nanpercentile(temp_spec, 75)
            std = np.std(temp_spec[(temp_spec > ht25) & (temp_spec < lt75)])
            # BAD PIXEL THRESHOLD --------------------------------------------------------
            bpthr = temp_spec > np.nanmedian(temp_spec) + 2e0 * std
            if True in bpthr:
                spec[bpthr] = np.nan
            # second
            temp_spec2 = spec / template
            ht25 = np.nanpercentile(temp_spec2, 25)
            lt75 = np.nanpercentile(temp_spec2, 75)
            selfin = np.isfinite(temp_spec2)
            std1 = np.nanstd(
                temp_spec2[selfin][
                    (temp_spec2[selfin] > ht25) & (temp_spec2[selfin] < lt75)
                ]
            )
            # BAD PIXEL THRESHOLD --------------------------------------------------------
            bpthr = temp_spec2[selfin] > np.nanmedian(temp_spec2) + 2e0 * std1
            spec_cut = spec.copy()
            data['SPECTRUM'][index] = spec_cut
            data['SPECTRUM_CLEAN'][index] = spec_cut
            # WAVELENGTH CALIBRATION -----------------------------------------------------
            disp_all = []
            if np.sum(np.isfinite(spec)) > (spec.size / 2):
                wavecalspec = spec
                # wavecalspec = spec.copy()
                if 'G430' in flttype:
                    wavecalspec = spec[:-1]
                finitespec = np.isfinite(wavecalspec)
                spec_norm = wavecalspec[finitespec] / np.max(
                    wavecalspec[finitespec]
                )
                phoenix_model = phoenix(set_wav)
                bin_spec = phoenix_model[1]
                mid = phoenix_model[0]
                bin_spec_norm = bin_spec / np.max(bin_spec)
                # select=spec_norm > 1e-1
                x = np.arange(len(wavecalspec))
                x_finite = x[finitespec]
                th_norm = tt / np.max(tt)
                f = itp.interp1d(
                    x_finite, spec_norm, bounds_error=False, fill_value=0
                )
                # f_x = f(x)
                mid_ang = mid * 10
                g = itp.interp1d(
                    wavett, th_norm, bounds_error=False, fill_value=0
                )
                cond_mid = np.where(
                    (mid_ang >= wavett[0]) & (mid_ang <= wavett[-1])
                )
                g_wav = g(mid_ang[cond_mid])
                # model = scipy.signal.medfilt(g_wav*bin_spec_norm[cond_mid], 5)
                # wave = np.arange(spec.size)*disper*1e-4 + shift
                if 'G750' in flttype:
                    x0 = (1.0 / 4.72, -1000, 1.0)
                else:
                    x0 = (1.0 / 2.72, -1000, 1.0)
                result = opt.minimize(
                    chisqfunc,
                    x0,
                    args=(g_wav, bin_spec_norm, cond_mid, f, mid_ang),
                    method='Nelder-Mead',
                )
                d_frc = result.x[0]
                d = 1.0 / result.x[0]
                dispersion_list.append(d)
                s = result.x[1]
                # sc = result.x[2]  # not used!
                calib_spec = f(s + (d_frc) * mid_ang[cond_mid])
                data['SPECTRUM'][index] = calib_spec * np.max(
                    wavecalspec[finitespec]
                )
                liref = itp.interp1d(
                    wavett, tt, bounds_error=False, fill_value=np.nan
                )
                phot2counts = liref(mid_ang[cond_mid])
                data['PHT2CNT'][index] = phot2counts
                data['WAVE'][index] = mid[cond_mid] * 0.001
                data['DISPERSION'][index] = d
                data['SHIFT'][index] = s
                err = data['SPECERR'][index]
                data['SPECERR'][index] = err[cond_mid]
                disp_all.append(np.median(dispersion_list))
                pass
            pass

    if 'G430' in flttype:
        for v in set(visits):
            select = (visits == v) & ~(data['IGNORED'])
            for index, valid in enumerate(select):
                spec_all = np.array(data['SPECTRUM'][index])
                wave_all = np.array(data['WAVE'][index])
                phot_all = np.array(data['PHT2CNT'][index])
                err_all = np.array(data['SPECERR'][index])
                cond_wavcut = np.where(wave_all > 0.45)
                data['SPECTRUM'][index] = spec_all[cond_wavcut]
                data['WAVE'][index] = wave_all[cond_wavcut]
                data['PHT2CNT'][index] = phot_all[cond_wavcut]
                data['SPECERR'][index] = err_all[cond_wavcut]
                pass
            pass
        pass
    allignore = data['IGNORED']
    allculprits = data['TRIAL']
    log.info(
        '>-- IGNORED: %s / %s', str(np.nansum(allignore)), str(len(allignore))
    )
    for index, ignore in enumerate(allignore):
        if ignore:
            log.info('>-- %s: %s', str(index), str(allculprits[index]))
        pass
    data.pop('EXP', None)
    for k, v in data.items():
        out['data'][k] = v
    calibrated = not np.all(data['IGNORED'])
    if calibrated:
        out['STATUS'].append(True)
    return calibrated


def find_target(target, hdu, verbose=False):
    '''
    query simbad to get proper motions
    http://simbad.u-strasbg.fr/simbad/tap/tapsearch.html
    '''
    service = vo.dal.TAPService("http://simbad.u-strasbg.fr/simbad/sim-tap")
    # GMR: seems this query was formatted in a very specific way, I might have broken it
    query = f"""
    SELECT basic.OID, ra, dec, main_id, pmra, pmdec
    FROM basic JOIN ident ON oidref = oid
    WHERE id = '{target}';
    """

    result = service.search(query)
    # Future check that simbad returned a value

    # set up astropy object
    coord = SkyCoord(
        ra=result['ra'][0] * astropy.units.deg,
        dec=result['dec'][0] * astropy.units.deg,
        distance=1 * astropy.units.pc,
        pm_ra_cosdec=result['pmra'][0] * astropy.units.mas / astropy.units.yr,
        pm_dec=result['pmdec'][0] * astropy.units.mas / astropy.units.yr,
        frame="icrs",
        obstime=Time("2000-1-1T00:00:00"),
    )

    # apply proper motion
    t = Time(hdu.header['DATE_OBS'], format='isot', scale='utc')
    coordpm = coord.apply_space_motion(new_obstime=t)

    # wcs coordinate translation
    try:
        wcs = WCS(hdu.header)
    except ValueError:
        # https://github.com/astropy/astropy/issues/10527
        hdu.header['NAXIS'] = 2
        wcs = WCS(hdu.header)

    pixcoord = wcs.wcs_world2pix([[coordpm.ra.value, coordpm.dec.value]], 0)

    if verbose:
        print("Simbad:", result)
        print("\nObs Date:", hdu.header['DATE_OBS'])
        print("NEW:", coordpm.ra, coordpm.dec)
        print("Pixels:", pixcoord[0])

    return pixcoord[0]


# ---------------------- ---------------------------------------------
# -------------------------- -----------------------------------------
def binnagem(t, nbins):
    '''binnagem ds'''
    tmax = t[-1]
    tmin = t[0]
    tbin = (tmax - tmin) * np.arange(nbins + 1) / nbins
    tbin = tbin + tmin
    lower = np.resize(tbin, len(tbin) - 1)
    tmid = lower + 0.5 * np.diff(tbin)
    higher = tmid + 0.5 * np.diff(tbin)
    binsize = np.diff(tbin)
    return tmid, lower, higher, binsize


# ---------------------- ---------------------------------------------
# -- SPITZER CALIBRATION -- ------------------------------------------
def spitzercal(clc, out):
    '''
    K. PEARSON: SPITZER data extraction
    '''
    calibrated = False
    dbs = os.path.join(dawgie.context.data_dbs, 'mast')

    data = {
        'LOC': [],
        'EXPLEN': [],
        'TIME': [],
        'FRAME': [],
        'NOISEPIXEL': [],
        'PHOT': [],
        'WX': [],
        'WY': [],
        'BG': [],
        'FAILED': [],
    }

    # c = 0
    # LOAD DATA --------------------------------------------------------------------------
    for loc in sorted(clc['LOC']):
        fullloc = os.path.join(dbs, loc)
        with pyfits.open(fullloc) as hdulist:
            for fits in hdulist:
                # if science data
                if (fits.size != 0) and (fits.header.get('exptype') == 'sci'):
                    start = fits.header.get('MJD_OBS') + 2400000.5

                    if fits.data.ndim == 2:
                        # full frame data - todo locate star in field
                        continue
                    if fits.data.ndim == 3:
                        dcube = fits.data.copy()

                    # convert from ADU to e/s
                    dcube /= float(fits.header.get('FLUXCONV', 0.1257))
                    dcube *= float(fits.header.get('GAIN', 3.7))

                    idur = fits.header.get('ATIMEEND') - fits.header.get(
                        'AINTBEG'
                    )
                    nimgs = dcube.shape[0]
                    exptime = idur / nimgs  # sec

                    # convert from Mjy/sr to DN/s then to e/s and finally e
                    dcube *= exptime

                    dt = idur / nimgs / (24 * 60 * 60)
                    dcube[np.isnan(dcube)] = 0
                    dcube[np.isinf(dcube)] = 0

                    # for each subframe
                    for i in range(nimgs):

                        try:
                            wx, wy, apers, bgs, npps = aper_phot(dcube[i])
                            data['TIME'].append(start + i * dt)
                            data['FRAME'].append(i)
                            data['WX'].append(wx)
                            data['WY'].append(wy)

                            # nd arrays
                            data['BG'].append(bgs)
                            data['PHOT'].append(apers)
                            data['NOISEPIXEL'].append(npps)

                            data['FAILED'].append(False)
                        except (ValueError, IndexError):
                            data['FAILED'].append(True)

    # transfer data
    for k, v in data.items():
        out['data'][k] = v
    calibrated = not np.all(data['FAILED'])
    if calibrated:
        out['STATUS'].append(True)
    return calibrated


def aper_phot(img):
    '''aper_phot ds'''
    # flux weighted centroid

    mc = np.unravel_index(np.argmax(img, axis=None), img.shape)
    xc, yc = mc[0], mc[1]
    # Geoff: I don't know why it thinks this is a problem. mesh_box returns 2 things
    xv, yv = mesh_box([xc, yc], 5)
    wx = np.sum(np.unique(xv) * img[yv, xv].sum(0)) / np.sum(img[yv, xv].sum(0))
    wy = np.sum(np.unique(yv) * img[yv, xv].sum(1)) / np.sum(img[yv, xv].sum(1))

    # loop through aper sizes
    apers = []
    bgs = []
    npps = []
    for r in np.arange(1.5, 4.1, 0.15):
        area, bg, npp = phot(img, wx, wy, r=r, dr=7)
        apers.append(area)
        bgs.append(bg)
        npps.append(npp)

    return wx, wy, apers, bgs, npps


def mesh_box(pos, box):
    '''mesh_box ds'''
    pos = [int(np.round(pos[0])), int(np.round(pos[1]))]
    x = np.arange(pos[0] - box, pos[0] + box + 1)
    y = np.arange(pos[1] - box, pos[1] + box + 1)
    xv, yv = np.meshgrid(x, y)
    return xv.astype(int), yv.astype(int)


def phot(data, xc, yc, r=5, dr=5):
    '''phot ds'''
    if dr > 0:
        bgflux = skybg_phot(data, xc, yc, r, dr)
    else:
        bgflux = 0
    positions = [(xc, yc)]
    data = data - bgflux
    data[data < 0] = 0

    apertures = CircularAperture(positions, r=r)
    phot_table = aperture_photometry(data, apertures, method='exact')
    subdata = (
        data[
            apertures.to_mask()[0]
            .bbox.iymin : apertures.to_mask()[0]
            .bbox.iymax,
            apertures.to_mask()[0]
            .bbox.ixmin : apertures.to_mask()[0]
            .bbox.ixmax,
        ]
        * apertures.to_mask()[0].data
    )
    npp = subdata.sum() ** 2 / np.sum(subdata**2)
    return float(phot_table['aperture_sum']), bgflux, npp


def skybg_phot(data, xc, yc, r=10, dr=5):
    '''skybg_phot ds'''
    # create a crude annulus to mask out bright background pixels
    xv, yv = mesh_box([xc, yc], np.round(r + dr))
    rv = ((xv - xc) ** 2 + (yv - yc) ** 2) ** 0.5
    mask = (rv > r) & (rv < (r + dr))
    cutoff = np.percentile(data[yv, xv][mask], 50)
    dat = np.copy(data)
    dat[dat > cutoff] = cutoff  # ignore bright pixels like stars
    return min(np.mean(dat[yv, xv][mask]), np.median(dat[yv, xv][mask]))


def linfit(x, y):
    '''linfit ds'''
    A = np.vstack([np.ones(len(x)), x]).T
    return np.linalg.lstsq(A, y, rcond=None)[0]  # b, m


def jwstcal_NIRISS(fin, clc, tim, tid, flttype, out, verbose=False):
    '''
    K. PEARSON: JWST NIRISS spectral extraction
    '''
    calibrated = False
    dbs = os.path.join(dawgie.context.data_dbs, 'mast')

    data = {
        'TIME': [],
        'SPEC': [],
        'WAVE': [],
        'LOC': [],
        'RAMP': [],
        'RAMP_OFFSET': [],
        'RAMP_NUM': [],
        'FAILED': [],
    }

    nframe = None
    dtframe = None
    ngroups = None

    for loc in sorted(clc['LOC']):
        fullloc = os.path.join(dbs, loc)
        with pyfits.open(fullloc) as hdulist:
            header0 = hdulist[0].header
            ftime = []
            exptime = []
            for fits in hdulist:
                if "primary" in fits.name.lower():
                    # keywords for NIRISS
                    start = fits.header.get("TIME-BJD")
                    ngroup = fits.header.get("ngroups", 1)
                    dtgroup = fits.header.get("tgroup") / (24 * 60 * 60.0)
                    nframe = fits.header.get("nframes", 1)
                    dtframe = fits.header.get("tframes")
                    data['TIME'].extend(
                        [start + i * dtgroup for i in range(ngroup)]
                    )
                    data['LOC'].append(loc)
                elif "sci" in fits.name.lower():
                    ngroup, nramp, height, width = fits.shape
                    slope = np.zeros((height, width))
                    offset = np.zeros((height, width))
                    # Future: optimize
                    for x in range(width):
                        for y in range(height):
                            pix = fits.data[:, :, y, x].flatten()
                            t = np.array(
                                [
                                    np.arange(1, nramp + 1)
                                    * dtgroup
                                    * 24
                                    * 60
                                    * 60
                                ]
                            ).flatten()
                            # Future: outlier rejection from neighboring pixels
                            b, m = linfit(t, pix)
                            slope[y, x] = m
                            offset[y, x] = b

                    data['RAMP'].append(slope)
                    data['RAMP_OFFSET'].append(offset)

                    # first order extraction
                    sub = slope[:96]
                    # find data above the background
                    mask = sub > np.percentile(sub, 65)
                    mask = binary_closing(mask)
                    mask = binary_erosion(mask)
                    mask = binary_dilation(mask)
                    mask = binary_fill_holes(mask)
                    labels, ngroups = label(mask)
                    inds, counts = np.unique(labels, return_counts=True)
                    sinds = np.argsort(counts)[::-1]
                    submask = (
                        labels == inds[sinds[1]]
                    )  # second biggest is data, first is bg

                    # approx wavecal
                    # https://jwst-docs.stsci.edu/near-infrared-imager-and-slitless-spectrograph/niriss-instrumentation/niriss-gr700xd-grism
                    minx = 18
                    maxx = 1487
                    x = np.array([112, 328, 545, 765, 988, 1213, 1441])
                    w = np.array([2.7, 2.4, 2.1, 1.8, 1.5, 1.2, 0.9])
                    fwave = itp.interp1d(
                        width * (x - minx) / (maxx - minx),
                        w,
                        fill_value="extrapolate",
                    )

                    # assumes ngroup = 1
                    # remake image by subtracting 0 read
                    dcube = np.zeros((nramp, height, width))
                    for k in range(nramp):
                        dcube[k] = fits.data[0, k] - offset
                        data['SPEC'].append(np.sum(dcube[k][:96] * submask, 0))
                        data['WAVE'].append(fwave(np.arange(width)))
                        data['RAMP_NUM'].append(k + 1)

    # transfer data
    for k, v in data.items():
        out['data'][k] = v
    calibrated = not np.all(data['FAILED']) or len(data['FAILED']) == 0
    if calibrated:
        out['STATUS'].append(True)

    # GMR: Fake use of inputs to satisfy pylint
    # 2/12/25 Geoff: nframe,dtframe,ngroups initializations added above
    #   to satisfy used-before-assignment error.
    #  if this part here is removed, then those three lines above can be removed too
    if verbose:
        _ = fin
        _ = tim
        _ = tid
        _ = flttype
        _ = header0
        _ = ftime
        _ = exptime
        _ = nframe
        _ = dtframe
        _ = ngroups
        pass
    # DELETE above when real use
    return calibrated
