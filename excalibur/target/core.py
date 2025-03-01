'''target core ds'''

# Heritage code shame:
# pylint: disable=too-many-arguments,too-many-branches,too-many-lines,too-many-locals,too-many-nested-blocks,too-many-positional-arguments,too-many-statements

# -- IMPORTS -- ------------------------------------------------------
import os
import re
import json
import time
import shutil

import requests
import tempfile
import subprocess
import logging

import dawgie
import dawgie.db

import excalibur.runtime.binding as rtbind

import excalibur.target.edit as trgedit
import excalibur.target.mast_api_utils as masttool

import numpy as np
import astropy.io.fits as pyfits
import urllib.request as urlrequest

from importlib import import_module as fetch  # avoid cicular dependencies

log = logging.getLogger(__name__)


# ------------- ------------------------------------------------------
# -- URLTRICK -- -----------------------------------------------------
class UrlTrick:
    '''
    # GMR: with statement generates an attribute __enter__ error
    # without the with statement it doesnt pass CI checks
    # pulling out dirty tricks
    '''

    def __init__(self, thisurl):
        '''__init__ ds'''
        self.thisurl = thisurl
        self.req = None
        return

    def __enter__(self):
        '''__enter__ ds'''
        self.req = urlrequest.urlopen(self.thisurl)
        return self.req.read()

    def __exit__(self, exc_type, exc_value, traceback):
        '''__exit__ ds'''
        if self.req:
            self.req.close()
        return

    pass


# ----------------- --------------------------------------------------
# -- TAP QUERY --  ---------------------------------------------------
def tap_query(base_url, query):
    '''table access protocol query'''
    # build url
    uri_full = base_url
    for k in query:
        if k != "format":
            uri_full += f"{k} {query[k]} "
        pass
    uri_full = uri_full[:-1] + f"&format={query.get('format', 'csv')}"
    uri_full = uri_full.replace(' ', '+')
    response = None
    with UrlTrick(uri_full) as test:
        response = test.decode('utf-8')
    return response


# ----------------- --------------------------------------------------
# -- SCRAPE IDS -- ---------------------------------------------------
def scrapeids(ds: dawgie.Dataset, out, web, gen_ids=True):
    '''
    Parses table from ipac exoplanetarchive
    '''
    targets = trgedit.targetlist.__doc__
    targets = targets.split('\n')
    targets = [t.strip() for t in targets if t.replace(' ', '')]
    tn = os.environ.get('TARGET_NAME', None)
    if tn is not None and tn != '':
        found_target_list = None
        for target in targets:
            if tn == target.split(':')[0].strip():
                found_target_list = target
            pass
        if found_target_list is None:
            # this is an ERROR.  the selected target should be in the target list
            # exit('ERROR: are you sure about that target?  it is not in the list')
            mssg = f'Obsolete target / Error in target name: {tn}'
            raise dawgie.NoValidOutputDataError(mssg)
        targets = [found_target_list]
        pass
    for target in targets:
        parsedstr = target.split(':')
        parsedstr = [t.strip() for t in parsedstr]
        out['starID'][parsedstr[0]] = {
            'planets': [],
            'PID': [],
            'aliases': [],
            'observatory': [],
            'datatable': [],
        }
        if parsedstr[1]:
            aliaslist = parsedstr[1].split(',')
            aliaslist = [a.strip() for a in aliaslist if a.strip()]
            out['starID'][parsedstr[0]]['aliases'].extend(aliaslist)
            pass
        if gen_ids:
            # pylint: disable=protected-access # because dawgie requires it
            dawgie.db.connect(
                fetch('excalibur.target').algorithms.create(),
                ds._bot(),
                parsedstr[0],
            ).load()
            pass
        pass
    # new additions 6/5/24:
    #   lower/upper limit flags (particularly for eccentricity)
    # new additions 4/14/23:
    #   pl_orblper (omega is used as prior in the transit fitting)
    #   pl_trandep (transit depth    is not used, but load it in, in case we want it later)
    #   pl_insol   (insolation       is not used, but load it in, in case we want it later)
    #   pl_trandur (transit duration is not used, but load it in, in case we want it later)
    #   pl_ratdor  (a/R*             is not used, but load it in, in case we want it later)
    #   pl_ratror  (Rp/R*            is not used, but load it in, in case we want it later)
    #   sy_dist    (same for distance; not used, but potentially useful)
    #  actually these additional references are only in the PSCompPars table, not the PS table
    # NO  sy_hmag_reflink (H magnitude reference can differ from the star param reference)
    # NO  sy_age_reflink  (age reference can differ from the star param reference)
    cols = "hostname,pl_letter,rowupdate,st_refname,pl_refname,sy_pnum,pl_orbper,pl_orbpererr1,pl_orbpererr2,pl_orbperlim,pl_orbsmax,pl_orbsmaxerr1,pl_orbsmaxerr2,pl_orbsmaxlim,pl_orbeccen,pl_orbeccenerr1,pl_orbeccenerr2,pl_orbeccenlim,pl_orbincl,pl_orbinclerr1,pl_orbinclerr2,pl_orbincllim,pl_bmassj,pl_bmassjerr1,pl_bmassjerr2,pl_bmassjlim,pl_radj,pl_radjerr1,pl_radjerr2,pl_radjlim,pl_dens,pl_denserr1,pl_denserr2,pl_denslim,pl_eqt,pl_eqterr1,pl_eqterr2,pl_eqtlim,pl_tranmid,pl_tranmiderr1,pl_tranmiderr2,pl_tranmidlim,pl_imppar,pl_impparerr1,pl_impparerr2,pl_impparlim,st_teff,st_tefferr1,st_tefferr2,st_tefflim,st_mass,st_masserr1,st_masserr2,st_masslim,st_rad,st_raderr1,st_raderr2,st_radlim,st_lum,st_lumerr1,st_lumerr2,st_lumlim,st_logg,st_loggerr1,st_loggerr2,st_logglim,st_dens,st_denserr1,st_denserr2,st_denslim,st_met,st_meterr1,st_meterr2,st_metlim,sy_jmag,sy_jmagerr1,sy_jmagerr2,sy_hmag,sy_hmagerr1,sy_hmagerr2,sy_kmag,sy_kmagerr1,sy_kmagerr2,sy_tmag,sy_tmagerr1,sy_tmagerr2,st_age,st_ageerr1,st_ageerr2,st_agelim,pl_orblper,pl_orblpererr1,pl_orblpererr2,pl_orblperlim,pl_trandep,pl_trandeperr1,pl_trandeperr2,pl_trandeplim,pl_insol,pl_insolerr1,pl_insolerr2,pl_insollim,pl_trandur,pl_trandurerr1,pl_trandurerr2,pl_trandurlim,pl_ratdor,pl_ratdorerr1,pl_ratdorerr2,pl_ratdorlim,pl_ratror,pl_ratrorerr1,pl_ratrorerr2,pl_ratrorlim,sy_dist,sy_disterr1,sy_disterr2,st_spectype"
    uri_ipac_query = {
        "select": cols,
        "from": 'ps',
        "where": "tran_flag = 1 and default_flag = 1",
        "order by": "pl_name",
        "format": "csv",
    }
    # First Exoplanet Archive query gets default row for each target
    uri_ipac_query['where'] = 'tran_flag = 1 and default_flag = 1'
    response = tap_query(web, uri_ipac_query)
    for line in response.split('\n'):
        out['nexsciDefaults'].append(line)

    # Second Exoplanet Archive query gets all rows, not just the defaults
    uri_ipac_query['where'] = 'tran_flag = 1'
    response = tap_query(web, uri_ipac_query)
    for line in response.split('\n'):
        out['nexsciFulltable'].append(line)

    out['STATUS'].append(True)
    ds.update()
    return


# ---------------- ---------------------------------------------------
# -- CREATE FILTERS -- -----------------------------------------------
def createfltrs(out):
    '''
    Create filter name
    '''
    created = False
    filters = [str(fn) for fn in rtbind.filter_names.values()]
    out['activefilters']['TOTAL'] = len(filters)
    out['activefilters']['NAMES'] = filters
    if filters:
        for flt in filters:
            out['activefilters'][flt] = {'ROOTNAME': [], 'LOC': [], 'TOTAL': []}
            pass
        created = True
        out['STATUS'].append(True)
        pass
    return created


# -------------------- -----------------------------------------------
# -- AUTOFILL -- -----------------------------------------------------
def autofillversion():
    '''
    1.1.1: manual trigger to autofill
    1.1.2: manual trigger for new spitzer priors
    1.1.3: manual trigger to fix new priors
    1.1.4: add stellar age integration
    1.1.5: added Hmag, L*; fixed FEH* problems
    1.2.0: new Exoplanet Archive tables and TAP protocol
    1.2.1: non-default parameters included; only two Exoplanet Archive queries; RHO* filled in by R*,M*
    2.0.1: requests are now using MAST api
    2.0.2: multiple attempts for the mast query with delay
    '''
    return dawgie.VERSION(2, 0, 2)


def autofill(ident, thistarget, out, allowed_filters, searchrad=0.2, ntrymax=4):
    '''Queries MAST for available data and parses NEXSCI data into tables'''

    out['starID'][thistarget] = ident['starIDs']['starID'][thistarget]
    target_id = [thistarget]

    # AUTOFILL WITH MAST QUERY ---------------------------------------
    solved = False
    # Using new MAST API instead of handmade urlrequest
    # 1 - solves for ra, dec
    # Cannot rely on target_name since proposers gets too inventive
    request = {
        'service': 'Mast.Name.Lookup',
        'params': {'input': thistarget, 'format': 'json'},
    }

    # Attempting with delay if server gives us the middle finger until ntrymax attempts
    def waitabit(myseconds, maxtry=ntrymax):
        def mydecorator(whatever):
            def wrapper(*args, **kwargs):
                ntry = 0
                while ntry < maxtry:
                    answer = whatever(*args, **kwargs)
                    if answer:
                        return answer
                    log.warning(
                        '--< Trying again, %s/%s with %s s delay >--',
                        ntry + 1,
                        maxtry,
                        myseconds * (10**ntry),
                    )
                    time.sleep(myseconds * (10**ntry))
                    ntry += 1
                    pass
                return ''

            return wrapper

        return mydecorator

    @waitabit(np.random.randint(low=1, high=10))
    def mastpoke(request, maxwaittime=42):
        hr, outstr = masttool.mast_query(request, maxwaittime=maxwaittime)
        if not outstr:
            log.warning('--< TARGET AUTOFILL: Empty MAST answer: %s >--', hr)
            pass
        return outstr

    outstr = mastpoke(request)
    # outstr = False  # temporarily skip the MAST stuff; just pull from Exoplanet Archive

    # assigning junk value to ra/dec avoids used-before-assignment error
    # but it would be better to just put the whole 'solved' part inside this if block
    ra = 123
    dec = 45

    if outstr:
        outjson = json.loads(outstr)
        if outjson['resolvedCoordinate']:
            ra = outjson['resolvedCoordinate'][0]['ra']
            dec = outjson['resolvedCoordinate'][0]['decl']
            solved = True
            pass
        pass

    # 2 - searches for observations
    # Cant say I dont comment my code anymore
    if solved:
        request = {
            'service': 'Mast.Caom.Cone',
            'params': {'ra': ra, 'dec': dec, 'radius': searchrad},
            'format': 'json',
            'pagesize': 250000,  # this seems to help HD 209458 (~221K obs)
            'removenullcolumns': True,
            'removecache': True,
        }
        targettable = []
        platformlist = []
        # GMR: mastpoke returns outstr only
        # errors (formerly _h) are logged and handled inside mastpoke
        # 42 seconds isn't enough time (for HAT-P-11 testrun)
        outstr = mastpoke(request, maxwaittime=420)
        if outstr:
            outjson = json.loads(
                outstr
            )  # this line fails without maskpoke loop

            # 3 - activefilters filtering on TELESCOPE INSTRUMENT FILTER
            if 'data' not in outjson:
                outjson['data'] = []  # special case for TOI-1338
            log.warning(
                '--< TARGET AUTOFILL: %s  #-of-obs %s >--',
                thistarget,
                len(outjson['data']),
            )
            for obs in outjson['data']:
                for f in allowed_filters:
                    if obs['obs_collection'] is not None:
                        pfcond = (
                            f.split('-')[0].upper() in obs['obs_collection']
                        )
                    else:
                        pfcond = False
                    if obs['instrument_name'] is not None:
                        nscond = f.split('-')[1] in obs['instrument_name']
                    else:
                        nscond = False
                    if (
                        f.split('-')[1] == 'sim'
                    ):  # skip any simulated instruments
                        flcond = False
                    elif (f.split('-')[0] not in ['Spitzer']) and (
                        obs['filters'] is not None
                    ):
                        flcond = f.split('-')[3] in obs['filters']
                    else:
                        sptzfl = None
                        if f.split('-')[3] in ['36']:
                            sptzfl = 'IRAC1'
                        if f.split('-')[3] in ['45']:
                            sptzfl = 'IRAC2'
                        if obs['filters'] is not None:
                            flcond = sptzfl in obs['filters']
                        else:
                            flcond = False
                    if pfcond and nscond and flcond:
                        targettable.append(obs)
                        platformlist.append(obs['obs_collection'])
                    pass
                pass
            log.warning(
                '--< TARGET AUTOFILL: %s  in targettable %s >--',
                thistarget,
                len(targettable),
            )
            pass
        else:
            log.warning(
                '--< TARGET AUTOFILL: MAST FAIL (despite mastpoke) %s >--',
                thistarget,
            )
            pass
        if not targettable:
            solved = False
            pass
        # Note: If we use data that are not known by MAST but are on disk
        # this will also return False

        # 4 - pid filtering
        pidlist = []
        aliaslist = []
        obslist = []
        for item, pf in zip(targettable, platformlist):
            if item['proposal_id'] not in pidlist:
                pidlist.append(item['proposal_id'])
                aliaslist.append(item['target_name'])
                obslist.append(pf)
                pass
            pass
        out['starID'][thistarget]['aliases'].extend(aliaslist)
        out['starID'][thistarget]['PID'].extend(pidlist)
        out['starID'][thistarget]['observatory'].extend(obslist)
        out['starID'][thistarget]['datatable'].extend(targettable)
        pass
    out['STATUS'].append(solved)

    # AUTOFILL WITH NEXSCI EXOPLANET TABLE, DEFAULTS ONLY ---------------------------
    merged = False
    target_id.extend(ident['starIDs']['starID'][thistarget]['aliases'])
    response = ident['starIDs']['nexsciDefaults']
    header = response[0].split(',')
    # list all keys contained in csv
    # using tuples of form (a, b) where a is the key name
    # and then b are the other types of columns for that key
    # that are included
    # GMR: They are continuously changing the format: creating translatekeys()
    matchlist = translatekeys(header)

    banlist = ['star', 'planet', 'update', 'ref', 'ref_st', 'ref_pl', 'np']
    plist = [
        'period',
        'period_uperr',
        'period_lowerr',
        'period_ref',
        'period_lim',
        'sma',
        'sma_uperr',
        'sma_lowerr',
        'sma_ref',
        'sma_lim',
        'ecc',
        'ecc_uperr',
        'ecc_lowerr',
        'ecc_ref',
        'ecc_lim',
        'inc',
        'inc_uperr',
        'inc_lowerr',
        'inc_ref',
        'inc_lim',
        'mass',
        'mass_uperr',
        'mass_lowerr',
        'mass_ref',
        'mass_lim',
        'rp',
        'rp_uperr',
        'rp_lowerr',
        'rp_ref',
        'rp_lim',
        'rho',
        'rho_uperr',
        'rho_lowerr',
        'rho_ref',
        'rho_lim',
        'teq',
        'teq_uperr',
        'teq_lowerr',
        'teq_ref',
        'teq_lim',
        't0',
        't0_uperr',
        't0_lowerr',
        't0_ref',
        't0_lim',
        'impact',
        'impact_uperr',
        'impact_lowerr',
        'impact_ref',
        'impact_lim',
        'omega',
        'omega_uperr',
        'omega_lowerr',
        'omega_ref',
        'omega_lim',
        'trandepth',
        'trandepth_uperr',
        'trandepth_lowerr',
        'trandepth_ref',
        'trandepth_lim',
        'trandur',
        'trandur_uperr',
        'trandur_lowerr',
        'trandur_ref',
        'trandur_lim',
        'insol',
        'insol_uperr',
        'insol_lowerr',
        'insol_ref',
        'insol_lim',
        'ars',
        'ars_uperr',
        'ars_lowerr',
        'ars_ref',
        'ars_lim',
        'rprs',
        'rprs_uperr',
        'rprs_lowerr',
        'rprs_ref',
        'rprs_lim',
    ]
    for line in response:
        elem = line.split(',')
        elem = clean_elems(
            elem
        )  # remove things such as bounding quotation marks
        trgt = elem[header.index('hostname')]
        if trgt in target_id:
            out['starID'][thistarget]['PLANETS #'] = elem[
                header.index('sy_pnum')
            ]
            thisplanet = elem[header.index('pl_letter')]
            if thisplanet not in out['starID'][thistarget].keys():
                # GMR: Init planet dict
                out['starID'][thistarget][thisplanet] = {}
                if (
                    'CREATED'
                    not in out['starID'][thistarget][thisplanet].keys()
                ):
                    out['starID'][thistarget][thisplanet]['CREATED'] = True
                    # GMR: Inits with empty strings
                    null_keys = ['date', 'date_ref']
                    for key in null_keys:
                        out['starID'][thistarget][thisplanet][key] = ['']
                        pass
                    # GMR: Inits with empty lists
                    blank_keys = [
                        'period_units',
                        'period_ref',
                        'period_lim',
                        'sma_units',
                        'sma_ref',
                        'sma_lim',
                        'ecc_units',
                        'ecc_ref',
                        'ecc_lim',
                        'inc_units',
                        'inc_ref',
                        'inc_lim',
                        'mass_units',
                        'mass_ref',
                        'mass_lim',
                        'rp_units',
                        'rp_ref',
                        'rp_lim',
                        'rho_units',
                        'rho_ref',
                        'rho_lim',
                        'teq_units',
                        'teq_ref',
                        'teq_lim',
                        't0_units',
                        't0_ref',
                        't0_lim',
                        'impact_units',
                        'impact_ref',
                        'impact_lim',
                        'omega_units',
                        'omega_ref',
                        'omega_lim',
                        'trandepth_units',
                        'trandepth_ref',
                        'trandepth_lim',
                        'trandur_units',
                        'trandur_ref',
                        'trandur_lim',
                        'insol_units',
                        'insol_ref',
                        'insol_lim',
                        'ars_units',
                        'ars_ref',
                        'ars_lim',
                        'rprs_units',
                        'rprs_ref',
                        'rprs_lim',
                    ]
                    for key in blank_keys:
                        out['starID'][thistarget][thisplanet][key] = []
                        pass
                    pass
                # GMR: Init stellar dict
                if 'CREATED' not in out['starID'][thistarget].keys():
                    out['starID'][thistarget]['CREATED'] = True
                    out['starID'][thistarget]['date'] = ['']
                    out['starID'][thistarget]['date_ref'] = ['']
                    out['starID'][thistarget]['R*_units'] = []
                    out['starID'][thistarget]['R*_ref'] = []
                    out['starID'][thistarget]['M*_units'] = []
                    out['starID'][thistarget]['M*_ref'] = []
                    out['starID'][thistarget]['T*_units'] = []
                    out['starID'][thistarget]['T*_ref'] = []
                    out['starID'][thistarget]['L*_units'] = []
                    out['starID'][thistarget]['L*_ref'] = []
                    out['starID'][thistarget]['LOGG*_units'] = []
                    out['starID'][thistarget]['LOGG*_ref'] = []
                    out['starID'][thistarget]['RHO*_units'] = []
                    out['starID'][thistarget]['RHO*_ref'] = []
                    out['starID'][thistarget]['FEH*_ref'] = []
                    out['starID'][thistarget]['FEH*_units'] = ['[dex]']
                    out['starID'][thistarget]['AGE*'] = []
                    out['starID'][thistarget]['AGE*_uperr'] = []
                    out['starID'][thistarget]['AGE*_lowerr'] = []
                    out['starID'][thistarget]['AGE*_units'] = []
                    out['starID'][thistarget]['AGE*_ref'] = []
                    out['starID'][thistarget]['Jmag_units'] = []
                    out['starID'][thistarget]['Jmag_ref'] = []
                    out['starID'][thistarget]['Hmag_units'] = []
                    out['starID'][thistarget]['Hmag_ref'] = []
                    out['starID'][thistarget]['Kmag_units'] = []
                    out['starID'][thistarget]['Kmag_ref'] = []
                    out['starID'][thistarget]['TESSmag_units'] = []
                    out['starID'][thistarget]['TESSmag_ref'] = []
                    out['starID'][thistarget]['dist_units'] = []
                    out['starID'][thistarget]['dist_ref'] = []
                    out['starID'][thistarget]['spTyp'] = []
                    out['starID'][thistarget]['spTyp_uperr'] = []
                    out['starID'][thistarget]['spTyp_lowerr'] = []
                    out['starID'][thistarget]['spTyp_units'] = []
                    out['starID'][thistarget]['spTyp_ref'] = []
                ref_pl = elem[header.index('pl_refname')]
                ref_pl = ref_pl.split('</a>')[0]
                ref_pl = ref_pl.split('target=ref>')[-1]
                ref_pl = ref_pl.strip()
                ref_st = elem[header.index('st_refname')]
                ref_st = ref_st.split('</a>')[0]
                ref_st = ref_st.split('target=ref>')[-1]
                ref_st = ref_st.strip()
                pass
            for idx, keymatch in enumerate(matchlist):
                if keymatch not in banlist:
                    # GMR: Planet dict
                    if keymatch in plist:
                        out_ref = out['starID'][thistarget][thisplanet]
                    # GMR: Stellar dict
                    else:
                        out_ref = out['starID'][thistarget]
                    # Init
                    if keymatch not in out_ref:
                        out_ref[keymatch] = []
                    out_ref[keymatch].append(elem[idx])
                    pass
                pass
            # GMR: Adding units to the planet dict
            out['starID'][thistarget][thisplanet]['period_units'].append(
                '[days]'
            )
            out['starID'][thistarget][thisplanet]['sma_units'].append('[AU]')
            out['starID'][thistarget][thisplanet]['ecc_units'].append('[]')
            out['starID'][thistarget][thisplanet]['inc_units'].append(
                '[degree]'
            )
            out['starID'][thistarget][thisplanet]['mass_units'].append(
                '[Jupiter mass]'
            )
            out['starID'][thistarget][thisplanet]['rp_units'].append(
                '[Jupiter radius]'
            )
            out['starID'][thistarget][thisplanet]['rho_units'].append(
                '[g.cm-3]'
            )
            out['starID'][thistarget][thisplanet]['teq_units'].append('[K]')
            out['starID'][thistarget][thisplanet]['t0_units'].append(
                '[Julian Days]'
            )
            out['starID'][thistarget][thisplanet]['impact_units'].append('[R*]')
            out['starID'][thistarget][thisplanet]['omega_units'].append(
                '[degree]'
            )
            out['starID'][thistarget][thisplanet]['trandepth_units'].append(
                '[%]'
            )
            out['starID'][thistarget][thisplanet]['trandur_units'].append(
                '[hour]'
            )
            out['starID'][thistarget][thisplanet]['insol_units'].append(
                '[Earth flux]'
            )
            out['starID'][thistarget][thisplanet]['ars_units'].append('[]')
            out['starID'][thistarget][thisplanet]['rprs_units'].append('[]')
            # GMR: Adding refs to the planet dict
            pkeys = [pk for pk in plist if '_' not in pk]
            for pk in pkeys:
                if out['starID'][thistarget][thisplanet][pk]:
                    addme = ref_pl
                else:
                    addme = ''
                out['starID'][thistarget][thisplanet][pk + '_ref'].append(addme)
                pass
            # GMR: Adding refs to the stellar dict
            skeys = [
                sk
                for sk in matchlist
                if (sk not in banlist) and (sk not in plist) and ('_' not in sk)
            ]
            for sk in skeys:
                if out['starID'][thistarget][sk]:
                    addme = ref_st
                else:
                    addme = ''
                # if not out['starID'][thistarget][sk+'_ref']:
                out['starID'][thistarget][sk + '_ref'].append(addme)
            # GMR: Adding units to the stellar dict
            out['starID'][thistarget]['R*_units'].append('[Rsun]')
            out['starID'][thistarget]['M*_units'].append('[Msun]')
            out['starID'][thistarget]['T*_units'].append('[K]')
            out['starID'][thistarget]['L*_units'].append('[Lsun]')
            out['starID'][thistarget]['LOGG*_units'].append('log10[cm.s-2]')
            out['starID'][thistarget]['RHO*_units'].append('[g.cm-3]')
            out['starID'][thistarget]['AGE*_units'].append('[Gyr]')
            out['starID'][thistarget]['Jmag_units'].append('[mag]')
            out['starID'][thistarget]['Hmag_units'].append('[mag]')
            out['starID'][thistarget]['Kmag_units'].append('[mag]')
            out['starID'][thistarget]['TESSmag_units'].append('[mag]')
            out['starID'][thistarget]['dist_units'].append('[pc]')
            out['starID'][thistarget]['spTyp_units'].append('')
            # spectral type doesn't have uncertainties. fill in here by hand
            out['starID'][thistarget]['spTyp_uperr'].append('')
            out['starID'][thistarget]['spTyp_lowerr'].append('')
            # out['starID'][thistarget]['spTyp_lim'].append('')
            merged = True
            pass
        pass

    # AUTOFILL WITH NEXSCI FULL TABLE, INCLUDING NON-DEFAULTS ----------------------------
    response = ident['starIDs']['nexsciFulltable']
    header = response[0].split(',')
    # print('header',header)
    matchkey = translatekeys(header)
    # print('matchkey',matchkey)
    for line in response:
        elem = line.split(',')
        elem = clean_elems(elem)
        if len(elem) != len(header):
            trgt = 'YOLO'
            thisplanet = 'CARPEDIEM'
            pass
        else:
            trgt = elem[header.index('hostname')]
            thisplanet = elem[header.index('pl_letter')]
            pass
        if trgt in target_id:
            ref_pl = elem[header.index('pl_refname')]
            ref_pl = ref_pl.split('</a>')[0]
            ref_pl = ref_pl.split('target=ref>')[-1]
            ref_pl = ref_pl.strip()
            ref_st = elem[header.index('st_refname')]
            ref_st = ref_st.split('</a>')[0]
            ref_st = ref_st.split('target=ref>')[-1]
            ref_st = ref_st.strip()
            numdate = elem[header.index('rowupdate')]
            strnum = ''
            for n in re.split('-|:| ', numdate):
                strnum = strnum + n.strip()
            numdate = float(strnum)
            for addme, key in zip(elem, matchkey):
                if key not in banlist:
                    refkey = key.split('_')[0]
                    if key in plist:
                        out['starID'][thistarget][thisplanet][key].append(addme)
                        if refkey == key:
                            out['starID'][thistarget][thisplanet][
                                refkey + '_ref'
                            ].append(ref_pl)
                        pass
                    else:
                        refkey = key.split('_')[0]
                        out['starID'][thistarget][key].append(addme)
                        if refkey == key:
                            out['starID'][thistarget][refkey + '_ref'].append(
                                ref_st
                            )
                            pass
                        # again, fill in spectral type have uncertainties by hand
                        if key == 'spTyp':
                            out['starID'][thistarget][key + '_uperr'].append('')
                            out['starID'][thistarget][key + '_lowerr'].append(
                                ''
                            )
                            # out['starID'][thistarget][key+'_lim'].append('')
                            pass
                        pass
                    pass
                pass
            merged = True
            pass
        pass

    # print('final out',out['starID'][thistarget])
    # print('final out',out['starID'][thistarget].keys())

    # FINALIZE OUTPUT ------------------------------------------------
    if merged:
        candidates = [
            'b',
            'c',
            'd',
            'e',
            'f',
            'g',
            'h',
            'i',
            'j',
            'k',
            'l',
            'm',
            'n',
            'o',
            'p',
            'q',
            'r',
            's',
            't',
            'u',
            'v',
            'w',
            'x',
            'y',
            'z',
        ]
        plist = []
        for letter in candidates:
            if letter in out['starID'][thistarget].keys():
                plist.append(letter)
                pass
            pass
        out['starID'][thistarget]['planets'] = plist
        out['candidates'].extend(candidates)
        out['starkeys'].extend(sorted(skeys))
        out['planetkeys'].extend(sorted(pkeys))
        out['exts'].extend(['_uperr', '_lowerr', '_units', '_ref'])
        out['STATUS'].append(True)
        pass
    # GMR: Some targets cannot be solved but we still need them for sims.
    # Changing this condition to OR
    status = solved or merged
    return status


def clean_elem(elem):
    '''remove formatting on a single element from TAP API'''
    if not elem:
        return elem
    if elem[0] == '"' and elem[-1] == '"':
        return elem[1:-1]
    return elem


def clean_elems(elems):
    '''removes formatting symbols and prepares values for saving'''
    return [clean_elem(elem) for elem in elems]


def translatekeys(header):
    '''GMR: Translate NEXSCI keywords into Excalibur ones'''
    matchlist = []
    for key in header:
        thiskey = key
        if 'pl_hostname' == thiskey:
            xclbrkey = 'star'
        elif 'hostname' == thiskey:
            xclbrkey = 'star'
        elif 'pl_letter' == thiskey:
            xclbrkey = 'planet'
        elif 'pl_reflink' == thiskey:
            xclbrkey = 'ref_pl'
        elif 'rowupdate' == thiskey:
            xclbrkey = 'update'
        elif 'pl_def_refname' == thiskey:
            xclbrkey = 'ref_pl'
        elif 'pl_refname' == thiskey:
            xclbrkey = 'ref_pl'
        elif 'pl_pnum' == thiskey:
            xclbrkey = 'np'
        elif 'sy_pnum' == thiskey:
            xclbrkey = 'np'
        elif 'pl_orbper' == thiskey:
            xclbrkey = 'period'
        elif 'pl_orbpererr1' == thiskey:
            xclbrkey = 'period_uperr'
        elif 'pl_orbpererr2' == thiskey:
            xclbrkey = 'period_lowerr'
        elif 'pl_orbperlim' == thiskey:
            xclbrkey = 'period_lim'
        elif 'pl_orbperreflink' == thiskey:
            xclbrkey = 'period_ref'
        elif 'pl_orbsmax' == thiskey:
            xclbrkey = 'sma'
        elif 'pl_smax' == thiskey:
            xclbrkey = 'sma'
        elif 'pl_orbsmaxerr1' == thiskey:
            xclbrkey = 'sma_uperr'
        elif 'pl_smaxerr1' == thiskey:
            xclbrkey = 'sma_uperr'
        elif 'pl_orbsmaxerr2' == thiskey:
            xclbrkey = 'sma_lowerr'
        elif 'pl_orbsmaxlim' == thiskey:
            xclbrkey = 'sma_lim'
        elif 'pl_smaxerr2' == thiskey:
            xclbrkey = 'sma_lowerr'
        elif 'pl_smaxlim' == thiskey:
            xclbrkey = 'sma_lim'
        elif 'pl_smaxreflink' == thiskey:
            xclbrkey = 'sma_ref'
        elif 'pl_orbeccen' == thiskey:
            xclbrkey = 'ecc'
        elif 'pl_eccen' == thiskey:
            xclbrkey = 'ecc'
        elif 'pl_orbeccenerr1' == thiskey:
            xclbrkey = 'ecc_uperr'
        elif 'pl_eccenerr1' == thiskey:
            xclbrkey = 'ecc_uperr'
        elif 'pl_orbeccenerr2' == thiskey:
            xclbrkey = 'ecc_lowerr'
        elif 'pl_orbeccenlim' == thiskey:
            xclbrkey = 'ecc_lim'
        elif 'pl_eccenerr2' == thiskey:
            xclbrkey = 'ecc_lowerr'
        elif 'pl_eccenlim' == thiskey:
            xclbrkey = 'ecc_lim'
        elif 'pl_eccenreflink' == thiskey:
            xclbrkey = 'ecc_ref'
        elif 'pl_orbincl' == thiskey:
            xclbrkey = 'inc'
        elif 'pl_orbinclerr1' == thiskey:
            xclbrkey = 'inc_uperr'
        elif 'pl_orbinclerr2' == thiskey:
            xclbrkey = 'inc_lowerr'
        elif 'pl_orbincllim' == thiskey:
            xclbrkey = 'inc_lim'
        elif 'pl_bmassj' == thiskey:
            xclbrkey = 'mass'
        elif 'pl_bmassjerr1' == thiskey:
            xclbrkey = 'mass_uperr'
        elif 'pl_bmassjerr2' == thiskey:
            xclbrkey = 'mass_lowerr'
        elif 'pl_bmassjlim' == thiskey:
            xclbrkey = 'mass_lim'
        elif 'pl_radj' == thiskey:
            xclbrkey = 'rp'
        elif 'pl_radjerr1' == thiskey:
            xclbrkey = 'rp_uperr'
        elif 'pl_radjerr2' == thiskey:
            xclbrkey = 'rp_lowerr'
        elif 'pl_radjlim' == thiskey:
            xclbrkey = 'rp_lim'
        elif 'pl_radreflink' == thiskey:
            xclbrkey = 'rp_ref'
        elif 'pl_dens' == thiskey:
            xclbrkey = 'rho'
        elif 'pl_denserr1' == thiskey:
            xclbrkey = 'rho_uperr'
        elif 'pl_denserr2' == thiskey:
            xclbrkey = 'rho_lowerr'
        elif 'pl_denslim' == thiskey:
            xclbrkey = 'rho_lim'
        elif 'pl_eqt' == thiskey:
            xclbrkey = 'teq'
        elif 'pl_eqterr1' == thiskey:
            xclbrkey = 'teq_uperr'
        elif 'pl_eqterr2' == thiskey:
            xclbrkey = 'teq_lowerr'
        elif 'pl_eqtlim' == thiskey:
            xclbrkey = 'teq_lim'
        elif 'pl_eqtreflink' == thiskey:
            xclbrkey = 'teq_ref'
        elif 'pl_tranmid' == thiskey:
            xclbrkey = 't0'
        elif 'pl_tranmiderr1' == thiskey:
            xclbrkey = 't0_uperr'
        elif 'pl_tranmiderr2' == thiskey:
            xclbrkey = 't0_lowerr'
        elif 'pl_tranmidlim' == thiskey:
            xclbrkey = 't0_lim'
        elif 'pl_imppar' == thiskey:
            xclbrkey = 'impact'
        elif 'pl_impparerr1' == thiskey:
            xclbrkey = 'impact_uperr'
        elif 'pl_impparerr2' == thiskey:
            xclbrkey = 'impact_lowerr'
        elif 'pl_impparlim' == thiskey:
            xclbrkey = 'impact_lim'
        elif 'pl_orblper' == thiskey:
            xclbrkey = 'omega'
        elif 'pl_orblpererr1' == thiskey:
            xclbrkey = 'omega_uperr'
        elif 'pl_orblpererr2' == thiskey:
            xclbrkey = 'omega_lowerr'
        elif 'pl_orblperlim' == thiskey:
            xclbrkey = 'omega_lim'
        elif 'pl_trandep' == thiskey:
            xclbrkey = 'trandepth'
        elif 'pl_trandeperr1' == thiskey:
            xclbrkey = 'trandepth_uperr'
        elif 'pl_trandeperr2' == thiskey:
            xclbrkey = 'trandepth_lowerr'
        elif 'pl_trandeplim' == thiskey:
            xclbrkey = 'trandepth_lim'
        elif 'pl_trandur' == thiskey:
            xclbrkey = 'trandur'
        elif 'pl_trandurerr1' == thiskey:
            xclbrkey = 'trandur_uperr'
        elif 'pl_trandurerr2' == thiskey:
            xclbrkey = 'trandur_lowerr'
        elif 'pl_trandurlim' == thiskey:
            xclbrkey = 'trandur_lim'
        elif 'pl_insol' == thiskey:
            xclbrkey = 'insol'
        elif 'pl_insolerr1' == thiskey:
            xclbrkey = 'insol_uperr'
        elif 'pl_insolerr2' == thiskey:
            xclbrkey = 'insol_lowerr'
        elif 'pl_insollim' == thiskey:
            xclbrkey = 'insol_lim'
        elif 'pl_ratdor' == thiskey:
            xclbrkey = 'ars'
        elif 'pl_ratdorerr1' == thiskey:
            xclbrkey = 'ars_uperr'
        elif 'pl_ratdorerr2' == thiskey:
            xclbrkey = 'ars_lowerr'
        elif 'pl_ratdorlim' == thiskey:
            xclbrkey = 'ars_lim'
        elif 'pl_ratror' == thiskey:
            xclbrkey = 'rprs'
        elif 'pl_ratrorerr1' == thiskey:
            xclbrkey = 'rprs_uperr'
        elif 'pl_ratrorerr2' == thiskey:
            xclbrkey = 'rprs_lowerr'
        elif 'pl_ratrorlim' == thiskey:
            xclbrkey = 'rprs_lim'
        elif 'st_teff' == thiskey:
            xclbrkey = 'T*'
        elif 'st_tefferr1' == thiskey:
            xclbrkey = 'T*_uperr'
        elif 'st_tefferr2' == thiskey:
            xclbrkey = 'T*_lowerr'
        elif 'st_tefflim' == thiskey:
            xclbrkey = 'T*_lim'
        elif 'st_teffreflink' == thiskey:
            xclbrkey = 'T*_ref'
        elif 'st_refname' == thiskey:
            xclbrkey = 'ref_st'
        elif 'st_mass' == thiskey:
            xclbrkey = 'M*'
        elif 'st_masserr1' == thiskey:
            xclbrkey = 'M*_uperr'
        elif 'st_masserr2' == thiskey:
            xclbrkey = 'M*_lowerr'
        elif 'st_masslim' == thiskey:
            xclbrkey = 'M*_lim'
        elif 'st_rad' == thiskey:
            xclbrkey = 'R*'
        elif 'st_raderr1' == thiskey:
            xclbrkey = 'R*_uperr'
        elif 'st_raderr2' == thiskey:
            xclbrkey = 'R*_lowerr'
        elif 'st_radlim' == thiskey:
            xclbrkey = 'R*_lim'
        elif 'st_radreflink' == thiskey:
            xclbrkey = 'R*_ref'
        elif 'st_lum' == thiskey:
            xclbrkey = 'L*'
        elif 'st_lumerr1' == thiskey:
            xclbrkey = 'L*_uperr'
        elif 'st_lumerr2' == thiskey:
            xclbrkey = 'L*_lowerr'
        elif 'st_lumlim' == thiskey:
            xclbrkey = 'L*_lim'
        elif 'st_logg' == thiskey:
            xclbrkey = 'LOGG*'
        elif 'st_loggerr1' == thiskey:
            xclbrkey = 'LOGG*_uperr'
        elif 'st_loggerr2' == thiskey:
            xclbrkey = 'LOGG*_lowerr'
        elif 'st_logglim' == thiskey:
            xclbrkey = 'LOGG*_lim'
        elif 'st_loggreflink' == thiskey:
            xclbrkey = 'LOGG*_ref'
        elif 'st_dens' == thiskey:
            xclbrkey = 'RHO*'
        elif 'st_denserr1' == thiskey:
            xclbrkey = 'RHO*_uperr'
        elif 'st_denserr2' == thiskey:
            xclbrkey = 'RHO*_lowerr'
        elif 'st_denslim' == thiskey:
            xclbrkey = 'RHO*_lim'
        elif 'st_metfe' == thiskey:
            xclbrkey = 'FEH*'
        elif 'st_met' == thiskey:
            xclbrkey = 'FEH*'
        elif 'st_metfeerr1' == thiskey:
            xclbrkey = 'FEH*_uperr'
        elif 'st_meterr1' == thiskey:
            xclbrkey = 'FEH*_uperr'
        elif 'st_metfeerr2' == thiskey:
            xclbrkey = 'FEH*_lowerr'
        elif 'st_meterr2' == thiskey:
            xclbrkey = 'FEH*_lowerr'
        elif 'st_metfelim' == thiskey:
            xclbrkey = 'FEH*_lim'
        elif 'st_metlim' == thiskey:
            xclbrkey = 'FEH*_lim'
        elif 'st_metratio' == thiskey:
            xclbrkey = 'FEH*_units'
        elif 'st_metreflink' == thiskey:
            xclbrkey = 'FEH*_ref'
        elif 'sy_jmag' == thiskey:
            xclbrkey = 'Jmag'
        elif 'sy_jmagerr1' == thiskey:
            xclbrkey = 'Jmag_uperr'
        elif 'sy_jmagerr2' == thiskey:
            xclbrkey = 'Jmag_lowerr'
        # elif 'sy_jmaglim' == thiskey: xclbrkey = 'Jmag_lim'
        elif 'sy_hmag' == thiskey:
            xclbrkey = 'Hmag'
        elif 'sy_hmagerr1' == thiskey:
            xclbrkey = 'Hmag_uperr'
        elif 'sy_hmagerr2' == thiskey:
            xclbrkey = 'Hmag_lowerr'
        # elif 'sy_hmaglim' == thiskey: xclbrkey = 'Hmag_lim'
        elif 'sy_kmag' == thiskey:
            xclbrkey = 'Kmag'
        elif 'sy_kmagerr1' == thiskey:
            xclbrkey = 'Kmag_uperr'
        elif 'sy_kmagerr2' == thiskey:
            xclbrkey = 'Kmag_lowerr'
        # elif 'sy_kmaglim' == thiskey: xclbrkey = 'Kmag_lim'
        elif 'sy_tmag' == thiskey:
            xclbrkey = 'TESSmag'
        elif 'sy_tmagerr1' == thiskey:
            xclbrkey = 'TESSmag_uperr'
        elif 'sy_tmagerr2' == thiskey:
            xclbrkey = 'TESSmag_lowerr'
        elif 'st_age' == thiskey:
            xclbrkey = 'AGE*'
        elif 'st_ageerr1' == thiskey:
            xclbrkey = 'AGE*_uperr'
        elif 'st_ageerr2' == thiskey:
            xclbrkey = 'AGE*_lowerr'
        elif 'st_agelim' == thiskey:
            xclbrkey = 'AGE*_lim'
        elif 'sy_dist' == thiskey:
            xclbrkey = 'dist'
        elif 'sy_disterr1' == thiskey:
            xclbrkey = 'dist_uperr'
        elif 'sy_disterr2' == thiskey:
            xclbrkey = 'dist_lowerr'
        # elif 'sy_distlim' == thiskey: xclbrkey = 'dist_lim'
        elif 'st_spectype' == thiskey:
            xclbrkey = 'spTyp'
        else:
            # print('HEY NO MATCH FOR THIS ARCHIVE KEYWORD!',thiskey)
            xclbrkey = None
        if xclbrkey is not None:
            matchlist.append(xclbrkey)
        pass
    if len(matchlist) != len(header):
        errstr = 'MISSING NEXSCI KEY MAPPING'
        log.warning('--< TARGET AUTOFILL: %s >--', errstr)
        pass
    return matchlist


# -------------- -----------------------------------------------------
# -- SCRAPE -- -------------------------------------------------------
def scrapeversion():
    '''
    2.0.0: GMR: Code changes to use the new MAST API
    2.0.1: GMR: See #539 bugfixes after PR #546
    '''
    return dawgie.VERSION(2, 0, 1)


# ------------ -------------------------------------------------------
# -- MAST -- ---------------------------------------------------------
def mastapi(tfl, out, dbs, download_url=None, hst_url=None, verbose=False):
    '''Uses MAST API tools to download data'''

    target = list(tfl['starID'].keys())[0]
    obstable = tfl['starID'][target]['datatable']
    obsids = [o['obsid'] for o in obstable]
    obsids = set(obsids)
    # obsids are generally repeated many times. this will save lots of time
    allsci = []
    allurl = []
    allmiss = []
    allraw = []
    for o in obsids:
        donmast = False
        request = {
            'service': 'Mast.Caom.Products',
            'params': {'obsid': o},
            'format': 'json',
        }
        # print('MAST download can take a while. use larger wait time')
        errmastq, datastr = masttool.mast_query(request, maxwaittime=1000)
        data = json.loads(datastr)
        if data['data']:
            donmast = True
        dtlvl = None
        clblvl = None
        obscol = ''
        thisurl = ''
        if donmast:
            if 'SPITZER' in data['data'][0].get("obs_collection", None).upper():
                # Relying on disk data for now.
                obscol = 'SPITZER'
                dtlvl = 'SEE WITH KYLE'
                clblvl = 666
                thisurl = 'https:do.not.dl.for.now'
                pass
            if 'JWST' in data['data'][0].get("obs_collection", None):
                obscol = 'JWST'
                dtlvl = 'CALINTS'
                clblvl = 2
                thisurl = download_url
                pass
            if 'HST' in data['data'][0].get("obs_collection", None):
                obscol = 'HST'
                # No comment on mast api missing ima files. See HST hack below.
                dtlvl = 'FLT'
                clblvl = 2
                thisurl = hst_url
                pass
            scidata = [
                x
                for x in data['data']
                if (x.get("productType", None) == 'SCIENCE')
                and (x.get("dataRights", None) == 'PUBLIC')
                and (x.get("calib_level", None) == clblvl)
                and (x.get("productSubGroupDescription", None) == dtlvl)
            ]
            # Associated L1b JWST data (UNCAL)
            if 'JWST' in obscol and scidata:
                sciids = [x.get('obsID', None) for x in scidata]
                rawdata = [
                    x
                    for x in data['data']
                    if (x.get("productType", None) == 'SCIENCE')
                    and (x.get("dataRights", None) == 'PUBLIC')
                    and (x.get("calib_level", None) == 1)
                    and (x.get("productSubGroupDescription", None) == 'UNCAL')
                    and (x.get('obsID') in sciids)
                ]
                allraw.extend(rawdata)
                # --<
                # Downloads JWST data only
                allsci.extend(scidata)
                allmiss.extend([obscol] * len(scidata))
                allurl.extend([thisurl] * len(scidata))
                pass
            # Downloads all missions
            # allsci.extend(scidata)
            # allmiss.extend([obscol]*len(scidata))
            # allurl.extend([thisurl]*len(scidata))
            # >--
            if verbose:
                log.warning('%s: %s: %s', obscol, o, len(scidata))
            pass
        else:
            log.warning('>-- No data in MAST query %s', errmastq)
        pass
    tempdir = tempfile.mkdtemp(
        dir=dawgie.context.data_stg, prefix=target.replace(' ', '') + '_'
    )
    thisobsids = [row['productFilename'].split('_')[-2] for row in allsci]

    for irow, row in enumerate(allsci):
        # HST: mast api hack
        if allmiss[irow] in ['HST']:
            thisobsid = row['productFilename'].split('_')[-2]
            # 6/13/24 note that GJ 1132 STIS is sometimes having trouble here
            # the MAST downloaded file is ldlm01m0q_flt.fits but
            # it's coming up as a ima.fits maybe?
            if row['project'] in ['CALSTIS']:
                thisobsid = thisobsid.upper() + '%2F' + thisobsid + '_flt.fits'
                pass
            else:
                # want ida504e9 --> IDA504E9Q%2Fida504e9q_ima.fits
                # and ida504e9q --> IDA504E9QQ%2Fida504e9q_ima.fits
                # (currently the second one gets a double 'qq' toward the end,
                # which fails)
                # special case for K2-3 and others with a couple weird 's' files
                if len(thisobsid) == 9 and thisobsid.endswith('s'):
                    thisobsid = (
                        thisobsid[:-1].upper()
                        + 'QQ%2F'
                        + thisobsid
                        + '_ima.fits'
                    )
                    pass
                # remove the second double-q (the lower-case one)
                elif (
                    len(thisobsid) == 9
                    and thisobsid.endswith('q')
                    and not thisobsid.startswith('ldl')
                ):
                    thisobsid = (
                        thisobsid.upper() + 'Q%2F' + thisobsid + '_ima.fits'
                    )
                    # thisobsid = thisobsid.replace('qq_ima','q_ima')
                    pass
                # special case for K2-3 and others with a couple weird 's' files
                elif thisobsid + 's' in thisobsids:
                    thisobsid = (
                        thisobsid.upper() + 'Q%2F' + thisobsid + 's_ima.fits'
                    )
                    pass
                else:
                    thisobsid = (
                        thisobsid.upper() + 'Q%2F' + thisobsid + 'q_ima.fits'
                    )
                    pass
            fileout = os.path.join(tempdir, os.path.basename(thisobsid))
            shellcom = "'".join(
                [
                    "curl -s -L -X GET ",
                    allurl[irow] + "product_name=" + thisobsid,
                    " --output ",
                    fileout,
                    "",
                ]
            )
            subprocess.run(shellcom, shell=True, check=False)
            pass
        # JWST
        elif allmiss[irow] in ['JWST']:
            payload = {"uri": row['dataURI']}
            resp = requests.get(allurl[irow], params=payload, timeout=42)
            fileout = os.path.join(
                tempdir, os.path.basename(row['productFilename'])
            )
            with open(fileout, 'wb') as flt:
                flt.write(resp.content)
                pass
            pass
        # SPITZER DO NOTHING FOR NOW DL IS TOO LONG
        else:
            pass
        if verbose:
            log.warning('>-- %s %s', os.path.getsize(fileout), fileout)
            pass
        pass
    if allraw:
        for irow, row in enumerate(allraw):
            # print('  downloading',irow,len(allraw))
            payload = {"uri": row['dataURI']}
            try:
                resp = requests.get(allurl[irow], params=payload, timeout=42)
                fileout = os.path.join(
                    tempdir, os.path.basename(row['productFilename'])
                )
                with open(fileout, 'wb') as flt:
                    flt.write(resp.content)
            except requests.exceptions.ReadTimeout:
                log.warning(
                    '--< TIMEDOUT on the allraw request loop for %s (%s/%s) >--',
                    target,
                    irow,
                    len(allraw),
                )
            pass
        pass
    locations = [tempdir]
    new = dbscp(target, locations, dbs, out)
    # Delete local copy /proj/sdp/data/stg
    # Disabling for now I wanna know what s happening there
    # --<
    # os.chmod(tempdir, int('0777', 8))
    shutil.rmtree(tempdir, True)
    # >--
    return new


# ---------- ---------------------------------------------------------
# -- DISK -- ---------------------------------------------------------
def disk(selfstart, out, diskloc, dbs):
    '''Query on disk data'''
    merge = False
    target_id = list(selfstart['starID'].keys())

    # This was originally designed because we wanted to ingest data from X or Y
    # that had crazy directory names.
    # Ditch it?
    # Will require standardized directory names, no Kepler but KEPLER
    # --<
    targets = trgedit.targetondisk.__doc__
    targets = targets.split('\n')
    targets = [t.strip() for t in targets if t.replace(' ', '')]
    locations = None
    for t in targets:
        parsedstr = t.split(':')
        parsedstr = [tt.strip() for tt in parsedstr]
        if parsedstr[0] in target_id:
            folders = parsedstr[1]
            folders = folders.split(',')
            folders = [f.strip() for f in folders]
            locations = [os.path.join(diskloc, fold) for fold in folders]
            pass
        pass
    if locations is None:
        log.warning('ADD data subdirectory name to list in target/edit.py!!')
        pass
    # >--
    lookforme = target_id[0]
    # Maybe get rid of this it may work with the original name
    # --<
    lookforme = lookforme.replace(' ', '')
    lookforme = lookforme.replace('-', '')
    lookforme = lookforme.upper()
    # >--
    stdlocations = [os.path.join(diskloc, lookforme)]
    if locations is None:
        locations = stdlocations
        pass
    # make sure that the data storage directory exists for this star
    for loc in locations:
        if not os.path.exists(loc):
            os.makedirs(loc)
            pass
        pass
    if locations is not None:
        merge = dbscp('on disk', locations, dbs, out, verbose=False)
        pass
    return merge


# ---------- ---------------------------------------------------------
# -- DBS COPY -- -----------------------------------------------------
def dbscp(target, locations, dbs, out, verbose=False):
    '''Format data into SV'''
    copied = False
    imalist = None
    for loc in locations:
        if imalist is None:
            imalist = [
                os.path.join(loc, imafits)
                for imafits in os.listdir(loc)
                if '.fits' in imafits
            ]
            pass
        else:
            imalist.extend(
                [
                    os.path.join(loc, imafits)
                    for imafits in os.listdir(loc)
                    if '.fits' in imafits
                ]
            )
            pass
        pass

    okfiles = []
    n_fails = 0
    for fitsfile in imalist:
        try:
            pyfits.open(fitsfile)
            okfiles.append(True)
            pass
        except OSError:
            log.warning('--< !!! Failed to read: %s ', fitsfile)
            okfiles.append(False)
            n_fails += 1
            pass
        pass
    log.warning(
        '--< %s: %i out of %i MAST files are unreadable >--',
        target,
        n_fails,
        len(imalist),
    )

    if n_fails == 0:
        for fitsfile in imalist:
            filedict = {
                'md5': None,
                'sha': None,
                'loc': None,
                'observatory': None,
                'instrument': None,
                'detector': None,
                'filter': None,
                'mode': None,
            }
            m = subprocess.check_output(['md5sum', '-b', fitsfile])
            s = subprocess.check_output(['sha1sum', '-b', fitsfile])
            m = m.decode('utf-8').split(' ')
            s = s.decode('utf-8').split(' ')
            md5 = m[0]
            sha = s[0]
            filedict['md5'] = md5
            filedict['sha'] = sha
            filedict['loc'] = fitsfile
            with pyfits.open(fitsfile) as pf:
                mainheader = pf[0].header
                # HST
                if 'SCI' in mainheader.get('FILETYPE', ''):
                    keys = [k for k in mainheader.keys() if k != '']
                    keys = list(set(keys))
                    if 'TELESCOP' in keys:
                        filedict['observatory'] = mainheader['TELESCOP']
                    if 'INSTRUME' in keys:
                        filedict['instrument'] = mainheader['INSTRUME']
                    if 'DETECTOR' in keys:
                        filedict['detector'] = mainheader['DETECTOR'].replace(
                            '-', '.'
                        )
                    if 'FILTER' in keys:
                        filedict['filter'] = mainheader['FILTER']
                    if ('SCAN_RAT' in keys) and (mainheader['SCAN_RAT'] > 0):
                        filedict['mode'] = 'SCAN'
                        pass
                    else:
                        filedict['mode'] = 'STARE'
                    if 'FOCUS' in keys and filedict['detector'] is None:
                        filedict['detector'] = mainheader['FOCUS']
                        filedict['mode'] = 'IMAGE'
                        pass
                    if filedict['instrument'] in ['STIS']:
                        filedict['filter'] = mainheader['OPT_ELEM']
                        pass
                    out['name'][mainheader['ROOTNAME']] = filedict
                    pass
                # Spitzer
                elif (
                    'spitzer' in mainheader.get('TELESCOP').lower()
                    and 'sci' in mainheader.get('EXPTYPE').lower()
                ):
                    filedict['observatory'] = mainheader.get('TELESCOP')
                    filedict['instrument'] = mainheader.get('INSTRUME')
                    filedict['mode'] = mainheader.get('READMODE')
                    filedict['detector'] = 'IR'
                    if mainheader.get('CHNLNUM') == 1:
                        filedict['filter'] = "36"
                    elif mainheader.get('CHNLNUM') == 2:
                        filedict['filter'] = "45"
                    else:
                        filedict['filter'] = mainheader.get('CHNLNUM')
                    out['name'][str(mainheader.get('DPID'))] = filedict
                    pass
                # JWST
                elif 'JWST' in mainheader.get('TELESCOP'):
                    filedict['observatory'] = mainheader.get('TELESCOP')
                    filedict['instrument'] = mainheader.get('INSTRUME')
                    filedict['detector'] = mainheader.get('DETECTOR')
                    filedict['filter'] = mainheader.get('FILTER')
                    if mainheader.get('PUPIL'):
                        filedict['mode'] = mainheader.get('PUPIL')
                    else:
                        filedict['mode'] = mainheader.get('GRATING', '')
                    key = mainheader.get("FILENAME").split('.')[0]
                    out['name'][key] = filedict
                    if verbose:
                        log.warning(
                            '--< %s: %s %s %s %s %s',
                            key,
                            filedict['observatory'],
                            filedict['instrument'],
                            filedict['detector'],
                            filedict['filter'],
                            filedict['mode'],
                        )
                        pass
                    pass
                pass
            mastout = os.path.join(dbs, md5 + '_' + sha)
            onmast = os.path.isfile(mastout)
            if not onmast:
                shutil.copyfile(fitsfile, mastout)
                os.chmod(mastout, int('0664', 8))
                pass
            pass
        copied = True
        pass
    else:
        copied = False
    out['STATUS'].append(copied)
    return copied


# -------------- -----------------------------------------------------
