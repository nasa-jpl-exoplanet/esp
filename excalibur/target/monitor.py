'''target.monitor ds'''

# -- IMPORTS -- ------------------------------------------------------
import email.message
import logging

import numpy
import smtplib
import math

log = logging.getLogger(__name__)
from collections import defaultdict

# ------------- ------------------------------------------------------

care_about = ['t0']


def _diff(vl):
    '''_diff ds'''
    if 1 < len(vl):
        try:
            d = numpy.nan
            a = float(vl[0])
            d = 0
            b = float(vl[1])
            d = (a - b) / b * 100
        except ValueError:
            pass
    else:
        d = 0
    return d


def _outlier(vl):
    '''Finds whether first element of vl is within 5 sigma of other elems'''
    # 1 or more element of vl is an empty string there '' and that triggers bad things
    if 1 < len(vl):
        # Cheap fix attempt with try except
        try:
            vl = [float(v) for v in vl]

            if math.isnan(vl[0]):
                is_outlier = False
            else:
                vl_prev = numpy.array(vl[1:])
                vl_prev = vl_prev[~numpy.isnan(vl_prev)]  # clear all nans
                if len(vl_prev) > 1:
                    mean = numpy.mean(vl_prev)
                    std = numpy.std(vl_prev)
                    is_outlier = abs(vl[0] - mean) > 5 * std
                else:
                    is_outlier = False
                pass
            pass
        except ValueError:
            is_outlier = False
    else:
        is_outlier = False  # only 1 or 0 elems; no outlier can exist
    return is_outlier


def alert(
    asp: {str: {str: {str: object}}}, known: [], table: []
) -> ([], [], []):
    '''alert ds'''
    changes, kwn, tab = [], [], []
    for target in asp:
        kwn.append(target)
        tab.append(asp[target]['target.variations_of.parameters']['last'])

        if target in known:
            index = known.index(target)
            for pk in tab[-1].keys():
                if pk in table[index]:
                    lp_isnan = numpy.isnan(tab[-1][pk])
                    tb_isnan = numpy.isnan(table[index][pk])

                    if any(
                        [
                            lp_isnan and not tb_isnan,
                            not lp_isnan and tb_isnan,
                            (
                                tab[-1][pk] != table[index][pk]
                                and not (lp_isnan and tb_isnan)
                            ),
                        ]
                    ):
                        pks = pk.split('_')
                        remark = (
                            str(target)
                            + '::'
                            + str(pks[0])
                            + ' '
                            + '_'.join(pks[1:])
                            + ' has transitiond from {0} to {1}'
                        )
                        remark = remark.format(
                            str(table[index][pk]), str(tab[-1][pk])
                        )
                        changes.append(remark)

    changes.sort()

    if changes:
        try:
            msg = email.message.EmailMessage()
            msg.set_content('\n'.join(changes))
            msg['Subject'] = 'Alert: target parameter changes detected'
            msg['From'] = 'do-not-reply@excalibur.jpl.nasa.gov'
            msg['To'] = 'sdp@jpl.nasa.gov'
            s = smtplib.SMTP('localhost')
            s.send_message(msg)
        except:  # fmt: skip # noqa: E722 # pylint: disable=bare-except # because we do not expect exceptions and do not want to crash
            log.exception('Could not send alert email')
        pass
    return changes, kwn, tab


def regress(
    planet: {}, rids: [], tl: {str: {str: {str: object}}}
) -> ({str: float}, {str: []}, []):
    '''regress ds'''
    for i, rid in enumerate(tl):
        if rid in rids:
            break
        rids.insert(i, rid)
        svv = list(tl[rid]['target.autofill.parameters']['starID'].values())
        svv = svv[0]
        for p in svv['planets']:
            for ca in care_about:
                k = '_'.join([p, ca])
                if k in planet:
                    planet[k].insert(i, svv[p][ca][0])
                else:
                    planet[k] = [svv[p][ca][0]]
                pass
            pass
        pass
    # last = dict([(pp, _diff (vl)) for pp,vl in planet.items()])
    # I believe this might have been a debug thingy, commenting the print statement
    # print (planet.items())
    last = {pp: _diff(vl) for pp, vl in planet.items()}
    # outliers = dict([(pp, _outlier (vl)) for pp,vl in planet.items()])
    outliers = {pp: _outlier(vl) for pp, vl in planet.items()}
    return last, outliers


def regress_for_frame_counts(
    data: {int: {str: int, str: int}}, tl: {str: {str: {str: object}}}
) -> ({str: float}, {str: []}, []):
    '''
    this functions regresses through the runids for target.scrape.databases.
    data currently contains:

    observatory-instrument-detector-filter-mode: # of occurences

    timeline iterates through descending order of runID. Therefore
    the newest runIDs go first, then code breaks once it hits a runID
    that has already been computed.
    #'''
    # print('in monitor.regress_through_runid')

    for rid, val in tl.items():
        # print(rid)
        # i think i need to loop through all rids because if there is a new thing
        # in the most recent runid, I actually need to go back and
        if rid in data:
            break
        data[rid] = defaultdict(int)
        # only one product so loops only once
        for _, d in val.items():
            frames_d = dict(d)['name']
            for _, frame_d in frames_d.items():
                if frame_d['observatory'] in ['HST', 'JWST']:
                    identifier = (
                        frame_d['observatory']
                        + '-'
                        + frame_d['instrument']
                        + '-'
                        + frame_d['detector']
                        + '-'
                        + frame_d['filter']
                        + '-'
                    )
                    if frame_d['mode'] is None:
                        identifier += 'STARE'
                    else:
                        identifier += frame_d['mode']
                    # print(identifier)
                    data[rid][identifier] += 1
    # print(data)

    # getting the most current frame counts and previous frames counts to compare
    sorted_keys = sorted(data.keys(), reverse=True)
    cur_frames = data[sorted_keys[0]]
    prev_frames = data[sorted_keys[1]]

    # print(f'current: {cur_frames}')
    # print(f'previous: {prev_frames}')

    status = 1
    for identifier in prev_frames:
        if identifier not in cur_frames:
            # checking if previously this key existed (ie, had some frames) but in current doesn't?
            status = -1
            break

    return data, status
