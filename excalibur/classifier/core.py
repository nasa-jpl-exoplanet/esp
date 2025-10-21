'''classifier core ds'''

# Heritage code shame:
# pylint: disable=too-many-branches,too-many-locals,too-many-nested-blocks,too-many-statements

# -- IMPORTS -- ------------------------------------------------------
import os
import dawgie
import numpy as np
import math
import logging

import excalibur
import excalibur.system.core as syscore

log = logging.getLogger(__name__)


# -------------------------------------------------------------------
def lc_resid_classification_version():
    '''lc_resid_classification_version ds'''
    return dawgie.VERSION(1, 0, 0)


def lc_resid_classification(transit_whitelight, ext, out):
    '''
    A. Arora
    Classify the shape of an exoplanet's light curve residual
    '''

    for planet in transit_whitelight['data']:
        if 'Spitzer' in ext:
            initial_dict = transit_whitelight['data'][planet][0]
        else:
            initial_dict = transit_whitelight['data'][planet]

        # note that Spitzer doesn't have 'postim', so if-then needed here
        try:
            sep = np.array(initial_dict['postsep'])
            whitelight = np.array(initial_dict['allwhite'])
            if 'postim' in initial_dict.keys():
                whitelight /= np.array(initial_dict['postim'])
            model = np.array(initial_dict['postlc'])
        except AttributeError:
            sep = []
            whitelight = []
            model = []
            for p in initial_dict:
                sep.extend(p['postsep'])
                whitelight.append(np.array(p['allwhite']))
                if 'postim' in p.keys():
                    whitelight[-1] /= np.array(p['postim'])
                model.append(p['postlc'])
                pass
            pass

        wl_flat = (
            np.concatenate(whitelight).ravel()
            if 'Spitzer' not in ext
            else whitelight
        )

        pred = 0  # fake bad torch result

        flags = {
            0: 'yellow',
            1: 'green',
        }

        flag_color = flags[pred]
        flag_descrip = ''
        if pred == 0:
            flag_descrip += 'Light curve residual has a bad shape.'
        else:
            flag_descrip += 'Light curve residual is satisfactory.'

        # avoid overwriting existing data
        try:
            out['data'][planet]['residual_shape'] = {}
        except KeyError:
            out['data'][planet] = {}
            out['data'][planet]['residual_shape'] = {}

        out['data'][planet]['residual_shape']['data'] = wl_flat
        out['data'][planet]['residual_shape']['model'] = model
        out['data'][planet]['residual_shape']['z'] = sep
        out['data'][planet]['residual_shape']['flag_color'] = flag_color
        out['data'][planet]['residual_shape']['flag_descrip'] = flag_descrip

        pass

    out['STATUS'].append(True)
    return True


def cpwlversion():
    '''cpwlversion ds'''
    return dawgie.VERSION(1, 0, 0)


# cpwl stands for count points whitelight.
def cpwl(transit_whitelight, priors, ext, out):
    '''
    K. Mccarthy
    Counts points in transit curve
    '''

    for planet in transit_whitelight['data']:

        if 'Spitzer' in ext:
            initial_dict = transit_whitelight['data'][planet][0]
        else:
            initial_dict = transit_whitelight['data'][planet]

        try:
            # test = transit_whitelight['data'][planet].keys()
            sep = initial_dict[
                'postsep'
            ]  # alternatively, could use orbital phase
            whitelight = np.array(initial_dict['allwhite'])
            pass
        except AttributeError:
            sep = []
            whitelight = []
            for p in initial_dict:
                sep.extend(p['postsep'])
                whitelight.append(p['allwhite'])
                pass
            pass

        wl_flat = (
            np.concatenate(whitelight).ravel().tolist()
            if 'Spitzer' not in ext
            else whitelight
        )

        # each point is in [sep, wl] format.
        points = np.column_stack((np.array(sep), np.array(wl_flat)))
        transit_points_tot = []  # between contact points 1 and 4
        transit_points_full = []  # between contact points 2 and 3

        # get variables required for z equations
        r_p = priors[planet]['rp']  # planet radius (Jupiter radii)
        r_s = priors['R*']  # star radius (solar radii)

        sma = priors[planet]['sma']  # semi major axis (in AU) of planet's orbit
        inc = priors[planet][
            'inc'
        ]  # inclination (in degrees) of planet's orbit

        # get constants for unit conversions
        ssc = syscore.ssconstants()
        b = (((math.cos(math.radians(inc))) * sma) / r_s) * ssc[
            'Rsun/AU'
        ]  # impact parameter
        radius_ratio = r_p / r_s * ssc['Rjup/Rsun']

        # use z equations to find postsep vals at the 4 contact points:
        # Note that z_1 = -1 * z_3, and z_2 = -1 * z_4.
        # Therefore, we can use abs(t) to find if a value t is between contact points 1 and 4 or between 2 and 3.
        z_3 = math.sqrt(
            ((1 - radius_ratio) ** 2) - (b**2)
        )  # postsep val at 3rd contact point
        z_4 = math.sqrt(
            ((1 + radius_ratio) ** 2) - (b**2)
        )  # postsep val at 4th contact point

        # count the number of points inside the transit
        for p in points:

            # between 1st and 4th contact point
            if abs(p[0]) <= z_4:
                transit_points_tot.append(p)

            # between 2nd and 3rd contact point
            if abs(p[0]) <= z_3:
                transit_points_full.append(p)

        # flag depending on number of transit points

        # color key
        flags = {0: 'green', 1: 'yellow', 2: 'red'}

        flag_val = 0
        flag_descrip = ""

        # No points between 1st and 4th contact.
        if len(transit_points_tot) == 0:
            flag_val = max(flag_val, 2)
            flag_descrip += (
                "Insufficient points between 1st and 4th contact points."
            )

        # Not many points between 1st and 4th contact.
        elif len(transit_points_tot) < 6:
            flag_val = max(flag_val, 1)
            flag_descrip += (
                "Insufficient points between 1st and 4th contact points."
            )

        else:
            flag_descrip += (
                "Sufficient points between 1st and 4th contact points."
            )

        # No points between 2nd and 3rd contact.
        if len(transit_points_full) == 0:

            # Not many points between 1st and 4th contact, either.
            if len(transit_points_tot) < 6:
                flag_val = max(flag_val, 2)
                flag_descrip += " No points between 2nd and 3rd contact points."

            # At least 6 points between 1st and 4th contact, so mark yellow instead of red.
            else:
                flag_val = max(flag_val, 1)
                flag_descrip += " No points between 2nd and 3rd contact points, but at least 6 points between 1st and 4th contact points."

        # Not many points between 2nd and 3rd contact.
        elif len(transit_points_full) < 5:
            flag_val = max(flag_val, 1)
            flag_descrip += (
                " Insufficient points between 2nd and 3rd contact points."
            )

        else:
            flag_descrip += (
                " Sufficient points between 2nd and 3rd contact points."
            )

        flag_color = flags[flag_val]

        # avoid overwriting existing data
        try:
            out['data'][planet]['count_points_wl'] = {}
        except KeyError:
            out['data'][planet] = {}
            out['data'][planet]['count_points_wl'] = {}

        out['data'][planet]['count_points_wl']['flag_color'] = flag_color
        out['data'][planet]['count_points_wl']['flag_descrip'] = flag_descrip
        out['data'][planet]['count_points_wl']['total'] = len(
            transit_points_tot
        )
        out['data'][planet]['count_points_wl']['full'] = len(
            transit_points_full
        )

        pass
    out['STATUS'].append(True)
    return True


def symwlversion():
    '''symwlversion ds'''
    return dawgie.VERSION(1, 0, 0)


# symwl stands for symmetry whitelight.
def symwl(transit_whitelight, priors, ext, out):
    '''
    K. McCarthy
    Check symmetry of light curve
    '''

    for planet in transit_whitelight['data']:

        if 'Spitzer' in ext:
            initial_dict = transit_whitelight['data'][planet][0]
        else:
            initial_dict = transit_whitelight['data'][planet]

        try:
            sep = initial_dict[
                'postsep'
            ]  # (alternatively, could use orbital phase)
            whitelight = np.array(initial_dict['allwhite'])
            pass
        except AttributeError:
            sep = []
            whitelight = []
            for p in initial_dict:
                sep.extend(p['postsep'])
                whitelight.append(p['allwhite'])
                pass
            pass

        wl_flat = (
            np.concatenate(whitelight).ravel().tolist()
            if 'Spitzer' not in ext
            else whitelight
        )

        # each point is in [sep, wl] format.
        points = np.column_stack((np.array(sep), np.array(wl_flat)))
        transit_points_before_1 = []  # before contact point 1
        transit_points_after_4 = []  # after contact point 4

        # get variables required for z equations
        r_p = priors[planet]['rp']  # planet radius (Jupiter radii)
        r_s = priors['R*']  # star radius (solar radii)

        sma = priors[planet]['sma']  # semi major axis (in AU) of planet's orbit
        inc = priors[planet][
            'inc'
        ]  # inclination (in degrees) of planet's orbit

        # get constants for unit conversions
        ssc = syscore.ssconstants()
        b = (((math.cos(math.radians(inc))) * sma) / r_s) * ssc[
            'Rsun/AU'
        ]  # impact parameter
        radius_ratio = r_p / r_s * ssc['Rjup/Rsun']

        # use z equations to find postsep vals at the 4 contact points:
        # Note that z_1 = -1 * z_3, and z_2 = -1 * z_4.
        # Therefore, we can use abs(t) to find if a value t is between contact points 1 and 4 or between 2 and 3.
        z_4 = math.sqrt(
            ((1 + radius_ratio) ** 2) - (b**2)
        )  # postsep val at 4th contact point
        z_1 = -1 * z_4

        # count the number of points outside of the transit on both sides
        for p in points:
            # outside 1st or 4th contact point
            if p[0] > z_4:
                transit_points_after_4.append(p)
            if p[0] < z_1:
                transit_points_before_1.append(p)

        # flag depending on number of pre- and post-transit points

        # color key
        flags = {0: 'green', 1: 'yellow', 2: 'red'}

        flag_val = 0
        flag_descrip = ""

        if (
            len(transit_points_before_1) == 0
            and len(transit_points_after_4) == 0
        ):
            flag_val = max(flag_val, 2)
            flag_descrip = "Insufficient points before 1st contact point and insufficient points after 4th contact point."
        elif len(transit_points_before_1) == 0:
            flag_val = max(flag_val, 1)
            flag_descrip = "Insufficient points before 1st contact point."
        elif len(transit_points_after_4) == 0:
            flag_val = max(flag_val, 1)
            flag_descrip = "Insufficient points after 4th contact point."
        else:
            flag_val = max(flag_val, 0)
            flag_descrip = "Sufficient points before 1st contact point and after 4th contact point."

        flag_color = flags[flag_val]

        # try/except to avoid overwriting any existing data in out['data'][planet].
        try:
            out['data'][planet]['symmetry_wl'] = {}
        except KeyError:
            out['data'][planet] = {}
            out['data'][planet]['symmetry_wl'] = {}

        out['data'][planet]['symmetry_wl']['flag_color'] = flag_color
        out['data'][planet]['symmetry_wl']['flag_descrip'] = flag_descrip
        out['data'][planet]['symmetry_wl']['left'] = len(
            transit_points_before_1
        )
        out['data'][planet]['symmetry_wl']['right'] = len(
            transit_points_after_4
        )

        pass

    out['STATUS'].append(True)
    return True


def rsdmversion():
    '''RSDMversion ds'''
    return dawgie.VERSION(1, 0, 0)


def rsdm(transit_spectrum, out):
    '''
    K. McCarthy
    Calculate RSDM (in transit.spectrum) and flag target accordingly
    '''

    for planet in transit_spectrum['data']:

        try:
            spectrum_data = transit_spectrum['data'][planet]
            pass

        except AttributeError:
            pass

        if 'LCFIT' in spectrum_data:
            spec_rsdpn = [
                np.nanstd(i['residuals']) / i['dnoise']
                for i in spectrum_data['LCFIT']
            ]
            avg_rsdpn = np.nanmean(spec_rsdpn)

            # flag depending on mean RSDM
            # color key
            flags = {0: 'green', 1: 'yellow', 2: 'red'}

            flag_val = 0
            flag_descrip = ""

            green_upper_bound = 8.804
            yellow_upper_bound = 13.818

            if avg_rsdpn >= yellow_upper_bound:
                flag_val = 2
                flag_descrip = "Mean RSDM above threshold."
            if avg_rsdpn >= green_upper_bound:
                flag_val = 1
                flag_descrip = "Mean RSDM above threshold."
            else:
                flag_val = 0
                flag_descrip = "Mean RSDM is below threshold."

            flag_color = flags[flag_val]

            # avoid overwriting existing planet data
            try:
                out['data'][planet]['rsdm'] = {}
            except KeyError:
                out['data'][planet] = {}
                out['data'][planet]['rsdm'] = {}

            out['data'][planet]['rsdm']['mean_rsdm'] = avg_rsdpn
            out['data'][planet]['rsdm']['flag_color'] = flag_color
            out['data'][planet]['rsdm']['flag_descrip'] = flag_descrip

    out['STATUS'].append(True)
    return True


def perc_rejected_version():
    '''perc_rejected_version ds'''
    return dawgie.VERSION(1, 0, 0)


def perc_rejected(transit_spectrum, out):
    '''
    K. McCarthy
    Calculate the percentage of spectral channels rejected in cumulative spectrum distribution (in transit.spectrum) and flag target accordingly
    '''

    for planet in transit_spectrum['data']:

        try:
            spectrum_data = transit_spectrum['data'][planet]
            pass

        except AttributeError:
            pass

        # % spectral channels rejected in cumulative spectrum distribution
        planet_spectrum = np.array(spectrum_data['ES'])
        vspectrum = np.array(planet_spectrum)
        vspectrum = vspectrum**2
        perc_rejected_value = (
            len([i for i in vspectrum if np.isnan(i)]) / len(vspectrum) * 100
        )

        # flag depending on percent rejected
        # color key
        flags = {0: 'green', 1: 'yellow', 2: 'red'}

        flag_val = 0
        flag_descrip = ""

        green_upper_bound = 31.844
        yellow_upper_bound = 52.166

        if perc_rejected_value >= yellow_upper_bound:
            flag_val = 2
            flag_descrip = "Percent Rejected in Cumulative Spectrum Distribution is very high."
        elif perc_rejected_value >= green_upper_bound:
            flag_val = 1
            flag_descrip = (
                "Percent Rejected in Cumulative Spectrum Distribution is high."
            )
        else:
            flag_val = 0
            flag_descrip = "Percent Rejected in Cumulative Spectrum Distribution is below threshold."

        flag_color = flags[flag_val]

        # avoid overwriting existing planet data
        try:
            out['data'][planet]['perc_rejected'] = {}
        except KeyError:
            out['data'][planet] = {}
            out['data'][planet]['perc_rejected'] = {}

        out['data'][planet]['perc_rejected'][
            'percent_rejected_value'
        ] = perc_rejected_value
        out['data'][planet]['perc_rejected']['flag_color'] = flag_color
        out['data'][planet]['perc_rejected']['flag_descrip'] = flag_descrip

    out['STATUS'].append(True)
    return True


def median_error_version():
    '''median_error_version ds'''
    return dawgie.VERSION(1, 0, 0)


def median_error(data_calibration, out):
    '''
    K. McCarthy
    Calculate median error in data.calibration and flag target accordingly
    '''

    try:
        # get error distribution (adapted from esp/excalibur/data/core.py)
        data = data_calibration['data']
        spec = np.array(
            [d for d, i in zip(data['SPECTRUM'], data['IGNORED']) if not i]
        )
        wave = np.array(
            [d for d, i in zip(data['WAVE'], data['IGNORED']) if not i]
        )
        errspec = np.array(
            [d for d, i in zip(data['SPECERR'], data['IGNORED']) if not i]
        )
        vrange = data['VRANGE']
        allerr = []
        for s, e, w in zip(spec, errspec, wave):
            select = (w > vrange[0]) & (w < vrange[1])
            allerr.extend(e[select] / np.sqrt(s[select]))
            pass
        allerr = np.array(allerr)
        select = np.isfinite(allerr)
        allerr = allerr[select]
        allerr = allerr[allerr > 0.9]

        median_err = np.median(allerr)

        # flag depending on error distribution
        # color key
        flags = {0: 'green', 1: 'yellow', 2: 'red'}

        flag_val = 0
        flag_descrip = ""

        green_upper_bound = 2.5

        if median_err >= green_upper_bound:
            flag_val = 1
            flag_descrip = (
                "Median error in data.calibration is above threshold."
            )
        else:
            flag_val = 0
            flag_descrip = (
                "Median error in data.calibration is below threshold."
            )

        flag_color = flags[flag_val]

        out['data']['median_error'] = {}
        out['data']['median_error']['median_error_value'] = median_err
        out['data']['median_error']['flag_color'] = flag_color
        out['data']['median_error']['flag_descrip'] = flag_descrip

    except KeyError:
        pass

    out['STATUS'].append(True)
    return True


# ---------------------------- ---------------------------------------
def savesv(aspects, fltrs):
    '''
    save the results as a csv file in /proj/data/spreadsheets
    '''
    svname = 'classifier.flags'

    # directory where the results are saved
    save_dir = excalibur.context['data_dir'] + '/spreadsheets/'
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    flag_types = [
        'count_points_wl',
        'symmetry_wl',
        'rsdm',
        'perc_rejected',
        'median_error',
        'residual_shape',
        'overall_flag',
    ]

    # file name where the results are saved
    out_file_name = svname.replace('.', '_') + '.csv'

    with open(save_dir + out_file_name, 'w', encoding='ascii') as outfile:

        # write the header row
        outfile.write('star,planet,filter,')

        for key in flag_types:
            outfile.write(key + ',')
        outfile.write('\n')

        # loop through each target, with one row per planet
        for trgt in aspects:
            # loop through each filters
            for fltr in fltrs:

                transit_alg_fltr_name = svname + '.transit-' + str(fltr)

                if transit_alg_fltr_name in aspects[trgt].keys():
                    alg_flag_data = aspects[trgt][transit_alg_fltr_name]

                    if alg_flag_data['STATUS'][-1]:
                        for k in alg_flag_data[
                            'data'
                        ].keys():  # k is either planet letter or median_error
                            if k == 'median_error':
                                # skip the weirdo one that's not associated with a specific planet
                                pass
                            else:
                                planet_letter = k

                                colors = {}
                                for flag_type in flag_types:
                                    colors[flag_type] = '-'
                                if (
                                    'median_error'
                                    in alg_flag_data['data'].keys()
                                ):
                                    # add in the weirdo one that's not associated with a specific planet
                                    colors['median_error'] = alg_flag_data[
                                        'data'
                                    ]['median_error']['flag_color']
                                for alg_flag in alg_flag_data['data'][k]:
                                    if alg_flag == 'overall_flag':
                                        alg_flag_color = alg_flag_data['data'][
                                            k
                                        ][alg_flag]
                                    else:
                                        alg_flag_color = alg_flag_data['data'][
                                            k
                                        ][alg_flag]['flag_color']

                                    colors[alg_flag] = alg_flag_color

                                    if alg_flag not in flag_types:
                                        log.warning(
                                            '--< UNEXPECTED FLAG TYPE >-- %s',
                                            alg_flag,
                                        )

                                outfile.write(trgt + ',')
                                outfile.write(planet_letter + ',')
                                outfile.write(fltr + ',')
                                for flag_type in flag_types:
                                    outfile.write(colors[flag_type] + ',')
                                outfile.write('\n')

    return
