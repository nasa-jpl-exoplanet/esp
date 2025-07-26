'''classifier algorithms ds'''

# Heritage code shame:
# pylint: disable=too-many-branches,too-many-locals,too-many-nested-blocks,too-many-statements

# -- IMPORTS -- ------------------------------------------------------
import logging

from collections import OrderedDict

import dawgie
import dawgie.context

import excalibur.transit as trn

# import excalibur.transit.core as trncore
import excalibur.transit.algorithms as trnalg

import excalibur.eclipse as ecl
import excalibur.eclipse.algorithms as eclalg

import excalibur.system as sys
import excalibur.system.algorithms as sysalg

# import excalibur.system.core as syscore

import excalibur.classifier.core as clscore
import excalibur.classifier.states as clsstates

# DATA MODULE
import excalibur.data as dat
import excalibur.data.algorithms as datalg

import excalibur.runtime.algorithms as rtalg
import excalibur.runtime.binding as rtbind

from excalibur.classifier.core import savesv
from importlib import import_module as fetch  # avoid cicular dependencies

from excalibur.util.checksv import checksv

log = logging.getLogger(__name__)

# -------------------------------------------------------------------
# -- ALGO RUN OPTIONS -----------------------------------------------
# FILTERS
fltrs = [str(fn) for fn in rtbind.filter_names.values()]


# -------------------------------------------------------------------
# -- ALGORITHMS -----------------------------------------------------
class SummarizeFlags(dawgie.Analyzer):
    '''K. MCCARTHY: Summarize flags across targets'''

    def __init__(self):
        '''__init__ ds'''
        self._version_ = dawgie.VERSION(1, 0, 0)
        self.__out = [clsstates.FlagSummarySV(ext) for ext in fltrs]
        return

    def name(self):
        '''name ds'''
        return 'summarize_flags'

    def feedback(self):
        '''feedback ds'''
        return []

    def state_vectors(self):
        '''state_vectors ds'''
        return self.__out

    def traits(self) -> [dawgie.SV_REF, dawgie.V_REF]:
        '''Input State Vectors'''

        return [
            *[
                dawgie.SV_REF(
                    fetch('excalibur.classifier').task,
                    Flags(),
                    Flags().state_vectors()[i],
                )
                for i in range(len(Flags().state_vectors()))
            ]
        ]

    def run(self, aspects: dawgie.Aspect):
        '''run ds'''
        flags_dict = OrderedDict()
        values_dict = OrderedDict()

        for idx, fltr in enumerate(fltrs):

            for trgt in aspects:

                state_vec_name = 'flags'

                transit_alg_fltr_name = (
                    'classifier.'
                    + str(state_vec_name)
                    + '.transit-'
                    + str(fltr)
                )
                # eclipse_alg_fltr_name = 'classifier.' + str(state_vec_name) + '.eclipse-' + str(fltr)

                if transit_alg_fltr_name in aspects[trgt].keys():

                    alg_flag_data = aspects[trgt][transit_alg_fltr_name]

                    if alg_flag_data['STATUS'][-1]:

                        for k in alg_flag_data['data'].keys():

                            # handle all planet-specific algorithms
                            if k != 'median_error':

                                for alg_flag in alg_flag_data['data'][k]:

                                    if alg_flag == 'overall_flag':
                                        alg_flag_color = alg_flag_data['data'][
                                            k
                                        ][alg_flag]

                                    else:
                                        alg_flag_color = alg_flag_data['data'][
                                            k
                                        ][alg_flag]['flag_color']

                                        # extract and save values
                                        for value_field in alg_flag_data[
                                            'data'
                                        ][k][alg_flag]:
                                            if value_field not in (
                                                "flag_color",
                                                "flag_descrip",
                                            ):
                                                try:
                                                    values_dict[idx][alg_flag][
                                                        value_field
                                                    ].append(
                                                        alg_flag_data['data'][
                                                            k
                                                        ][alg_flag][value_field]
                                                    )
                                                except KeyError:
                                                    try:
                                                        values_dict[idx][
                                                            alg_flag
                                                        ][value_field] = [
                                                            alg_flag_data[
                                                                'data'
                                                            ][k][alg_flag][
                                                                value_field
                                                            ]
                                                        ]
                                                    except KeyError:
                                                        try:
                                                            values_dict[idx][
                                                                alg_flag
                                                            ] = {}
                                                            values_dict[idx][
                                                                alg_flag
                                                            ][value_field] = [
                                                                alg_flag_data[
                                                                    'data'
                                                                ][k][alg_flag][
                                                                    value_field
                                                                ]
                                                            ]
                                                        except KeyError:
                                                            values_dict[idx] = (
                                                                {}
                                                            )
                                                            values_dict[idx][
                                                                alg_flag
                                                            ] = {}
                                                            values_dict[idx][
                                                                alg_flag
                                                            ][value_field] = [
                                                                alg_flag_data[
                                                                    'data'
                                                                ][k][alg_flag][
                                                                    value_field
                                                                ]
                                                            ]

                                    try:
                                        count = flags_dict[idx][alg_flag][
                                            alg_flag_color
                                        ][1]
                                        targets = flags_dict[idx][alg_flag][
                                            alg_flag_color
                                        ][0]
                                        targets.append(trgt)
                                        flags_dict[idx][alg_flag][
                                            alg_flag_color
                                        ] = (targets, count + 1)
                                    except KeyError:
                                        try:
                                            flags_dict[idx][alg_flag][
                                                alg_flag_color
                                            ] = ([trgt], 1)
                                        except KeyError:
                                            try:
                                                flags_dict[idx][alg_flag] = {}
                                                flags_dict[idx][alg_flag][
                                                    alg_flag_color
                                                ] = ([trgt], 1)
                                            except KeyError:
                                                # orderedDict so that display in summarize_flags is in the same consistent order
                                                flags_dict[idx] = OrderedDict()
                                                flags_dict[idx][alg_flag] = {}
                                                flags_dict[idx][alg_flag][
                                                    alg_flag_color
                                                ] = ([trgt], 1)

                            # handle all non-planet-specific algorithms
                            else:

                                if alg_flag_data['STATUS'][-1]:

                                    alg_flag_color = alg_flag_data['data'][k][
                                        'flag_color'
                                    ]

                                    # extract and save values
                                    for value_field in alg_flag_data['data'][k]:
                                        if value_field not in (
                                            "flag_color",
                                            "flag_descrip",
                                        ):
                                            try:
                                                values_dict[idx][k][
                                                    value_field
                                                ].append(
                                                    alg_flag_data['data'][k][
                                                        value_field
                                                    ]
                                                )
                                            except KeyError:
                                                try:
                                                    values_dict[idx][k][
                                                        value_field
                                                    ] = [
                                                        alg_flag_data['data'][
                                                            k
                                                        ][value_field]
                                                    ]
                                                except KeyError:
                                                    try:
                                                        values_dict[idx][k] = {}
                                                        values_dict[idx][k][
                                                            value_field
                                                        ] = [
                                                            alg_flag_data[
                                                                'data'
                                                            ][k][value_field]
                                                        ]
                                                    except KeyError:
                                                        values_dict[idx] = {}
                                                        values_dict[idx][k] = {}
                                                        values_dict[idx][k][
                                                            value_field
                                                        ] = [
                                                            alg_flag_data[
                                                                'data'
                                                            ][k][value_field]
                                                        ]

                                    try:
                                        count = flags_dict[idx][k][
                                            alg_flag_color
                                        ][1]
                                        targets = flags_dict[idx][k][
                                            alg_flag_color
                                        ][0]
                                        targets.append(trgt)
                                        flags_dict[idx][k][alg_flag_color] = (
                                            targets,
                                            count + 1,
                                        )
                                    except KeyError:
                                        try:
                                            flags_dict[idx][k][
                                                alg_flag_color
                                            ] = ([trgt], 1)
                                        except KeyError:
                                            try:
                                                flags_dict[idx][k] = {}
                                                flags_dict[idx][k][
                                                    alg_flag_color
                                                ] = ([trgt], 1)
                                            except KeyError:
                                                # orderedDict so that display in summarize_flags is in the same consistent order
                                                flags_dict[idx] = OrderedDict()
                                                flags_dict[idx][k] = {}
                                                flags_dict[idx][k][
                                                    alg_flag_color
                                                ] = ([trgt], 1)

            # Create output dictionary for classifier_flags
            try:
                self.__out[idx]['data']['classifier_flags'] = {}
            except KeyError:
                try:
                    self.__out[idx]['data'] = {}
                except KeyError:
                    self.__out[idx] = {}
                    self.__out[idx]['data'] = {}
                self.__out[idx]['data']['classifier_vals'] = {}

            if idx in flags_dict:
                for k in flags_dict[idx]:
                    self.__out[idx]['data']['classifier_flags'][k] = flags_dict[
                        idx
                    ][k]

            # Create output dictionary for classifier_vals
            try:
                self.__out[idx]['data']['classifier_vals'] = {}
            except KeyError:
                try:
                    self.__out[idx]['data'] = {}
                except KeyError:
                    self.__out[idx] = {}
                    self.__out[idx]['data'] = {}
                self.__out[idx]['data']['classifier_vals'] = {}

            if idx in values_dict:
                for k in values_dict[idx]:
                    self.__out[idx]['data']['classifier_vals'][k] = values_dict[
                        idx
                    ][k]

            self.__out[idx]['STATUS'].append(True)

        # save classifier-flag results as .csv file (in /proj/data/spreadsheets/)
        savesv(aspects, fltrs)

        aspects.ds().update()

        return


class Flags(dawgie.Algorithm):
    '''K. MCCARTHY: Perform data quality assessment checks and flag target accordingly'''

    def __init__(self):
        '''__init__ ds'''
        self._version_ = dawgie.VERSION(1, 0, 0)
        self.__rt = rtalg.Autofill()
        self.__out = [clsstates.FlagsSV('transit-' + ext) for ext in fltrs]
        self.__out.extend(
            [clsstates.FlagsSV('eclipse-' + ext) for ext in fltrs]
        )

        # storing input state vectors in dictionary to avoid pylint "too-many-instance-attributes" warning
        self.__state_vecs = {
            'finalize': sysalg.Finalize(),
            'spectrum': trnalg.Spectrum(),
            'eclspectrum': eclalg.Spectrum(),
            'whitelight': trnalg.WhiteLight(),
            'eclwhitelight': eclalg.WhiteLight(),
            'data_calib': datalg.Calibration(),
        }

        return

    def name(self):
        '''name ds'''
        return 'flags'

    def previous(self):
        '''previous ds'''
        return [
            dawgie.ALG_REF(sys.task, self.__state_vecs['finalize']),
            dawgie.ALG_REF(trn.task, self.__state_vecs['spectrum']),
            dawgie.ALG_REF(ecl.task, self.__state_vecs['eclspectrum']),
            dawgie.ALG_REF(trn.task, self.__state_vecs['whitelight']),
            dawgie.ALG_REF(ecl.task, self.__state_vecs['eclwhitelight']),
            dawgie.ALG_REF(dat.task, self.__state_vecs['data_calib']),
        ] + self.__rt.refs_for_validity()

    def state_vectors(self):
        '''state_vectors ds'''
        return self.__out

    def run(self, ds, ps):
        '''run ds'''

        # stop here if it is not a runtime target
        if not self.__rt.is_valid():
            log.info(
                '--< CLASSFIER.%s: not a valid target >--', self.name().upper()
            )

        else:
            svupdate = []
            vfin, sfin = checksv(
                self.__state_vecs['finalize'].sv_as_dict()['parameters']
            )

            for ext in fltrs:

                # transit.spectrum
                vsp, ssp = checksv(
                    self.__state_vecs['spectrum'].sv_as_dict()[ext]
                )
                e_vsp, e_ssp = checksv(
                    self.__state_vecs['eclspectrum'].sv_as_dict()[ext]
                )

                # transit.whitelight
                vwl, swl = checksv(
                    self.__state_vecs['whitelight'].sv_as_dict()[ext]
                )
                e_vwl, e_swl = checksv(
                    self.__state_vecs['eclwhitelight'].sv_as_dict()[ext]
                )

                # data.calibration
                vdc, sdc = checksv(
                    self.__state_vecs['data_calib'].sv_as_dict()[ext]
                )

                # ======================  COUNT_POINTS_WL  ====================== #####
                metric_name = "Point Count"

                # if transit.whitelight exists
                if vwl and vfin:
                    log.info('--< IN-TRANSIT POINT COUNT: %s >--', ext)

                    status = clscore.cpwl(
                        self.__state_vecs['whitelight'].sv_as_dict()[ext],
                        self.__state_vecs['finalize'].sv_as_dict()[
                            'parameters'
                        ]['priors'],
                        ext,
                        self.__out[fltrs.index(ext)],
                    )

                    if status:
                        svupdate.append(self.__out[fltrs.index(ext)])
                    pass
                else:
                    errstr = [m for m in [swl, sfin] if m is not None]
                    self._failure(errstr[0], metric_name)
                    pass

                # if eclipse.whitelight exists
                if e_vwl and vfin:
                    log.info('--< IN-ECLIPSE POINT COUNT: %s >--', ext)

                    status = clscore.cpwl(
                        self.__state_vecs['eclwhitelight'].sv_as_dict()[ext],
                        self.__state_vecs['finalize'].sv_as_dict()[
                            'parameters'
                        ]['priors'],
                        ext,
                        self.__out[len(fltrs) + fltrs.index(ext)],
                    )

                    if status:
                        svupdate.append(
                            self.__out[len(fltrs) + fltrs.index(ext)]
                        )
                    pass
                else:
                    errstr = [m for m in [e_swl, vfin] if m is not None]
                    self._failure(errstr[0], metric_name)
                    pass
                pass
                # ==================================================== #####

                # ======================  SYMMETRY_WL  ====================== #####
                metric_name = "Light Curve Symmetry"

                # if transit.whitelight exists
                if vwl and vfin:
                    log.info('--< TRANSIT LIGHT CURVE SYMMETRY: %s >--', ext)
                    status = clscore.symwl(
                        self.__state_vecs['whitelight'].sv_as_dict()[ext],
                        self.__state_vecs['finalize'].sv_as_dict()[
                            'parameters'
                        ]['priors'],
                        ext,
                        self.__out[fltrs.index(ext)],
                    )

                    if status:
                        svupdate.append(self.__out[fltrs.index(ext)])

                else:
                    errstr = [m for m in [swl, sfin] if m is not None]
                    self._failure(errstr[0], metric_name)
                    pass

                # if eclipse.whitelight exists
                if e_vwl and vfin:
                    log.info('--< ECLIPSE LIGHT CURVE SYMMETRY: %s >--', ext)

                    status = clscore.symwl(
                        self.__state_vecs['eclwhitelight'].sv_as_dict()[ext],
                        self.__state_vecs['finalize'].sv_as_dict()[
                            'parameters'
                        ]['priors'],
                        ext,
                        self.__out[len(fltrs) + fltrs.index(ext)],
                    )

                    if status:
                        svupdate.append(
                            self.__out[len(fltrs) + fltrs.index(ext)]
                        )
                    pass
                else:
                    errstr = [m for m in [e_swl, vfin] if m is not None]
                    self._failure(errstr[0], metric_name)
                    pass
                pass
                # ==================================================== #####

                # ======================  RSDM  ====================== #####
                if 'Spitzer' not in ext:
                    metric_name = "RSDM"

                    # if transit.spectrum exists
                    if vsp and vfin:  # Q need to check for vfin here?
                        log.info('--< IN-TRANSIT RSDM: %s >--', ext)

                        status = clscore.rsdm(
                            self.__state_vecs['spectrum'].sv_as_dict()[ext],
                            self.__out[fltrs.index(ext)],
                        )

                        if status:
                            svupdate.append(self.__out[fltrs.index(ext)])
                        pass
                    else:
                        errstr = [m for m in [ssp, sfin] if m is not None]
                        self._failure(errstr[0], metric_name)
                        pass

                    # if eclipse.spectrum exists
                    if e_vsp and vfin:
                        log.info('--< IN-ECLIPSE RSDM: %s >--', ext)

                        status = clscore.rsdm(
                            self.__state_vecs['eclspectrum'].sv_as_dict()[ext],
                            self.__out[len(fltrs) + fltrs.index(ext)],
                        )

                        if status:
                            svupdate.append(
                                self.__out[len(fltrs) + fltrs.index(ext)]
                            )
                        pass
                    else:
                        errstr = [m for m in [e_ssp, sfin] if m is not None]
                        self._failure(errstr[0], metric_name)
                        pass
                    pass
                pass
                # ==================================================== #####

                # ======================  PERC_REJECTED  ====================== #####
                if 'Spitzer' not in ext:
                    metric_name = "Percent Rejected"

                    # if transit.spectrum exists
                    if vsp and vfin:
                        log.info(
                            '--< IN-TRANSIT PERCENT REJECTED: %s >--', ext
                        )

                        status = clscore.perc_rejected(
                            self.__state_vecs['spectrum'].sv_as_dict()[ext],
                            self.__out[fltrs.index(ext)],
                        )

                        if status:
                            svupdate.append(self.__out[fltrs.index(ext)])
                        pass
                    else:
                        errstr = [m for m in [ssp, sfin] if m is not None]
                        self._failure(errstr[0], metric_name)
                        pass

                    # if eclipse.spectrum exists
                    if e_vsp and vfin:
                        log.info(
                            '--< IN-ECLIPSE PERCENT REJECTED: %s >--', ext
                        )

                        status = clscore.perc_rejected(
                            self.__state_vecs['eclspectrum'].sv_as_dict()[ext],
                            self.__out[len(fltrs) + fltrs.index(ext)],
                        )

                        if status:
                            svupdate.append(
                                self.__out[len(fltrs) + fltrs.index(ext)]
                            )
                        pass
                    else:
                        errstr = [m for m in [e_ssp, sfin] if m is not None]
                        self._failure(errstr[0], metric_name)
                        pass
                    pass
                pass
                # ==================================================== #####

                # ======================  MEDIAN_ERROR  ====================== #####
                if 'Spitzer' not in ext:
                    metric_name = "data.calibration Median Error"

                    if vdc and vfin:
                        log.info(
                            '--< DATA.CALIBRATION MEDIAN ERROR: %s >--', ext
                        )

                        status = clscore.median_error(
                            self.__state_vecs['data_calib'].sv_as_dict()[ext],
                            self.__out[fltrs.index(ext)],
                        )

                        if status:
                            svupdate.append(self.__out[fltrs.index(ext)])
                        pass
                    else:
                        errstr = [m for m in [sdc, sfin] if m is not None]
                        self._failure(errstr[0], metric_name)
                        pass
                    pass
                pass
                # ==================================================== #####

                # ===============  LC_RESIDUAL_CLASSIFICATION  ============== #####

                metric_name = "Light Curve Residual Shape"

                # if transit.whitelight exists
                if vwl and vfin:
                    log.info('--< IN-TRANSIT RESIDUAL SHAPE: %s >--', ext)

                    status = clscore.lc_resid_classification(
                        self.__state_vecs['whitelight'].sv_as_dict()[ext],
                        ext,
                        self.__out[fltrs.index(ext)],
                    )

                    if status:
                        svupdate.append(self.__out[fltrs.index(ext)])
                    pass
                else:
                    errstr = [m for m in [swl, sfin] if m is not None]
                    self._failure(errstr[0], metric_name)
                    pass

                # if eclipse.spectrum exists
                if e_vsp and vfin:
                    log.info('--< IN-ECLIPSE RESIDUAL SHAPE: %s >--', ext)

                    status = clscore.lc_resid_classification(
                        self.__state_vecs['eclwhitelight'].sv_as_dict()[ext],
                        ext,
                        self.__out[len(fltrs) + fltrs.index(ext)],
                    )

                    if status:
                        svupdate.append(
                            self.__out[len(fltrs) + fltrs.index(ext)]
                        )
                    pass
                else:
                    errstr = [m for m in [e_ssp, sfin] if m is not None]
                    self._failure(errstr[0], metric_name)
                    pass
                pass
                # =========================================================== #####

                # ======================  CALCULATE OVERALL PLANET FLAGS  ====================== #####
                flag_vals = {'red': 2, 'yellow': 1, 'green': 0}

                flag_colors = ['green', 'yellow', 'red']

                for planet in self.__out[fltrs.index(ext)]['data']:  # transit

                    if (
                        planet != "median_error"
                    ):  # median error is not planet-specific

                        planet_flag_val = 0

                        for alg in self.__out[fltrs.index(ext)]['data'][planet]:
                            flag_color = self.__out[fltrs.index(ext)]['data'][
                                planet
                            ][alg]['flag_color']
                            planet_flag_val = max(
                                flag_vals[flag_color], planet_flag_val
                            )

                        # factor in the flags that are not planet-specific (e.g. median_error)
                        try:
                            flag_color = self.__out[fltrs.index(ext)]['data'][
                                'median_error'
                            ]['flag_color']
                            planet_flag_val = max(
                                flag_vals[flag_color], planet_flag_val
                            )
                        except KeyError:
                            pass

                        planet_flag_color = flag_colors[planet_flag_val]

                        self.__out[fltrs.index(ext)]['data'][planet][
                            'overall_flag'
                        ] = planet_flag_color

                for planet in self.__out[len(fltrs) + fltrs.index(ext)][
                    'data'
                ]:  # eclipse

                    if (
                        planet != "median_error"
                    ):  # median error is not planet-specific

                        planet_flag_val = 0

                        for alg in self.__out[len(fltrs) + fltrs.index(ext)][
                            'data'
                        ][planet]:
                            flag_color = self.__out[
                                len(fltrs) + fltrs.index(ext)
                            ]['data'][planet][alg]['flag_color']
                            planet_flag_val = max(
                                flag_vals[flag_color], planet_flag_val
                            )

                        # factor in the flags that are not planet-specific (e.g. median_error)
                        try:
                            flag_color = self.__out[
                                len(fltrs) + fltrs.index(ext)
                            ]['data']['median_error']['flag_color']
                            planet_flag_val = max(
                                flag_vals[flag_color], planet_flag_val
                            )
                        except KeyError:
                            pass

                        planet_flag_color = flag_colors[planet_flag_val]

                        self.__out[len(fltrs) + fltrs.index(ext)]['data'][
                            planet
                        ]['overall_flag'] = planet_flag_color
                # ============================================================================ #####

            if self.__out:
                # print(self.__out)
                ds.update()
            else:
                log.warning(
                    '--< NO OUTPUT CREATED FOR CLASSIFIER.%s >--',
                    self.name().upper(),
                )
        return

    @staticmethod
    def _failure(errstr, metric_name):
        '''_failure ds'''
        #  don't print 'FAILED' for all the empty datasets
        # log.warning('--< FAILED %s CHECK: %s >--', str(metric_name).upper(), errstr)
        if errstr == 'garbagio':
            print('jenkins is a ' + errstr)
        errstr = metric_name
        return

    pass
