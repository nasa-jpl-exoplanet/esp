'''ancillary algoithms dc'''

# Heritage code shame:
# pylint: disable=too-many-branches

# -- IMPORTS -- ------------------------------------------------------
import dawgie
import dawgie.context

from collections import defaultdict

import logging

import excalibur.ancillary.core as anccore
import excalibur.ancillary.states as ancstates

import excalibur.runtime.algorithms as rtalg

import excalibur.system as sys
import excalibur.system.algorithms as sysalg

from excalibur.target.targetlists import get_target_lists

from excalibur.ancillary.core import savesv

from excalibur.util.checksv import checksv

from importlib import import_module as fetch  # avoid cicular dependencies

log = logging.getLogger(__name__)


# ---------------------- ---------------------------------------------
# -- ALGORITHMS -- ---------------------------------------------------
class Estimate(dawgie.Algorithm):
    '''estimate ds'''

    def __init__(self):
        '''__init__ ds'''
        self._version_ = anccore.estimateversion()
        self._type = 'ancillary'
        self.__rt = rtalg.Autofill()
        self.__fin = sysalg.Finalize()
        self.__out = [ancstates.EstimateSV('parameters')]
        return

    def name(self):
        '''name ds'''
        return 'estimate'

    def previous(self):
        '''previous ds'''
        return [
            dawgie.ALG_REF(sys.task, self.__fin)
        ] + self.__rt.refs_for_validity()

    def state_vectors(self):
        '''state_vectors ds'''
        return self.__out

    def run(self, ds, ps):
        '''run ds'''

        # stop here if it is not a runtime target
        if not self.__rt.is_valid():
            log.info(
                '--< ANCILLARY.%s: not a valid target >--', self.name().upper()
            )

        else:
            update = False
            vfin, _ = checksv(self.__fin.sv_as_dict()['parameters'])
            if vfin:
                update = anccore.estimate(
                    self.__fin.sv_as_dict()['parameters'], self.__out[0]
                )
            if update:
                ds.update()
            else:
                raise dawgie.NoValidOutputDataError(
                    f'No output created for ANCILLARY.{self.name()}'
                )

        return


class Population(dawgie.Analyzer):
    '''population ds'''

    def __init__(self):
        '''__init__ ds'''
        self._version_ = dawgie.VERSION(1, 0, 3)
        self.__est = Estimate()
        self.__out = ancstates.PopulationSV('statistics')
        return

    def previous(self):
        '''Input State Vectors: ancillary.estimate'''
        return [dawgie.ALG_REF(sys.task, self.__est)]

    def feedback(self):
        '''feedback ds'''
        return []

    def name(self):
        '''name ds'''
        return 'population'

    def traits(self) -> [dawgie.SV_REF, dawgie.V_REF]:
        '''traits ds'''
        return [
            dawgie.SV_REF(
                fetch('excalibur.ancillary').task,
                Estimate(),
                Estimate().state_vectors()[0],
            )
        ]

    def state_vectors(self):
        '''state_vectors ds'''
        return [self.__out]

    def run(self, aspects: dawgie.Aspect):
        '''run ds'''

        targetlists = get_target_lists()

        aspecttargets = []
        for a in aspects:
            aspecttargets.append(a)

        # group together values by attribute
        svname = 'ancillary.estimate.parameters'

        st_attrs = defaultdict(list)
        pl_attrs = defaultdict(list)
        # include a second set of attributes, for comparison within each histogram
        st_attrs_roudier62 = defaultdict(list)
        pl_attrs_roudier62 = defaultdict(list)

        for trgt in targetlists['active']:
            if trgt not in aspecttargets:
                log.warning(
                    '--< ANCILLARY.POPULATION: target missing: %s >--', trgt
                )

        # Only consider the 'active' stars, not aliases/misspellings/dropped/etc
        # for trgt in targetlists['active']:
        # for trgt in filter(lambda tgt: 'STATUS' in aspects[tgt][svname], targetlists['active']):
        for trgt in filter(
            lambda tgt: tgt in aspecttargets, targetlists['active']
        ):
            anc_data = aspects[trgt][svname]

            # verify SV succeeded for target
            if anc_data['STATUS'][-1] or 'planets' in anc_data['data']:
                # get stellar attributes
                #  only include the basic data, no extensions
                for key in anc_data['data'].keys():
                    if (
                        (not key == 'planets')
                        and (key not in anc_data['data']['planets'])
                        and (not any(ext in key for ext in anccore.SV_EXTS))
                    ):
                        st_attrs[key].append(anc_data['data'][key])
                # get planetary attributes
                for pl in anc_data['data']['planets']:
                    pl_keys = [
                        i
                        for i in anc_data['data'][pl].keys()
                        if not any(ext in i for ext in anccore.SV_EXTS)
                    ]
                    for key in pl_keys:
                        pl_attrs[key].append(anc_data['data'][pl][key])

        # Loop through a second group of targets.  (this subset will be overplotted in the histos)
        # for trgt in filter(lambda tgt: 'STATUS' in aspects[tgt][svname], targetlists['roudier62']):
        for trgt in filter(
            lambda tgt: tgt in aspecttargets, targetlists['roudier62']
        ):
            anc_data = aspects[trgt][svname]

            # verify SV succeeded for target
            if anc_data['STATUS'][-1] or 'planets' in anc_data['data']:
                # get stellar attributes
                #  only include the basic data, no extensions
                for key in anc_data['data'].keys():
                    if (
                        (not key == 'planets')
                        and (key not in anc_data['data']['planets'])
                        and (not any(ext in key for ext in anccore.SV_EXTS))
                    ):
                        st_attrs_roudier62[key].append(anc_data['data'][key])
                # get planetary attributes
                for pl in anc_data['data']['planets']:
                    pl_keys = [
                        i
                        for i in anc_data['data'][pl].keys()
                        if not any(ext in i for ext in anccore.SV_EXTS)
                    ]
                    for key in pl_keys:
                        pl_attrs_roudier62[key].append(
                            anc_data['data'][pl][key]
                        )

        # Add to SV
        self.__out['data']['st_attrs'] = st_attrs
        self.__out['data']['pl_attrs'] = pl_attrs
        self.__out['data']['st_attrs_roudier62'] = st_attrs_roudier62
        self.__out['data']['pl_attrs_roudier62'] = pl_attrs_roudier62
        self.__out['STATUS'].append(True)
        aspects.ds().update()

        # save ancillary-estimate results as a .csv file (in /proj/data/spreadsheets/)
        savesv(aspects, targetlists)

        return

    pass


# ---------------- ---------------------------------------------------
