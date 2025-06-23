'''system algorithms ds'''

# Heritage code shame:
# pylint: disable=too-many-branches,too-many-locals

# -- IMPORTS -- ------------------------------------------------------
import logging

import dawgie

from collections import defaultdict

import excalibur
import excalibur.runtime.algorithms as rtalg

import excalibur.system.core as syscore
import excalibur.system.states as sysstates
import excalibur.system.overwriter as sysoverwriter

from excalibur.util.checksv import checksv

import excalibur.target as trg
import excalibur.target.algorithms as trgalg

from excalibur.system.consistency import consistency_checks
from excalibur.target.targetlists import get_target_lists
from excalibur.system.core import savesv

from importlib import import_module as fetch  # avoid cicular dependencies

log = logging.getLogger(__name__)


# ------------- ------------------------------------------------------
# -- ALGORITHMS -- ---------------------------------------------------
class Validate(dawgie.Algorithm):
    '''Pulls out formatted info from NEXSCI and checks what is missing'''

    def __init__(self):
        '''__init__ ds'''
        # self._version_ = dawgie.VERSION(1,1,4)  # typos fixed in core/ssconstants
        self._version_ = dawgie.VERSION(
            1, 2, 0
        )  # new bestValue() parameter selection
        self.__autofill = trgalg.Autofill()
        self.__rt = rtalg.Autofill()
        self.__out = sysstates.PriorsSV('parameters')
        return

    def name(self):
        '''Database name for subtask extension'''
        return 'validate'

    def previous(self):
        '''Input State Vectors: target.autofill'''
        return [
            dawgie.ALG_REF(trg.task, self.__autofill)
        ] + self.__rt.refs_for_validity()

    def state_vectors(self):
        '''Output State Vectors: system.validate'''
        return [self.__out]

    def run(self, ds, ps):
        '''Top level algorithm call'''

        # stop here if it is not a runtime target
        if not self.__rt.is_valid():
            log.warning(
                '--< SYSTEM.%s: not a valid target >--', self.name().upper()
            )

        else:
            autofill = self.__autofill.sv_as_dict()['parameters']
            runtime = self.__rt.sv_as_dict()['status']

            runtime_params = syscore.SystemParams(
                maximizeSelfConsistency=runtime[
                    'system_validate_maximizeSelfConsistency'
                ],
                selectMostRecent=runtime['system_validate_selectMostRecent'],
            )

            update = False
            valid, errstring = checksv(autofill)
            if valid:
                update = self._validate(autofill, runtime_params, self.__out)
            else:
                self._failure(errstring)

            if update:
                _ = excalibur.lagger()
                ds.update()
                pass
            elif valid:
                raise dawgie.NoValidOutputDataError(
                    f'No output created for SYSTEM.{self.name()}'
                )
        return

    @staticmethod
    def _validate(autofill, runtime_params, out):
        '''Core code call'''
        afilled = syscore.buildsp(autofill, runtime_params, out)
        return afilled

    @staticmethod
    def _failure(errstr):
        '''Failure log'''
        log.warning('--< SYSTEM VALIDATE: %s >--', errstr)
        return

    pass


class Finalize(dawgie.Algorithm):
    '''Generates SV used as input in other algorithm calls,
    throws out incomplete systems'''

    def __init__(self):
        self._version_ = dawgie.VERSION(1, 1, 4)
        self.__rt = rtalg.Autofill()
        self.__val = Validate()
        self.__out = sysstates.PriorsSV('parameters')
        return

    def name(self):
        '''Database name for subtask extension'''
        return 'finalize'

    def previous(self):
        '''Input State Vectors: system.validate'''
        return [
            dawgie.ALG_REF(fetch('excalibur.system').task, self.__val)
        ] + self.__rt.refs_for_validity()

    def state_vectors(self):
        '''Output State Vectors: system.finalize'''
        return [self.__out]

    def run(self, ds, ps):
        '''Top level algorithm call'''

        # FIXMEE: this code needs repaired by moving out to config
        target = repr(self).split('.')[1]

        # stop here if it is not a runtime target
        if not self.__rt.is_valid():
            log.warning(
                '--< SYSTEM.%s: %s not a valid target >--',
                target,
                self.name().upper(),
            )

        else:
            update = False
            val = self.__val.sv_as_dict()['parameters']

            valid, errstring = checksv(val)
            if valid:
                overwrite = sysoverwriter.ppar()
                for key in val:
                    self.__out[key] = val.copy()[key]
                if target.startswith('test'):
                    overwritetargetname = target[:-3]
                else:
                    overwritetargetname = target
                if overwritetargetname in overwrite:
                    update = self._priority(overwrite[overwritetargetname],
                                            self.__out)
                    if not update:
                        log.warning('>-- STILL MISSING DICT INFO')
                        log.warning('>-- ADD MORE KEYS TO SYSTEM/OVERWRITER')
                    pass
                elif not self.__out['PP'][-1]:
                    update = True
                else:
                    log.warning('>-- MISSING DICT INFO: %s --<', target)
                    log.warning('>-- ADD KEY TO SYSTEM/OVERWRITER --<')

                # consistency checks
                inconsistencies = consistency_checks(
                    self.__out['priors'], self.__out['ignore']
                )
                for inconsistency in inconsistencies:
                    self.__out['autofill'].append(
                        'inconsistent:' + inconsistency
                    )

                # log warnings moved to the very end (previously were before forcepar)
                # 6/16/24 target name added to log, otherwise can't tell which one has the error
                log.warning(
                    '>-- FORCE PARAMETER: %s %s',
                    target,
                    str(self.__out['PP'][-1]),
                )
                log.warning(
                    '>-- MISSING MANDATORY PARAMETERS: %s %s',
                    target,
                    str(self.__out['needed']),
                )
                log.warning(
                    '>-- MISSING PLANET PARAMETERS: %s %s',
                    target,
                    str(self.__out['pneeded']),
                )
                log.warning(
                    '>-- PLANETS IGNORED: %s %s',
                    target,
                    str(self.__out['ignore']),
                )
                log.warning(
                    '>-- INCONSISTENCIES: %s %s', target, str(inconsistencies)
                )
                log.warning(
                    '>-- AUTOFILL: %s %s', target, str(self.__out['autofill'])
                )
                pass
            else:
                self._failure(errstring)

            if update:
                _ = excalibur.lagger()
                ds.update()
                pass
            elif valid:
                raise dawgie.NoValidOutputDataError(
                    f'No output created for SYSTEM.{self.name()}'
                )
        return

    @staticmethod
    def _priority(overwrite, out):
        '''Core code call'''
        ffill = syscore.forcepar(overwrite, out)
        return ffill

    @staticmethod
    def _failure(errstr):
        '''Failure log'''
        log.warning('--< SYSTEM FINALIZE: %s >--', errstr)
        return

    pass


# ---------------- ---------------------------------------------------


class Population(dawgie.Analyzer):
    '''population ds'''

    def __init__(self):
        '''__init__ ds'''
        self._version_ = dawgie.VERSION(1, 0, 3)
        self.__fin = Finalize()
        self.__out = sysstates.PopulationSV('statistics')
        return

    def previous(self):
        '''Input State Vectors: system.finalize'''
        return [dawgie.ALG_REF(fetch('excalibur.system').task, self.__fin)]

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
                fetch('excalibur.system').task,
                Finalize(),
                Finalize().state_vectors()[0],
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

        # cross-check the 'active' list against __all__
        for target in targetlists['active']:
            if target not in aspecttargets:
                log.warning(
                    '--< SYSTEM.POPULATION: target missing: %s >--', target
                )
        # this prints out all the excluded targets.  don't bother
        # for target in aspecttargets:
        #    if target not in targetlists['active']:
        #        log.warning('--< SYSTEM.POPULATION: add to active list: %s >--', target)

        # group together values by attribute
        svname = 'system.finalize.parameters'

        st_attrs = defaultdict(list)
        pl_attrs = defaultdict(list)
        st_attrs_roudier62 = defaultdict(list)
        pl_attrs_roudier62 = defaultdict(list)

        # for trgt in aspects:
        # for trgt in targetlists['active']:
        # for trgt in filter(lambda tgt: 'STATUS' in aspects[tgt][svname], targetlists['active']):
        for trgt in filter(
            lambda tgt: tgt in aspecttargets, targetlists['active']
        ):

            system_data = aspects[trgt][svname]

            # verify SV succeeded for target
            if system_data['STATUS'][-1]:
                # get stellar attributes
                st_keys = system_data['starmdt']  # mandatory params
                pl_keys = system_data['planetmdt']

                for key in st_keys:
                    st_attrs[key].append(system_data['priors'][key])

                # get planetary attributes
                for planet_letter in system_data['priors']['planets']:
                    for key in pl_keys:
                        pl_attrs[key].append(
                            system_data['priors'][planet_letter][key]
                        )

        ntarg = 0
        nplanet = 0
        #  *** second loop for second (overlapping) histogram ***
        # for trgt in targetlists['roudier62']:
        # for trgt in filter(lambda tgt: 'STATUS' in aspects[tgt][svname], targetlists['roudier62']):
        for trgt in filter(
            lambda tgt: tgt in aspecttargets, targetlists['roudier62']
        ):
            system_data = aspects[trgt][svname]
            ntarg += 1

            # verify SV succeeded for target
            if system_data['STATUS'][-1]:
                st_keys = system_data['starmdt']  # mandatory params
                # st_keys.extend(system_data['starnonmdt'])  # non-mandatory params
                pl_keys = system_data['planetmdt']
                # pl_keys.extend(system_data['planetnonmdt'])

                # get stellar attributes
                for key in st_keys:
                    st_attrs_roudier62[key].append(system_data['priors'][key])
                # get planetary attributes
                for planet_letter in system_data['priors']['planets']:
                    nplanet += 1
                    for key in pl_keys:
                        pl_attrs_roudier62[key].append(
                            system_data['priors'][planet_letter][key]
                        )
        # print('system population: # of stars  ',ntarg)
        # print('system population: # of planets',nplanet)
        # Add to SV
        self.__out['data']['st_attrs'] = st_attrs
        self.__out['data']['st_attrs_roudier62'] = st_attrs_roudier62
        self.__out['data']['pl_attrs'] = pl_attrs
        self.__out['data']['pl_attrs_roudier62'] = pl_attrs_roudier62
        self.__out['STATUS'].append(True)
        aspects.ds().update()

        # save system-finalize results as .csv file (in /proj/data/spreadsheets/)
        savesv(aspects, targetlists)

        return

    pass


# ---------------- ---------------------------------------------------
