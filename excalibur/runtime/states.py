'''Runtime configuration products'''

# Heritage code shame:
#  no-member is for "Instance of HiLoValue has no _hi"
# pylint: disable=no-member,method-hidden,

import dawgie
import excalibur
import logging

log = logging.getLogger(__name__)


class BoolValue(dawgie.Value):
    '''helper value for boolean type'''

    def __bool__(self):
        '''allows class to be treated like boolean using its __state'''
        return self.__state

    def __init__(self, state: bool = False):
        '''init the boolean'''
        self.__state = bool(state)
        self._version_ = dawgie.VERSION(1, 0, 0)
        return

    def __str__(self):
        '''define the string format of this class'''
        return str(self.__state)

    def features(self):
        '''contains no features'''
        return []

    def new(self, state=None):
        '''hide explicit requirement for dawgie'''
        return BoolValue(state if state is not None else self.__state)

    pass


class HiLoValue(dawgie.Value):
    '''helper value for hi-lo type'''

    def __init__(self, hi: float = 1, lo: float = 0):
        '''init the hi-lo'''
        self.hi = hi
        self.lo = lo
        self._version_ = dawgie.VERSION(1, 0, 0)
        return

    def __str__(self):
        '''define the string format of this class'''
        # fails. HiLoValue object has no attribute '_HiLoValue__state'
        # return str(self.__state)
        return str(self.__getstate__())
        #    it now prints this:
        # "{'hi': 1.5, 'lo': 0.75, '_version_seal_': VERSION(design=1, impl=0, bugfix=0)}"

    def features(self):
        '''contains no features'''
        return []

    def new(self, hilo=None):
        '''hide explicit requirement for dawgie'''
        return HiLoValue(
            *((float(hilo.hi), float(hilo.lo)) if hilo else (1, 0))
        )

    def hi(self):
        return self._hi

    def lo(self):
        return self._lo


class MoleculeValue(dawgie.Value):
    '''helper value for molecule-list type'''

    def __init__(self, molecules: list = []):
        '''init the molecule list'''
        self.molecules = list(molecules)
        self._version_ = dawgie.VERSION(1, 0, 0)
        return

    def __str__(self):
        '''define the string format of this class'''
        return str(self.__getstate__())

    def features(self):
        '''contains no features'''
        return []

    def new(self, state=None):
        '''hide explicit requirement for dawgie'''
        # print('  SHOULD RETURN THIS',state.molecule)
        moleculeVals = MoleculeValue(
            state.molecule if state is not None else self.__state)

        if not isinstance(moleculeVals.molecules, list):
            log.error('ERROR: molecules should be a list!')
        else:
            # has to be pickleable; convert each molecule to a string
            print('NEW() before fix',moleculeVals)

            cleanVals = []
            for item in moleculeVals.molecules:
                print('  item fore',item,type(item))
                cleanVals.append(str(item))
                print('  item aftr',item,type(str(item)))
            moleculeVals.molecules = cleanVals

            print('NEW() after fix ',moleculeVals)

        return moleculeVals

    def molecules(self):
        return self._molecules


class CompositeSV(dawgie.StateVector):
    '''State representation of the configuration file'''

    def __init__(self, constituents: [dawgie.StateVector]):
        '''init the state vector with empty values'''
        self._version_ = dawgie.VERSION(1, 0, 0)
        for constituent in constituents:
            self[constituent.name()] = constituent
        return

    def name(self):
        '''database name'''
        return 'composite'

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''Show the configutation information'''
        for key in sorted(self):
            self[key].view(caller, visitor)
        return

    pass


class ControlsSV(dawgie.StateVector, dawgie.Value):
    '''State representation of the parameter controls'''

    def __init__(self):
        '''init the state vector with empty values'''
        self._version_ = dawgie.VERSION(1, 0, 0)
        self['system_validate_maximizeSelfConsistency'] = BoolValue()
        self['system_validate_selectMostRecent'] = BoolValue()
        self['ariel_simspectrum_thorngrenMassMetals'] = BoolValue()
        self['ariel_simspectrum_includeMetallicityDispersion'] = BoolValue()
        self['ariel_simspectrum_randomCloudProperties'] = BoolValue()
        self['ariel_simspectrum_tier'] = excalibur.ValueScalar()
        self['ariel_simspectrum_randomseed'] = excalibur.ValueScalar()
        self['ariel_simspectrum_SNRadjustment'] = excalibur.ValueScalar()
        self['ariel_simspectrum_metallicityDispersion'] = (
            excalibur.ValueScalar()
        )
        self['ariel_simspectrum_CtoOaverage'] = excalibur.ValueScalar()
        self['ariel_simspectrum_CtoOdispersion'] = excalibur.ValueScalar()
        self['cerberus_atmos_sliceSampler'] = BoolValue()
        self['cerberus_atmos_fitCloudParameters'] = BoolValue()
        self['cerberus_atmos_fitT'] = BoolValue()
        self['cerberus_atmos_fitCtoO'] = BoolValue()
        self['cerberus_atmos_fitNtoO'] = BoolValue()
        self['cerberus_crbmodel_isothermal'] = BoolValue()
        self['cerberus_crbmodel_lbroadening'] = BoolValue()
        self['cerberus_crbmodel_lshifting'] = BoolValue()
        self['cerberus_crbmodel_nlevels'] = excalibur.ValueScalar()
        self['cerberus_crbmodel_solrad'] = excalibur.ValueScalar()
        self['cerberus_crbmodel_Hsmax'] = excalibur.ValueScalar()
        self['cerberus_crbmodel_fitmolecules'] = MoleculeValue()
        self['cerberus_crbmodel_HITEMPmolecules'] = MoleculeValue()
        self['cerberus_crbmodel_HITRANmolecules'] = MoleculeValue()
        self['cerberus_crbmodel_EXOMOLmolecules'] = MoleculeValue()
        self['cerberus_atmos_bounds_Teq'] = HiLoValue()
        self['cerberus_atmos_bounds_abundances'] = HiLoValue()
        self['cerberus_atmos_bounds_CTP'] = HiLoValue()
        self['cerberus_atmos_bounds_HLoc'] = HiLoValue()
        self['cerberus_atmos_bounds_HScale'] = HiLoValue()
        self['cerberus_atmos_bounds_HThick'] = HiLoValue()
        self['cerberus_plotters_cornerBins'] = excalibur.ValueScalar()
        self['cerberus_results_randomseed'] = excalibur.ValueScalar()
        self['cerberus_results_nrandomwalkers'] = excalibur.ValueScalar()
        self['selftest_Nrepeats'] = excalibur.ValueScalar()
        return

    def features(self):
        '''contains no features'''
        return []

    def name(self):
        '''database name'''
        return 'controls'

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''Show the configutation information'''
        visitor.add_declaration_inline('', div='<div><hr>')
        table = visitor.add_table(
            # ['Switch', 'State'],
            ['Algorithm', 'Parameter', 'Value'],
            len(self) + 1,
            'Processing Control Parameters',
            # 'Processing Control Switches',
        )
        #  let's drop this alphabetical sorting and organize more chronologically
        # for row, key in enumerate(sorted(self)):
        for row, key in enumerate(self):
            isplitter = key.rfind('_')
            algorithm = key[:isplitter]
            param = key[isplitter + 1 :]
            table.get_cell(row + 1, 0).add_primitive(algorithm)
            table.get_cell(row + 1, 1).add_primitive(param)
            if isinstance(self[key], excalibur.ValueScalar):
                table.get_cell(row + 1, 2).add_primitive(self[key].value())
            else:
                table.get_cell(row + 1, 2).add_primitive(self[key])
            # table.get_cell(row + 1, 0).add_primitive(key)
            # table.get_cell(row + 1, 1).add_primitive(
            #     'on' if self[key] else 'off'
            # )
        visitor.add_declaration_inline('', div='</div>')
        return

    pass


class FilterSV(dawgie.StateVector, dawgie.Value):
    '''State representation of the filters to be included/excluded'''

    def __init__(self):
        '''init the state vector with empty values'''
        self._version_ = dawgie.VERSION(1, 0, 0)
        self['includes'] = excalibur.ValuesList()
        self['excludes'] = excalibur.ValuesList()
        return

    def features(self):
        '''contains no features'''
        return []

    def name(self):
        '''datebase name'''
        return 'filters'

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''Show the configutation information'''
        visitor.add_declaration_inline('', div='<div><hr>')
        table_len = max(len(self['excludes']), len(self['includes']))
        if table_len == 1:
            visitor.add_declaration_inline('No target filters set', tag='b')
        else:
            self['excludes'].sort()
            self['includes'].sort()
            table = visitor.add_table(
                ['Exclude', 'Include'],
                table_len + 1,
                'Mission/Instrument Filters',
            )
            for row in range(table_len):
                for col, filt in enumerate(['excludes', 'includes']):
                    if col < len(self[filt]):
                        content = (
                            self[filt][row] if row < len(self[filt]) else ' '
                        )
                        table.get_cell(row + 1, col).add_primitive(content)
        visitor.add_declaration_inline('', div='</div>')
        return

    pass


class PymcSV(dawgie.StateVector, dawgie.Value):
    '''State representation of the PYMC parameters'''

    def __init__(self, name: str = 'undefined'):
        '''init the state vector with empty values'''
        self.__name = name
        self._version_ = dawgie.VERSION(1, 0, 0)
        self['default'] = excalibur.ValueScalar()
        self['overrides'] = excalibur.ValuesDict()
        return

    def features(self):
        '''contains no features'''
        return []

    def name(self):
        '''database name'''
        return f'pymc-{self.__name}'

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''Show the configutation information'''
        visitor.add_declaration_inline('', div='<div><hr>')
        if (self.__name).endswith('chainlen'):
            paramname = 'Chain Length'
            algorithmname = (self.__name)[:-8]
        elif (self.__name).endswith('chains'):
            paramname = '# of Chains'
            algorithmname = (self.__name)[:-6]
        else:
            paramname = 'Value'
            algorithmname = self.__name
        visitor.add_declaration_inline(
            f'PYMC in {algorithmname}: default {paramname} = '
            f'{self["default"].value()}',
            tag='b',
        )
        if self['overrides']:
            table = visitor.add_table(
                ['Target', paramname],
                len(self['overrides']) + 1,
                'Overrides:',
            )
            for row, tn in enumerate(sorted(self['overrides'])):
                table.get_cell(row + 1, 0).add_primitive(tn)
                table.get_cell(row + 1, 1).add_primitive(self['overrides'][tn])
        else:
            visitor.add_primitive(' Overrides:  None')
        visitor.add_declaration_inline('', div='</div>')
        return

    pass


class StatusSV(dawgie.StateVector):
    '''State representation of how the AE should view this target'''

    def __init__(self):
        '''init the state vector with empty values'''
        self._version_ = dawgie.VERSION(1, 0, 0)
        self['allowed_filter_names'] = excalibur.ValuesList()
        self['ariel_simspectrum_includeMetallicityDispersion'] = BoolValue()
        self['ariel_simspectrum_randomCloudProperties'] = BoolValue()
        self['ariel_simspectrum_thorngrenMassMetals'] = BoolValue()
        self['ariel_simspectrum_tier'] = excalibur.ValueScalar()
        self['ariel_simspectrum_randomseed'] = excalibur.ValueScalar()
        self['ariel_simspectrum_SNRadjustment'] = excalibur.ValueScalar()
        self['ariel_simspectrum_metallicityDispersion'] = (
            excalibur.ValueScalar()
        )
        self['ariel_simspectrum_CtoOaverage'] = excalibur.ValueScalar()
        self['ariel_simspectrum_CtoOdispersion'] = excalibur.ValueScalar()
        self['cerberus_atmos_fitCloudParameters'] = BoolValue()
        self['cerberus_atmos_fitNtoO'] = BoolValue()
        self['cerberus_atmos_fitCtoO'] = BoolValue()
        self['cerberus_atmos_fitT'] = BoolValue()
        self['cerberus_crbmodel_lbroadening'] = BoolValue()
        self['cerberus_crbmodel_lshifting'] = BoolValue()
        self['cerberus_crbmodel_isothermal'] = BoolValue()
        self['cerberus_crbmodel_nlevels'] = excalibur.ValueScalar()
        self['cerberus_crbmodel_solrad'] = excalibur.ValueScalar()
        self['cerberus_crbmodel_Hsmax'] = excalibur.ValueScalar()
        self['cerberus_crbmodel_fitmolecules'] = MoleculeValue()
        self['cerberus_crbmodel_HITEMPmolecules'] = MoleculeValue()
        self['cerberus_crbmodel_HITRANmolecules'] = MoleculeValue()
        self['cerberus_crbmodel_EXOMOLmolecules'] = MoleculeValue()
        self['cerberus_atmos_bounds_Teq'] = HiLoValue()
        self['cerberus_atmos_bounds_abundances'] = HiLoValue()
        self['cerberus_atmos_bounds_CTP'] = HiLoValue()
        self['cerberus_atmos_bounds_HLoc'] = HiLoValue()
        self['cerberus_atmos_bounds_HScale'] = HiLoValue()
        self['cerberus_atmos_bounds_HThick'] = HiLoValue()
        self['cerberus_plotters_cornerBins'] = excalibur.ValueScalar()
        self['cerberus_results_randomseed'] = excalibur.ValueScalar()
        self['cerberus_results_nrandomwalkers'] = excalibur.ValueScalar()
        self['cerberus_chains'] = excalibur.ValueScalar()
        self['cerberus_steps'] = excalibur.ValueScalar()
        self['cerberus_atmos_sliceSampler'] = BoolValue()
        self['isValidTarget'] = BoolValue()
        self['runTarget'] = BoolValue(True)
        self['spectrum_chains'] = excalibur.ValueScalar()
        self['spectrum_steps'] = excalibur.ValueScalar()
        self['system_validate_selectMostRecent'] = BoolValue()
        self['system_validate_maximizeSelfConsistency'] = BoolValue()
        self['selftest_Nrepeats'] = excalibur.ValueScalar()

    def name(self):
        '''database name'''
        return 'status'

    def proceed(self, ext: str = None):
        '''determine if those that care should proceed'''
        allowed = ext in self['allowed_filter_names'] if ext else True
        run = self['runTarget']
        valid = self['isValidTarget']
        if not all([allowed, run, valid]):
            msg = (
                'Determined that should not process this target for ext:\n'
                f'  validTarget:    {valid}\n'
                f'  runnable:       {run}\n'
                f'  ext is allowed: {allowed}'
            )
            log.info(msg)
            raise dawgie.NoValidInputDataError(msg)

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''Show the state for this target'''
        visitor.add_declaration_inline(
            'Briefly: ', div='<div><hr><h3>', tag='span'
        )
        try:
            self.proceed()
            visitor.add_declaration_inline(
                'process', style='color:green', tag='span'
            )
        except dawgie.NoValidInputDataError:
            visitor.add_declaration_inline(
                'DO NOT process', style='color:red', tag='span'
            )
        visitor.add_declaration_inline(' this target', tag='span')
        visitor.add_declaration_inline('', div='</h3></div>')
        visitor.add_declaration_inline('', div='<div><hr>')
        visitor.add_declaration_inline('PYMC Sampling Parameters', tag='b')
        table = visitor.add_table(
            ['Algorithm', 'Sampler', '# of chains', 'Chain length'], 2
        )
        table.get_cell(0, 0).add_primitive('cerberus')
        table.get_cell(0, 1).add_primitive('metropolis-hastings')
        table.get_cell(0, 2).add_primitive(self['cerberus_chains'].value())
        table.get_cell(0, 3).add_primitive(self['cerberus_steps'].value())
        table.get_cell(1, 0).add_primitive('spectrum')
        table.get_cell(1, 1).add_primitive('metropolis-hastings')
        table.get_cell(1, 2).add_primitive(self['spectrum_chains'].value())
        table.get_cell(1, 3).add_primitive(self['spectrum_steps'].value())
        visitor.add_declaration_inline('', div='</div>')
        visitor.add_declaration_inline('', div='<div><hr>')
        visitor.add_declaration_inline(
            'Processing Control Parameters',
            tag='b',
            # 'Control switches/parameters and their state', tag='b'
        )
        switches = [
            'runTarget',
            'isValidTarget',
            'system_validate_selectMostRecent',
            'system_validate_maximizeSelfConsistency',
            'ariel_simspectrum_includeMetallicityDispersion',
            'ariel_simspectrum_randomCloudProperties',
            'ariel_simspectrum_thorngrenMassMetals',
            'ariel_simspectrum_tier',
            'ariel_simspectrum_randomseed',
            'ariel_simspectrum_SNRadjustment',
            'ariel_simspectrum_metallicityDispersion',
            'ariel_simspectrum_CtoOaverage',
            'ariel_simspectrum_CtoOdispersion',
            'cerberus_atmos_sliceSampler',
            'cerberus_atmos_fitT',
            'cerberus_atmos_fitCtoO',
            'cerberus_atmos_fitNtoO',
            'cerberus_atmos_fitCloudParameters',
            'cerberus_crbmodel_lbroadening',
            'cerberus_crbmodel_lshifting',
            'cerberus_crbmodel_isothermal',
            'cerberus_crbmodel_nlevels',
            'cerberus_crbmodel_solrad',
            'cerberus_crbmodel_Hsmax',
            'cerberus_crbmodel_fitmolecules',
            'cerberus_crbmodel_HITEMPmolecules',
            'cerberus_crbmodel_HITRANmolecules',
            'cerberus_crbmodel_EXOMOLmolecules',
            'cerberus_atmos_bounds_Teq',
            'cerberus_atmos_bounds_abundances',
            'cerberus_atmos_bounds_CTP',
            'cerberus_atmos_bounds_HLoc',
            'cerberus_atmos_bounds_HScale',
            'cerberus_atmos_bounds_HThick',
            'cerberus_plotters_cornerBins',
            'cerberus_results_nrandomwalkers',
            'cerberus_results_randomseed',
            'selftest_Nrepeats',
        ]
        table = visitor.add_table(['Switch', 'State'], len(switches))
        for row, switch in enumerate(switches):
            table.get_cell(row, 0).add_primitive(switch)
            # table.get_cell(row, 1).add_primitive(
            #    'on' if self[switch] else 'off'
            if isinstance(self[switch], excalibur.ValueScalar):
                table.get_cell(row, 1).add_primitive(self[switch].value())
            else:
                table.get_cell(row, 1).add_primitive(self[switch])
        visitor.add_declaration_inline('', div='</div>')
        visitor.add_declaration_inline('', div='<div><hr><ul>')
        visitor.add_declaration_inline(
            'Mission/Instrument Filters to Process', tag='b'
        )
        for name in self['allowed_filter_names']:
            visitor.add_declaration_inline(name, tag='li')
        visitor.add_declaration_inline('', div='</ul></div>')
        return

    pass


class TargetsSV(dawgie.StateVector, dawgie.Value):
    '''State representation of the targets to sequester'''

    def __init__(self, name: str = 'undefined'):
        '''init the state vector with empty values'''
        self._name = name
        self._version_ = dawgie.VERSION(1, 0, 0)
        self['targets'] = excalibur.ValuesList()
        return

    def features(self):
        '''contains no features'''
        return []

    def name(self):
        '''database name'''
        return self._name

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''Show the configuration information'''
        visitor.add_declaration_inline('', div='<div><hr>')
        if self._name == 'run_only':
            title = 'Run only these targets:'
            if not self['targets']:
                title = 'Run ALL targets'
        else:
            title = 'Never run these targets:'
            if not self['targets']:
                title += ' NONE'

        visitor.add_declaration_inline(title, tag='b')

        if self['targets']:
            table = visitor.add_table(
                ['Target', 'Why'], len(self['targets']) + 1, ''
            )
            for row, tn in enumerate(
                sorted(self['targets'], key=lambda t: t[0])
            ):
                table.get_cell(row + 1, 0).add_primitive(tn[0])
                table.get_cell(row + 1, 1).add_primitive(tn[1])
        visitor.add_declaration_inline('', div='</div>')
        return

    pass
