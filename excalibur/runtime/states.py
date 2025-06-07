'''Runtime configuration products'''

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
        self['target_autofill_selectMostRecent'] = BoolValue()
        self['target_autofill_maximizeSelfConsistency'] = BoolValue()
        self['ariel_simulate_spectra_includeMetallicityDispersion'] = (
            BoolValue()
        )
        self['cerberus_atmos_sliceSampler'] = BoolValue()
        self['cerberus_atmos_fitCloudParameters'] = BoolValue()
        self['cerberus_atmos_fitT'] = BoolValue()
        self['cerberus_atmos_fitCtoO'] = BoolValue()
        self['cerberus_atmos_fitNtoO'] = BoolValue()
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
            ['Switch', 'State'],
            len(self) + 1,
            'Processing Control Switches and Other Parameters',
        )
        for row, key in enumerate(sorted(self)):
            table.get_cell(row + 1, 0).add_primitive(key)
            table.get_cell(row + 1, 1).add_primitive(
                'on' if self[key] else 'off'
            )
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
                'Mission-Platform-Instrument-Mode Filters',
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
        visitor.add_declaration_inline(
            f'PYMC for {self.__name} default chain '
            f'length: {self["default"].value()}',
            tag='b',
        )
        if self['overrides']:
            table = visitor.add_table(
                ['Target', 'Chainlength'],
                len(self['overrides']) + 1,
                f'Overrides for {self.__name}',
            )
            for row, tn in enumerate(sorted(self['overrides'])):
                table.get_cell(row + 1, 0).add_primitive(tn)
                table.get_cell(row + 1, 1).add_primitive(self['overrides'][tn])
        else:
            visitor.add_primitive(' with no overrides')
        visitor.add_declaration_inline('', div='</div>')
        return

    pass


class StatusSV(dawgie.StateVector):
    '''State representation of how the AE should view this target'''

    def __init__(self):
        '''init the state vector with empty values'''
        self._version_ = dawgie.VERSION(1, 0, 0)
        self['allowed_filter_names'] = excalibur.ValuesList()
        self['ariel_simulate_spectra_includeMetallicityDispersion'] = (
            BoolValue()
        )
        self['cerberus_atmos_fitCloudParameters'] = BoolValue()
        self['cerberus_atmos_fitNtoO'] = BoolValue()
        self['cerberus_atmos_fitCtoO'] = BoolValue()
        self['cerberus_atmos_fitT'] = BoolValue()
        self['cerberus_chains'] = excalibur.ValueScalar()
        self['cerberus_steps'] = excalibur.ValueScalar()
        self['cerberus_atmos_sliceSampler'] = BoolValue()
        self['isValidTarget'] = BoolValue()
        self['runTarget'] = BoolValue(True)
        self['spectrum_chains'] = excalibur.ValueScalar()
        self['spectrum_steps'] = excalibur.ValueScalar()
        self['target_autofill_selectMostRecent'] = BoolValue()
        self['target_autofill_maximizeSelfConsistency'] = BoolValue()

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
        visitor.add_declaration_inline('Sampling parameters for PYMC', tag='b')
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
            'Control switches/parameters and their state', tag='b'
        )
        switches = [
            'runTarget',
            'isValidTarget',
            'target_autofill_selectMostRecent',
            'target_autofill_maximizeSelfConsistency',
            'ariel_simulate_spectra_includeMetallicityDispersion',
            'cerberus_atmos_sliceSampler',
            'cerberus_atmos_fitT',
            'cerberus_atmos_fitCtoO',
            'cerberus_atmos_fitNtoO',
            'cerberus_atmos_fitCloudParameters',
        ]
        table = visitor.add_table(['Switch', 'State'], len(switches))
        for row, switch in enumerate(switches):
            table.get_cell(row, 0).add_primitive(switch)
            table.get_cell(row, 1).add_primitive(
                'on' if self[switch] else 'off'
            )
        visitor.add_declaration_inline('', div='</div>')
        visitor.add_declaration_inline('', div='<div><hr><ul>')
        visitor.add_declaration_inline(
            'Mission-Platform-Instrument-Mode filters to process', tag='b'
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
                title = 'Run ALL targets.'
        else:
            title = 'Never run these targets:'
            if not self['targets']:
                title += ' NONE'

        visitor.add_declaration_inline(title, tag='b')

        if self['targets']:
            table = visitor.add_table(
                ['Target', 'Why'], len(self['targets']) + 1, 'Target table'
            )
            for row, tn in enumerate(
                sorted(self['targets'], key=lambda t: t[0])
            ):
                table.get_cell(row + 1, 0).add_primitive(tn[0])
                table.get_cell(row + 1, 1).add_primitive(tn[1])
        visitor.add_declaration_inline('', div='</div>')
        return

    pass
