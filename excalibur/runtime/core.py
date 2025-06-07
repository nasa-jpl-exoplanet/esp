'''science functionality separated from dawgie'''

import logging

import os
import re

from . import binding

log = logging.getLogger(__name__)

ENV_NAME = 'EXCALIBUR_LEVER_AND_KNOB_SETTINGS'


def _sequester2sv(sequester_type, sv, targets):
    for tn in sequester_type.target:
        if tn.isRegex:
            regex = re.compile(tn.value())
            matching_targets = filter(
                lambda t, rex=regex: rex.match(t), targets
            )
            if matching_targets:
                for mt in matching_targets:
                    sv['targets'].append((mt, tn.because))
            else:
                log.warning(
                    'Sequester regex %s produced no matches', tn.value()
                )
        elif tn.value() in targets:
            sv['targets'].append((tn.value(), tn.because))
        else:
            log.warning('Sequester target %s is not known', tn.value())


def isolate(sv: {}, table: {str: {}}, tn: str) -> None:
    '''isolate target specific state from the global table'''
    if table['filters']['includes']:
        allowed_names = table['filters']['includes']
    else:
        allowed_names = binding.filter_names.itervalues()
    # make a copy of unique names ditching the old types along the way
    allowed_names = set(allowed_names)
    for exclude in table['filters']['excludes']:
        allowed_names.discard(exclude)
    sv['allowed_filter_names'].extend(allowed_names)
    for key in [
        'ariel_simulate_spectra_includeMetallicityDispersion',
        'cerberus_atmos_fitCloudParameters',
        'cerberus_atmos_fitNtoO',
        'cerberus_atmos_fitCtoO',
        'cerberus_atmos_fitT',
        'cerberus_atmos_sliceSampler',
        'target_autofill_selectMostRecent',
        'target_autofill_maximizeSelfConsistency',
    ]:
        sv[key] = table['controls'][key].new()
    pymc = table['pymc-cerberuschainlen']
    default = pymc['default'].value()
    sv['cerberus_steps'] = sv['cerberus_steps'].new(
        pymc['overrides'].get(tn, default)
    )
    # asdf
    #    pymc = table['pymc-cerberuschains']
    #    default = pymc['default'].value()
    #    sv['cerberus_chains'] = sv['cerberus_chains'].new(
    #        pymc['overrides'].get(tn, default)
    #    )
    sv['isValidTarget'] = sv['isValidTarget'].new(
        tn
        not in [
            targetandreason[0]
            for targetandreason in table['sequester']['targets']
        ]
    )
    if table['run_only']['targets']:
        sv['runTarget'] = sv['runTarget'].new(
            tn
            in [
                targetandreason[0]
                for targetandreason in table['run_only']['targets']
            ]
        )
    pymc = table['pymc-spectrumchainlen']
    default = pymc['default'].value()
    sv['spectrum_steps'] = sv['spectrum_steps'].new(
        pymc['overrides'].get(tn, default)
    )


# asdf
#    pymc = table['pymc-spectrumchains']
#    default = pymc['default'].value()
#    sv['spectrum_chains'] = sv['spectrum_chains'].new(
#        pymc['overrides'].get(tn, default)
#    )


def load(sv_dict: {str: {}}, targets) -> None:
    '''load the configuation file into state vectors

    1. read the filename from environment variable
    2. load the file and treat it as XML
    3. have pybgen binding module parse the XML
    4. move data into state vectors
    '''
    fn = os.environ[ENV_NAME]
    with open(fn, 'rt', encoding='utf-8') as file:
        xml = file.read()
    settings = binding.CreateFromDocument(xml)
    controls = sv_dict['controls']
    for knob in controls:
        controls[knob] = controls[knob].new(getattr(settings.controls, knob))
    sv_dict['filters']['excludes'].extend(
        [str(s) for s in settings.filters.exclude]
    )
    sv_dict['filters']['includes'].extend(
        [str(s) for s in settings.filters.include]
    )
    #    for pymc in ['cerberus-chains', 'cerberus-chainlen',
    #                 'spectrum-chains', 'spectrum-chainlen']:
    # asdf
    #    for pymc in ['cerberus-chainlen', 'spectrum-chainlen']:
    for pymc in ['cerberuschainlen', 'spectrumchainlen']:
        cf = getattr(settings.pymc, pymc)
        sv = sv_dict[f'pymc-{pymc}']
        sv['default'] = sv['default'].new(cf.default)
        for override in cf.target:
            sv['overrides'][override.name] = override.steps
    _sequester2sv(settings.run_only, sv_dict['run_only'], targets)
    _sequester2sv(settings.sequester, sv_dict['sequester'], targets)
