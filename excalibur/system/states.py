'''System Database Products View'''

# Heritage code shame:
# pylint: disable=too-many-locals

# -- IMPORTS -- ------------------------------------------------------
import dawgie

import excalibur
from excalibur.util.plotters import distrplot
from excalibur.util.svs import ExcaliburSV


# ------------- ------------------------------------------------------
# -- SV -- -----------------------------------------------------------
class PriorsSV(ExcaliburSV):
    '''General format for system State Vector view'''

    def __init__(self, name):
        ExcaliburSV.__init__(self, name, dawgie.VERSION(1, 1, 4))
        self['PP'] = excalibur.ValuesList()
        self['priors'] = excalibur.ValuesDict()
        self['pignore'] = excalibur.ValuesDict()
        self['ignore'] = excalibur.ValuesList()
        self['autofill'] = excalibur.ValuesList()
        self['needed'] = excalibur.ValuesList()
        self['pneeded'] = excalibur.ValuesList()
        self['starmdt'] = excalibur.ValuesList()
        self['starnonmdt'] = excalibur.ValuesList()
        self['planetmdt'] = excalibur.ValuesList()
        self['planetnonmdt'] = excalibur.ValuesList()
        self['starkeys'] = excalibur.ValuesList()
        self['planetkeys'] = excalibur.ValuesList()
        self['exts'] = excalibur.ValuesList()
        self['PP'].append(False)
        return

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''view ds'''
        if self['STATUS'][-1]:
            vlabels = [
                'FORCE PARAMETER',
                'MISSING MANDATORY PARAMETERS',
                'MISSING PLANET PARAMETERS',
                'PLANETS IGNORED',
                'AUTOFILL',
            ]
            hlabels = ['/', 'VALUE']
            table = visitor.add_table(clabels=hlabels, rows=len(vlabels))
            table.get_cell(0, 0).add_primitive(vlabels[0])
            table.get_cell(0, 1).add_primitive(self['PP'][-1])
            table.get_cell(1, 0).add_primitive(vlabels[1])
            table.get_cell(1, 1).add_primitive(self['needed'])
            table.get_cell(2, 0).add_primitive(vlabels[2])
            table.get_cell(2, 1).add_primitive(self['pneeded'])
            table.get_cell(3, 0).add_primitive(vlabels[3])
            table.get_cell(3, 1).add_primitive(self['ignore'])
            table.get_cell(4, 0).add_primitive(vlabels[4])
            table.get_cell(4, 1).add_primitive(self['autofill'])

            allstar = []
            for key in self['starmdt']:
                listkeys = [key]
                listkeys.extend([key + x for x in self['exts']])
                allstar.append(listkeys)
                pass
            pkeys = self['planetmdt']
            allplanet = []
            for key in pkeys:
                listkeys = [key]
                listkeys.extend([key + x for x in self['exts']])
                allplanet.append(listkeys)
                pass
            labels = ['STAR', 'UPPER ERR', 'LOWER ERR', 'UNITS', 'REF']
            table = visitor.add_table(clabels=labels, rows=len(allstar))
            for starlabels in allstar:
                i = allstar.index(starlabels)
                for starlabel in starlabels:
                    table.get_cell(
                        i, starlabels.index(starlabel)
                    ).add_primitive(starlabel)
                    table.get_cell(
                        i, starlabels.index(starlabel)
                    ).add_primitive(self['priors'][starlabel])
                    pass
                pass
            for planet_letter in self['priors']['planets']:
                labels = [
                    'PLANET ' + planet_letter,
                    'UPPER ERR',
                    'LOWER ERR',
                    'UNITS',
                    'REF',
                ]
                table = visitor.add_table(clabels=labels, rows=len(allplanet))
                for starlabels in allplanet:
                    i = allplanet.index(starlabels)
                    for starlabel in starlabels:
                        table.get_cell(
                            i, starlabels.index(starlabel)
                        ).add_primitive(starlabel)
                        table.get_cell(
                            i, starlabels.index(starlabel)
                        ).add_primitive(self['priors'][c][starlabel])

        return


# -------- -----------------------------------------------------------
class PopulationSV(ExcaliburSV):
    '''PopulationSV ds'''

    def __init__(self, name):
        ExcaliburSV.__init__(self, name, dawgie.VERSION(2, 0, 0))

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''view ds'''
        to_process = [
            (
                '----------------------Stellar Population Distributions----------------------',
                self['data']['st_attrs'],
                self['data']['st_attrs_roudier62'],
                False,
            ),
            (
                '---------------------Planetary Population Distributions---------------------',
                self['data']['pl_attrs'],
                self['data']['pl_attrs_roudier62'],
                True,
            ),
        ]
        # for title, attrs, attrs_roudier62,is_planet in to_process:
        for title, attrs, attrs_roudier62, _ in to_process:
            visitor.add_primitive(title)
            for key in attrs:
                distrplot(
                    key, attrs[key], attrs_roudier62[key], visitor, 'no units?'
                )

        return


# -------------------------------------------------------------------
