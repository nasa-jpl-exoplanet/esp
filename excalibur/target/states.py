'''Target Database Products View'''

# Heritage code shame:
# pylint: disable=too-many-branches,too-many-locals,too-many-nested-blocks,too-many-statements

# -- IMPORTS -- ------------------------------------------------------
import bokeh.embed
import bokeh.plotting  # the awesome plotting engine
import dawgie
import pandas as pd

import excalibur


# ------------- ------------------------------------------------------
# -- TARGET -- -------------------------------------------------------
class TargetSV(dawgie.StateVector):
    '''target view'''

    def __init__(self, name):
        '''__init__ ds'''
        self._version_ = dawgie.VERSION(1, 1, 1)
        self.__name = name
        self['STATUS'] = excalibur.ValuesList()
        self['starID'] = excalibur.ValuesDict()
        self['nexsciDefaults'] = excalibur.ValuesList()
        self['nexsciFulltable'] = excalibur.ValuesList()
        self['candidates'] = excalibur.ValuesList()
        self['starkeys'] = excalibur.ValuesList()
        self['planetkeys'] = excalibur.ValuesList()
        self['exts'] = excalibur.ValuesList()
        self['STATUS'].append(False)
        return

    def name(self):
        '''name ds'''
        return self.__name

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''view ds'''
        if self['STATUS'][-1]:
            targetlist = list(self['starID'].keys())
            targetlist = sorted(targetlist)
            ntarget = len(targetlist)
            labels = ['Star', 'Aliases', 'Planet', 'Proposal']
            table = visitor.add_table(clabels=labels, rows=ntarget)
            for target in targetlist:
                starinfo = self['starID'][target]
                i = targetlist.index(target)
                table.get_cell(i, 0).add_primitive(target)
                for alias in starinfo['aliases']:
                    table.get_cell(i, 1).add_primitive(alias)
                    pass
                for planet in starinfo['planets']:
                    table.get_cell(i, 2).add_primitive(planet)
                    pass
                for pid in starinfo['PID']:
                    table.get_cell(i, 3).add_primitive(pid)
                    pass
                pass
            starinfo = self['starID'][targetlist[0]]
            if len(self['STATUS']) > 2:
                allstar = []
                skeys = self['starkeys']
                exts = self['exts']
                for key in skeys:
                    listkeys = [key]
                    listkeys.extend([key + x for x in exts])
                    allstar.append(listkeys)
                    pass
                allplanet = []
                pkeys = self['planetkeys']
                for key in pkeys:
                    listkeys = [key]
                    listkeys.extend([key + x for x in exts])
                    allplanet.append(listkeys)
                    pass
                labels = [
                    targetlist[0],
                    'UPPER ERR',
                    'LOWER ERR',
                    'UNITS',
                    'REF',
                ]
                table = visitor.add_table(clabels=labels, rows=len(allstar))
                for starlabels in allstar:
                    i = allstar.index(starlabels)
                    for starlabel in starlabels:
                        j = starlabels.index(starlabel)
                        table.get_cell(i, j).add_primitive(starlabel)
                        elem = starinfo[starlabel][0]
                        if elem:
                            table.get_cell(i, j).add_primitive(elem)
                            pass
                        else:
                            table.get_cell(i, j).add_primitive('NA')
                            pass
                        pass
                    pass
                for c in self['starID'][targetlist[0]]['planets']:
                    labels = [
                        'PLANET ' + c,
                        'UPPER ERR',
                        'LOWER ERR',
                        'UNITS',
                        'REF',
                    ]
                    table = visitor.add_table(
                        clabels=labels, rows=len(allplanet)
                    )
                    for starlabels in allplanet:
                        i = allplanet.index(starlabels)
                        for starlabel in starlabels:
                            j = starlabels.index(starlabel)
                            table.get_cell(i, j).add_primitive(starlabel)
                            elem = starinfo[c][starlabel][0]
                            if elem:
                                table.get_cell(i, j).add_primitive(elem)
                            else:
                                table.get_cell(i, j).add_primitive('NA')

        return


# ------------ -------------------------------------------------------
# -- FILTER -- -------------------------------------------------------
class FilterSV(dawgie.StateVector):
    '''filter view'''

    def __init__(self, name):
        '''__init__ ds'''
        self._version_ = dawgie.VERSION(1, 1, 1)
        self.__name = name
        self['STATUS'] = excalibur.ValuesList()
        self['PROCESS'] = excalibur.ValuesDict()
        self['activefilters'] = excalibur.ValuesDict()
        self['STATUS'].append(False)
        return

    def name(self):
        '''name ds'''
        return self.__name

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''view ds'''
        if self['STATUS'][-1]:
            if len(self['STATUS']) < 3:
                labels = ['Active Filters']
                nf = self['activefilters']['TOTAL']
                table = visitor.add_table(clabels=labels, rows=nf)
                for name in self['activefilters']['NAMES']:
                    i = self['activefilters']['NAMES'].index(name)
                    table.get_cell(i, 0).add_primitive(name)
                    pass
                pass
            if len(self['STATUS']) > 2:
                ignorekeys = ['NAMES', 'TOTAL']
                actflts = [
                    f
                    for f in self['activefilters'].keys()
                    if f not in ignorekeys
                ]
                labels = ['Filter', 'Frames collected']
                table = visitor.add_table(clabels=labels, rows=len(actflts))
                for flt in actflts:
                    i = actflts.index(flt)
                    number = len(self['activefilters'][flt]['TOTAL'])
                    table.get_cell(i, 0).add_primitive(flt)
                    table.get_cell(i, 1).add_primitive(number)
        return


# ------------ -------------------------------------------------------
# -- DATABASE -- -----------------------------------------------------
class DatabaseSV(dawgie.StateVector):
    '''target.scrape view'''

    def __init__(self, name):
        self._version_ = dawgie.VERSION(1, 1, 1)
        self.__name = name
        self['STATUS'] = excalibur.ValuesList()
        self['name'] = excalibur.ValuesDict()
        self['STATUS'].append(False)
        return

    def name(self):
        '''__init__ ds'''
        return self.__name

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''view ds'''
        if self['STATUS'][-1]:
            ordlab = ['observatory', 'instrument', 'detector', 'filter', 'mode']
            table = visitor.add_table(clabels=ordlab, rows=1)
            for label in ordlab:
                vlist = [self['name'][n][label] for n in self['name'].keys()]
                i = ordlab.index(label)
                if vlist is not None:
                    for v in set(vlist):
                        total = len(vlist)
                        nb = vlist.count(v)
                        if v is not None:
                            percent = 1e2 * vlist.count(v) / total
                            out = (
                                str(v)
                                + ': '
                                + str(int(nb))
                                + ' ('
                                + str(round(percent))
                                + '%)'
                            )
                            table.get_cell(0, i).add_primitive(out)

        return


class ScrapeValidationSV(dawgie.StateVector):
    '''SoCustomSV ds'''

    def __init__(self, name):
        '''__init__ ds'''
        self._version_ = dawgie.VERSION(1, 1, 1)

        # data's structure
        # {runID: {'jwst': int, 'hst': int}}
        self['data'] = excalibur.ValuesList()
        self['data'].append({})

        # quality's structure
        # {rid: }
        # -1 = bad, 0 = iffy, 1 = good. only -1 and 1 used so far
        self['quality'] = excalibur.ValuesList()
        self['quality'].append({})

        self['STATUS'] = excalibur.ValuesList()
        self['STATUS'].append(False)
        self.__name = name
        return

    def name(self):
        '''name ds'''
        return self.__name

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''view ds'''
        df = pd.DataFrame(self['data'][0]).T
        df = df.fillna(value=0)

        p = bokeh.plotting.figure(
            title="# of Frames vs RunID",
            x_axis_label="RunIDs",
            width=800,
            height=400,
        )
        colors = bokeh.palettes.magma(len(df.columns))
        legend_items = []

        for i, col in enumerate(df.columns):
            color = colors[i]

            source = bokeh.models.ColumnDataSource(
                data={
                    'runid': df.index,
                    'count': df[col],
                    'colname': [col] * len(df),
                }
            )

            line = p.line(
                'runid', 'count', source=source, line_width=2, color=color
            )
            dots = p.scatter(
                'runid', 'count', source=source, size=5, color=color
            )

            hover = bokeh.models.HoverTool(
                renderers=[dots],
                tooltips=[
                    ('Instrument', '@colname'),
                    ('RunID', '@runid'),
                    ('Count', '@count'),
                ],
            )
            p.add_tools(hover)
            legend_items.append(
                bokeh.models.LegendItem(label=col, renderers=[line, dots])
            )

        legend = bokeh.models.Legend(items=legend_items, location="center")
        p.add_layout(legend, 'right')

        p.legend.label_text_font_size = '10pt'
        p.legend.spacing = 2
        p.legend.label_standoff = 5

        js, div = bokeh.embed.components(p)
        visitor.add_declaration(None, div=div, js=js)

        raw = self['quality'][0]

        # show just the last 15 statuses
        x = list(raw.keys())[-15:]
        y = list(raw.values())[-15:]
        source = bokeh.models.ColumnDataSource(data={'runid': x, 'status': y})

        p = bokeh.plotting.figure(
            title="Status vs RunID (1: Good, -1: Bad)",
            x_axis_label="RunIDs",
            y_axis_label="Status",
            width=800,
            height=400,
        )

        line = p.line(
            'runid', 'status', source=source, line_width=2, color="orange"
        )
        dots = p.circle(
            'runid', 'status', source=source, size=5, color="orange"
        )

        legend = bokeh.models.Legend(
            items=[
                bokeh.models.LegendItem(label="Status", renderers=[line, dots])
            ]
        )
        p.add_layout(legend, 'above')
        js, div = bokeh.embed.components(p)
        visitor.add_declaration(None, div=div, js=js)
        return


# -------------- -----------------------------------------------------
