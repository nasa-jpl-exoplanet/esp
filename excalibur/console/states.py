'''state vectors to hold the processed metric data'''

import bokeh.embed
import bokeh.layouts
import bokeh.models
import bokeh.plotting
import dawgie
import excalibur
import logging
import math
import numpy
import warnings

log = logging.getLogger(__name__)


class Accomplished(dawgie.StateVector):
    def __init__(self):
        self._version_ = dawgie.VERSION(1, 0, 0)
        self['algs'] = excalibur.ValuesDict()
        self['cube'] = excalibur.ValuesDict()
        self['maxes'] = excalibur.ValuesList()
        self['rids'] = excalibur.ValuesDict()
        self['tns'] = excalibur.ValuesDict()

    def fill(self, mds):
        # lots of temp varaibles to build the data cube in a bokeh friendly way
        # pylint: disable=too-many-locals
        algs = set()
        tns = set()
        rids = set()
        ak = []
        rk = []
        rx = {}
        tk = []
        for md in mds:
            alg = '.'.join([md.task, md.alg_name])
            algs.add(alg)
            ak.append(alg)
            rids.add(md.run_id)
            rk.append(md.run_id)
            n = '.'.join([md.target, md.task, md.alg_name])
            rx[n] = max(md.run_id, rx.get(n, -1))
            tns.add(md.target)
            tk.append(md.target)
        algs = dict((s, i) for i, s in enumerate(sorted(algs)))
        rids = dict((s, i) for i, s in enumerate(sorted(rids)))
        tns = dict((s, i) for i, s in enumerate(sorted(tns)))
        a = []
        r = []
        t = []
        maxes = []
        for i in range(len(ak)):  # pylint: disable=consider-using-enumerate
            n = tk[i] + '.' + ak[i]
            maxes.append(rk[i] == rx[n])
            a.append(algs[ak[i]])
            r.append(rids[rk[i]])
            t.append(tns[tk[i]])
        self['algs'].clear()
        self['algs'].update(algs)
        self['cube'].clear()
        self['cube'].update(
            {'alg': a, 'rid': r, 'tid': t, 'algn': ak, 'rn': rk, 'tn': tk}
        )
        self['maxes'].clear()
        self['maxes'].extend(maxes)
        self['rids'].clear()
        self['rids'].update(rids)
        self['tns'].clear()
        self['tns'].update(tns)

    def features(self):
        '''contains no features'''
        return []

    def name(self):
        '''database name'''
        return "accomplishments"

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''Show the configutation information'''
        fig = plot_cube(
            self['cube'], self['algs'], self['rids'], self['tns'], self['maxes']
        )
        js, div = bokeh.embed.components(fig)
        visitor.add_declaration(None, div=div, js=js)


class CpuAndMem(dawgie.StateVector):
    '''State representation of metric data'''

    def __init__(self):
        '''init the state vector by building it from condensed form'''
        self._version_ = dawgie.VERSION(1, 0, 0)
        self['rids'] = excalibur.ValuesDict()
        self['talgs'] = excalibur.ValuesDict()
        self['tns'] = excalibur.ValuesDict()
        self['x_rid'] = excalibur.ValuesDict()
        self['x_tn'] = excalibur.ValuesDict()

    def fill(self, mds):
        self['rids'].clear()
        self['talgs'].clear()
        self['tns'].clear()
        self['x_rid'].clear()
        self['x_rid'].update(
            {
                'talg': [],
                'x': [],
                'ym': [],
                'ymn': [],
                'ymx': [],
                'yt': [],
                'ytn': [],
                'ytx': [],
            }
        )
        self['x_tn'].clear()
        self['x_tn'].update(
            {
                'talg': [],
                'x': [],
                'ym': [],
                'ymn': [],
                'ymx': [],
                'yt': [],
                'ytn': [],
                'ytx': [],
            }
        )
        rids = set()
        table = {}
        tns = set()
        for md in mds:
            name = '.'.join([md.task, md.alg_name])
            sv = {
                'db_memory': md.sv['db_memory'].value(),
                'task_memory': md.sv['task_memory'].value(),
                'task_wall': md.sv['task_wall'].value(),
            }
            if name not in table:
                table[name] = {'byrid': {}, 'bytn': {}}
            if md.run_id not in table[name]['byrid']:
                table[name]['byrid'][md.run_id] = []
            if md.target not in table[name]['bytn']:
                table[name]['bytn'][md.target] = []
            rids.add(md.run_id)
            tns.add(md.target)
            table[name]['byrid'][md.run_id].append(sv)
            table[name]['bytn'][md.target].append(sv)
        self['rids'].update(dict((n, i) for i, n in enumerate(sorted(rids))))
        self['talgs'].update(dict((n, i) for i, n in enumerate(sorted(table))))
        self['tns'].update(dict((n, i) for i, n in enumerate(sorted(tns))))
        for talg in self['talgs']:
            unravel(
                self['talgs'][talg],
                self['x_rid'],
                table[talg]['byrid'],
                self['rids'],
            )
            unravel(
                self['talgs'][talg],
                self['x_tn'],
                table[talg]['bytn'],
                self['tns'],
            )

    def features(self):
        '''contains no features'''
        return []

    def name(self):
        '''database name'''
        return "resource usage"

    def plot_resources(self):
        talg, idx = next(iter(self['talgs'].items()))
        rid = bokeh.models.ColumnDataSource(data=self['x_rid'])
        tn = bokeh.models.ColumnDataSource(data=self['x_tn'])
        rf = bokeh.models.BooleanFilter(
            [n == idx for n in self['x_rid']['talg']]
        )
        rv = bokeh.models.CDSView(filter=rf)
        tf = bokeh.models.BooleanFilter(
            [n == idx for n in self['x_rid']['talg']]
        )
        tv = bokeh.models.CDSView(filter=tf)
        ll = bokeh.plotting.figure(
            height=400,
            tools='pan,ywheel_zoom,box_zoom,xwheel_zoom,save,reset',
            width=600,
            x_axis_location='below',
        )
        ll.segment('x', 'ytn', 'x', 'ytx', color='cyan', source=rid, view=rv)
        ll.scatter('x', 'yt', marker='*', color='blue', source=rid, view=rv)
        ll.xaxis.axis_label = 'Run ID'
        ll.xaxis.formatter = bokeh.models.CustomJSTickFormatter(code=f'''
var labels = {list(self['rids'])};
return labels[tick];''')
        ll.xaxis.major_label_orientation = math.pi / 4
        ll.yaxis.axis_label = 'Wall Clock [d H:M:S]'
        ll.yaxis.formatter = bokeh.models.CustomJSTickFormatter(code='''
function convertSecondsToDHMS(totalSeconds) {
  totalSeconds = Number(totalSeconds);
  const days = Math.floor(totalSeconds / (24 * 3600));
  let remainingSeconds = totalSeconds % (24 * 3600);
  const hours = Math.floor(remainingSeconds / 3600);
  remainingSeconds %= 3600;
  const minutes = Math.floor(remainingSeconds / 60);
  const seconds = remainingSeconds % 60;
  let result = "";
  result += `${days} `;
  result += `${hours.toString().padStart(2, '0')}:`;
  result += `${minutes.toString().padStart(2, '0')}:`;
  result += `${seconds.toString().padStart(2, '0')}`;
  return result;
}
return convertSecondsToDHMS(tick);''')
        lr = bokeh.plotting.figure(
            height=400,
            tools='pan,ywheel_zoom,box_zoom,xwheel_zoom,save,reset',
            width=600,
            x_axis_location='below',
            y_axis_location='right',
        )
        lr.segment('x', 'ytn', 'x', 'ytx', color='cyan', source=tn, view=tv)
        lr.scatter('x', 'yt', marker='*', color='blue', source=tn, view=tv)
        lr.xaxis.axis_label = 'Target'
        lr.xaxis.formatter = bokeh.models.CustomJSTickFormatter(code=f'''
var labels = {list(self['tns'])};
return labels[tick];''')
        lr.xaxis.major_label_orientation = math.pi / 4
        lr.yaxis.axis_label = 'Wall Clock [d H:M:S]'
        lr.yaxis.formatter = bokeh.models.CustomJSTickFormatter(code='''
function convertSecondsToDHMS(totalSeconds) {
  totalSeconds = Number(totalSeconds);
  const days = Math.floor(totalSeconds / (24 * 3600));
  let remainingSeconds = totalSeconds % (24 * 3600);
  const hours = Math.floor(remainingSeconds / 3600);
  remainingSeconds %= 3600;
  const minutes = Math.floor(remainingSeconds / 60);
  const seconds = remainingSeconds % 60;
  let result = "";
  result += `${days} `;
  result += `${hours.toString().padStart(2, '0')}:`;
  result += `${minutes.toString().padStart(2, '0')}:`;
  result += `${seconds.toString().padStart(2, '0')}`;
  return result;
}
return convertSecondsToDHMS(tick);''')
        ul = bokeh.plotting.figure(
            height=400,
            tools='pan,ywheel_zoom,box_zoom,xwheel_zoom,save,reset',
            width=600,
            x_axis_location='above',
        )
        ul.segment('x', 'ymn', 'x', 'ymx', color='cyan', source=rid, view=rv)
        ul.scatter('x', 'ym', marker='*', color='blue', source=rid, view=rv)
        ul.xaxis.axis_label = 'Run ID'
        ul.xaxis.formatter = bokeh.models.CustomJSTickFormatter(code=f'''
var labels = {list(self['rids'])};
return labels[tick];''')
        ul.xaxis.major_label_orientation = math.pi / 4
        ul.yaxis.axis_label = 'Memory [B]'
        ul.yaxis.formatter = bokeh.models.CustomJSTickFormatter(code='''
function format_engineering(value) {
    if (value === 0) return "0";
    let sign = value < 0 ? "-" : "";
    let abs_value = Math.abs(value);
    let exponent = Math.floor(Math.log10(abs_value) / 3) * 3;
    let mantissa = abs_value / Math.pow(10, exponent);
    return sign + mantissa.toFixed(2) + "e" + exponent;
}
return format_engineering(tick);''')
        ur = bokeh.plotting.figure(
            height=400,
            tools='pan,ywheel_zoom,box_zoom,xwheel_zoom,save,reset',
            width=600,
            x_axis_location='above',
            y_axis_location='right',
        )
        ur.segment('x', 'ymn', 'x', 'ymx', color='cyan', source=tn, view=tv)
        ur.scatter('x', 'ym', marker='*', color='blue', source=tn, view=tv)
        ur.xaxis.axis_label = 'Target'
        ur.xaxis.formatter = bokeh.models.CustomJSTickFormatter(code=f'''
            var labels = {list(self['tns'])};
            return labels[tick];''')
        ur.xaxis.major_label_orientation = math.pi / 4
        ur.yaxis.axis_label = 'Memory [B]'
        ur.yaxis.formatter = bokeh.models.CustomJSTickFormatter(code='''
function format_engineering(value) {
    if (value === 0) return "0";
    let sign = value < 0 ? "-" : "";
    let abs_value = Math.abs(value);
    let exponent = Math.floor(Math.log10(abs_value) / 3) * 3;
    let mantissa = abs_value / Math.pow(10, exponent);
    return sign + mantissa.toFixed(2) + "e" + exponent;
}
return format_engineering(tick);''')
        sel = bokeh.models.Select(
            title='Task.Algorithm:',
            value=talg,
            options=list(self['talgs']),
        )
        sel.js_on_change(
            'value',
            bokeh.models.CustomJS(
                args={
                    'select': sel,
                    'xr': self['x_rid']['talg'],
                    'xt': self['x_tn']['talg'],
                    'rf': rf,
                    'tf': tf,
                    'talgs': self['talgs'],
                },
                code="""
    const N = talgs[select.value]
    for (let i = 0 ; i < xr.length ; i++) {
       rf.booleans[i] = xr[i] === N
    }
    for (let i = 0 ; i < xt.length ; i++) {
       tf.booleans[i] = xt[i] === N
    }
    rf.change.emit()
    tf.change.emit()
        """,
            ),
        )

        ur.y_range = ul.y_range
        lr.y_range = ll.y_range
        ll.x_range = ul.x_range
        lr.x_range = ur.x_range
        fig = bokeh.layouts.gridplot(
            [
                [sel],
                [ul, ur],
                [ll, lr],
            ],
            toolbar_location='above',
        )
        return fig

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''Show the configutation information'''
        visitor.add_declaration_inline('', div='<div><hr>')
        visitor.add_declaration_inline('Dark blue dots are the median')
        visitor.add_declaration_inline('Thin cyan lines are the min/max')
        visitor.add_declaration_inline(
            'When x-axis is run id, then min/max is over targets and run id when x-axis is target'
        )
        visitor.add_declaration_inline('', div='</div>')
        js, div = bokeh.embed.components(self.plot_resources())
        visitor.add_declaration(None, div=div, js=js)


def find_items(base, value):
    return [v == value for v in base]


def plot_cube(data, algs, rids, tns, max_rids):
    # plotting is a mess because have to keep handles on so much stuff
    # pylint: disable=too-many-locals,use-dict-literal
    an = next(iter(algs))
    tn = next(iter(tns))
    ff = bokeh.models.BooleanFilter(find_items(data['algn'], an))
    fv = bokeh.models.CDSView(filter=ff)
    sf = bokeh.models.BooleanFilter(find_items(data['tn'], tn))
    sv = bokeh.models.CDSView(filter=sf)
    tf0 = bokeh.models.BooleanFilter(
        [
            a and not b and not c
            for a, b, c in zip(max_rids, ff.booleans, sf.booleans)
        ]
    )
    tf1 = bokeh.models.BooleanFilter(
        [a and b for a, b in zip(max_rids, ff.booleans)]
    )
    tf2 = bokeh.models.BooleanFilter(
        [a and b for a, b in zip(max_rids, sf.booleans)]
    )
    tv0 = bokeh.models.CDSView(filter=tf0)
    tv1 = bokeh.models.CDSView(filter=tf1)
    tv2 = bokeh.models.CDSView(filter=tf2)
    cds = bokeh.models.ColumnDataSource(data=data)
    top = bokeh.plotting.figure(
        active_drag='pan',
        active_scroll='wheel_zoom',
        height=400,
        title='Top View',
        toolbar_location='above',
        tools='hover,pan,reset,wheel_zoom',
        tooltips=[("ID", "@rid.@tn.@algn")],
        width=800,
        y_axis_location='right',
    )
    t0 = top.scatter(
        x='tid', y='alg', color='blue', marker='*', source=cds, view=tv0
    )
    t1 = top.scatter(
        x='tid', y='alg', color='green', marker='*', source=cds, view=tv1
    )
    t2 = top.scatter(
        x='tid', y='alg', color='red', marker='*', source=cds, view=tv2
    )
    top.select_one(bokeh.models.HoverTool).renderers = [t0, t1, t2]
    top.xaxis.formatter = bokeh.models.CustomJSTickFormatter(code=f'''
var labels = {sorted(tns)};
return labels[tick];''')
    top.xaxis.major_label_orientation = math.pi / 2
    top.yaxis.formatter = bokeh.models.CustomJSTickFormatter(code=f'''
var labels = {sorted(algs)};
return labels[tick];''')
    front = bokeh.plotting.figure(
        active_drag='pan',
        active_scroll='xwheel_zoom',
        height=400,
        title=f'Front View of plane {an}',
        title_location='below',
        tools='pan,xwheel_zoom',
        width=800,
        x_axis_location='above',
    )
    front.toolbar_location = None
    front.scatter(
        x='tid', y='rid', color='green', marker='x', source=cds, view=fv
    )
    front.yaxis.axis_label = 'Run ID'
    front.xaxis.formatter = bokeh.models.CustomJSTickFormatter(code=f'''
var labels = {sorted(rids)};
return labels[tick];''')
    front.xaxis.major_label_text_font_size = "0pt"
    front.x_range = top.x_range
    side = bokeh.plotting.figure(
        active_drag='pan',
        active_scroll='ywheel_zoom',
        height=400,
        title=f'Side View of Plane {tn}',
        tools='pan,ywheel_zoom',
        width=600,
    )
    side.toolbar_location = None
    side.scatter(x='rid', y='alg', color='red', marker='x', source=cds, view=sv)
    side.xaxis.axis_label = 'Run ID'
    side.xaxis.formatter = bokeh.models.CustomJSTickFormatter(code=f'''
var labels = {sorted(rids)};
return labels[tick];''')
    side.yaxis.major_label_text_font_size = "0pt"
    side.y_range = top.y_range
    htcb = bokeh.models.CustomJS(
        args=dict(
            cds=cds,
            front=front,
            side=side,
            ff=ff,
            sf=sf,
            tf0=tf0,
            tf1=tf1,
            tf2=tf2,
            mxr=max_rids,
        ),
        code='''
    const {indices} = cb_data.index
    if (indices.length > 0) {
       let an = cds.data["algn"][indices[0]]
       let tn = cds.data["tn"][indices[0]]
       front.title.text = "Front View of Plane " + an
       side.title.text = "Side View of Plane " + tn
       for (let i = 0; i < ff.booleans.length; i++) {
          ff.booleans[i] = cds.data["algn"][i] === an
          sf.booleans[i] = cds.data["tn"][i] === tn
          tf0.booleans[i] = mxr[i] && !ff.booleans[i] && !sf.booleans[i]
          tf1.booleans[i] = mxr[i] && ff.booleans[i]
          tf2.booleans[i] = mxr[i] && sf.booleans[i]
       }
       ff.change.emit()
       sf.change.emit()
       tf0.change.emit()
       tf1.change.emit()
       tf2.change.emit()
    }
    ''',
    )
    ht0 = bokeh.models.HoverTool(
        callback=htcb,
        renderers=[t0],
        tooltips=[],
    )
    ht1 = bokeh.models.HoverTool(
        callback=htcb,
        renderers=[t1],
        tooltips=[],
    )
    ht2 = bokeh.models.HoverTool(
        callback=htcb,
        renderers=[t2],
        tooltips=[],
    )
    top.add_tools(ht0, ht1, ht2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        fig = bokeh.layouts.gridplot(
            [
                [top, side],
                [
                    front,
                ],
            ],
            toolbar_location='above',
        )
    return fig


def unravel(n, src, by, idx):
    for x, ys in by.items():
        dm = numpy.empty(len(ys))
        tm = numpy.empty(len(ys))
        tw = numpy.empty(len(ys))
        for i, y in enumerate(ys):
            dm[i] = y['db_memory']
            tm[i] = y['task_memory']
            tw[i] = y['task_wall']
        dm[dm < 0] = numpy.nan
        tm[tm < 0] = numpy.nan
        tw[tw < 0] = numpy.nan
        m = dm + tm
        src['talg'].append(n)
        src['x'].append(idx[x])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            src['ym'].append(numpy.nanmedian(m))
            src['ymn'].append(numpy.nanmin(m))
            src['ymx'].append(numpy.nanmax(m))
            src['yt'].append(numpy.nanmedian(tw))
            src['ytn'].append(numpy.nanmin(tw))
            src['ytx'].append(numpy.nanmax(tw))
