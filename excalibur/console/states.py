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

log = logging.getLogger(__name__)


class Accomplished(dawgie.StateVector):
    def __init__(self, mds):
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
        algs = excalibur.ValuesDict((s, i) for i, s in enumerate(sorted(algs)))
        rids = excalibur.ValuesDict((s, i) for i, s in enumerate(sorted(rids)))
        tns = excalibur.ValuesDict((s, i) for i, s in enumerate(sorted(tns)))
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
        self['algs'] = algs
        self['cube'] = excalibur.ValuesDict(
            alg=a, rid=r, tid=t, algn=ak, rn=rk, tn=tk
        )
        self['maxes'] = maxes
        self['rids'] = rids
        self['tns'] = tns

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

    def __init__(self, name, mds):
        '''init the state vector by building it from condensed form'''
        self._version_ = dawgie.VERSION(1, 0, 0)
        self._name = name
        self['byrid'] = excalibur.ValuesDict()
        self['bytn'] = excalibur.ValuesDict()
        for md in mds:
            sv = {
                'db_memory': md.sv['db_memory'].value(),
                'task_memory': md.sv['task_memory'].value(),
                'task_wall': md.sv['task_wall'].value(),
            }
            if md.run_id not in self['byrid']:
                self['byrid'][md.run_id] = []
            self['byrid'][md.run_id].append(sv)
            if md.target not in self['bytn']:
                self['bytn'][md.target] = []
            self['bytn'][md.target].append(sv)

    def features(self):
        '''contains no features'''
        return []

    def name(self):
        '''database name'''
        return "resource usage"

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''Show the configutation information'''
        visitor.add_declaration_inline('', div='<div><hr>')
        visitor.add_declaration_inline('Dark blue dots are the median')
        visitor.add_declaration_inline('Thin cyan lines are the min/max')
        visitor.add_declaration_inline(
            'When x-axis is run id, then min/max is over targets and run id when x-axis is target'
        )
        visitor.add_declaration_inline('', div='</div>')
        c1r1 = plot_resource(
            self['byrid'],
            sorted(self['byrid'].keys()),
            self._name,
            'Run ID',
            True,
        )
        c1r2 = plot_resource(
            self['byrid'],
            sorted(self['byrid'].keys()),
            self._name,
            'Run ID',
            False,
        )
        c1r2.x_range = c1r1.x_range
        c2r1 = plot_resource(
            self['bytn'],
            sorted(self['bytn'].keys()),
            self._name,
            'Target Name',
            True,
        )
        c2r2 = plot_resource(
            self['bytn'],
            sorted(self['bytn'].keys()),
            self._name,
            'Target Name',
            False,
        )
        c2r2.x_range = c2r1.x_range
        c2r1.y_range = c1r1.y_range
        c2r2.y_range = c1r2.y_range
        f = bokeh.layouts.gridplot(
            [[c1r1, c2r1], [c1r2, c2r2]], toolbar_location='right'
        )
        js, div = bokeh.embed.components(f)
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
    top.xaxis.formatter = bokeh.models.CustomJSTickFormatter(
        code=f'''
var labels = {sorted(tns)};
return labels[tick];'''
    )
    top.xaxis.major_label_orientation = math.pi / 2
    top.yaxis.formatter = bokeh.models.CustomJSTickFormatter(
        code=f'''
var labels = {sorted(algs)};
return labels[tick];'''
    )
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
    front.xaxis.formatter = bokeh.models.CustomJSTickFormatter(
        code=f'''
var labels = {sorted(rids)};
return labels[tick];'''
    )
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
    side.xaxis.formatter = bokeh.models.CustomJSTickFormatter(
        code=f'''
var labels = {sorted(rids)};
return labels[tick];'''
    )
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


def plot_resource(data_table, keys, tan, xlabel, mem):
    res = 'Memory' if mem else 'Processing Time'
    svs = [data_table[k] for k in keys]
    fig = bokeh.plotting.figure(
        height=400,
        title=f'{res}: AE component {tan}',
        tools='pan,ywheel_zoom,box_zoom,xwheel_zoom,save,reset',
        width=600,
    )
    x = list(range(len(keys)))
    y = []
    yn = []
    yx = []
    for svl in svs:
        w = [
            (
                sv['db_memory' if mem else 'task_wall']
                + sv['task_memory' if mem else 'task_wall']
            )
            * (1024 if mem else 1)
            for sv in svl
        ]
        y.append(numpy.nanmedian([numpy.nan if a < 0 else a for a in w]))
        yn.append(min(w))
        yx.append(max(w))
    fig.segment(x, yn, x, yx, color='cyan')
    fig.scatter(x, y, marker='*', color='blue')
    fig.xaxis.axis_label = xlabel
    fig.xaxis.formatter = bokeh.models.CustomJSTickFormatter(
        code=f'''
var labels = {keys};
return labels[tick];'''
    )
    fig.xaxis.major_label_orientation = math.pi / 4
    fig.yaxis.axis_label = 'Memory [B]' if mem else 'Wall Clock [days H:M:S]'
    if mem:
        fig.yaxis.formatter = bokeh.models.CustomJSTickFormatter(
            code='''
function format_engineering(value) {
    if (value === 0) return "0";
    let sign = value < 0 ? "-" : "";
    let abs_value = Math.abs(value);
    let exponent = Math.floor(Math.log10(abs_value) / 3) * 3;
    let mantissa = abs_value / Math.pow(10, exponent);
    return sign + mantissa.toFixed(2) + "e" + exponent;
}
return format_engineering(tick);'''
        )
    else:
        fig.yaxis.formatter = bokeh.models.CustomJSTickFormatter(
            code='''
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
return convertSecondsToDHMS(tick);'''
        )
    return fig
