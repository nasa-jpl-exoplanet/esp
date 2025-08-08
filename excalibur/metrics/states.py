'''state vectors to hold the processed metric data'''

import bokeh.embed
import bokeh.layouts
import bokeh.models
import bokeh.plotting
import dawgie
import excalibur
import logging

log = logging.getLogger(__name__)


class CpuAndMem(dawgie.StateVector):
    '''State representation of metric data'''

    def __init__(self, name, mds):
        '''init the state vector by building it from condensed form'''
        self._version_ = dawgie.VERSION(1, 0, 0)
        self._name = name
        self['byrid'] = excalibur.ValuesDict()
        self['bytn'] = excalibur.ValuesDict()
        for md in mds:
            if md.run_id not in self['byrid']:
                sel['byrid'][md.run_id] = []
            self['byrid'][md.run_id].append(md.sv)
            if md.target not in self['bytn']:
                self['bytn'][md.target] = []
            self['bytn'][md.target].append(md.sv)

    def features(self):
        '''contains no features'''
        return []

    def name(self):
        '''database name'''
        return self._name.replace('.', '::')

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


def plot_resource(data_table, keys, tan, xlabel, mem):
    res = 'Memory' if mem else 'Processing Time'
    svs = [data_table[k] for k in keys]
    fig = bokeh.plotting.figure(
        height=400,
        title=f'{res}: AE component {tan}',
        tools='pan,ywheel_zoom,box_zoom,xwheel_zoom,save,reset',
        width=600,
    )
    x = [i for i in range(len(keys))]
    y = []
    yn = []
    yx = []
    for svl in svs:
        w = [
            (
                sv['db_memory' if mem else 'task_wall'].value()
                + sv['task_memory' if mem else 'task_wall'].value()
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
