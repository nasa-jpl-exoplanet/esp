'''Phasecurve Database Products View'''

# Heritage code shame:
# pylint: disable=too-many-locals,too-many-branches

# -- IMPORTS -- ------------------------------------------------------

import dawgie
import excalibur

import matplotlib.pyplot as plt
from excalibur.util.plotters import plot_normalized_byvisit, save_plot_toscreen
from excalibur.util.svs import ExcaliburSV


# ------------- ------------------------------------------------------
# -- SV -- -----------------------------------------------------------
class NormSV(ExcaliburSV):
    '''phasecurve.normalization view'''

    def __init__(self, name):
        '''__init__ ds'''
        ExcaliburSV.__init__(self, name, dawgie.VERSION(1, 1, 0))

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''view ds'''
        if self['STATUS'][-1]:
            for p in self['data'].keys():
                if 'vignore' in self['data'][p]:
                    for v, m in zip(
                        self['data'][p]['vignore'], self['data'][p]['trial']
                    ):
                        strignore = str(int(v)) + ' ' + m
                        visitor.add_declaration('VISIT IGNORED: ' + strignore)
                    pass
                if 'vrange' in self['data'][p]:
                    vrange = self['data'][p]['vrange']
                    plot_normalized_byvisit(self['data'][p], vrange, visitor)
            pass
        pass


class WhiteLightSV(ExcaliburSV):
    '''phasecurve.whitelight view'''

    def __init__(self, name):
        '''__init__ ds'''
        ExcaliburSV.__init__(self, name, dawgie.VERSION(1, 1, 1))

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''view ds'''
        if self['STATUS'][-1]:
            for p in self['data'].keys():

                if 'HST' in self.name():

                    visits = self['data'][p]['visits']
                    phase = self['data'][p]['phase']
                    allwhite = self['data'][p]['allwhite']
                    postim = self['data'][p]['postim']
                    postphase = self['data'][p]['postphase']
                    postlc = self['data'][p]['postlc']
                    postflatphase = self['data'][p]['postflatphase']
                    myfig = plt.figure(figsize=(10, 6))
                    plt.title(p)
                    for index, v in enumerate(visits):
                        plt.plot(phase[index], allwhite[index], 'k+')
                        plt.plot(
                            postphase[index],
                            allwhite[index] / postim[index],
                            'o',
                            label=str(v),
                        )
                        pass
                    if len(visits) > 14:
                        ncol = 2
                    else:
                        ncol = 1
                    plt.plot(postflatphase, postlc, '^', label='M')
                    plt.xlabel('Orbital Phase')
                    plt.ylabel('Normalized Post White Light Curve')
                    plt.legend(
                        bbox_to_anchor=(1 + 0.1 * (ncol - 0.5), 0.5),
                        loc=5,
                        ncol=ncol,
                        mode='expand',
                        numpoints=1,
                        borderaxespad=0.0,
                        frameon=False,
                    )
                    plt.tight_layout(rect=[0, 0, (1 - 0.1 * ncol), 1])
                    save_plot_toscreen(myfig, visitor)
                elif 'Spitzer' in self.name():
                    # for each event
                    for i in range(len(self['data'][p])):
                        # plots are saved into sv
                        visitor.add_image(
                            '...', ' ', self['data'][p][i]['plot_bestfit']
                        )
                        visitor.add_image(
                            '...', ' ', self['data'][p][i]['plot_residual_fft']
                        )
                        visitor.add_image(
                            '...', ' ', self['data'][p][i]['plot_posterior']
                        )
                        visitor.add_image(
                            '...', ' ', self['data'][p][i]['plot_pixelmap']
                        )
                        # another centroid timeseries plot?
        return


# -------- -----------------------------------------------------------


class FlaresSV(ExcaliburSV):
    '''phasecurve.flaredetection view'''

    def __init__(self, name):
        '''__init__ ds'''
        ExcaliburSV.__init__(self, name, dawgie.VERSION(1, 2, 0))

    @staticmethod
    def _format_value(value):
        if value is None or value == '':
            return 'N/A'
        if isinstance(value, float):
            return f'{value:.6g}'
        return str(value)

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''view ds'''
        if self['STATUS'][-1]:
            metadata = (
                self['data']['metadata'] if 'metadata' in self['data'] else {}
            )
            if metadata:
                visitor.add_declaration(
                    f'Flare detection results for '
                    f'{metadata.get("target", "UNKNOWN_TARGET")} '
                    f'({metadata.get("filter", self.name())})'
                )
                visitor.add_declaration(
                    'Distance [pc]: '
                    f'{self._format_value(metadata.get("distance_pc"))}; '
                    'Flux density [mJy]: '
                    f'{self._format_value(metadata.get("flux_density_mjy"))}; '
                    'Quiescent luminosity [W]: '
                    f'{self._format_value(metadata.get("quiescent_luminosity_w"))}; '
                    'C_bol: '
                    f'{self._format_value(metadata.get("c_bol"))}'
                )
                if metadata.get('results_dir'):
                    visitor.add_declaration(
                        f'Results directory: {metadata["results_dir"]}'
                    )

            if 'frequency_summary' in self['data']:
                summary = self['data']['frequency_summary']
                visitor.add_declaration(
                    'Flare frequency [1/day]: '
                    f'{self._format_value(summary.get("flare_frequency_per_day"))}; '
                    'Flare frequency [1/hour]: '
                    f'{self._format_value(summary.get("flare_frequency_per_hour"))}; '
                    'Total observed time [days]: '
                    f'{self._format_value(summary.get("total_observed_time_days"))}; '
                    'Unique flare intervals: '
                    f'{self._format_value(summary.get("unique_flare_intervals"))}; '
                    'Detected flare rows: '
                    f'{self._format_value(summary.get("detected_flare_rows"))}'
                )

            visit_rows = []
            flare_rows = []
            for p in self['data'].keys():
                if p in ('target', 'filter', 'metadata', 'frequency_summary'):
                    continue
                for ivisit in range(len(self['data'][p])):
                    visit_data = self['data'][p][ivisit]
                    visit_rows.append(
                        {
                            'planet': p,
                            'visit': visit_data.get('visit', ivisit),
                            'n_flares': visit_data.get('n_flares', 0),
                            'status': (
                                'error' if visit_data.get('error') else 'ok'
                            ),
                            'error': visit_data.get('error', ''),
                        }
                    )
                    visitor.add_declaration(
                        f'{p} visit {visit_data.get("visit", ivisit)}: '
                        f'{visit_data.get("n_flares", 0)} flare(s)'
                    )
                    if 'plot_lightcurve' in visit_data:
                        visitor.add_image(
                            '...',
                            (
                                f'Example detrended light curve: {p} '
                                f'visit {visit_data.get("visit", ivisit)}'
                            ),
                            visit_data['plot_lightcurve'],
                        )
                    for flare_index, flare_data in enumerate(
                        visit_data.get('flares', [])
                    ):
                        flare_rows.append(
                            {
                                'planet': p,
                                'visit': visit_data.get('visit', ivisit),
                                'flare': flare_index,
                                'start': flare_data.get('start', ''),
                                'stop': flare_data.get('stop', ''),
                                'tpeak': flare_data.get('tpeak', ''),
                                'observed_duration_minutes': flare_data.get(
                                    'observed_duration_minutes',
                                    '',
                                ),
                                'fwhm_minutes': flare_data.get(
                                    'fwhm_minutes',
                                    '',
                                ),
                                'ampl': flare_data.get('ampl', ''),
                                'ED_seconds': flare_data.get('ED_seconds', ''),
                                'peak_flare_luminosity_w': flare_data.get(
                                    'peak_flare_luminosity_w',
                                    '',
                                ),
                                'E_band_ergs': flare_data.get(
                                    'E_band_ergs',
                                    '',
                                ),
                                'E_bol_ergs': flare_data.get(
                                    'E_bol_ergs',
                                    '',
                                ),
                            }
                        )
                        if 'plot_fit' in flare_data:
                            visitor.add_image(
                                '...',
                                (
                                    f'Example flare fit: {p} '
                                    f'visit {visit_data.get("visit", ivisit)} '
                                    f'flare {flare_index}'
                                ),
                                flare_data['plot_fit'],
                            )
                        if 'plot_posterior' in flare_data:
                            visitor.add_image(
                                '...',
                                (
                                    f'Example flare posterior: {p} '
                                    f'visit {visit_data.get("visit", ivisit)} '
                                    f'flare {flare_index}'
                                ),
                                flare_data['plot_posterior'],
                            )

            if visit_rows:
                visit_table = visitor.add_table(
                    clabels=['planet', 'visit', 'n_flares', 'status', 'error'],
                    rows=len(visit_rows),
                )
                for row_index, visit_row in enumerate(visit_rows):
                    for col_index, key in enumerate(
                        ['planet', 'visit', 'n_flares', 'status', 'error']
                    ):
                        visit_table.get_cell(
                            row_index, col_index
                        ).add_primitive(
                            self._format_value(visit_row.get(key, ''))
                        )

            if flare_rows:
                flare_columns = [
                    'planet',
                    'visit',
                    'flare',
                    'start',
                    'stop',
                    'tpeak',
                    'observed_duration_minutes',
                    'fwhm_minutes',
                    'ampl',
                    'ED_seconds',
                    'peak_flare_luminosity_w',
                    'E_band_ergs',
                    'E_bol_ergs',
                ]
                flare_table = visitor.add_table(
                    clabels=flare_columns,
                    rows=len(flare_rows),
                )
                for row_index, flare_row in enumerate(flare_rows):
                    for col_index, key in enumerate(flare_columns):
                        flare_table.get_cell(
                            row_index, col_index
                        ).add_primitive(
                            self._format_value(flare_row.get(key, ''))
                        )
        return


# -------- -----------------------------------------------------------
