'''Classifier Database Products View'''

# Heritage code shame:
# pylint: disable=too-many-branches,too-many-locals,too-many-nested-blocks,too-many-statements

# -- IMPORTS --------------------------------------------------------
import dawgie
import excalibur
from excalibur.util.plotters import save_plot_toscreen
from excalibur.util.svs import ExcaliburSV

import matplotlib.pyplot as plt
import numpy as np
import logging

log = logging.getLogger(__name__)


# -- SV -------------------------------------------------------------
class PredictSV(ExcaliburSV):
    '''PredictSV ds'''

    def __init__(self, name):
        '''__init__ ds'''
        ExcaliburSV.__init__(self, name, dawgie.VERSION(1, 1, 0))

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''view ds'''
        if self['STATUS'][-1]:
            for p in self['data'].keys():
                visitor.add_declaration('PLANET: ' + p)
                visitor.add_declaration(
                    'PREDICTION: ' + str(self['data'][p]['prediction'])
                )
                #                 allwhite = self['data'][p]['allwhite']
                #                 postlight curve residual cl = self['data'][p]['postlc']
                #                 postsep = self['data'][p]['postsep']
                #                 myfig = plt.figure(figsize=(10, 6))
                #                 plt.title(p)
                #                 plt.plot(postsep, allwhite, 'o')
                #                 plt.plot(postsep, postlc, '^', label='M')
                #                 plt.xlabel('Star / Planet separation [R$^*$]')
                #                 plt.ylabel('Normalized Post White Light Curve')
                #                 plt.legend(bbox_to_anchor=(1 + 0.1*(0.5), 0.5),
                #                            loc=5, ncol=1, mode='expand', numpoints=1,
                #                            borderaxespad=0., frameon=False)
                #                 plt.tight_layout(rect=[0,0,(1 - 0.1),1])
                #                 save_plot_toscreen(myfig, visitor)
            pass
        pass

    pass


# Summarize flags for each target
class FlagsSV(ExcaliburSV):
    '''Flags_SV ds'''

    def __init__(self, name):
        '''__init__ ds'''
        ExcaliburSV.__init__(self, name, dawgie.VERSION(1, 0, 0))

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''view ds'''

        flag_algs_info = {
            'count_points_wl': {
                'title': 'Points in Full and Total Transit',
                'field_descriptions': {
                    'full': 'Number of points in full transit',
                    'total': 'Number of points in total transit',
                },
            },
            'symmetry_wl': {
                'title': 'Light Curve Symmetry',
                'field_descriptions': {
                    'left': 'Number of points before transit',
                    'right': 'Number of points after transit',
                },
            },
            'rsdm': {
                'title': 'Residual Standard Deviation Mean (RSDM)',
                'field_descriptions': {'mean_rsdm': 'Average RSDM'},
            },
            'perc_rejected': {
                'title': 'Percent Rejected',
                'field_descriptions': {
                    'percent_rejected_value': 'Percentage of spectral channels rejected in cumulative spectrum distribution'
                },
            },
            'median_error': {
                'title': 'Median Error',
                'field_descriptions': {'median_error_value': 'Median error'},
            },
            'residual_shape': {
                'title': 'Light Curve Residual Shape',
                'field_descriptions': {'data': '', 'model': '', 'z': ''},
            },
        }

        if self['STATUS'][-1]:

            for k in self['data'].keys():

                if k not in flag_algs_info:  # if it's a planet

                    visitor.add_declaration("_____")

                    visitor.add_declaration("PLANET: " + str(k))

                    visitor.add_declaration(
                        "Overall Flag: " + str(self['data'][k]['overall_flag'])
                    )

                    for alg in self['data'][k].keys():
                        if alg != 'overall_flag':
                            visitor.add_declaration(
                                str(flag_algs_info[alg]['title']).upper()
                            )
                            visitor.add_declaration(
                                "Flag: "
                                + str(self['data'][k][alg]['flag_color'])
                            )
                            visitor.add_declaration(
                                "Flag Description: "
                                + str(self['data'][k][alg]['flag_descrip'])
                            )

                            if alg == 'residual_shape':
                                model = self['data'][k][alg]['model']
                                data = self['data'][k][alg]['data']
                                z = self['data'][k][alg]['z']

                                # Graph 1: data, model
                                fig, ax = plt.subplots(1, 2, figsize=(10, 4))
                                ax[0].set_title(
                                    'Planet ' + str(k) + ': Whitelight Curve'
                                )
                                ax[0].plot(z, data, '.', label='data')
                                ax[0].plot(z, model, '.', label='model')
                                ax[0].legend()
                                ax[0].set_xlabel(str('z'))
                                ax[0].set_ylabel(
                                    str('Normallized white light curve')
                                )

                                # Graph 2: data - model
                                ax[1].set_title('Residual Curve')
                                ax[1].axvline(1, alpha=0.25)
                                ax[1].axvline(-1, alpha=0.25)
                                ax[1].axhline(0, alpha=0.25)
                                ax[1].plot(z, data - model, '.', color='black')
                                ax[1].set_xlabel(str('z'))
                                ax[1].set_ylabel(str('data - model'))
                                y_lim = max(
                                    data - model, abs(min(data - model))
                                )
                                x_lim = max(z, abs(min(z)))
                                ax[1].set_ylim(y_lim * -1.5, y_lim * 1.5)
                                ax[1].set_xlim(x_lim * -1.25, x_lim * 1.25)
                                plt.tight_layout(pad=3.0)
                                save_plot_toscreen(fig, visitor)
                            else:
                                for field in self['data'][k][alg].keys():
                                    if field not in (
                                        "flag_color",
                                        "flag_descrip",
                                    ):
                                        visitor.add_declaration(
                                            str(
                                                flag_algs_info[alg][
                                                    'field_descriptions'
                                                ][field]
                                            )
                                            + ": "
                                            + str(self['data'][k][alg][field])
                                        )

                else:  # if it's a metric that's not planet-specific
                    visitor.add_declaration("_____")

                    visitor.add_declaration(
                        str(flag_algs_info[k]['title']).upper()
                    )

                    visitor.add_declaration(
                        "Flag: " + str(self['data'][k]['flag_color'])
                    )
                    visitor.add_declaration(
                        "Flag Description: "
                        + str(self['data'][k]['flag_descrip'])
                    )

                    for field in self['data'][k].keys():
                        if field not in ("flag_color", "flag_descrip"):
                            visitor.add_declaration(
                                str(
                                    flag_algs_info[k]['field_descriptions'][
                                        field
                                    ]
                                )
                                + ": "
                                + str(self['data'][k][field])
                            )

        pass

    pass


# Summarize flags across ALL targets
class FlagSummarySV(ExcaliburSV):
    '''Flag_Summary_SV ds'''

    def __init__(self, name):
        '''__init__ ds'''
        ExcaliburSV.__init__(self, name, dawgie.VERSION(1, 0, 0))

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''view ds'''

        if 'classifier_flags' in self['data']:

            # first table displays all of the data quality metrics with the count of green, yellow, and red flags
            # second table displays all of the target names associated with yellow/red flags
            flag_colors = ['green', 'yellow', 'red']
            vlabels = list(self['data']['classifier_flags'].keys())
            hlabels_count = ['Data Quality Metric', 'Green', 'Yellow', 'Red']
            hlabels_trgts = ['Data Quality Metric', 'Yellow', 'Red']

            flag_count_table = visitor.add_table(
                clabels=hlabels_count, rows=len(vlabels)
            )
            flag_trgts_table = visitor.add_table(
                clabels=hlabels_trgts, rows=len(vlabels)
            )

            # label rows with algorithm names
            for row, ele in enumerate(vlabels):
                flag_count_table.get_cell(row, 0).add_primitive(vlabels[row])
                flag_trgts_table.get_cell(row, 0).add_primitive(vlabels[row])

                # to work around pylint "unused-variable" warning and "enumerate" requirement
                if ele:
                    pass

            # go through all the classifier algorithms
            for a, vlabel in enumerate(vlabels):

                alg_data = self['data']['classifier_flags'][vlabel]

                # to work around pylint "unused-variable" warning and "enumerate" requirement
                if vlabel:
                    pass

                # go through all potential flag colors
                for i, c in enumerate(flag_colors):

                    # to work around pylint "unused-variable" warning and "enumerate" requirement
                    if c:
                        pass

                    color = flag_colors[i]
                    table_column = (
                        i + 1
                    )  # add 1 to account for 'Data Quality Metric' header.

                    if color in alg_data:
                        flag_count = alg_data[color][1]
                        flag_count_table.get_cell(
                            a, table_column
                        ).add_primitive(flag_count)
                        if color in ('red', 'yellow'):
                            flag_trgts = alg_data[color][0]
                            flag_trgts_table.get_cell(
                                a, table_column - 1
                            ).add_primitive(flag_trgts)
                    else:
                        flag_count_table.get_cell(
                            a, table_column
                        ).add_primitive(0)
                        if color in ('red', 'yellow'):
                            flag_trgts_table.get_cell(
                                a, table_column - 1
                            ).add_primitive(
                                'No targets with ' + str(color) + ' flags.'
                            )
            pass

        # note that 'gold' below refers to a yellow flag. the 'gold' matplotlib color to represent yellow.
        flag_algs_info = {
            'count_points_wl': {
                'suptitle': 'Points in Transit',
                'subplot_titles': {
                    'full': 'Number of points in full transit',
                    'total': 'Number of points in total transit',
                    'x_axis_label': 'Number of points',
                },
                'thresh_vals': {
                    'full': {'yellow': 5, 'red': 0},
                    'total': {'yellow': 6, 'red': 0},
                },
                'xscale': 'log',
            },
            'symmetry_wl': {
                'suptitle': 'Light Curve Symmetry',
                'subplot_titles': {
                    'left': 'Number of points before transit',
                    'right': 'Number of points after transit',
                    'x_axis_label': 'Number of points',
                },
                'thresh_vals': {'yellow': 0, 'red': '-inf'},
                'xscale': 'log',
            },
            'rsdm': {
                'suptitle': 'Residual Standard Deviation Metric (RSDM)',
                'subplot_titles': {'mean_rsdm': 'Average RSDM'},
                'thresh_vals': {'yellow': 8.804, 'red': 13.818},
                'xscale': 'log',
            },
            'perc_rejected': {
                'suptitle': 'Spectral Channels Rejected',
                'subplot_titles': {
                    'percent_rejected_value': 'Percent Rejected'
                },
                'thresh_vals': {'yellow': 31.844, 'red': 52.166},
                'xscale': 'linear',
            },
            'median_error': {
                'suptitle': 'Median Error',
                'subplot_titles': {'median_error_value': 'Median error'},
                'thresh_vals': {'yellow': 2.5, 'red': 'inf'},
                'xscale': 'linear',
            },
        }

        if 'classifier_vals' in self['data']:

            for metric in (
                metric
                for metric in self['data']['classifier_vals']
                if metric != 'residual_shape'
            ):

                metric_info = self['data']['classifier_vals'][metric]

                suptitle = flag_algs_info[metric]['suptitle']

                myfig = plt.figure()

                for count, subplot in enumerate(metric_info):

                    subplot_title = flag_algs_info[metric]['subplot_titles'][
                        subplot
                    ]
                    points_to_plot = metric_info[subplot]

                    plt.subplot(1, len(metric_info), (count + 1))

                    edgecolor = '#000'
                    color = '#000'
                    if 'thresh_vals' not in flag_algs_info[metric]:
                        edgecolor = '#000a36'
                        color = '#5c6280'

                    if flag_algs_info[metric]['xscale'] == 'log':
                        b = 25
                        _, bins = np.histogram(points_to_plot, bins=b)
                        logbins = np.logspace(0, np.log10(bins[-1]), len(bins))
                        # log.info(logbins)
                        plt.hist(
                            points_to_plot,
                            bins=sorted(logbins),
                            edgecolor=edgecolor,
                            color=color,
                            alpha=0.8,
                            zorder=3,
                        )
                        plt.xscale('log')
                    else:  # xscale is set to linear by default.
                        plt.hist(
                            points_to_plot,
                            edgecolor=edgecolor,
                            color=color,
                            alpha=0.8,
                            zorder=3,
                        )

                    if (
                        'x_axis_label'
                        in flag_algs_info[metric]['subplot_titles']
                    ):
                        plt.xlabel(
                            flag_algs_info[metric]['subplot_titles'][
                                'x_axis_label'
                            ]
                        )
                        plt.title(subplot_title)
                    else:
                        plt.xlabel(subplot_title)

                    plt.suptitle(str(suptitle) + ' Across Targets')

                    # get limits of graph's x-axis
                    ax = plt.gca()
                    xmin, xmax = ax.get_xlim()

                    # plot threshold lines
                    if 'thresh_vals' in flag_algs_info[metric]:

                        if subplot in flag_algs_info[metric]['thresh_vals']:
                            red_val = flag_algs_info[metric]['thresh_vals'][
                                subplot
                            ]['red']
                            yellow_val = flag_algs_info[metric]['thresh_vals'][
                                subplot
                            ]['yellow']

                        else:
                            red_val = flag_algs_info[metric]['thresh_vals'][
                                'red'
                            ]
                            yellow_val = flag_algs_info[metric]['thresh_vals'][
                                'yellow'
                            ]

                        if red_val == '-inf':
                            plt.axvline(
                                yellow_val,
                                color='#EAAA00',
                                linestyle='dashed',
                                linewidth=1,
                                zorder=2,
                            )
                            xmin, xmax = ax.get_xlim()
                            plt.axvspan(
                                xmin,
                                yellow_val,
                                color='#F6CD00',
                                alpha=0.1,
                                zorder=0,
                            )  # yellow
                            xmin, xmax = ax.get_xlim()
                            plt.axvspan(
                                yellow_val,
                                xmax,
                                color='#00DA5B',
                                alpha=0.1,
                                zorder=0,
                            )  # green

                        elif red_val == 'inf':
                            plt.axvspan(
                                xmin,
                                yellow_val,
                                color='#00DA5B',
                                alpha=0.1,
                                zorder=0,
                            )  # green
                            xmin, xmax = ax.get_xlim()
                            plt.axvspan(
                                yellow_val,
                                xmax,
                                color='#F6CD00',
                                alpha=0.1,
                                zorder=0,
                            )  # yellow

                        elif yellow_val < red_val:
                            plt.axvspan(
                                xmin,
                                yellow_val,
                                color='#00DA5B',
                                alpha=0.1,
                                zorder=0,
                            )  # green
                            plt.axvspan(
                                yellow_val,
                                red_val,
                                color='#F6CD00',
                                alpha=0.1,
                                zorder=0,
                            )  # yellow
                            xmin, xmax = ax.get_xlim()
                            plt.axvspan(
                                red_val,
                                xmax,
                                color='#C72125',
                                alpha=0.1,
                                zorder=0,
                            )  # red
                            plt.axvline(
                                red_val,
                                color='#C72125',
                                linestyle='dashed',
                                linewidth=1,
                                zorder=2,
                            )

                        elif yellow_val > red_val:
                            plt.axvspan(
                                yellow_val,
                                xmax,
                                color='#00DA5B',
                                alpha=0.1,
                                zorder=0,
                            )  # green
                            plt.axvspan(
                                red_val,
                                yellow_val,
                                color='#F6CD00',
                                alpha=0.1,
                                zorder=0,
                            )  # yellow
                            plt.axvline(
                                red_val,
                                color='#C72125',
                                linestyle='dashed',
                                linewidth=1,
                                zorder=2,
                            )
                            xmin, xmax = ax.get_xlim()
                            plt.axvspan(
                                xmin,
                                red_val,
                                color='#C72125',
                                alpha=0.1,
                                zorder=0,
                            )  # red

                        plt.axvline(
                            yellow_val,
                            color='#EAAA00',
                            linestyle='dashed',
                            linewidth=1,
                            zorder=2,
                        )

                    if metric == 'median_error':
                        plt.ylabel('Number of targets')
                    else:
                        plt.ylabel('Number of planets')

                # display the figure
                plt.tight_layout()
                save_plot_toscreen(myfig, visitor)

        pass
