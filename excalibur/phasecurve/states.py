'''Phasecurve Database Products View'''

# Heritage code shame:
# pylint: disable=too-many-locals

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
