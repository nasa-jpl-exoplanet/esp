'''hwo Database Products View'''

# Heritage code shame:
# pylint: disable=too-many-nested-blocks,too-many-branches,too-many-locals

# -- IMPORTS -- ------------------------------------------------------
import dawgie

import excalibur
from excalibur.util.plotters import save_plot_toscreen
from excalibur.util.svs import ExcaliburSV

import os

import matplotlib.image as img
import matplotlib.pyplot as plt


# ------------- ------------------------------------------------------
# -- SV -- -----------------------------------------------------------
class SimSpectrumSV(ExcaliburSV):
    '''General format for hwo State Vector view'''

    def __init__(self, name):
        ExcaliburSV.__init__(self, name, dawgie.VERSION(1, 1, 4))

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''view ds'''
        if self['STATUS'][-1]:

            target = self['data']['target']

            for planet_letter in self['data']['planets']:
                for model in self['data']['models']:
                    # problem: for multiplanet systems, individual models may fail
                    # so there might be some missing planets here
                    # but there should be all models present, if any models present
                    if planet_letter in self['data'].keys():
                        if model in self['data'][planet_letter].keys():
                            for dictkey in self['data'][planet_letter][
                                model
                            ].keys():
                                if dictkey.startswith('plot_'):
                                    if dictkey == 'plot_depthprobed':
                                        title = '------ atmospheric depths probed for '
                                    else:
                                        title = (
                                            '------ simulated HWO spectrum for '
                                        )
                                    visitor.add_image(
                                        '...',
                                        title
                                        + target
                                        + ' '
                                        + planet_letter
                                        + '  MODEL:'
                                        + model
                                        + ' ------',
                                        self['data'][planet_letter][model][
                                            dictkey
                                        ],
                                    )

        return


# -------- -----------------------------------------------------------
