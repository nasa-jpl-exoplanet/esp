'''testcerb Database Products View'''

# Heritage code shame:
# pylint: disable=too-many-nested-blocks

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
    '''General format for testcerb.simspectrum State Vector view'''

    def __init__(self, name):
        ExcaliburSV.__init__(self, name, dawgie.VERSION(1, 1, 4))

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''view ds'''
        if self['STATUS'][-1]:

            target = self['data']['target']

            for planet_letter in self['data']['planets']:
                for model in self['data']['models']:
                    if planet_letter in self['data'].keys():
                        if model in self['data'][planet_letter].keys():
                            visitor.add_image(
                                '...',
                                '------ simulated Testcerb spectrum for '
                                + target
                                + ' '
                                + planet_letter
                                + '  MODEL:'
                                + model
                                + ' ------',
                                self['data'][planet_letter][model][
                                    'plot_simspectrum'
                                ],
                            )

        return


# -------- -----------------------------------------------------------
