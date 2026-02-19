'''Gemli Database Products View'''

# -- IMPORTS -- ------------------------------------------------------

import dawgie

import excalibur
from excalibur.util.plotters import save_plot_toscreen
from excalibur.util.svs import ExcaliburSV

import matplotlib.image as img
import matplotlib.pyplot as plt

import os

# ------------- ------------------------------------------------------
# -- SV -- -----------------------------------------------------------


class AtmosSv(ExcaliburSV):
    '''gemli.atmos view'''

    def __init__(self, name):
        '''__init__ ds'''
        ExcaliburSV.__init__(self, name, dawgie.VERSION(1, 1, 0))

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''view ds'''
        if self['STATUS'][-1]:
            myfig = plt.figure()
            gemlilogo = img.imread(
                os.path.join(
                    excalibur.context['data_dir'], 'CERBERUS/cerberus.png'
                )
            )
            plt.imshow(gemlilogo)
            plt.axis('off')
            save_plot_toscreen(myfig, visitor, headertext='GEMLI ')
        return


# -------- -----------------------------------------------------------
class ResSv(ExcaliburSV):
    '''gemli.results view'''

    def __init__(self, name):
        '''__init__ ds'''
        ExcaliburSV.__init__(self, name, dawgie.VERSION(1, 0, 0))
        self['target'] = excalibur.ValuesList()
        self['planets'] = excalibur.ValuesList()

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''view ds'''
        if self['STATUS'][-1]:
            for target, planet_letter in zip(self['target'], self['planets']):
                for savedresult in self['data'][planet_letter].keys():
                    if 'plot' in savedresult:
                        if savedresult.startswith('plot_spectrum'):
                            plotlabel = 'best-fit spectrum'
                        elif savedresult.startswith('plot_corner'):
                            plotlabel = 'corner plot'
                        elif savedresult.startswith('plot_vsprior'):
                            plotlabel = 'improvement past prior'
                        elif savedresult.startswith('plot_walkerevol'):
                            plotlabel = 'walker evolution'
                        else:
                            plotlabel = 'unknown plottype plot'
                        if savedresult.endswith('PHOTOCHEM'):
                            plotlabel = plotlabel + ' : DISEQ MODEL'
                        else:
                            plotlabel = plotlabel + ' : TEQ MODEL'
                        textlabel = (
                            '--------- '
                            + plotlabel
                            + ' for '
                            + target
                            + ' '
                            + planet_letter
                            + ' ---------'
                        )
                        visitor.add_image(
                            '...',
                            textlabel,
                            self['data'][planet_letter][savedresult],
                        )
        return


# -------- -----------------------------------------------------------
class AnalysisSv(ExcaliburSV):
    '''PopulationSV ds'''

    def __init__(self, name):
        ExcaliburSV.__init__(self, name, dawgie.VERSION(1, 0, 0))

    def view(self, caller: excalibur.Identity, visitor: dawgie.Visitor) -> None:
        '''view ds'''
        if self['STATUS'][-1]:
            for savedresult in self['data'].keys():
                if 'plot' in savedresult:
                    if savedresult in (
                        'plot_massVmetals',
                        'plot_mass_v_metals',
                    ):
                        plotlabel = 'Planet Mass vs Metallicity'
                    elif savedresult == 'plot_fitT':
                        plotlabel = 'T_eff'
                    elif savedresult == 'plot_fitMetal':
                        plotlabel = 'Metallicity'
                    elif savedresult == 'plot_fitCO':
                        plotlabel = 'C/O'
                    elif savedresult == 'plot_fitNO':
                        plotlabel = 'N/O'
                    else:
                        plotlabel = 'unknown plottype plot'
                    # the plot titles are different for real data vs simulated
                    # use __name to decide it it's a comparison against truth
                    if savedresult in [
                        'plot_fitT',
                        'plot_fitMetal',
                        'plot_fitCO',
                        'plot_fitNO',
                    ]:
                        if 'sim' in self.name():
                            plotlabel = (
                                plotlabel + ' : retrieved vs input values'
                            )
                        else:
                            plotlabel = (
                                plotlabel
                                + ' : retrieved values and uncertainties'
                            )

                    textlabel = '--------- ' + plotlabel + ' ---------'
                    visitor.add_image(
                        '...', textlabel, self['data'][savedresult]
                    )
        return


# -------------------------------------------------------------------
