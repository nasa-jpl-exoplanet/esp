'''Gemli Database Products View'''

# this might be temporary - classes are currently same as in cerb but need edits
# Heritage code shame:
# pylint: disable=duplicate-code


# -- IMPORTS -- ------------------------------------------------------

import dawgie

import excalibur
from excalibur.util.svs import ExcaliburSV

# ------------- ------------------------------------------------------
# -- SV -- -----------------------------------------------------------


# -------- -----------------------------------------------------------
class MLfitSv(ExcaliburSV):
    '''gemli.mlfit view'''

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
                        if savedresult.startswith('plot_MLfitvstruth'):
                            plotlabel = 'machine learning fit vs truth'
                        elif savedresult.startswith('plot_spectrum'):
                            plotlabel = 'cerberus-fit spectrum'
                        elif savedresult.startswith('plot_corner'):
                            plotlabel = 'cerberus corner plot'
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
