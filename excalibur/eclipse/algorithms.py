'''eclipse algorithms dc'''

# -- IMPORTS -- ------------------------------------------------------
import dawgie
import dawgie.context

import logging

import excalibur.system as sys
import excalibur.transit.algorithms as trnalg
import excalibur.transit.core as trncore

from importlib import import_module as fetch  # avoid cicular dependencies

log = logging.getLogger(__name__)


# ---------------------- ---------------------------------------------
# -- ALGORITHMS -- ---------------------------------------------------
# ECLIPSE CLASSES INHERIT FROM TRANSIT CLASSES
class Normalization(trnalg.Normalization):
    '''Normalize to out of transit, inherits from transit.normalization'''

    def __init__(self):
        '''__init__ ds'''
        trnalg.Normalization.__init__(self)
        self._version_ = trncore.normversion()
        self._type = 'eclipse'
        return

    pass


class WhiteLight(trnalg.WhiteLight):
    '''Create White Light Curves, inherits from transit.whitelight'''

    def __init__(self):
        '''__init__ ds'''
        trnalg.WhiteLight.__init__(self, nrm=Normalization())
        self._version_ = trncore.wlversion()
        self._type = 'eclipse'
        return

    def previous(self):
        '''Input State Vectors: eclipse.normalization, system.finalize'''
        return (
            [
                dawgie.ALG_REF(fetch('excalibur.eclipse').task, self._nrm),
                dawgie.ALG_REF(sys.task, self.__fin),
            ]
            + self.__rt.trigger('spectrum')
            + self.__rt.refs_for_proceed()
        )

    pass


class Spectrum(trnalg.Spectrum):
    '''Create emission spectrum, inherits from transit.spectrum'''

    def __init__(self):
        '''__init__ ds'''
        trnalg.Spectrum.__init__(self, nrm=Normalization(), wht=WhiteLight())
        self._version_ = trncore.spectrumversion()
        self._type = 'eclipse'
        return

    def previous(self):
        '''Input State Vectors: system.finalize, eclipse.normalization,
        eclipse.whitelight'''
        return (
            [
                dawgie.ALG_REF(sys.task, self.__fin),
                dawgie.ALG_REF(fetch('excalibur.eclipse').task, self._nrm),
                dawgie.ALG_REF(fetch('excalibur.eclipse').task, self._wht),
            ]
            + self.__rt.trigger('spectrum')
            + self.__rt.refs_for_proceed()
        )

    pass


# ---------------- ---------------------------------------------------
