# ./excalibur/runtime/binding.py
# -*- coding: utf-8 -*-
# PyXB bindings for NM:e92452c8d3e28a9e27abfc9994d2007779e7f4c9
# Generated 2025-10-14 10:04:56.479796 by PyXB version 1.3.3 using Python 3.12.3.final.0
# Namespace AbsentNamespace0

from __future__ import unicode_literals
import pyxb
import pyxb.binding
import pyxb.binding.saxer
import io
import pyxb.utils.utility
import pyxb.utils.domutils
import sys
import pyxb.utils.sal as _six
# Unique identifier for bindings created at the same time
_GenerationUID = pyxb.utils.utility.UniqueIdentifier('urn:uuid:e8dc5ee1-a91f-11f0-9cc9-0cc47aaa0c06')

# Version of PyXB used to generate the bindings
_PyXBVersion = '1.3.3'
# Generated bindings are not compatible across PyXB versions
if pyxb.__version__ != _PyXBVersion:
    raise pyxb.PyXBVersionError(_PyXBVersion)

# A holder for module-level binding classes so we can access them from
# inside class definitions where property names may conflict.
_module_typeBindings = pyxb.utils.utility.Object()

# Import bindings for namespaces imported into schema
import pyxb.binding.datatypes

# NOTE: All namespace declarations are reserved within the binding
Namespace = pyxb.namespace.CreateAbsentNamespace()
Namespace.configureCategories(['typeBinding', 'elementBinding'])

def CreateFromDocument (xml_text, fallback_namespace=None, location_base=None, default_namespace=None):
    """Parse the given XML and use the document element to create a
    Python instance.

    @param xml_text An XML document.  This should be data (Python 2
    str or Python 3 bytes), or a text (Python 2 unicode or Python 3
    str) in the L{pyxb._InputEncoding} encoding.

    @keyword fallback_namespace An absent L{pyxb.Namespace} instance
    to use for unqualified names when there is no default namespace in
    scope.  If unspecified or C{None}, the namespace of the module
    containing this function will be used, if it is an absent
    namespace.

    @keyword location_base: An object to be recorded as the base of all
    L{pyxb.utils.utility.Location} instances associated with events and
    objects handled by the parser.  You might pass the URI from which
    the document was obtained.

    @keyword default_namespace An alias for @c fallback_namespace used
    in PyXB 1.1.4 through 1.2.6.  It behaved like a default namespace
    only for absent namespaces.
    """

    if pyxb.XMLStyle_saxer != pyxb._XMLStyle:
        dom = pyxb.utils.domutils.StringToDOM(xml_text)
        return CreateFromDOM(dom.documentElement)
    if fallback_namespace is None:
        fallback_namespace = default_namespace
    if fallback_namespace is None:
        fallback_namespace = Namespace.fallbackNamespace()
    saxer = pyxb.binding.saxer.make_parser(fallback_namespace=fallback_namespace, location_base=location_base)
    handler = saxer.getContentHandler()
    xmld = xml_text
    if isinstance(xmld, _six.text_type):
        xmld = xmld.encode(pyxb._InputEncoding)
    saxer.parse(io.BytesIO(xmld))
    instance = handler.rootObject()
    return instance

def CreateFromDOM (node, fallback_namespace=None, default_namespace=None):
    """Create a Python instance from the given DOM node.
    The node tag must correspond to an element declaration in this module.

    @deprecated: Forcing use of DOM interface is unnecessary; use L{CreateFromDocument}."""
    if fallback_namespace is None:
        fallback_namespace = default_namespace
    if fallback_namespace is None:
        fallback_namespace = Namespace.fallbackNamespace()
    return pyxb.binding.basis.element.AnyCreateFromDOM(node, fallback_namespace)


# Atomic simple type: HITEMPmolecule
class HITEMPmolecule (pyxb.binding.datatypes.string, pyxb.binding.basis.enumeration_mixin):

    """An atomic simple type."""

    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'HITEMPmolecule')
    _XSDLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 38, 2)
    _Documentation = None
HITEMPmolecule._CF_enumeration = pyxb.binding.facets.CF_enumeration(value_datatype=HITEMPmolecule, enum_prefix=None)
HITEMPmolecule.C2H2 = HITEMPmolecule._CF_enumeration.addEnumeration(unicode_value='C2H2', tag='C2H2')
HITEMPmolecule.CH3OH = HITEMPmolecule._CF_enumeration.addEnumeration(unicode_value='CH3OH', tag='CH3OH')
HITEMPmolecule.CH4 = HITEMPmolecule._CF_enumeration.addEnumeration(unicode_value='CH4', tag='CH4')
HITEMPmolecule.CO = HITEMPmolecule._CF_enumeration.addEnumeration(unicode_value='CO', tag='CO')
HITEMPmolecule.CO2 = HITEMPmolecule._CF_enumeration.addEnumeration(unicode_value='CO2', tag='CO2')
HITEMPmolecule.H2O = HITEMPmolecule._CF_enumeration.addEnumeration(unicode_value='H2O', tag='H2O')
HITEMPmolecule.HNO3 = HITEMPmolecule._CF_enumeration.addEnumeration(unicode_value='HNO3', tag='HNO3')
HITEMPmolecule.N2 = HITEMPmolecule._CF_enumeration.addEnumeration(unicode_value='N2', tag='N2')
HITEMPmolecule.N2O = HITEMPmolecule._CF_enumeration.addEnumeration(unicode_value='N2O', tag='N2O')
HITEMPmolecule.NO = HITEMPmolecule._CF_enumeration.addEnumeration(unicode_value='NO', tag='NO')
HITEMPmolecule.NO2 = HITEMPmolecule._CF_enumeration.addEnumeration(unicode_value='NO2', tag='NO2')
HITEMPmolecule.O2 = HITEMPmolecule._CF_enumeration.addEnumeration(unicode_value='O2', tag='O2')
HITEMPmolecule.O3 = HITEMPmolecule._CF_enumeration.addEnumeration(unicode_value='O3', tag='O3')
HITEMPmolecule.OH = HITEMPmolecule._CF_enumeration.addEnumeration(unicode_value='OH', tag='OH')
HITEMPmolecule._InitializeFacetMap(HITEMPmolecule._CF_enumeration)
Namespace.addCategoryObject('typeBinding', 'HITEMPmolecule', HITEMPmolecule)
_module_typeBindings.HITEMPmolecule = HITEMPmolecule

# Atomic simple type: HITRANmolecule
class HITRANmolecule (pyxb.binding.datatypes.string, pyxb.binding.basis.enumeration_mixin):

    """An atomic simple type."""

    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'HITRANmolecule')
    _XSDLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 56, 2)
    _Documentation = None
HITRANmolecule._CF_enumeration = pyxb.binding.facets.CF_enumeration(value_datatype=HITRANmolecule, enum_prefix=None)
HITRANmolecule.H2_H = HITRANmolecule._CF_enumeration.addEnumeration(unicode_value='H2-H', tag='H2_H')
HITRANmolecule.H2_H2 = HITRANmolecule._CF_enumeration.addEnumeration(unicode_value='H2-H2', tag='H2_H2')
HITRANmolecule.H2_He = HITRANmolecule._CF_enumeration.addEnumeration(unicode_value='H2-He', tag='H2_He')
HITRANmolecule.He_H = HITRANmolecule._CF_enumeration.addEnumeration(unicode_value='He-H', tag='He_H')
HITRANmolecule._InitializeFacetMap(HITRANmolecule._CF_enumeration)
Namespace.addCategoryObject('typeBinding', 'HITRANmolecule', HITRANmolecule)
_module_typeBindings.HITRANmolecule = HITRANmolecule

# Atomic simple type: EXOMOLmolecule
class EXOMOLmolecule (pyxb.binding.datatypes.string, pyxb.binding.basis.enumeration_mixin):

    """An atomic simple type."""

    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'EXOMOLmolecule')
    _XSDLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 64, 2)
    _Documentation = None
EXOMOLmolecule._CF_enumeration = pyxb.binding.facets.CF_enumeration(value_datatype=EXOMOLmolecule, enum_prefix=None)
EXOMOLmolecule.ALCL = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='ALCL', tag='ALCL')
EXOMOLmolecule.ALH = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='ALH', tag='ALH')
EXOMOLmolecule.ALO = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='ALO', tag='ALO')
EXOMOLmolecule.BEH = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='BEH', tag='BEH')
EXOMOLmolecule.C2 = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='C2', tag='C2')
EXOMOLmolecule.C2H2 = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='C2H2', tag='C2H2')
EXOMOLmolecule.C2H4 = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='C2H4', tag='C2H4')
EXOMOLmolecule.C2H6 = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='C2H6', tag='C2H6')
EXOMOLmolecule.C3 = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='C3', tag='C3')
EXOMOLmolecule.C3H8 = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='C3H8', tag='C3H8')
EXOMOLmolecule.CAO = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='CAO', tag='CAO')
EXOMOLmolecule.CAOH = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='CAOH', tag='CAOH')
EXOMOLmolecule.CH3 = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='CH3', tag='CH3')
EXOMOLmolecule.CH3CL = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='CH3CL', tag='CH3CL')
EXOMOLmolecule.CH3CHO = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='CH3CHO', tag='CH3CHO')
EXOMOLmolecule.CH4 = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='CH4', tag='CH4')
EXOMOLmolecule.CN = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='CN', tag='CN')
EXOMOLmolecule.CO = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='CO', tag='CO')
EXOMOLmolecule.CO2 = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='CO2', tag='CO2')
EXOMOLmolecule.CP = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='CP', tag='CP')
EXOMOLmolecule.CS = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='CS', tag='CS')
EXOMOLmolecule.H2 = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='H2', tag='H2')
EXOMOLmolecule.H2CO = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='H2CO', tag='H2CO')
EXOMOLmolecule.H2O = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='H2O', tag='H2O')
EXOMOLmolecule.H2S = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='H2S', tag='H2S')
EXOMOLmolecule.HBO = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='HBO', tag='HBO')
EXOMOLmolecule.HCL = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='HCL', tag='HCL')
EXOMOLmolecule.HCN = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='HCN', tag='HCN')
EXOMOLmolecule.KCL = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='KCL', tag='KCL')
EXOMOLmolecule.KOH = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='KOH', tag='KOH')
EXOMOLmolecule.LIH = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='LIH', tag='LIH')
EXOMOLmolecule.LIOH = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='LIOH', tag='LIOH')
EXOMOLmolecule.MGH = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='MGH', tag='MGH')
EXOMOLmolecule.MGO = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='MGO', tag='MGO')
EXOMOLmolecule.NACL = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='NACL', tag='NACL')
EXOMOLmolecule.NAH = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='NAH', tag='NAH')
EXOMOLmolecule.NAOH = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='NAOH', tag='NAOH')
EXOMOLmolecule.NH3 = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='NH3', tag='NH3')
EXOMOLmolecule.NO = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='NO', tag='NO')
EXOMOLmolecule.NS = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='NS', tag='NS')
EXOMOLmolecule.PH3 = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='PH3', tag='PH3')
EXOMOLmolecule.PN = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='PN', tag='PN')
EXOMOLmolecule.PO = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='PO', tag='PO')
EXOMOLmolecule.PS = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='PS', tag='PS')
EXOMOLmolecule.SIH4 = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='SIH4', tag='SIH4')
EXOMOLmolecule.SIN = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='SIN', tag='SIN')
EXOMOLmolecule.SIO2 = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='SIO2', tag='SIO2')
EXOMOLmolecule.SIS = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='SIS', tag='SIS')
EXOMOLmolecule.SO = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='SO', tag='SO')
EXOMOLmolecule.SO2 = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='SO2', tag='SO2')
EXOMOLmolecule.SO3 = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='SO3', tag='SO3')
EXOMOLmolecule.TIO = EXOMOLmolecule._CF_enumeration.addEnumeration(unicode_value='TIO', tag='TIO')
EXOMOLmolecule._InitializeFacetMap(EXOMOLmolecule._CF_enumeration)
Namespace.addCategoryObject('typeBinding', 'EXOMOLmolecule', EXOMOLmolecule)
_module_typeBindings.EXOMOLmolecule = EXOMOLmolecule

# Atomic simple type: filter_names
class filter_names (pyxb.binding.datatypes.normalizedString, pyxb.binding.basis.enumeration_mixin):

    """An atomic simple type."""

    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'filter_names')
    _XSDLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 360, 2)
    _Documentation = None
filter_names._CF_enumeration = pyxb.binding.facets.CF_enumeration(value_datatype=filter_names, enum_prefix=None)
filter_names.Ariel_sim = filter_names._CF_enumeration.addEnumeration(unicode_value='Ariel-sim', tag='Ariel_sim')
filter_names.HST_STIS_CCD_G430L_STARE = filter_names._CF_enumeration.addEnumeration(unicode_value='HST-STIS-CCD-G430L-STARE', tag='HST_STIS_CCD_G430L_STARE')
filter_names.HST_STIS_CCD_G750L_STARE = filter_names._CF_enumeration.addEnumeration(unicode_value='HST-STIS-CCD-G750L-STARE', tag='HST_STIS_CCD_G750L_STARE')
filter_names.HST_WFC3_IR_G102_SCAN = filter_names._CF_enumeration.addEnumeration(unicode_value='HST-WFC3-IR-G102-SCAN', tag='HST_WFC3_IR_G102_SCAN')
filter_names.HST_WFC3_IR_G141_SCAN = filter_names._CF_enumeration.addEnumeration(unicode_value='HST-WFC3-IR-G141-SCAN', tag='HST_WFC3_IR_G141_SCAN')
filter_names.JWST_NIRCAM_IMAGE_F210M_WLP8 = filter_names._CF_enumeration.addEnumeration(unicode_value='JWST-NIRCAM-IMAGE-F210M-WLP8', tag='JWST_NIRCAM_IMAGE_F210M_WLP8')
filter_names.JWST_NIRCAM_NRCALONG_F322W2_GRISMR = filter_names._CF_enumeration.addEnumeration(unicode_value='JWST-NIRCAM-NRCALONG-F322W2-GRISMR', tag='JWST_NIRCAM_NRCALONG_F322W2_GRISMR')
filter_names.JWST_NIRISS_NIS_CLEAR_GR700XD = filter_names._CF_enumeration.addEnumeration(unicode_value='JWST-NIRISS-NIS-CLEAR-GR700XD', tag='JWST_NIRISS_NIS_CLEAR_GR700XD')
filter_names.JWST_NIRSPEC_NRS_CLEAR_PRISM = filter_names._CF_enumeration.addEnumeration(unicode_value='JWST-NIRSPEC-NRS-CLEAR-PRISM', tag='JWST_NIRSPEC_NRS_CLEAR_PRISM')
filter_names.JWST_NIRSPEC_NRS_F290LP_G395H = filter_names._CF_enumeration.addEnumeration(unicode_value='JWST-NIRSPEC-NRS-F290LP-G395H', tag='JWST_NIRSPEC_NRS_F290LP_G395H')
filter_names.JWST_NIRSPEC_NRS_F290LP_G395M = filter_names._CF_enumeration.addEnumeration(unicode_value='JWST-NIRSPEC-NRS-F290LP-G395M', tag='JWST_NIRSPEC_NRS_F290LP_G395M')
filter_names.Spitzer_IRAC_IR_36_SUB = filter_names._CF_enumeration.addEnumeration(unicode_value='Spitzer-IRAC-IR-36-SUB', tag='Spitzer_IRAC_IR_36_SUB')
filter_names.Spitzer_IRAC_IR_45_SUB = filter_names._CF_enumeration.addEnumeration(unicode_value='Spitzer-IRAC-IR-45-SUB', tag='Spitzer_IRAC_IR_45_SUB')
filter_names._InitializeFacetMap(filter_names._CF_enumeration)
Namespace.addCategoryObject('typeBinding', 'filter_names', filter_names)
_module_typeBindings.filter_names = filter_names

# Union simple type: fitMolecule
# superclasses pyxb.binding.datatypes.anySimpleType
class fitMolecule (pyxb.binding.basis.STD_union):

    """Simple type that is a union of HITEMPmolecule, HITRANmolecule, EXOMOLmolecule."""

    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'fitMolecule')
    _XSDLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 135, 2)
    _Documentation = None

    _MemberTypes = ( HITEMPmolecule, HITRANmolecule, EXOMOLmolecule, )
fitMolecule._CF_pattern = pyxb.binding.facets.CF_pattern()
fitMolecule._CF_enumeration = pyxb.binding.facets.CF_enumeration(value_datatype=fitMolecule)
fitMolecule.C2H2 = 'C2H2'                         # originally HITEMPmolecule.C2H2
fitMolecule.CH3OH = 'CH3OH'                       # originally HITEMPmolecule.CH3OH
fitMolecule.CH4 = 'CH4'                           # originally HITEMPmolecule.CH4
fitMolecule.CO = 'CO'                             # originally HITEMPmolecule.CO
fitMolecule.CO2 = 'CO2'                           # originally HITEMPmolecule.CO2
fitMolecule.H2O = 'H2O'                           # originally HITEMPmolecule.H2O
fitMolecule.HNO3 = 'HNO3'                         # originally HITEMPmolecule.HNO3
fitMolecule.N2 = 'N2'                             # originally HITEMPmolecule.N2
fitMolecule.N2O = 'N2O'                           # originally HITEMPmolecule.N2O
fitMolecule.NO = 'NO'                             # originally HITEMPmolecule.NO
fitMolecule.NO2 = 'NO2'                           # originally HITEMPmolecule.NO2
fitMolecule.O2 = 'O2'                             # originally HITEMPmolecule.O2
fitMolecule.O3 = 'O3'                             # originally HITEMPmolecule.O3
fitMolecule.OH = 'OH'                             # originally HITEMPmolecule.OH
fitMolecule.H2_H = 'H2-H'                         # originally HITRANmolecule.H2_H
fitMolecule.H2_H2 = 'H2-H2'                       # originally HITRANmolecule.H2_H2
fitMolecule.H2_He = 'H2-He'                       # originally HITRANmolecule.H2_He
fitMolecule.He_H = 'He-H'                         # originally HITRANmolecule.He_H
fitMolecule.ALCL = 'ALCL'                         # originally EXOMOLmolecule.ALCL
fitMolecule.ALH = 'ALH'                           # originally EXOMOLmolecule.ALH
fitMolecule.ALO = 'ALO'                           # originally EXOMOLmolecule.ALO
fitMolecule.BEH = 'BEH'                           # originally EXOMOLmolecule.BEH
fitMolecule.C2 = 'C2'                             # originally EXOMOLmolecule.C2
fitMolecule.C2H2 = 'C2H2'                         # originally EXOMOLmolecule.C2H2
fitMolecule.C2H4 = 'C2H4'                         # originally EXOMOLmolecule.C2H4
fitMolecule.C2H6 = 'C2H6'                         # originally EXOMOLmolecule.C2H6
fitMolecule.C3 = 'C3'                             # originally EXOMOLmolecule.C3
fitMolecule.C3H8 = 'C3H8'                         # originally EXOMOLmolecule.C3H8
fitMolecule.CAO = 'CAO'                           # originally EXOMOLmolecule.CAO
fitMolecule.CAOH = 'CAOH'                         # originally EXOMOLmolecule.CAOH
fitMolecule.CH3 = 'CH3'                           # originally EXOMOLmolecule.CH3
fitMolecule.CH3CL = 'CH3CL'                       # originally EXOMOLmolecule.CH3CL
fitMolecule.CH3CHO = 'CH3CHO'                     # originally EXOMOLmolecule.CH3CHO
fitMolecule.CH4 = 'CH4'                           # originally EXOMOLmolecule.CH4
fitMolecule.CN = 'CN'                             # originally EXOMOLmolecule.CN
fitMolecule.CO = 'CO'                             # originally EXOMOLmolecule.CO
fitMolecule.CO2 = 'CO2'                           # originally EXOMOLmolecule.CO2
fitMolecule.CP = 'CP'                             # originally EXOMOLmolecule.CP
fitMolecule.CS = 'CS'                             # originally EXOMOLmolecule.CS
fitMolecule.H2 = 'H2'                             # originally EXOMOLmolecule.H2
fitMolecule.H2CO = 'H2CO'                         # originally EXOMOLmolecule.H2CO
fitMolecule.H2O = 'H2O'                           # originally EXOMOLmolecule.H2O
fitMolecule.H2S = 'H2S'                           # originally EXOMOLmolecule.H2S
fitMolecule.HBO = 'HBO'                           # originally EXOMOLmolecule.HBO
fitMolecule.HCL = 'HCL'                           # originally EXOMOLmolecule.HCL
fitMolecule.HCN = 'HCN'                           # originally EXOMOLmolecule.HCN
fitMolecule.KCL = 'KCL'                           # originally EXOMOLmolecule.KCL
fitMolecule.KOH = 'KOH'                           # originally EXOMOLmolecule.KOH
fitMolecule.LIH = 'LIH'                           # originally EXOMOLmolecule.LIH
fitMolecule.LIOH = 'LIOH'                         # originally EXOMOLmolecule.LIOH
fitMolecule.MGH = 'MGH'                           # originally EXOMOLmolecule.MGH
fitMolecule.MGO = 'MGO'                           # originally EXOMOLmolecule.MGO
fitMolecule.NACL = 'NACL'                         # originally EXOMOLmolecule.NACL
fitMolecule.NAH = 'NAH'                           # originally EXOMOLmolecule.NAH
fitMolecule.NAOH = 'NAOH'                         # originally EXOMOLmolecule.NAOH
fitMolecule.NH3 = 'NH3'                           # originally EXOMOLmolecule.NH3
fitMolecule.NO = 'NO'                             # originally EXOMOLmolecule.NO
fitMolecule.NS = 'NS'                             # originally EXOMOLmolecule.NS
fitMolecule.PH3 = 'PH3'                           # originally EXOMOLmolecule.PH3
fitMolecule.PN = 'PN'                             # originally EXOMOLmolecule.PN
fitMolecule.PO = 'PO'                             # originally EXOMOLmolecule.PO
fitMolecule.PS = 'PS'                             # originally EXOMOLmolecule.PS
fitMolecule.SIH4 = 'SIH4'                         # originally EXOMOLmolecule.SIH4
fitMolecule.SIN = 'SIN'                           # originally EXOMOLmolecule.SIN
fitMolecule.SIO2 = 'SIO2'                         # originally EXOMOLmolecule.SIO2
fitMolecule.SIS = 'SIS'                           # originally EXOMOLmolecule.SIS
fitMolecule.SO = 'SO'                             # originally EXOMOLmolecule.SO
fitMolecule.SO2 = 'SO2'                           # originally EXOMOLmolecule.SO2
fitMolecule.SO3 = 'SO3'                           # originally EXOMOLmolecule.SO3
fitMolecule.TIO = 'TIO'                           # originally EXOMOLmolecule.TIO
fitMolecule._InitializeFacetMap(fitMolecule._CF_pattern,
   fitMolecule._CF_enumeration)
Namespace.addCategoryObject('typeBinding', 'fitMolecule', fitMolecule)
_module_typeBindings.fitMolecule = fitMolecule

# Complex type hilo with content type EMPTY
class hilo (pyxb.binding.basis.complexTypeDefinition):
    """
      """
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_EMPTY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'hilo')
    _XSDLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 29, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Attribute lo uses Python identifier lo
    __lo = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'lo'), 'lo', '__AbsentNamespace0_hilo_lo', pyxb.binding.datatypes.float, required=True)
    __lo._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 34, 4)
    __lo._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 34, 4)
    
    lo = property(__lo.value, __lo.set, None, None)

    
    # Attribute hi uses Python identifier hi
    __hi = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'hi'), 'hi', '__AbsentNamespace0_hilo_hi', pyxb.binding.datatypes.float, required=True)
    __hi._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 35, 4)
    __hi._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 35, 4)
    
    hi = property(__hi.value, __hi.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __lo.name() : __lo,
        __hi.name() : __hi
    })
_module_typeBindings.hilo = hilo
Namespace.addCategoryObject('typeBinding', 'hilo', hilo)


# Complex type HITEMPmoleculeList with content type ELEMENT_ONLY
class HITEMPmoleculeList (pyxb.binding.basis.complexTypeDefinition):
    """Complex type HITEMPmoleculeList with content type ELEMENT_ONLY"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'HITEMPmoleculeList')
    _XSDLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 139, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element molecule uses Python identifier molecule
    __molecule = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'molecule'), 'molecule', '__AbsentNamespace0_HITEMPmoleculeList_molecule', True, pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 141, 6), )

    
    molecule = property(__molecule.value, __molecule.set, None, None)

    _ElementMap.update({
        __molecule.name() : __molecule
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.HITEMPmoleculeList = HITEMPmoleculeList
Namespace.addCategoryObject('typeBinding', 'HITEMPmoleculeList', HITEMPmoleculeList)


# Complex type HITRANmoleculeList with content type ELEMENT_ONLY
class HITRANmoleculeList (pyxb.binding.basis.complexTypeDefinition):
    """Complex type HITRANmoleculeList with content type ELEMENT_ONLY"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'HITRANmoleculeList')
    _XSDLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 144, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element molecule uses Python identifier molecule
    __molecule = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'molecule'), 'molecule', '__AbsentNamespace0_HITRANmoleculeList_molecule', True, pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 146, 6), )

    
    molecule = property(__molecule.value, __molecule.set, None, None)

    _ElementMap.update({
        __molecule.name() : __molecule
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.HITRANmoleculeList = HITRANmoleculeList
Namespace.addCategoryObject('typeBinding', 'HITRANmoleculeList', HITRANmoleculeList)


# Complex type EXOMOLmoleculeList with content type ELEMENT_ONLY
class EXOMOLmoleculeList (pyxb.binding.basis.complexTypeDefinition):
    """Complex type EXOMOLmoleculeList with content type ELEMENT_ONLY"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'EXOMOLmoleculeList')
    _XSDLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 149, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element molecule uses Python identifier molecule
    __molecule = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'molecule'), 'molecule', '__AbsentNamespace0_EXOMOLmoleculeList_molecule', True, pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 151, 6), )

    
    molecule = property(__molecule.value, __molecule.set, None, None)

    _ElementMap.update({
        __molecule.name() : __molecule
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.EXOMOLmoleculeList = EXOMOLmoleculeList
Namespace.addCategoryObject('typeBinding', 'EXOMOLmoleculeList', EXOMOLmoleculeList)


# Complex type fitMoleculeList with content type ELEMENT_ONLY
class fitMoleculeList (pyxb.binding.basis.complexTypeDefinition):
    """Complex type fitMoleculeList with content type ELEMENT_ONLY"""
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'fitMoleculeList')
    _XSDLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 154, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element molecule uses Python identifier molecule
    __molecule = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'molecule'), 'molecule', '__AbsentNamespace0_fitMoleculeList_molecule', True, pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 156, 6), )

    
    molecule = property(__molecule.value, __molecule.set, None, None)

    _ElementMap.update({
        __molecule.name() : __molecule
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.fitMoleculeList = fitMoleculeList
Namespace.addCategoryObject('typeBinding', 'fitMoleculeList', fitMoleculeList)


# Complex type control_type with content type ELEMENT_ONLY
class control_type (pyxb.binding.basis.complexTypeDefinition):
    """
        These attributes allow specific control within the AE. The attribute name
        is the task.algorithm.varName being set. The value of true means it will
        be exercised while false means it will not.
      """
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'control_type')
    _XSDLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 161, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element cerberus.crbmodel.fitmolecules uses Python identifier cerberus_crbmodel_fitmolecules
    __cerberus_crbmodel_fitmolecules = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'cerberus.crbmodel.fitmolecules'), 'cerberus_crbmodel_fitmolecules', '__AbsentNamespace0_control_type_cerberus_crbmodel_fitmolecules', False, pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 170, 6), )

    
    cerberus_crbmodel_fitmolecules = property(__cerberus_crbmodel_fitmolecules.value, __cerberus_crbmodel_fitmolecules.set, None, None)

    
    # Element cerberus.crbmodel.HITEMPmolecules uses Python identifier cerberus_crbmodel_HITEMPmolecules
    __cerberus_crbmodel_HITEMPmolecules = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'cerberus.crbmodel.HITEMPmolecules'), 'cerberus_crbmodel_HITEMPmolecules', '__AbsentNamespace0_control_type_cerberus_crbmodel_HITEMPmolecules', False, pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 172, 6), )

    
    cerberus_crbmodel_HITEMPmolecules = property(__cerberus_crbmodel_HITEMPmolecules.value, __cerberus_crbmodel_HITEMPmolecules.set, None, None)

    
    # Element cerberus.crbmodel.HITRANmolecules uses Python identifier cerberus_crbmodel_HITRANmolecules
    __cerberus_crbmodel_HITRANmolecules = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'cerberus.crbmodel.HITRANmolecules'), 'cerberus_crbmodel_HITRANmolecules', '__AbsentNamespace0_control_type_cerberus_crbmodel_HITRANmolecules', False, pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 174, 6), )

    
    cerberus_crbmodel_HITRANmolecules = property(__cerberus_crbmodel_HITRANmolecules.value, __cerberus_crbmodel_HITRANmolecules.set, None, None)

    
    # Element cerberus.crbmodel.EXOMOLmolecules uses Python identifier cerberus_crbmodel_EXOMOLmolecules
    __cerberus_crbmodel_EXOMOLmolecules = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'cerberus.crbmodel.EXOMOLmolecules'), 'cerberus_crbmodel_EXOMOLmolecules', '__AbsentNamespace0_control_type_cerberus_crbmodel_EXOMOLmolecules', False, pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 176, 6), )

    
    cerberus_crbmodel_EXOMOLmolecules = property(__cerberus_crbmodel_EXOMOLmolecules.value, __cerberus_crbmodel_EXOMOLmolecules.set, None, None)

    
    # Element cerberus.atmos.bounds.Teq uses Python identifier cerberus_atmos_bounds_Teq
    __cerberus_atmos_bounds_Teq = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.bounds.Teq'), 'cerberus_atmos_bounds_Teq', '__AbsentNamespace0_control_type_cerberus_atmos_bounds_Teq', False, pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 178, 6), )

    
    cerberus_atmos_bounds_Teq = property(__cerberus_atmos_bounds_Teq.value, __cerberus_atmos_bounds_Teq.set, None, None)

    
    # Element cerberus.atmos.bounds.abundances uses Python identifier cerberus_atmos_bounds_abundances
    __cerberus_atmos_bounds_abundances = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.bounds.abundances'), 'cerberus_atmos_bounds_abundances', '__AbsentNamespace0_control_type_cerberus_atmos_bounds_abundances', False, pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 180, 6), )

    
    cerberus_atmos_bounds_abundances = property(__cerberus_atmos_bounds_abundances.value, __cerberus_atmos_bounds_abundances.set, None, None)

    
    # Element cerberus.atmos.bounds.CTP uses Python identifier cerberus_atmos_bounds_CTP
    __cerberus_atmos_bounds_CTP = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.bounds.CTP'), 'cerberus_atmos_bounds_CTP', '__AbsentNamespace0_control_type_cerberus_atmos_bounds_CTP', False, pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 182, 6), )

    
    cerberus_atmos_bounds_CTP = property(__cerberus_atmos_bounds_CTP.value, __cerberus_atmos_bounds_CTP.set, None, None)

    
    # Element cerberus.atmos.bounds.HLoc uses Python identifier cerberus_atmos_bounds_HLoc
    __cerberus_atmos_bounds_HLoc = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.bounds.HLoc'), 'cerberus_atmos_bounds_HLoc', '__AbsentNamespace0_control_type_cerberus_atmos_bounds_HLoc', False, pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 184, 6), )

    
    cerberus_atmos_bounds_HLoc = property(__cerberus_atmos_bounds_HLoc.value, __cerberus_atmos_bounds_HLoc.set, None, None)

    
    # Element cerberus.atmos.bounds.HScale uses Python identifier cerberus_atmos_bounds_HScale
    __cerberus_atmos_bounds_HScale = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.bounds.HScale'), 'cerberus_atmos_bounds_HScale', '__AbsentNamespace0_control_type_cerberus_atmos_bounds_HScale', False, pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 186, 6), )

    
    cerberus_atmos_bounds_HScale = property(__cerberus_atmos_bounds_HScale.value, __cerberus_atmos_bounds_HScale.set, None, None)

    
    # Element cerberus.atmos.bounds.HThick uses Python identifier cerberus_atmos_bounds_HThick
    __cerberus_atmos_bounds_HThick = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.bounds.HThick'), 'cerberus_atmos_bounds_HThick', '__AbsentNamespace0_control_type_cerberus_atmos_bounds_HThick', False, pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 188, 6), )

    
    cerberus_atmos_bounds_HThick = property(__cerberus_atmos_bounds_HThick.value, __cerberus_atmos_bounds_HThick.set, None, None)

    
    # Attribute cerberus.atmos.fitCloudParameters uses Python identifier cerberus_atmos_fitCloudParameters
    __cerberus_atmos_fitCloudParameters = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.fitCloudParameters'), 'cerberus_atmos_fitCloudParameters', '__AbsentNamespace0_control_type_cerberus_atmos_fitCloudParameters', pyxb.binding.datatypes.boolean, required=True)
    __cerberus_atmos_fitCloudParameters._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 191, 4)
    __cerberus_atmos_fitCloudParameters._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 191, 4)
    
    cerberus_atmos_fitCloudParameters = property(__cerberus_atmos_fitCloudParameters.value, __cerberus_atmos_fitCloudParameters.set, None, None)

    
    # Attribute cerberus.atmos.fitNtoO uses Python identifier cerberus_atmos_fitNtoO
    __cerberus_atmos_fitNtoO = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.fitNtoO'), 'cerberus_atmos_fitNtoO', '__AbsentNamespace0_control_type_cerberus_atmos_fitNtoO', pyxb.binding.datatypes.boolean, required=True)
    __cerberus_atmos_fitNtoO._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 193, 4)
    __cerberus_atmos_fitNtoO._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 193, 4)
    
    cerberus_atmos_fitNtoO = property(__cerberus_atmos_fitNtoO.value, __cerberus_atmos_fitNtoO.set, None, None)

    
    # Attribute cerberus.atmos.fitCtoO uses Python identifier cerberus_atmos_fitCtoO
    __cerberus_atmos_fitCtoO = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.fitCtoO'), 'cerberus_atmos_fitCtoO', '__AbsentNamespace0_control_type_cerberus_atmos_fitCtoO', pyxb.binding.datatypes.boolean, required=True)
    __cerberus_atmos_fitCtoO._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 195, 4)
    __cerberus_atmos_fitCtoO._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 195, 4)
    
    cerberus_atmos_fitCtoO = property(__cerberus_atmos_fitCtoO.value, __cerberus_atmos_fitCtoO.set, None, None)

    
    # Attribute cerberus.atmos.fitT uses Python identifier cerberus_atmos_fitT
    __cerberus_atmos_fitT = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.fitT'), 'cerberus_atmos_fitT', '__AbsentNamespace0_control_type_cerberus_atmos_fitT', pyxb.binding.datatypes.boolean, required=True)
    __cerberus_atmos_fitT._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 197, 4)
    __cerberus_atmos_fitT._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 197, 4)
    
    cerberus_atmos_fitT = property(__cerberus_atmos_fitT.value, __cerberus_atmos_fitT.set, None, None)

    
    # Attribute cerberus.atmos.sliceSampler uses Python identifier cerberus_atmos_sliceSampler
    __cerberus_atmos_sliceSampler = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.sliceSampler'), 'cerberus_atmos_sliceSampler', '__AbsentNamespace0_control_type_cerberus_atmos_sliceSampler', pyxb.binding.datatypes.boolean, required=True)
    __cerberus_atmos_sliceSampler._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 199, 4)
    __cerberus_atmos_sliceSampler._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 199, 4)
    
    cerberus_atmos_sliceSampler = property(__cerberus_atmos_sliceSampler.value, __cerberus_atmos_sliceSampler.set, None, None)

    
    # Attribute cerberus.crbmodel.lbroadening uses Python identifier cerberus_crbmodel_lbroadening
    __cerberus_crbmodel_lbroadening = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'cerberus.crbmodel.lbroadening'), 'cerberus_crbmodel_lbroadening', '__AbsentNamespace0_control_type_cerberus_crbmodel_lbroadening', pyxb.binding.datatypes.boolean, required=True)
    __cerberus_crbmodel_lbroadening._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 201, 4)
    __cerberus_crbmodel_lbroadening._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 201, 4)
    
    cerberus_crbmodel_lbroadening = property(__cerberus_crbmodel_lbroadening.value, __cerberus_crbmodel_lbroadening.set, None, None)

    
    # Attribute cerberus.crbmodel.lshifting uses Python identifier cerberus_crbmodel_lshifting
    __cerberus_crbmodel_lshifting = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'cerberus.crbmodel.lshifting'), 'cerberus_crbmodel_lshifting', '__AbsentNamespace0_control_type_cerberus_crbmodel_lshifting', pyxb.binding.datatypes.boolean, required=True)
    __cerberus_crbmodel_lshifting._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 203, 4)
    __cerberus_crbmodel_lshifting._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 203, 4)
    
    cerberus_crbmodel_lshifting = property(__cerberus_crbmodel_lshifting.value, __cerberus_crbmodel_lshifting.set, None, None)

    
    # Attribute cerberus.crbmodel.isothermal uses Python identifier cerberus_crbmodel_isothermal
    __cerberus_crbmodel_isothermal = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'cerberus.crbmodel.isothermal'), 'cerberus_crbmodel_isothermal', '__AbsentNamespace0_control_type_cerberus_crbmodel_isothermal', pyxb.binding.datatypes.boolean, required=True)
    __cerberus_crbmodel_isothermal._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 205, 4)
    __cerberus_crbmodel_isothermal._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 205, 4)
    
    cerberus_crbmodel_isothermal = property(__cerberus_crbmodel_isothermal.value, __cerberus_crbmodel_isothermal.set, None, None)

    
    # Attribute cerberus.crbmodel.nlevels uses Python identifier cerberus_crbmodel_nlevels
    __cerberus_crbmodel_nlevels = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'cerberus.crbmodel.nlevels'), 'cerberus_crbmodel_nlevels', '__AbsentNamespace0_control_type_cerberus_crbmodel_nlevels', pyxb.binding.datatypes.integer, required=True)
    __cerberus_crbmodel_nlevels._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 207, 4)
    __cerberus_crbmodel_nlevels._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 207, 4)
    
    cerberus_crbmodel_nlevels = property(__cerberus_crbmodel_nlevels.value, __cerberus_crbmodel_nlevels.set, None, None)

    
    # Attribute cerberus.crbmodel.solrad uses Python identifier cerberus_crbmodel_solrad
    __cerberus_crbmodel_solrad = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'cerberus.crbmodel.solrad'), 'cerberus_crbmodel_solrad', '__AbsentNamespace0_control_type_cerberus_crbmodel_solrad', pyxb.binding.datatypes.float, required=True)
    __cerberus_crbmodel_solrad._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 209, 4)
    __cerberus_crbmodel_solrad._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 209, 4)
    
    cerberus_crbmodel_solrad = property(__cerberus_crbmodel_solrad.value, __cerberus_crbmodel_solrad.set, None, None)

    
    # Attribute cerberus.crbmodel.Hsmax uses Python identifier cerberus_crbmodel_Hsmax
    __cerberus_crbmodel_Hsmax = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'cerberus.crbmodel.Hsmax'), 'cerberus_crbmodel_Hsmax', '__AbsentNamespace0_control_type_cerberus_crbmodel_Hsmax', pyxb.binding.datatypes.integer, required=True)
    __cerberus_crbmodel_Hsmax._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 211, 4)
    __cerberus_crbmodel_Hsmax._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 211, 4)
    
    cerberus_crbmodel_Hsmax = property(__cerberus_crbmodel_Hsmax.value, __cerberus_crbmodel_Hsmax.set, None, None)

    
    # Attribute cerberus.plotters.cornerBins uses Python identifier cerberus_plotters_cornerBins
    __cerberus_plotters_cornerBins = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'cerberus.plotters.cornerBins'), 'cerberus_plotters_cornerBins', '__AbsentNamespace0_control_type_cerberus_plotters_cornerBins', pyxb.binding.datatypes.integer, required=True)
    __cerberus_plotters_cornerBins._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 213, 4)
    __cerberus_plotters_cornerBins._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 213, 4)
    
    cerberus_plotters_cornerBins = property(__cerberus_plotters_cornerBins.value, __cerberus_plotters_cornerBins.set, None, None)

    
    # Attribute cerberus.results.nrandomwalkers uses Python identifier cerberus_results_nrandomwalkers
    __cerberus_results_nrandomwalkers = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'cerberus.results.nrandomwalkers'), 'cerberus_results_nrandomwalkers', '__AbsentNamespace0_control_type_cerberus_results_nrandomwalkers', pyxb.binding.datatypes.integer, required=True)
    __cerberus_results_nrandomwalkers._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 215, 4)
    __cerberus_results_nrandomwalkers._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 215, 4)
    
    cerberus_results_nrandomwalkers = property(__cerberus_results_nrandomwalkers.value, __cerberus_results_nrandomwalkers.set, None, None)

    
    # Attribute cerberus.results.randomseed uses Python identifier cerberus_results_randomseed
    __cerberus_results_randomseed = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'cerberus.results.randomseed'), 'cerberus_results_randomseed', '__AbsentNamespace0_control_type_cerberus_results_randomseed', pyxb.binding.datatypes.integer, required=True)
    __cerberus_results_randomseed._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 217, 4)
    __cerberus_results_randomseed._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 217, 4)
    
    cerberus_results_randomseed = property(__cerberus_results_randomseed.value, __cerberus_results_randomseed.set, None, None)

    
    # Attribute system.validate.selectMostRecent uses Python identifier system_validate_selectMostRecent
    __system_validate_selectMostRecent = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'system.validate.selectMostRecent'), 'system_validate_selectMostRecent', '__AbsentNamespace0_control_type_system_validate_selectMostRecent', pyxb.binding.datatypes.boolean, required=True)
    __system_validate_selectMostRecent._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 219, 4)
    __system_validate_selectMostRecent._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 219, 4)
    
    system_validate_selectMostRecent = property(__system_validate_selectMostRecent.value, __system_validate_selectMostRecent.set, None, None)

    
    # Attribute system.validate.maximizeSelfConsistency uses Python identifier system_validate_maximizeSelfConsistency
    __system_validate_maximizeSelfConsistency = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'system.validate.maximizeSelfConsistency'), 'system_validate_maximizeSelfConsistency', '__AbsentNamespace0_control_type_system_validate_maximizeSelfConsistency', pyxb.binding.datatypes.boolean, required=True)
    __system_validate_maximizeSelfConsistency._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 221, 4)
    __system_validate_maximizeSelfConsistency._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 221, 4)
    
    system_validate_maximizeSelfConsistency = property(__system_validate_maximizeSelfConsistency.value, __system_validate_maximizeSelfConsistency.set, None, None)

    
    # Attribute ariel.simspectrum.includeMetallicityDispersion uses Python identifier ariel_simspectrum_includeMetallicityDispersion
    __ariel_simspectrum_includeMetallicityDispersion = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'ariel.simspectrum.includeMetallicityDispersion'), 'ariel_simspectrum_includeMetallicityDispersion', '__AbsentNamespace0_control_type_ariel_simspectrum_includeMetallicityDispersion', pyxb.binding.datatypes.boolean, required=True)
    __ariel_simspectrum_includeMetallicityDispersion._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 223, 4)
    __ariel_simspectrum_includeMetallicityDispersion._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 223, 4)
    
    ariel_simspectrum_includeMetallicityDispersion = property(__ariel_simspectrum_includeMetallicityDispersion.value, __ariel_simspectrum_includeMetallicityDispersion.set, None, None)

    
    # Attribute ariel.simspectrum.randomCloudProperties uses Python identifier ariel_simspectrum_randomCloudProperties
    __ariel_simspectrum_randomCloudProperties = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'ariel.simspectrum.randomCloudProperties'), 'ariel_simspectrum_randomCloudProperties', '__AbsentNamespace0_control_type_ariel_simspectrum_randomCloudProperties', pyxb.binding.datatypes.boolean, required=True)
    __ariel_simspectrum_randomCloudProperties._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 225, 4)
    __ariel_simspectrum_randomCloudProperties._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 225, 4)
    
    ariel_simspectrum_randomCloudProperties = property(__ariel_simspectrum_randomCloudProperties.value, __ariel_simspectrum_randomCloudProperties.set, None, None)

    
    # Attribute ariel.simspectrum.SNRadjustment uses Python identifier ariel_simspectrum_SNRadjustment
    __ariel_simspectrum_SNRadjustment = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'ariel.simspectrum.SNRadjustment'), 'ariel_simspectrum_SNRadjustment', '__AbsentNamespace0_control_type_ariel_simspectrum_SNRadjustment', pyxb.binding.datatypes.float, required=True)
    __ariel_simspectrum_SNRadjustment._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 227, 4)
    __ariel_simspectrum_SNRadjustment._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 227, 4)
    
    ariel_simspectrum_SNRadjustment = property(__ariel_simspectrum_SNRadjustment.value, __ariel_simspectrum_SNRadjustment.set, None, None)

    
    # Attribute ariel.simspectrum.thorngrenMassMetals uses Python identifier ariel_simspectrum_thorngrenMassMetals
    __ariel_simspectrum_thorngrenMassMetals = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'ariel.simspectrum.thorngrenMassMetals'), 'ariel_simspectrum_thorngrenMassMetals', '__AbsentNamespace0_control_type_ariel_simspectrum_thorngrenMassMetals', pyxb.binding.datatypes.boolean, required=True)
    __ariel_simspectrum_thorngrenMassMetals._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 229, 4)
    __ariel_simspectrum_thorngrenMassMetals._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 229, 4)
    
    ariel_simspectrum_thorngrenMassMetals = property(__ariel_simspectrum_thorngrenMassMetals.value, __ariel_simspectrum_thorngrenMassMetals.set, None, None)

    
    # Attribute ariel.simspectrum.chachanMassMetals uses Python identifier ariel_simspectrum_chachanMassMetals
    __ariel_simspectrum_chachanMassMetals = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'ariel.simspectrum.chachanMassMetals'), 'ariel_simspectrum_chachanMassMetals', '__AbsentNamespace0_control_type_ariel_simspectrum_chachanMassMetals', pyxb.binding.datatypes.boolean, required=True)
    __ariel_simspectrum_chachanMassMetals._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 231, 4)
    __ariel_simspectrum_chachanMassMetals._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 231, 4)
    
    ariel_simspectrum_chachanMassMetals = property(__ariel_simspectrum_chachanMassMetals.value, __ariel_simspectrum_chachanMassMetals.set, None, None)

    
    # Attribute ariel.simspectrum.tier uses Python identifier ariel_simspectrum_tier
    __ariel_simspectrum_tier = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'ariel.simspectrum.tier'), 'ariel_simspectrum_tier', '__AbsentNamespace0_control_type_ariel_simspectrum_tier', pyxb.binding.datatypes.integer, required=True)
    __ariel_simspectrum_tier._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 233, 4)
    __ariel_simspectrum_tier._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 233, 4)
    
    ariel_simspectrum_tier = property(__ariel_simspectrum_tier.value, __ariel_simspectrum_tier.set, None, None)

    
    # Attribute ariel.simspectrum.randomseed uses Python identifier ariel_simspectrum_randomseed
    __ariel_simspectrum_randomseed = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'ariel.simspectrum.randomseed'), 'ariel_simspectrum_randomseed', '__AbsentNamespace0_control_type_ariel_simspectrum_randomseed', pyxb.binding.datatypes.integer, required=True)
    __ariel_simspectrum_randomseed._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 235, 4)
    __ariel_simspectrum_randomseed._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 235, 4)
    
    ariel_simspectrum_randomseed = property(__ariel_simspectrum_randomseed.value, __ariel_simspectrum_randomseed.set, None, None)

    
    # Attribute ariel.simspectrum.metallicityDispersion uses Python identifier ariel_simspectrum_metallicityDispersion
    __ariel_simspectrum_metallicityDispersion = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'ariel.simspectrum.metallicityDispersion'), 'ariel_simspectrum_metallicityDispersion', '__AbsentNamespace0_control_type_ariel_simspectrum_metallicityDispersion', pyxb.binding.datatypes.float, required=True)
    __ariel_simspectrum_metallicityDispersion._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 237, 4)
    __ariel_simspectrum_metallicityDispersion._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 237, 4)
    
    ariel_simspectrum_metallicityDispersion = property(__ariel_simspectrum_metallicityDispersion.value, __ariel_simspectrum_metallicityDispersion.set, None, None)

    
    # Attribute ariel.simspectrum.CtoOdaSilva uses Python identifier ariel_simspectrum_CtoOdaSilva
    __ariel_simspectrum_CtoOdaSilva = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'ariel.simspectrum.CtoOdaSilva'), 'ariel_simspectrum_CtoOdaSilva', '__AbsentNamespace0_control_type_ariel_simspectrum_CtoOdaSilva', pyxb.binding.datatypes.boolean, required=True)
    __ariel_simspectrum_CtoOdaSilva._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 239, 4)
    __ariel_simspectrum_CtoOdaSilva._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 239, 4)
    
    ariel_simspectrum_CtoOdaSilva = property(__ariel_simspectrum_CtoOdaSilva.value, __ariel_simspectrum_CtoOdaSilva.set, None, None)

    
    # Attribute ariel.simspectrum.CtoOaverage uses Python identifier ariel_simspectrum_CtoOaverage
    __ariel_simspectrum_CtoOaverage = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'ariel.simspectrum.CtoOaverage'), 'ariel_simspectrum_CtoOaverage', '__AbsentNamespace0_control_type_ariel_simspectrum_CtoOaverage', pyxb.binding.datatypes.float, required=True)
    __ariel_simspectrum_CtoOaverage._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 241, 4)
    __ariel_simspectrum_CtoOaverage._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 241, 4)
    
    ariel_simspectrum_CtoOaverage = property(__ariel_simspectrum_CtoOaverage.value, __ariel_simspectrum_CtoOaverage.set, None, None)

    
    # Attribute ariel.simspectrum.CtoOdispersion uses Python identifier ariel_simspectrum_CtoOdispersion
    __ariel_simspectrum_CtoOdispersion = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'ariel.simspectrum.CtoOdispersion'), 'ariel_simspectrum_CtoOdispersion', '__AbsentNamespace0_control_type_ariel_simspectrum_CtoOdispersion', pyxb.binding.datatypes.float, required=True)
    __ariel_simspectrum_CtoOdispersion._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 243, 4)
    __ariel_simspectrum_CtoOdispersion._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 243, 4)
    
    ariel_simspectrum_CtoOdispersion = property(__ariel_simspectrum_CtoOdispersion.value, __ariel_simspectrum_CtoOdispersion.set, None, None)

    
    # Attribute transit.pymc.sliceSampler uses Python identifier transit_pymc_sliceSampler
    __transit_pymc_sliceSampler = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'transit.pymc.sliceSampler'), 'transit_pymc_sliceSampler', '__AbsentNamespace0_control_type_transit_pymc_sliceSampler', pyxb.binding.datatypes.boolean, required=True)
    __transit_pymc_sliceSampler._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 245, 4)
    __transit_pymc_sliceSampler._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 245, 4)
    
    transit_pymc_sliceSampler = property(__transit_pymc_sliceSampler.value, __transit_pymc_sliceSampler.set, None, None)

    
    # Attribute selftest.Nrepeats uses Python identifier selftest_Nrepeats
    __selftest_Nrepeats = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'selftest.Nrepeats'), 'selftest_Nrepeats', '__AbsentNamespace0_control_type_selftest_Nrepeats', pyxb.binding.datatypes.integer, required=True)
    __selftest_Nrepeats._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 247, 4)
    __selftest_Nrepeats._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 247, 4)
    
    selftest_Nrepeats = property(__selftest_Nrepeats.value, __selftest_Nrepeats.set, None, None)

    _ElementMap.update({
        __cerberus_crbmodel_fitmolecules.name() : __cerberus_crbmodel_fitmolecules,
        __cerberus_crbmodel_HITEMPmolecules.name() : __cerberus_crbmodel_HITEMPmolecules,
        __cerberus_crbmodel_HITRANmolecules.name() : __cerberus_crbmodel_HITRANmolecules,
        __cerberus_crbmodel_EXOMOLmolecules.name() : __cerberus_crbmodel_EXOMOLmolecules,
        __cerberus_atmos_bounds_Teq.name() : __cerberus_atmos_bounds_Teq,
        __cerberus_atmos_bounds_abundances.name() : __cerberus_atmos_bounds_abundances,
        __cerberus_atmos_bounds_CTP.name() : __cerberus_atmos_bounds_CTP,
        __cerberus_atmos_bounds_HLoc.name() : __cerberus_atmos_bounds_HLoc,
        __cerberus_atmos_bounds_HScale.name() : __cerberus_atmos_bounds_HScale,
        __cerberus_atmos_bounds_HThick.name() : __cerberus_atmos_bounds_HThick
    })
    _AttributeMap.update({
        __cerberus_atmos_fitCloudParameters.name() : __cerberus_atmos_fitCloudParameters,
        __cerberus_atmos_fitNtoO.name() : __cerberus_atmos_fitNtoO,
        __cerberus_atmos_fitCtoO.name() : __cerberus_atmos_fitCtoO,
        __cerberus_atmos_fitT.name() : __cerberus_atmos_fitT,
        __cerberus_atmos_sliceSampler.name() : __cerberus_atmos_sliceSampler,
        __cerberus_crbmodel_lbroadening.name() : __cerberus_crbmodel_lbroadening,
        __cerberus_crbmodel_lshifting.name() : __cerberus_crbmodel_lshifting,
        __cerberus_crbmodel_isothermal.name() : __cerberus_crbmodel_isothermal,
        __cerberus_crbmodel_nlevels.name() : __cerberus_crbmodel_nlevels,
        __cerberus_crbmodel_solrad.name() : __cerberus_crbmodel_solrad,
        __cerberus_crbmodel_Hsmax.name() : __cerberus_crbmodel_Hsmax,
        __cerberus_plotters_cornerBins.name() : __cerberus_plotters_cornerBins,
        __cerberus_results_nrandomwalkers.name() : __cerberus_results_nrandomwalkers,
        __cerberus_results_randomseed.name() : __cerberus_results_randomseed,
        __system_validate_selectMostRecent.name() : __system_validate_selectMostRecent,
        __system_validate_maximizeSelfConsistency.name() : __system_validate_maximizeSelfConsistency,
        __ariel_simspectrum_includeMetallicityDispersion.name() : __ariel_simspectrum_includeMetallicityDispersion,
        __ariel_simspectrum_randomCloudProperties.name() : __ariel_simspectrum_randomCloudProperties,
        __ariel_simspectrum_SNRadjustment.name() : __ariel_simspectrum_SNRadjustment,
        __ariel_simspectrum_thorngrenMassMetals.name() : __ariel_simspectrum_thorngrenMassMetals,
        __ariel_simspectrum_chachanMassMetals.name() : __ariel_simspectrum_chachanMassMetals,
        __ariel_simspectrum_tier.name() : __ariel_simspectrum_tier,
        __ariel_simspectrum_randomseed.name() : __ariel_simspectrum_randomseed,
        __ariel_simspectrum_metallicityDispersion.name() : __ariel_simspectrum_metallicityDispersion,
        __ariel_simspectrum_CtoOdaSilva.name() : __ariel_simspectrum_CtoOdaSilva,
        __ariel_simspectrum_CtoOaverage.name() : __ariel_simspectrum_CtoOaverage,
        __ariel_simspectrum_CtoOdispersion.name() : __ariel_simspectrum_CtoOdispersion,
        __transit_pymc_sliceSampler.name() : __transit_pymc_sliceSampler,
        __selftest_Nrepeats.name() : __selftest_Nrepeats
    })
_module_typeBindings.control_type = control_type
Namespace.addCategoryObject('typeBinding', 'control_type', control_type)


# Complex type filter_type with content type ELEMENT_ONLY
class filter_type (pyxb.binding.basis.complexTypeDefinition):
    """
        The filters are used to determine which generic state vectors should be
        computed and which should be ignored.

        If there are no includes, then all filters will be processed unless
        explicitly named with an exclude.

        If there are includes, then only those in the list will be processed
        unless explicity named with an exclude.
      """
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'filter_type')
    _XSDLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 251, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element exclude uses Python identifier exclude
    __exclude = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'exclude'), 'exclude', '__AbsentNamespace0_filter_type_exclude', True, pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 265, 6), )

    
    exclude = property(__exclude.value, __exclude.set, None, None)

    
    # Element include uses Python identifier include
    __include = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'include'), 'include', '__AbsentNamespace0_filter_type_include', True, pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 267, 6), )

    
    include = property(__include.value, __include.set, None, None)

    _ElementMap.update({
        __exclude.name() : __exclude,
        __include.name() : __include
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.filter_type = filter_type
Namespace.addCategoryObject('typeBinding', 'filter_type', filter_type)


# Complex type lever_type with content type ELEMENT_ONLY
class lever_type (pyxb.binding.basis.complexTypeDefinition):
    """
        The levers of power (controls) available to the scientist
        are limited to:
          
        controls: these are used by the AE to make decisions on
                  what processing should be done.
        filters: controls what filters should be processed and/or ignored.
        pymc: controls specific to the PYMC processing.
        restrict_to: when empty do all targets, otherwise just those listed
        sequester: never do these targets dispite restrict_to
      """
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'lever_type')
    _XSDLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 272, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element controls uses Python identifier controls
    __controls = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'controls'), 'controls', '__AbsentNamespace0_lever_type_controls', False, pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 287, 6), )

    
    controls = property(__controls.value, __controls.set, None, None)

    
    # Element filters uses Python identifier filters
    __filters = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'filters'), 'filters', '__AbsentNamespace0_lever_type_filters', False, pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 288, 6), )

    
    filters = property(__filters.value, __filters.set, None, None)

    
    # Element pymc uses Python identifier pymc
    __pymc = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'pymc'), 'pymc', '__AbsentNamespace0_lever_type_pymc', False, pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 289, 6), )

    
    pymc = property(__pymc.value, __pymc.set, None, None)

    
    # Element run_only uses Python identifier run_only
    __run_only = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'run_only'), 'run_only', '__AbsentNamespace0_lever_type_run_only', False, pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 290, 6), )

    
    run_only = property(__run_only.value, __run_only.set, None, None)

    
    # Element sequester uses Python identifier sequester
    __sequester = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'sequester'), 'sequester', '__AbsentNamespace0_lever_type_sequester', False, pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 291, 6), )

    
    sequester = property(__sequester.value, __sequester.set, None, None)

    
    # Attribute index uses Python identifier index
    __index = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'index'), 'index', '__AbsentNamespace0_lever_type_index', pyxb.binding.datatypes.normalizedString, unicode_default='registry')
    __index._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 293, 4)
    __index._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 293, 4)
    
    index = property(__index.value, __index.set, None, None)

    _ElementMap.update({
        __controls.name() : __controls,
        __filters.name() : __filters,
        __pymc.name() : __pymc,
        __run_only.name() : __run_only,
        __sequester.name() : __sequester
    })
    _AttributeMap.update({
        __index.name() : __index
    })
_module_typeBindings.lever_type = lever_type
Namespace.addCategoryObject('typeBinding', 'lever_type', lever_type)


# Complex type pymc_count_type with content type ELEMENT_ONLY
class pymc_count_type (pyxb.binding.basis.complexTypeDefinition):
    """
        Specify how many iterations for MCMC. The default is used for
        all unspecified targets. Targets can also be individually specified.
      """
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'pymc_count_type')
    _XSDLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 296, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element target uses Python identifier target
    __target = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'target'), 'target', '__AbsentNamespace0_pymc_count_type_target', True, pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 304, 6), )

    
    target = property(__target.value, __target.set, None, None)

    
    # Attribute default uses Python identifier default
    __default = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'default'), 'default', '__AbsentNamespace0_pymc_count_type_default', pyxb.binding.datatypes.positiveInteger, required=True)
    __default._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 307, 4)
    __default._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 307, 4)
    
    default = property(__default.value, __default.set, None, None)

    _ElementMap.update({
        __target.name() : __target
    })
    _AttributeMap.update({
        __default.name() : __default
    })
_module_typeBindings.pymc_count_type = pymc_count_type
Namespace.addCategoryObject('typeBinding', 'pymc_count_type', pymc_count_type)


# Complex type pymc_type with content type ELEMENT_ONLY
class pymc_type (pyxb.binding.basis.complexTypeDefinition):
    """
        Control the max number of iterations for PYMC.
      """
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'pymc_type')
    _XSDLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 310, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element cerberuschains uses Python identifier cerberuschains
    __cerberuschains = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'cerberuschains'), 'cerberuschains', '__AbsentNamespace0_pymc_type_cerberuschains', False, pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 317, 6), )

    
    cerberuschains = property(__cerberuschains.value, __cerberuschains.set, None, None)

    
    # Element cerberuschainlen uses Python identifier cerberuschainlen
    __cerberuschainlen = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'cerberuschainlen'), 'cerberuschainlen', '__AbsentNamespace0_pymc_type_cerberuschainlen', False, pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 318, 6), )

    
    cerberuschainlen = property(__cerberuschainlen.value, __cerberuschainlen.set, None, None)

    
    # Element spectrumchains uses Python identifier spectrumchains
    __spectrumchains = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'spectrumchains'), 'spectrumchains', '__AbsentNamespace0_pymc_type_spectrumchains', False, pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 319, 6), )

    
    spectrumchains = property(__spectrumchains.value, __spectrumchains.set, None, None)

    
    # Element spectrumchainlen uses Python identifier spectrumchainlen
    __spectrumchainlen = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'spectrumchainlen'), 'spectrumchainlen', '__AbsentNamespace0_pymc_type_spectrumchainlen', False, pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 320, 6), )

    
    spectrumchainlen = property(__spectrumchainlen.value, __spectrumchainlen.set, None, None)

    _ElementMap.update({
        __cerberuschains.name() : __cerberuschains,
        __cerberuschainlen.name() : __cerberuschainlen,
        __spectrumchains.name() : __spectrumchains,
        __spectrumchainlen.name() : __spectrumchainlen
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.pymc_type = pymc_type
Namespace.addCategoryObject('typeBinding', 'pymc_type', pymc_type)


# Complex type sequester_type with content type ELEMENT_ONLY
class sequester_type (pyxb.binding.basis.complexTypeDefinition):
    """
        Sequester the targets to the island of no processing.
      """
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'sequester_type')
    _XSDLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 324, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element target uses Python identifier target
    __target = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'target'), 'target', '__AbsentNamespace0_sequester_type_target', True, pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 331, 6), )

    
    target = property(__target.value, __target.set, None, None)

    _ElementMap.update({
        __target.name() : __target
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.sequester_type = sequester_type
Namespace.addCategoryObject('typeBinding', 'sequester_type', sequester_type)


# Complex type target_type with content type SIMPLE
class target_type (pyxb.binding.basis.complexTypeDefinition):
    """
        The name of the target to not process and the reason why.
      """
    _TypeDefinition = pyxb.binding.datatypes.normalizedString
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_SIMPLE
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'target_type')
    _XSDLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 336, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.normalizedString
    
    # Attribute because uses Python identifier because
    __because = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'because'), 'because', '__AbsentNamespace0_target_type_because', pyxb.binding.datatypes.normalizedString, required=True)
    __because._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 344, 8)
    __because._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 344, 8)
    
    because = property(__because.value, __because.set, None, None)

    
    # Attribute isRegex uses Python identifier isRegex
    __isRegex = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'isRegex'), 'isRegex', '__AbsentNamespace0_target_type_isRegex', pyxb.binding.datatypes.boolean, unicode_default='false')
    __isRegex._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 345, 8)
    __isRegex._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 345, 8)
    
    isRegex = property(__isRegex.value, __isRegex.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __because.name() : __because,
        __isRegex.name() : __isRegex
    })
_module_typeBindings.target_type = target_type
Namespace.addCategoryObject('typeBinding', 'target_type', target_type)


# Complex type target_override_type with content type EMPTY
class target_override_type (pyxb.binding.basis.complexTypeDefinition):
    """
        Allows a specific target to be named and a given number of iterations.
      """
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_EMPTY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'target_override_type')
    _XSDLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 350, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Attribute name uses Python identifier name
    __name = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'name'), 'name', '__AbsentNamespace0_target_override_type_name', pyxb.binding.datatypes.normalizedString, required=True)
    __name._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 356, 4)
    __name._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 356, 4)
    
    name = property(__name.value, __name.set, None, None)

    
    # Attribute steps uses Python identifier steps
    __steps = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'steps'), 'steps', '__AbsentNamespace0_target_override_type_steps', pyxb.binding.datatypes.positiveInteger, required=True)
    __steps._DeclarationLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 357, 4)
    __steps._UseLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 357, 4)
    
    steps = property(__steps.value, __steps.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __name.name() : __name,
        __steps.name() : __steps
    })
_module_typeBindings.target_override_type = target_override_type
Namespace.addCategoryObject('typeBinding', 'target_override_type', target_override_type)


# Complex type [anonymous] with content type ELEMENT_ONLY
class CTD_ANON (lever_type):
    """
          This terrible construct is so that xjc can autodetect this as the
          root node for processing. Many things would be better but this is
          the most workable solution especially if the making of the binding
          code is automated in the pom. The only other real solution is to
          modify one of the classes generated by hand.
        """
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_ELEMENT_ONLY
    _Abstract = False
    _ExpandedName = None
    _XSDLocation = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 379, 4)
    _ElementMap = lever_type._ElementMap.copy()
    _AttributeMap = lever_type._AttributeMap.copy()
    # Base type is lever_type
    
    # Element controls (controls) inherited from lever_type
    
    # Element filters (filters) inherited from lever_type
    
    # Element pymc (pymc) inherited from lever_type
    
    # Element run_only (run_only) inherited from lever_type
    
    # Element sequester (sequester) inherited from lever_type
    
    # Attribute index inherited from lever_type
    _ElementMap.update({
        
    })
    _AttributeMap.update({
        
    })
_module_typeBindings.CTD_ANON = CTD_ANON


levers = pyxb.binding.basis.element(pyxb.namespace.ExpandedName(Namespace, 'levers'), CTD_ANON, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 378, 2))
Namespace.addCategoryObject('elementBinding', levers.name().localName(), levers)



HITEMPmoleculeList._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'molecule'), HITEMPmolecule, scope=HITEMPmoleculeList, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 141, 6)))

def _BuildAutomaton ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton
    del _BuildAutomaton
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 141, 6))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(HITEMPmoleculeList._UseForTag(pyxb.namespace.ExpandedName(None, 'molecule')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 141, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=None)
HITEMPmoleculeList._Automaton = _BuildAutomaton()




HITRANmoleculeList._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'molecule'), HITRANmolecule, scope=HITRANmoleculeList, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 146, 6)))

def _BuildAutomaton_ ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_
    del _BuildAutomaton_
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 146, 6))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(HITRANmoleculeList._UseForTag(pyxb.namespace.ExpandedName(None, 'molecule')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 146, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=None)
HITRANmoleculeList._Automaton = _BuildAutomaton_()




EXOMOLmoleculeList._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'molecule'), EXOMOLmolecule, scope=EXOMOLmoleculeList, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 151, 6)))

def _BuildAutomaton_2 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_2
    del _BuildAutomaton_2
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 151, 6))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(EXOMOLmoleculeList._UseForTag(pyxb.namespace.ExpandedName(None, 'molecule')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 151, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=None)
EXOMOLmoleculeList._Automaton = _BuildAutomaton_2()




fitMoleculeList._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'molecule'), fitMolecule, scope=fitMoleculeList, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 156, 6)))

def _BuildAutomaton_3 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_3
    del _BuildAutomaton_3
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(fitMoleculeList._UseForTag(pyxb.namespace.ExpandedName(None, 'molecule')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 156, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
         ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
fitMoleculeList._Automaton = _BuildAutomaton_3()




control_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'cerberus.crbmodel.fitmolecules'), fitMoleculeList, scope=control_type, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 170, 6)))

control_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'cerberus.crbmodel.HITEMPmolecules'), HITEMPmoleculeList, scope=control_type, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 172, 6)))

control_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'cerberus.crbmodel.HITRANmolecules'), HITRANmoleculeList, scope=control_type, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 174, 6)))

control_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'cerberus.crbmodel.EXOMOLmolecules'), EXOMOLmoleculeList, scope=control_type, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 176, 6)))

control_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.bounds.Teq'), hilo, scope=control_type, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 178, 6)))

control_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.bounds.abundances'), hilo, scope=control_type, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 180, 6)))

control_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.bounds.CTP'), hilo, scope=control_type, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 182, 6)))

control_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.bounds.HLoc'), hilo, scope=control_type, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 184, 6)))

control_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.bounds.HScale'), hilo, scope=control_type, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 186, 6)))

control_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.bounds.HThick'), hilo, scope=control_type, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 188, 6)))

def _BuildAutomaton_4 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_4
    del _BuildAutomaton_4
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = None
    symbol = pyxb.binding.content.ElementUse(control_type._UseForTag(pyxb.namespace.ExpandedName(None, 'cerberus.crbmodel.fitmolecules')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 170, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(control_type._UseForTag(pyxb.namespace.ExpandedName(None, 'cerberus.crbmodel.HITEMPmolecules')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 172, 6))
    st_1 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_1)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(control_type._UseForTag(pyxb.namespace.ExpandedName(None, 'cerberus.crbmodel.HITRANmolecules')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 174, 6))
    st_2 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_2)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(control_type._UseForTag(pyxb.namespace.ExpandedName(None, 'cerberus.crbmodel.EXOMOLmolecules')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 176, 6))
    st_3 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_3)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(control_type._UseForTag(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.bounds.Teq')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 178, 6))
    st_4 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_4)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(control_type._UseForTag(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.bounds.abundances')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 180, 6))
    st_5 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_5)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(control_type._UseForTag(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.bounds.CTP')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 182, 6))
    st_6 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_6)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(control_type._UseForTag(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.bounds.HLoc')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 184, 6))
    st_7 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_7)
    final_update = None
    symbol = pyxb.binding.content.ElementUse(control_type._UseForTag(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.bounds.HScale')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 186, 6))
    st_8 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_8)
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(control_type._UseForTag(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.bounds.HThick')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 188, 6))
    st_9 = fac.State(symbol, is_initial=False, final_update=final_update, is_unordered_catenation=False)
    states.append(st_9)
    transitions = []
    transitions.append(fac.Transition(st_1, [
         ]))
    st_0._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_2, [
         ]))
    st_1._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_3, [
         ]))
    st_2._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_4, [
         ]))
    st_3._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_5, [
         ]))
    st_4._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_6, [
         ]))
    st_5._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_7, [
         ]))
    st_6._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_8, [
         ]))
    st_7._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_9, [
         ]))
    st_8._set_transitionSet(transitions)
    transitions = []
    st_9._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
control_type._Automaton = _BuildAutomaton_4()




filter_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'exclude'), filter_names, scope=filter_type, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 265, 6)))

filter_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'include'), filter_names, scope=filter_type, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 267, 6)))

def _BuildAutomaton_5 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_5
    del _BuildAutomaton_5
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 265, 6))
    counters.add(cc_0)
    cc_1 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 267, 6))
    counters.add(cc_1)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(filter_type._UseForTag(pyxb.namespace.ExpandedName(None, 'exclude')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 265, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_1, False))
    symbol = pyxb.binding.content.ElementUse(filter_type._UseForTag(pyxb.namespace.ExpandedName(None, 'include')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 267, 6))
    st_1 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_1)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    transitions.append(fac.Transition(st_1, [
        fac.UpdateInstruction(cc_0, False) ]))
    st_0._set_transitionSet(transitions)
    transitions = []
    transitions.append(fac.Transition(st_1, [
        fac.UpdateInstruction(cc_1, True) ]))
    st_1._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=None)
filter_type._Automaton = _BuildAutomaton_5()




lever_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'controls'), control_type, scope=lever_type, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 287, 6)))

lever_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'filters'), filter_type, scope=lever_type, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 288, 6)))

lever_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'pymc'), pymc_type, scope=lever_type, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 289, 6)))

lever_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'run_only'), sequester_type, scope=lever_type, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 290, 6)))

lever_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'sequester'), sequester_type, scope=lever_type, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 291, 6)))

def _BuildAutomaton_7 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_7
    del _BuildAutomaton_7
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(lever_type._UseForTag(pyxb.namespace.ExpandedName(None, 'controls')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 287, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_8 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_8
    del _BuildAutomaton_8
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(lever_type._UseForTag(pyxb.namespace.ExpandedName(None, 'filters')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 288, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_9 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_9
    del _BuildAutomaton_9
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(lever_type._UseForTag(pyxb.namespace.ExpandedName(None, 'pymc')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 289, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_10 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_10
    del _BuildAutomaton_10
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(lever_type._UseForTag(pyxb.namespace.ExpandedName(None, 'run_only')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 290, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_11 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_11
    del _BuildAutomaton_11
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(lever_type._UseForTag(pyxb.namespace.ExpandedName(None, 'sequester')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 291, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_6 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_6
    del _BuildAutomaton_6
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    sub_automata = []
    sub_automata.append(_BuildAutomaton_7())
    sub_automata.append(_BuildAutomaton_8())
    sub_automata.append(_BuildAutomaton_9())
    sub_automata.append(_BuildAutomaton_10())
    sub_automata.append(_BuildAutomaton_11())
    final_update = set()
    symbol = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 286, 4)
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=True)
    st_0._set_subAutomata(*sub_automata)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
lever_type._Automaton = _BuildAutomaton_6()




pymc_count_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'target'), target_override_type, scope=pymc_count_type, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 304, 6)))

def _BuildAutomaton_12 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_12
    del _BuildAutomaton_12
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 304, 6))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(pymc_count_type._UseForTag(pyxb.namespace.ExpandedName(None, 'target')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 304, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=None)
pymc_count_type._Automaton = _BuildAutomaton_12()




pymc_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'cerberuschains'), pymc_count_type, scope=pymc_type, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 317, 6)))

pymc_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'cerberuschainlen'), pymc_count_type, scope=pymc_type, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 318, 6)))

pymc_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'spectrumchains'), pymc_count_type, scope=pymc_type, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 319, 6)))

pymc_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'spectrumchainlen'), pymc_count_type, scope=pymc_type, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 320, 6)))

def _BuildAutomaton_14 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_14
    del _BuildAutomaton_14
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(pymc_type._UseForTag(pyxb.namespace.ExpandedName(None, 'cerberuschains')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 317, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_15 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_15
    del _BuildAutomaton_15
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(pymc_type._UseForTag(pyxb.namespace.ExpandedName(None, 'cerberuschainlen')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 318, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_16 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_16
    del _BuildAutomaton_16
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(pymc_type._UseForTag(pyxb.namespace.ExpandedName(None, 'spectrumchains')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 319, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_17 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_17
    del _BuildAutomaton_17
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(pymc_type._UseForTag(pyxb.namespace.ExpandedName(None, 'spectrumchainlen')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 320, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_13 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_13
    del _BuildAutomaton_13
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    sub_automata = []
    sub_automata.append(_BuildAutomaton_14())
    sub_automata.append(_BuildAutomaton_15())
    sub_automata.append(_BuildAutomaton_16())
    sub_automata.append(_BuildAutomaton_17())
    final_update = set()
    symbol = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 316, 4)
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=True)
    st_0._set_subAutomata(*sub_automata)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
pymc_type._Automaton = _BuildAutomaton_13()




sequester_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'target'), target_type, scope=sequester_type, location=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 331, 6)))

def _BuildAutomaton_18 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_18
    del _BuildAutomaton_18
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 331, 6))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(sequester_type._UseForTag(pyxb.namespace.ExpandedName(None, 'target')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 331, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=None)
sequester_type._Automaton = _BuildAutomaton_18()




def _BuildAutomaton_20 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_20
    del _BuildAutomaton_20
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'controls')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 287, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_21 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_21
    del _BuildAutomaton_21
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'filters')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 288, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_22 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_22
    del _BuildAutomaton_22
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'pymc')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 289, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_23 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_23
    del _BuildAutomaton_23
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'run_only')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 290, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_24 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_24
    del _BuildAutomaton_24
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'sequester')), pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 291, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_19 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_19
    del _BuildAutomaton_19
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    sub_automata = []
    sub_automata.append(_BuildAutomaton_20())
    sub_automata.append(_BuildAutomaton_21())
    sub_automata.append(_BuildAutomaton_22())
    sub_automata.append(_BuildAutomaton_23())
    sub_automata.append(_BuildAutomaton_24())
    final_update = set()
    symbol = pyxb.utils.utility.Location('/home/bryden/gitrepos/esp/excalibur/runtime/levers.xsd', 286, 4)
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=True)
    st_0._set_subAutomata(*sub_automata)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
CTD_ANON._Automaton = _BuildAutomaton_19()

