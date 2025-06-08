# ./excalibur/runtime/binding.py
# -*- coding: utf-8 -*-
# PyXB bindings for NM:e92452c8d3e28a9e27abfc9994d2007779e7f4c9
# Generated 2025-02-19 09:41:02.601057 by PyXB version 1.3.3 using Python 3.12.3.final.0
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
_GenerationUID = pyxb.utils.utility.UniqueIdentifier('urn:uuid:b0172147-eee8-11ef-86dd-d493900190e3')

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


# Atomic simple type: filter_names
class filter_names (pyxb.binding.datatypes.normalizedString, pyxb.binding.basis.enumeration_mixin):

    """An atomic simple type."""

    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'filter_names')
    _XSDLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 158, 2)
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

# Complex type control_type with content type EMPTY
class control_type (pyxb.binding.basis.complexTypeDefinition):
    """
        These attributes allow specific control within the AE. The attribute name
        is the task.algorithm.varName being set. The value of true means it will
        be exercised while false means it will not.
      """
    _TypeDefinition = None
    _ContentTypeTag = pyxb.binding.basis.complexTypeDefinition._CT_EMPTY
    _Abstract = False
    _ExpandedName = pyxb.namespace.ExpandedName(Namespace, 'control_type')
    _XSDLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 29, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Attribute cerberus.atmos.fitCloudParameters uses Python identifier cerberus_atmos_fitCloudParameters
    __cerberus_atmos_fitCloudParameters = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.fitCloudParameters'), 'cerberus_atmos_fitCloudParameters', '__AbsentNamespace0_control_type_cerberus_atmos_fitCloudParameters', pyxb.binding.datatypes.boolean, required=True)
    __cerberus_atmos_fitCloudParameters._DeclarationLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 37, 4)
    __cerberus_atmos_fitCloudParameters._UseLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 37, 4)
    
    cerberus_atmos_fitCloudParameters = property(__cerberus_atmos_fitCloudParameters.value, __cerberus_atmos_fitCloudParameters.set, None, None)

    
    # Attribute cerberus.atmos.fitNtoO uses Python identifier cerberus_atmos_fitNtoO
    __cerberus_atmos_fitNtoO = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.fitNtoO'), 'cerberus_atmos_fitNtoO', '__AbsentNamespace0_control_type_cerberus_atmos_fitNtoO', pyxb.binding.datatypes.boolean, required=True)
    __cerberus_atmos_fitNtoO._DeclarationLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 39, 4)
    __cerberus_atmos_fitNtoO._UseLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 39, 4)
    
    cerberus_atmos_fitNtoO = property(__cerberus_atmos_fitNtoO.value, __cerberus_atmos_fitNtoO.set, None, None)

    
    # Attribute cerberus.atmos.fitCtoO uses Python identifier cerberus_atmos_fitCtoO
    __cerberus_atmos_fitCtoO = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.fitCtoO'), 'cerberus_atmos_fitCtoO', '__AbsentNamespace0_control_type_cerberus_atmos_fitCtoO', pyxb.binding.datatypes.boolean, required=True)
    __cerberus_atmos_fitCtoO._DeclarationLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 41, 4)
    __cerberus_atmos_fitCtoO._UseLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 41, 4)
    
    cerberus_atmos_fitCtoO = property(__cerberus_atmos_fitCtoO.value, __cerberus_atmos_fitCtoO.set, None, None)

    
    # Attribute cerberus.atmos.fitT uses Python identifier cerberus_atmos_fitT
    __cerberus_atmos_fitT = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.fitT'), 'cerberus_atmos_fitT', '__AbsentNamespace0_control_type_cerberus_atmos_fitT', pyxb.binding.datatypes.boolean, required=True)
    __cerberus_atmos_fitT._DeclarationLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 43, 4)
    __cerberus_atmos_fitT._UseLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 43, 4)
    
    cerberus_atmos_fitT = property(__cerberus_atmos_fitT.value, __cerberus_atmos_fitT.set, None, None)


    # Attribute target.autofill.selectMostRecent uses Python identifier target_autofill_selectMostRecent
    __target_autofill_selectMostRecent = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'target.autofill.selectMostRecent'), 'target_autofill_selectMostRecent', '__AbsentNamespace0_control_type_target_autofill_selectMostRecent', pyxb.binding.datatypes.boolean, required=True)
    __target_autofill_selectMostRecent._DeclarationLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 45, 4)
    __target_autofill_selectMostRecent._UseLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 45, 4)
    
    target_autofill_selectMostRecent = property(__target_autofill_selectMostRecent.value, __target_autofill_selectMostRecent.set, None, None)

    
    # Attribute ariel.simspectra.includeMetallicityDispersion uses Python identifier ariel_simspectrum_includeMetallicityDispersion
    __ariel_simspectrum_includeMetallicityDispersion = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'ariel.simspectrum.includeMetallicityDispersion'), 'ariel_simspectrum_includeMetallicityDispersion', '__AbsentNamespace0_control_type_ariel_simspectrum_includeMetallicityDispersion', pyxb.binding.datatypes.boolean, required=True)
    __ariel_simspectrum_includeMetallicityDispersion._DeclarationLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 47, 4)
    __ariel_simspectrum_includeMetallicityDispersion._UseLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 47, 4)
    
    ariel_simspectrum_includeMetallicityDispersion = property(__ariel_simspectrum_includeMetallicityDispersion.value, __ariel_simspectrum_includeMetallicityDispersion.set, None, None)

    # Attribute cerberus.atmos.sliceSampler uses Python identifier cerberus_atmos_sliceSampler
    __cerberus_atmos_sliceSampler = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.sliceSampler'), 'cerberus_atmos_sliceSampler', '__AbsentNamespace0_control_type_cerberus_atmos_sliceSampler', pyxb.binding.datatypes.boolean, required=True)
    __cerberus_atmos_sliceSampler._DeclarationLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 49, 4)
    __cerberus_atmos_sliceSampler._UseLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 49, 4)

    cerberus_atmos_sliceSampler = property(__cerberus_atmos_sliceSampler.value, __cerberus_atmos_sliceSampler.set, None, None)

    # Attribute target.autofill.maximizeSelfConsistency uses Python identifier target_autofill_maximizeSelfConsistency
    __target_autofill_maximizeSelfConsistency = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'target.autofill.maximizeSelfConsistency'), 'target_autofill_maximizeSelfConsistency', '__AbsentNamespace0_control_type_target_autofill_maximizeSelfConsistency', pyxb.binding.datatypes.boolean, required=True)
    __target_autofill_maximizeSelfConsistency._DeclarationLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 51, 4)
    __target_autofill_maximizeSelfConsistency._UseLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 51, 4)

    target_autofill_maximizeSelfConsistency = property(__target_autofill_maximizeSelfConsistency.value, __target_autofill_maximizeSelfConsistency.set, None, None)

    # Attribute ariel.simspectrum.randomCloudProperties uses Python identifier ariel_simspectrum_randomCloudProperties
    __ariel_simspectrum_randomCloudProperties = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'ariel.simspectrum.randomCloudProperties'), 'ariel_simspectrum_randomCloudProperties', '__AbsentNamespace0_control_type_ariel_simspectrum_randomCloudProperties', pyxb.binding.datatypes.boolean, required=True)
    __ariel_simspectrum_randomCloudProperties._DeclarationLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 53, 4)
    __ariel_simspectrum_randomCloudProperties._UseLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 53, 4)

    ariel_simspectrum_randomCloudProperties = property(__ariel_simspectrum_randomCloudProperties.value, __ariel_simspectrum_randomCloudProperties.set, None, None)

    # Attribute ariel.simspectrum.thorngrenMassMetals uses Python identifier ariel_simspectrum_thorngrenMassMetals
    __ariel_simspectrum_thorngrenMassMetals = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'ariel.simspectrum.thorngrenMassMetals'), 'ariel_simspectrum_thorngrenMassMetals', '__AbsentNamespace0_control_type_ariel_simspectrum_thorngrenMassMetals', pyxb.binding.datatypes.boolean, required=True)
    __ariel_simspectrum_thorngrenMassMetals._DeclarationLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 55, 4)
    __ariel_simspectrum_thorngrenMassMetals._UseLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 55, 4)

    ariel_simspectrum_thorngrenMassMetals = property(__ariel_simspectrum_thorngrenMassMetals.value, __ariel_simspectrum_thorngrenMassMetals.set, None, None)

        # Attribute cerberus.atmos.crbmodel.lbroadening uses Python identifier cerberus_atmos_crbmodel_lbroadening
    __cerberus_atmos_crbmodel_lbroadening = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.crbmodel.lbroadening'), 'cerberus_atmos_crbmodel_lbroadening', '__AbsentNamespace0_control_type_cerberus_atmos_crbmodel_lbroadening', pyxb.binding.datatypes.boolean, required=True)
    __cerberus_atmos_crbmodel_lbroadening._DeclarationLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 57, 4)
    __cerberus_atmos_crbmodel_lbroadening._UseLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 57, 4)

    cerberus_atmos_crbmodel_lbroadening = property(__cerberus_atmos_crbmodel_lbroadening.value, __cerberus_atmos_crbmodel_lbroadening.set, None, None)

    # Attribute cerberus.atmos.crbmodel.lshifting uses Python identifier cerberus_atmos_crbmodel_lshifting
    __cerberus_atmos_crbmodel_lshifting = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.crbmodel.lshifting'), 'cerberus_atmos_crbmodel_lshifting', '__AbsentNamespace0_control_type_cerberus_atmos_crbmodel_lshifting', pyxb.binding.datatypes.boolean, required=True)
    __cerberus_atmos_crbmodel_lshifting._DeclarationLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 59, 4)
    __cerberus_atmos_crbmodel_lshifting._UseLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 59, 4)

    cerberus_atmos_crbmodel_lshifting = property(__cerberus_atmos_crbmodel_lshifting.value, __cerberus_atmos_crbmodel_lshifting.set, None, None)

    # Attribute cerberus.atmos.crbmodel.isothermal uses Python identifier cerberus_atmos_crbmodel_isothermal
    __cerberus_atmos_crbmodel_isothermal = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'cerberus.atmos.crbmodel.isothermal'), 'cerberus_atmos_crbmodel_isothermal', '__AbsentNamespace0_control_type_cerberus_atmos_crbmodel_isothermal', pyxb.binding.datatypes.boolean, required=True)
    __cerberus_atmos_crbmodel_isothermal._DeclarationLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 61, 4)
    __cerberus_atmos_crbmodel_isothermal._UseLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 61, 4)

    cerberus_atmos_crbmodel_isothermal = property(__cerberus_atmos_crbmodel_isothermal.value, __cerberus_atmos_crbmodel_isothermal.set, None, None)

    # Attribute ariel.simspectrum.tier uses Python identifier ariel_simspectrum_tier
    __ariel_simspectrum_tier = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'ariel.simspectrum.tier'), 'ariel_simspectrum_tier', '__AbsentNamespace0_control_type_ariel_simspectrum_tier', pyxb.binding.datatypes.integer, required=True)
    __ariel_simspectrum_tier._DeclarationLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 63, 4)
    __ariel_simspectrum_tier._UseLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 63, 4)

    ariel_simspectrum_tier = property(__ariel_simspectrum_tier.value, __ariel_simspectrum_tier.set, None, None)

    _ElementMap.update({
        
    })
    _AttributeMap.update({
        __cerberus_atmos_fitCloudParameters.name() : __cerberus_atmos_fitCloudParameters,
        __cerberus_atmos_fitNtoO.name() : __cerberus_atmos_fitNtoO,
        __cerberus_atmos_fitCtoO.name() : __cerberus_atmos_fitCtoO,
        __cerberus_atmos_fitT.name() : __cerberus_atmos_fitT,
        __cerberus_atmos_sliceSampler.name() : __cerberus_atmos_sliceSampler,
        __cerberus_atmos_crbmodel_lbroadening.name() : __cerberus_atmos_crbmodel_lbroadening,
        __cerberus_atmos_crbmodel_lshifting.name() : __cerberus_atmos_crbmodel_lshifting,
        __cerberus_atmos_crbmodel_isothermal.name() : __cerberus_atmos_crbmodel_isothermal,
        __target_autofill_selectMostRecent.name() : __target_autofill_selectMostRecent,
        __target_autofill_maximizeSelfConsistency.name() : __target_autofill_maximizeSelfConsistency,
        __ariel_simspectrum_includeMetallicityDispersion.name() : __ariel_simspectrum_includeMetallicityDispersion,
        __ariel_simspectrum_randomCloudProperties.name() : __ariel_simspectrum_randomCloudProperties,
        __ariel_simspectrum_thorngrenMassMetals.name() : __ariel_simspectrum_thorngrenMassMetals,
        __ariel_simspectrum_tier.name() : __ariel_simspectrum_tier
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
    _XSDLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 51, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType

    # Element exclude uses Python identifier exclude
    __exclude = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'exclude'), 'exclude', '__AbsentNamespace0_filter_type_exclude', True, pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 65, 6), )
    
    exclude = property(__exclude.value, __exclude.set, None, None)

    
    # Element include uses Python identifier include
    __include = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'include'), 'include', '__AbsentNamespace0_filter_type_include', True, pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 67, 6), )

    
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
    _XSDLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 72, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element controls uses Python identifier controls
    __controls = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'controls'), 'controls', '__AbsentNamespace0_lever_type_controls', False, pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 87, 6), )

    
    controls = property(__controls.value, __controls.set, None, None)

    
    # Element filters uses Python identifier filters
    __filters = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'filters'), 'filters', '__AbsentNamespace0_lever_type_filters', False, pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 88, 6), )

    
    filters = property(__filters.value, __filters.set, None, None)

    
    # Element pymc uses Python identifier pymc
    __pymc = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'pymc'), 'pymc', '__AbsentNamespace0_lever_type_pymc', False, pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 89, 6), )

    
    pymc = property(__pymc.value, __pymc.set, None, None)

    
    # Element run_only uses Python identifier run_only
    __run_only = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'run_only'), 'run_only', '__AbsentNamespace0_lever_type_run_only', False, pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 90, 6), )

    
    run_only = property(__run_only.value, __run_only.set, None, None)

    
    # Element sequester uses Python identifier sequester
    __sequester = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'sequester'), 'sequester', '__AbsentNamespace0_lever_type_sequester', False, pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 91, 6), )

    
    sequester = property(__sequester.value, __sequester.set, None, None)

    
    # Attribute index uses Python identifier index
    __index = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'index'), 'index', '__AbsentNamespace0_lever_type_index', pyxb.binding.datatypes.normalizedString, unicode_default='registry')
    __index._DeclarationLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 93, 4)
    __index._UseLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 93, 4)
    
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
    _XSDLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 96, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element target uses Python identifier target
    __target = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'target'), 'target', '__AbsentNamespace0_pymc_count_type_target', True, pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 104, 6), )

    
    target = property(__target.value, __target.set, None, None)

    
    # Attribute default uses Python identifier default
    __default = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'default'), 'default', '__AbsentNamespace0_pymc_count_type_default', pyxb.binding.datatypes.positiveInteger, required=True)
    __default._DeclarationLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 107, 4)
    __default._UseLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 107, 4)
    
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
    _XSDLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 110, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element cerberuschainlen uses Python identifier cerberuschainlen
    __cerberuschainlen = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'cerberuschainlen'), 'cerberuschainlen', '__AbsentNamespace0_pymc_type_cerberuschainlen', False, pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 117, 6), )

    cerberuschainlen = property(__cerberuschainlen.value, __cerberuschainlen.set, None, None)

    # Element spectrumchainlen uses Python identifier spectrumchainlen
    __spectrumchainlen = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'spectrumchainlen'), 'spectrumchainlen', '__AbsentNamespace0_pymc_type_spectrumchainlen', False, pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 118, 6), )
    
    spectrumchainlen = property(__spectrumchainlen.value, __spectrumchainlen.set, None, None)

    # Element cerberuschains uses Python identifier cerberuschains
    __cerberuschains = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'cerberuschains'), 'cerberuschains', '__AbsentNamespace0_pymc_type_cerberuschains', False, pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 119, 6), )

    cerberuschains = property(__cerberuschains.value, __cerberuschains.set, None, None)

    # Element spectrumchains uses Python identifier spectrumchains
    __spectrumchains = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'spectrumchains'), 'spectrumchains', '__AbsentNamespace0_pymc_type_spectrumchains', False, pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 120, 6), )

    spectrumchains = property(__spectrumchains.value, __spectrumchains.set, None, None)

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
    _XSDLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 122, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Element target uses Python identifier target
    __target = pyxb.binding.content.ElementDeclaration(pyxb.namespace.ExpandedName(None, 'target'), 'target', '__AbsentNamespace0_sequester_type_target', True, pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 129, 6), )

    
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
    _XSDLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 134, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.normalizedString
    
    # Attribute because uses Python identifier because
    __because = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'because'), 'because', '__AbsentNamespace0_target_type_because', pyxb.binding.datatypes.normalizedString, required=True)
    __because._DeclarationLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 142, 8)
    __because._UseLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 142, 8)
    
    because = property(__because.value, __because.set, None, None)

    
    # Attribute isRegex uses Python identifier isRegex
    __isRegex = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'isRegex'), 'isRegex', '__AbsentNamespace0_target_type_isRegex', pyxb.binding.datatypes.boolean, unicode_default='false')
    __isRegex._DeclarationLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 143, 8)
    __isRegex._UseLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 143, 8)
    
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
    _XSDLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 148, 2)
    _ElementMap = {}
    _AttributeMap = {}
    # Base type is pyxb.binding.datatypes.anyType
    
    # Attribute name uses Python identifier name
    __name = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'name'), 'name', '__AbsentNamespace0_target_override_type_name', pyxb.binding.datatypes.normalizedString, required=True)
    __name._DeclarationLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 154, 4)
    __name._UseLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 154, 4)
    
    name = property(__name.value, __name.set, None, None)

    
    # Attribute steps uses Python identifier steps
    __steps = pyxb.binding.content.AttributeUse(pyxb.namespace.ExpandedName(None, 'steps'), 'steps', '__AbsentNamespace0_target_override_type_steps', pyxb.binding.datatypes.positiveInteger, required=True)
    __steps._DeclarationLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 155, 4)
    __steps._UseLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 155, 4)
    
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
    _XSDLocation = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 177, 4)
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


levers = pyxb.binding.basis.element(pyxb.namespace.ExpandedName(Namespace, 'levers'), CTD_ANON, location=pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 176, 2))
Namespace.addCategoryObject('elementBinding', levers.name().localName(), levers)



filter_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'exclude'), filter_names, scope=filter_type, location=pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 65, 6)))

filter_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'include'), filter_names, scope=filter_type, location=pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 67, 6)))

def _BuildAutomaton ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton
    del _BuildAutomaton
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 65, 6))
    counters.add(cc_0)
    cc_1 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 67, 6))
    counters.add(cc_1)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(filter_type._UseForTag(pyxb.namespace.ExpandedName(None, 'exclude')), pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 65, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_1, False))
    symbol = pyxb.binding.content.ElementUse(filter_type._UseForTag(pyxb.namespace.ExpandedName(None, 'include')), pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 67, 6))
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
filter_type._Automaton = _BuildAutomaton()




lever_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'controls'), control_type, scope=lever_type, location=pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 87, 6)))

lever_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'filters'), filter_type, scope=lever_type, location=pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 88, 6)))

lever_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'pymc'), pymc_type, scope=lever_type, location=pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 89, 6)))

lever_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'run_only'), sequester_type, scope=lever_type, location=pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 90, 6)))

lever_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'sequester'), sequester_type, scope=lever_type, location=pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 91, 6)))

def _BuildAutomaton_2 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_2
    del _BuildAutomaton_2
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(lever_type._UseForTag(pyxb.namespace.ExpandedName(None, 'controls')), pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 87, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_3 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_3
    del _BuildAutomaton_3
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(lever_type._UseForTag(pyxb.namespace.ExpandedName(None, 'filters')), pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 88, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_4 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_4
    del _BuildAutomaton_4
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(lever_type._UseForTag(pyxb.namespace.ExpandedName(None, 'pymc')), pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 89, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_5 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_5
    del _BuildAutomaton_5
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(lever_type._UseForTag(pyxb.namespace.ExpandedName(None, 'run_only')), pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 90, 6))
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
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(lever_type._UseForTag(pyxb.namespace.ExpandedName(None, 'sequester')), pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 91, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_ ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_
    del _BuildAutomaton_
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    sub_automata = []
    sub_automata.append(_BuildAutomaton_2())
    sub_automata.append(_BuildAutomaton_3())
    sub_automata.append(_BuildAutomaton_4())
    sub_automata.append(_BuildAutomaton_5())
    sub_automata.append(_BuildAutomaton_6())
    final_update = set()
    symbol = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 86, 4)
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=True)
    st_0._set_subAutomata(*sub_automata)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
lever_type._Automaton = _BuildAutomaton_()




pymc_count_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'target'), target_override_type, scope=pymc_count_type, location=pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 104, 6)))

def _BuildAutomaton_7 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_7
    del _BuildAutomaton_7
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 104, 6))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(pymc_count_type._UseForTag(pyxb.namespace.ExpandedName(None, 'target')), pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 104, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=None)
pymc_count_type._Automaton = _BuildAutomaton_7()




pymc_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'cerberuschainlen'), pymc_count_type, scope=pymc_type, location=pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 117, 6)))

pymc_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'spectrumchainlen'), pymc_count_type, scope=pymc_type, location=pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 118, 6)))

pymc_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'cerberuschains'), pymc_count_type, scope=pymc_type, location=pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 119, 6)))

pymc_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'spectrumchains'), pymc_count_type, scope=pymc_type, location=pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 120, 6)))

def _BuildAutomaton_11 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_11
    del _BuildAutomaton_11
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(pymc_type._UseForTag(pyxb.namespace.ExpandedName(None, 'cerberuschains')), pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 117, 6))
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
    symbol = pyxb.binding.content.ElementUse(pymc_type._UseForTag(pyxb.namespace.ExpandedName(None, 'cerberuschainlen')), pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 117, 6))
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
    symbol = pyxb.binding.content.ElementUse(pymc_type._UseForTag(pyxb.namespace.ExpandedName(None, 'spectrumchainlen')), pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 118, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_12 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_12
    del _BuildAutomaton_12
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(pymc_type._UseForTag(pyxb.namespace.ExpandedName(None, 'spectrumchains')), pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 118, 6))
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
    sub_automata = []
    sub_automata.append(_BuildAutomaton_9())
    sub_automata.append(_BuildAutomaton_10())
    sub_automata.append(_BuildAutomaton_11())
    sub_automata.append(_BuildAutomaton_12())
    final_update = set()
    symbol = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 116, 4)
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=True)
    st_0._set_subAutomata(*sub_automata)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
pymc_type._Automaton = _BuildAutomaton_8()




sequester_type._AddElement(pyxb.binding.basis.element(pyxb.namespace.ExpandedName(None, 'target'), target_type, scope=sequester_type, location=pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 129, 6)))

def _BuildAutomaton_11 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_11
    del _BuildAutomaton_11
    import pyxb.utils.fac as fac

    counters = set()
    cc_0 = fac.CounterCondition(min=0, max=None, metadata=pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 129, 6))
    counters.add(cc_0)
    states = []
    final_update = set()
    final_update.add(fac.UpdateInstruction(cc_0, False))
    symbol = pyxb.binding.content.ElementUse(sequester_type._UseForTag(pyxb.namespace.ExpandedName(None, 'target')), pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 129, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    transitions.append(fac.Transition(st_0, [
        fac.UpdateInstruction(cc_0, True) ]))
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, True, containing_state=None)
sequester_type._Automaton = _BuildAutomaton_11()




def _BuildAutomaton_13 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_13
    del _BuildAutomaton_13
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'controls')), pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 87, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_14 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_14
    del _BuildAutomaton_14
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    final_update = set()
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'filters')), pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 88, 6))
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
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'pymc')), pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 89, 6))
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
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'run_only')), pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 90, 6))
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
    symbol = pyxb.binding.content.ElementUse(CTD_ANON._UseForTag(pyxb.namespace.ExpandedName(None, 'sequester')), pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 91, 6))
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=False)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=st_0)

def _BuildAutomaton_12 ():
    # Remove this helper function from the namespace after it is invoked
    global _BuildAutomaton_12
    del _BuildAutomaton_12
    import pyxb.utils.fac as fac

    counters = set()
    states = []
    sub_automata = []
    sub_automata.append(_BuildAutomaton_13())
    sub_automata.append(_BuildAutomaton_14())
    sub_automata.append(_BuildAutomaton_15())
    sub_automata.append(_BuildAutomaton_16())
    sub_automata.append(_BuildAutomaton_17())
    final_update = set()
    symbol = pyxb.utils.utility.Location('/home/niessner/Projects/esp/excalibur/runtime/levers.xsd', 86, 4)
    st_0 = fac.State(symbol, is_initial=True, final_update=final_update, is_unordered_catenation=True)
    st_0._set_subAutomata(*sub_automata)
    states.append(st_0)
    transitions = []
    st_0._set_transitionSet(transitions)
    return fac.Automaton(states, counters, False, containing_state=None)
CTD_ANON._Automaton = _BuildAutomaton_12()

