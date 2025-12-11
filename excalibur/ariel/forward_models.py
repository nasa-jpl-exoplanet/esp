'''ariel forward_models ds'''

# Heritage code shame:
# pylint: disable=too-many-arguments,too-many-locals,too-many-positional-arguments

# import os
# import excalibur
import excalibur.system.core as syscore
from excalibur.cerberus.core import hazelib

# from excalibur.cerberus.forward_model import crbmodel
from excalibur.cerberus.forward_model import crbFM


# ----------------------------------------------------------------------------------------------
def make_cerberus_atmos(
    runtime_params,
    wavelength_um,
    model_params,
    xslib,
    planet_letter,
    chemistry='TEC',
    mixratios=None,
):
    '''
    Create a simulated spectrum using the code that's better than the other ones
    '''
    # print('modelparams', model_params)

    # EQUILIBRIUM TEMPERATURE
    Teq = model_params['Teq']

    # CLOUD/HAZE PARAMETERS
    ctp = model_params['CTP']
    hazescale = model_params['HScale']
    hazeloc = model_params['HLoc']
    hazethick = model_params['HThick']

    # ABUNDANCES
    if mixratios:
        # use the given mixing ratios
        tceqdict = None
    else:
        # equilibrium chemistry
        tceqdict = {}
        tceqdict['XtoH'] = model_params['metallicity']
        tceqdict['CtoO'] = model_params['C/O']
        tceqdict['NtoO'] = 0
        # print('cloudfree forward model input chem =', tceqdict)

    ssc = syscore.ssconstants(mks=True)
    rp0 = model_params['Rp'] * ssc['Rjup']  # MK
    # orbp = {
    #    'R*': model_params['R*'],
    #    planet_letter: {'logg': model_params['logg']},
    # }
    # model_params doesn't have same structure as system.finalize info
    #  there's no planet_letter dictionary, so one mod needed for planet logg
    model_params[planet_letter] = {'logg': model_params['logg']}

    crbhzlib = {'PROFILE': []}
    # hazedir = os.path.join(excalibur.context['data_dir'], 'CERBERUS/HAZE')
    # print('hazedir: ',hazedir)
    hazelib(crbhzlib)
    # print('haze lib',crbhzlib)

    # CERBERUS FORWARD MODEL
    fmc = crbFM().crbmodel(
        float(Teq),
        float(ctp),
        hazescale=float(hazescale),
        hazeloc=float(hazeloc),
        hazethick=float(hazethick),
        hzlib=crbhzlib,
        chemistry=chemistry,
        cheq=tceqdict,
        mixratio=mixratios,
        rp0=rp0,
        xsecs=xslib['data'][planet_letter]['XSECS'],
        qtgrid=xslib['data'][planet_letter]['QTGRID'],
        wgrid=wavelength_um,
        planet=planet_letter,
        orbp=model_params,
        knownspecies=runtime_params.knownspecies,
        cialist=runtime_params.cialist,
        xmollist=runtime_params.xmollist,
        lbroadening=runtime_params.lbroadening,
        lshifting=runtime_params.lshifting,
        nlevels=runtime_params.nlevels,
        Hsmax=runtime_params.Hsmax,
        solrad=runtime_params.solrad,
        break_down_by_molecule=True,
    )

    return (
        fmc.spectrum,
        fmc.breakdown_by_molecule,
        fmc.pressureGrid,
        fmc.opticalDepthProfiles,
    )
