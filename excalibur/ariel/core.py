'''ariel core ds'''
# -- IMPORTS -- ------------------------------------------------------
import logging; log = logging.getLogger(__name__)

import excalibur
import excalibur.system.core as syscore
import excalibur.util.cerberus as crbutil

from excalibur.ariel.metallicity import \
    massMetalRelation, massMetalRelationDisp, randomCtoO
from excalibur.ariel.arielInstrumentModel import load_ariel_instrument
from excalibur.ariel.arielObservingPlan import make_tier_table
from excalibur.ariel.forwardModels import makeTaurexAtmos, makeCerberusAtmos
from excalibur.cerberus.core import myxsecs

import os
import io
import sys
# import pickle
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as cst

# ---------------------------- ---------------------------------------
# -- SV VALIDITY -- --------------------------------------------------
def checksv(sv):
    '''
    Tests for empty SV shell
    '''
    if sv['STATUS'][-1]:
        valid = True
        errstring = None
    else:
        valid = False
        errstring = sv.name()+' IS EMPTY'
    return valid, errstring
# ----------------- --------------------------------------------------
# -- SIMULATE ARIEL SPECTRA ------------------------------------------
def simulate_spectra(target, system_dict, out):
    '''
    Simulate Ariel spectra, adding noise based on the Ariel instrument model
    Mulitple spectra are now calculated, allowing a choice within cerberus.atmos fitting
    1) both Taurex and Cerberus atmosphere models
    2) with or without clouds
    3) two models for metallicity/mmw (mmw=2.3 or FINESSE mass-metallicity relation)
    4) TEC vs DISEQ models [NOT IMPLEMENTED YET!]
    '''

    # ** two key parameters for adjusting the instrument noise model **
    # no longer necessary; Ariel-rad already does these two steps
    # noise_factor = 1./0.7
    # noise_floor_ppm = 50.

    observing_plan = make_tier_table()

    system_params = system_dict['priors']

    # maybe save/read the xslib file from disk.
    # for debugging at least, since it's a lot faster that way
    # xslibSaveDir = os.path.join(excalibur.context['data_dir'], 'bryden/')

    # specify which models should be calculated (use these as keys within data)
    atmosModels = ['cerberus']
    # atmosModels = ['cerberus', 'cerberusNoclouds']
    #  these 8 models were calculated in 2023 runs; less models used in 2024
    # atmosModels = ['cerberus', 'cerberusNoclouds',
    #               'cerberuslowmmw', 'cerberuslowmmwNoclouds',
    #               'taurex', 'taurexNoclouds',
    #               'taurexlowmmw', 'taurexlowmmwNoclouds']
    # atmosModels = ['taurexNoclouds']
    # atmosModels = ['taurex', 'taurexNoclouds',
    #               'taurexlowmmw', 'taurexlowmmwNoclouds']
    # atmosModels = ['cerberusNoclouds']
    # atmosModels = ['cerberusNoclouds','taurexNoclouds']
    # atmosModels = ['cerberusNoclouds','cerberus',
    #               'cerberuslowmmw', 'cerberuslowmmwNoclouds']
    # atmosModels = ['cerberusNoclouds', 'taurexNoclouds', 'cerberus']
    # atmosModels = ['cerberuslowmmw']
    # atmosModels = ['cerberuslowmmwNoclouds']
    # atmosModels = ['cerberusNoclouds',
    #               'cerberusNocloudsR1000',
    #               'cerberusNocloudsR3000',
    #               'cerberusNocloudsR6000',
    #               'cerberusNocloudsR10000']
    #               'cerberusNocloudsR10000',
    #               'cerberusNocloudsR15000']
    out['data']['models'] = atmosModels
    # save target,planet names, for plotting (in states.py)
    out['data']['target'] = target
    out['data']['planets'] = system_params['planets']

    # save the cross-section library from one model to the next (no need to re-calculate)
    # but be careful, it has to be defined for each planet letter
    xslib = {'data':{}, 'STATUS':[False]}

    completed_at_least_one_planet = False
    for planetLetter in system_params['planets']:
        # set the random seed as a function of target name
        #  (otherwise all spectra have identical noise realizations!)
        #  (also they would have the same C/O and metallicity offset!)
        intFromTarget = 1
        for char in target+' '+planetLetter:
            intFromTarget = (123 * intFromTarget + ord(char)) % 1000000
        np.random.seed(intFromTarget)

        # load in the wavelength bins and the noise model
        # there is a separate SNR file for each planet
        ariel_instrument = load_ariel_instrument(target+' '+planetLetter)

        if ariel_instrument:

            # asdf : LATER : add in uncertainty scatter to these model parameters
            model_params = {
                'T*':system_params['T*'],
                'R*':system_params['R*'],
                'Rp':system_params[planetLetter]['rp'],
                'Mp':system_params[planetLetter]['mass'],
                'logg':system_params[planetLetter]['logg'],
                'sma':system_params[planetLetter]['sma'],
                'Teq':system_params[planetLetter]['teq']}

            # Calculate the atmosphere scale height
            #  cerberus wants it, to normalize the spectrum
            #  (some of this code is copied from transit/core)
            sscmks = syscore.ssconstants(mks=True)
            #  NOTE : sma is system, but should be model value.  EDIT LATER?
            eqtemp = model_params['T*']*np.sqrt(model_params['R*']*sscmks['Rsun/AU']/
                                                (2.*model_params['sma']))
            pgrid = np.exp(np.arange(np.log(10.)-15., np.log(10.)+15./100, 15./99))
            pressure = pgrid[::-1]
            # Assume solar metallicity here but then below use each model's metallicity
            mixratio, fH2, fHe = crbutil.crbce(pressure, eqtemp)
            # X2Hr=cheq['XtoH'])
            # assume solar C/O and N/O for now
            # C2Or=cheq['CtoO'], N2Or=cheq['NtoO'])
            mmwsolar, fH2, fHe = crbutil.getmmw(mixratio, protosolar=False, fH2=fH2, fHe=fHe)
            mmw = mmwsolar * cst.m_p  # [kg]
            Hs = cst.Boltzmann * eqtemp / (mmw*1e-2*(10.**float(model_params['logg'])))  # [m]
            HoverRmax = Hs / (model_params['Rp']*sscmks['Rjup'])
            Hssolar = Hs / (model_params['R*']*sscmks['Rsun'])  # this is used for plot scaling

            # skip non-converging atmospheres!!
            # print()
            # print(target,planetLetter,'  scale height / planet radius (solar metallicity) =',HoverRmax)
            # print()
            # this limit corresponds to the maximum scale height allowing for ~1e11 pressure range
            # (our actual range is typically 10 bar to 1 microbar, requiring 0.06 max H/R)
            if HoverRmax > 0.04:
                # print('SKIP:',target,planetLetter,'  scale height / planet radius =',HoverRmax)
                log.warning('--< SKIP UNBOUND ATMOS: %s %s ; scale height / planet radius = %s >--',target,planetLetter,HoverRmax)
            else:

                # planet metallicity is from an assumed mass-metallicity relation
                #  with scatter included
                # planet metallicity should be defined relative to the stellar metallicity
                metallicity_star_dex = system_params['FEH*']
                M_p = model_params['Mp']
                includeMetallicityDispersion = True
                if includeMetallicityDispersion:
                    metallicity_planet_dex = massMetalRelationDisp(metallicity_star_dex, M_p,
                                                                   thorngren=True)
                else:
                    metallicity_planet_dex = massMetalRelation(metallicity_star_dex, M_p,
                                                               thorngren=True)
                CtoO_planet = randomCtoO()  # this is linear C/O, not dex

                # Load the instrument model and rescale based on #-of-transits
                uncertainties = ariel_instrument['noise']
                if target+' '+planetLetter in observing_plan:
                    visits = observing_plan[target+' '+planetLetter]['number of visits']
                    tier = observing_plan[target+' '+planetLetter]['tier']
                else:
                    # default to a single transit observation, if it's not in the Ariel list
                    # print(target+' '+planetLetter,'not found in the Ariel observing plan')
                    errstr = target+' '+planetLetter+' not found in observing plan'
                    log.warning('--< ARIEL SIM_SPECTRUM: %s >--', errstr)
                    visits = '1'
                    tier = '1'
                # print('# of visits:',visits,'  tier',tier,'  ',target+' '+planetLetter)
                uncertainties /= np.sqrt(float(visits))

                # ________LOOP OVER ALL SELECTED MODELS_______
                for atmosModel in atmosModels:
                    # print()
                    # print('starting Atmospheric Model:',atmosModel)

                    # ABUNDANCES
                    if 'lowmmw' in atmosModel:
                        # print(' - using a low mmw')
                        model_params['metallicity'] = 0.  # dex
                        model_params['C/O'] = 0.     # [C/O] (relative to solar)
                    else:
                        model_params['metallicity*'] = metallicity_star_dex
                        model_params['metallicity'] = metallicity_star_dex + metallicity_planet_dex

                        # planet C/O ratio is assumed to be solar
                        #  (0.54951 is the default in ACEChemistry, so it actually has no effect)
                        # actually, let's consider a distribution of C/O, as done for FINESSE
                        model_params['C/O'] = np.log10(CtoO_planet/0.54951)

                    if 'cerberus' in atmosModel:
                        # CLOUD PARAMETERS
                        if 'Noclouds' in atmosModel:
                            model_params['CTP'] = 3.       # cloud deck is very deep - 1000 bars
                            model_params['HScale'] = -10.  # small number means essentially no haze
                            model_params['HLoc'] = 0.
                            model_params['HThick'] = 0.
                        else:
                            # CLOUD TOP PRESSURE
                            # model_params['CTP'] = -2.  # 0.01 bar = 10 mbar
                            # HAZE SCAT CROSS SECTION SCALE FACTOR
                            # model_params['HScale'] = 0.  # some haze (1x times the nominal Jupiter model)
                            # model_params['HLoc'] = 0.
                            # model_params['HThick'] = 0.

                            # median values from Estrela et al 2022, Table 2 (TEC column)
                            model_params['CTP'] = -1.52     # 0.03 bar
                            model_params['HScale'] = -2.10  # 0.008
                            model_params['HLoc'] = -2.30    # 0.005 bar location (meaningless)
                            model_params['HThick'] = 9.76   # very wide vertical thickness

                        if 'R' in atmosModel:
                            iwave1 = atmosModel.index('R')
                            # print('atmosModel',atmosModel)
                            Nhires = atmosModel[iwave1+1:]
                            # Nhires = int(atmosModel[iwave1+1:])
                            # print('Nhires',Nhires)
                            wavelength_um = np.logspace(-0.5,1.0,int(Nhires))
                        else:
                            wavelength_um = ariel_instrument['wavelength']

                        if not xslib['STATUS'][-1]:
                            tempspc = {'data':{planetLetter:{'WB':wavelength_um}}}
                            # print('CALCulating cross-sections START')
                            _ = myxsecs(tempspc, xslib)
                            # print('CALCulating cross-sections DONE')
                        else:
                            # make sure that it exists for this planet letter
                            if planetLetter in xslib['data']:
                                # print('XSLIB ALREADY CALCULATED',planetLetter)
                                pass
                            else:
                                # print('XSLIB TRANSFERRED TO NEW PLANETLETTER')
                                existingPlanetLetter = list(xslib['data'].keys())[-1]
                                # print('existingPlanetLetter',existingPlanetLetter)
                                xslib['data'][planetLetter] = xslib['data'][existingPlanetLetter]
                        # print('xslib check',xslib.keys())
                        # print('xslib status',xslib['STATUS'])
                        # print('xslib data keys:',xslib['data'].keys())
                        # print('xslib data keys:',list(xslib['data'].keys())[-1])

                        cerbModel, cerbModel_by_molecule = makeCerberusAtmos(
                            wavelength_um, model_params, xslib, planetLetter, Hsmax=20)

                        fluxDepth = cerbModel
                        fluxDepth_by_molecule = cerbModel_by_molecule

                    elif 'taurex' in atmosModel:
                        if 'Noclouds' in atmosModel:
                            taurexModel = makeTaurexAtmos(model_params, clouds=False)
                        else:
                            taurexModel = makeTaurexAtmos(model_params)

                        wn,fluxDepth = taurexModel[:2]
                        fluxDepth_by_molecule = {}
                        wavelength_um = 1e4 / wn

                    else:
                        sys.exit('ERROR: unknown model')

                    bads = np.where(np.isnan(fluxDepth))
                    if len(bads[0]) > 0:
                        print('   bad fluxDepths:',len(bads[0]),'out of',len(fluxDepth))

                    # REBIN to the Ariel spectral resolution
                    wavelength_um_rebin = []
                    fluxDepth_rebin = []
                    fluxDepth_by_molecule_rebin = {}
                    molecules = fluxDepth_by_molecule.keys()
                    for molecule in molecules:
                        fluxDepth_by_molecule_rebin[molecule] = []
                    for wavelow,wavehigh in zip(ariel_instrument['wavelow'],
                                                ariel_instrument['wavehigh']):
                        thisWaveBin = np.where((wavelength_um >= wavelow) &
                                               (wavelength_um <= wavehigh))
                        fluxDepth_rebin.append(np.average(fluxDepth[thisWaveBin]))
                        wavelength_um_rebin.append(np.average(wavelength_um[thisWaveBin]))
                        for molecule in molecules:
                            # print('molecule:',molecule)
                            fluxDepth_by_molecule_rebin[molecule].append(np.average(
                                fluxDepth_by_molecule[molecule][thisWaveBin]))

                    wavelength_um_rebin = np.array(wavelength_um_rebin)
                    fluxDepth_rebin = np.array(fluxDepth_rebin)
                    for molecule in molecules:
                        fluxDepth_by_molecule_rebin[molecule] = np.array(
                            fluxDepth_by_molecule_rebin[molecule])
                    # print ('Nan check on raw',np.isnan(np.max(fluxDepth)))
                    # print ('Nan check on rebin',np.isnan(np.max(fluxDepth_rebin)))

                    # ADD OBSERVATIONAL NOISE TO THE TRUE SPECTRUM
                    fluxDepth_observed = fluxDepth_rebin + np.random.normal(scale=uncertainties)

                    # SAVE THE RESULTS
                    if planetLetter not in out['data'].keys():
                        out['data'][planetLetter] = {}

                    # careful - ES and ESerr are supposed to be radii, not transit depth
                    #  need to take a sqrt of them
                    # careful2 - watch out for sqrt of negative numbers
                    signedSqrt = np.sign(fluxDepth_observed) * np.sqrt(np.abs(fluxDepth_observed))
                    # move WB location down a level; it should be independent of atmosModel
                    out['data'][planetLetter]['WB'] = wavelength_um_rebin
                    out['data'][planetLetter][atmosModel] = {
                        # 'WB':wavelength_um_rebin,
                        'ES':signedSqrt,
                        'ESerr':0.5 * uncertainties / signedSqrt}
                    # 'ES':np.sqrt(fluxDepth_observed),
                    # 'ESerr':0.5 * uncertainties / np.sqrt(fluxDepth_observed)}

                    # cerberus also wants the scale height, to normalize the spectrum
                    #  keep Hs as it's own param (not inside of system_ or model_param)
                    #  it has to be this way to match the formatting for regular spectra itk

                    # redo the chemsitry/mmw calculation for this metallicity
                    # print('pressure',pressure)
                    # print('eqtemp',eqtemp)
                    # print('metallicity [X/H]dex:',model_params['metallicity'])
                    mixratio, fH2, fHe = crbutil.crbce(pressure, eqtemp,
                                                       X2Hr=model_params['metallicity'])
                    # print(' mixratio',mixratio)
                    # print(' fh2,fhe',fH2,fHe)
                    # assume solar C/O and N/O for now
                    # C2Or=cheq['CtoO'], N2Or=cheq['NtoO'])
                    mmwnow, fH2, fHe = crbutil.getmmw(mixratio, protosolar=False, fH2=fH2, fHe=fHe)
                    # print(' mmwnow,mmwsolar',mmwnow,mmwsolar)
                    out['data'][planetLetter][atmosModel]['Hs'] = Hssolar * mmwsolar / mmwnow
                    # print('Hs calculation',Hssolar,mmwsolar,mmwnow)
                    # save the true spectrum (both raw and binned)
                    out['data'][planetLetter][atmosModel]['true_spectrum'] = {
                        'fluxDepth':fluxDepth_rebin,
                        'wavelength_um':wavelength_um_rebin,
                        'fluxDepth_norebin':fluxDepth,
                        'wavelength_norebin':wavelength_um}

                    # also save the Tier level and the number of visits; add these to the plot
                    out['data'][planetLetter][atmosModel]['tier'] = tier
                    out['data'][planetLetter][atmosModel]['visits'] = visits

                    # save the parameters used to create the spectrum. some could be useful
                    # should save both observed value and truth with scatter added in
                    #  'system_params' = the info in system() task
                    #  'model_params' = what is actually used to create forward model
                    out['data'][planetLetter][atmosModel]['system_params'] = {
                        'R*':system_params['R*'],
                        'T*':system_params['T*'],
                        'Rp':system_params[planetLetter]['rp'],
                        'Teq':system_params[planetLetter]['teq'],
                        'Mp':system_params[planetLetter]['mass']}

                    out['data'][planetLetter][atmosModel]['model_params'] = model_params

                    # convert to percentage depth (just for plotting, not for the saved results)
                    fluxDepth = 100 * fluxDepth
                    fluxDepth_rebin = 100 * fluxDepth_rebin
                    fluxDepth_observed = 100 * fluxDepth_observed
                    # careful - uncertainties are reused, so don't change them permanently
                    uncertainties_percent = 100 * uncertainties
                    molecules = fluxDepth_by_molecule_rebin.keys()
                    for molecule in molecules:
                        fluxDepth_by_molecule_rebin[molecule] = \
                            fluxDepth_by_molecule_rebin[molecule] * 100

                    # PLOT THE SPECTRA
                    myfig, ax = plt.subplots(figsize=(8,4))
                    plt.title('Ariel simulation : '+target+' '+planetLetter+' : Tier-'+tier+' '+visits+' visits',
                              fontsize=16)
                    plt.xlabel(str('Wavelength [$\\mu m$]'), fontsize=14)
                    plt.ylabel(str('$(R_p/R_*)^2$ [%]'), fontsize=14)

                    # plot the true (model) spectrum - raw
                    plt.plot(wavelength_um, fluxDepth,
                             # color='palegoldenrod',ls='-',lw=0.1,
                             color='palegoldenrod',ls='-',lw=1,
                             zorder=1, label='truth raw')
                    # plot the true (model) spectrum - binned
                    plt.plot(wavelength_um_rebin, fluxDepth_rebin,
                             color='k',ls='-',lw=1,
                             zorder=3, label='truth binned')
                    for imole,molecule in enumerate(molecules):
                        colorlist = ['red','orange',
                                     'palegreen','lightseagreen',
                                     'blueviolet','fuchsia']
                        stylelist = ['--','--','--','--','--','--',
                                     ':',':',':',':',':',':']
                        plt.plot(wavelength_um_rebin, fluxDepth_by_molecule_rebin[molecule],
                                 color=colorlist[imole % len(colorlist)],
                                 ls=stylelist[imole % len(stylelist)],
                                 lw=1, zorder=2, label=molecule)
                    yrange = plt.ylim()
                    # plot the simulated data points
                    #  maybe leave off the simulated points?
                    #  it's too messy when plotting by molecule
                    # but keep it for taurex, which isn't broken down by molecule
                    if 'taurex' in atmosModel:
                        plt.scatter(wavelength_um_rebin, fluxDepth_observed,
                                    marker='o',s=20, color='None',edgecolor='k',
                                    zorder=4, label='simulated data')
                        plt.errorbar(wavelength_um_rebin, fluxDepth_observed,
                                     yerr=uncertainties_percent,
                                     linestyle='None',lw=0.2, color='grey', zorder=2)
                    plt.ylim(yrange)
                    plt.xlim(0.,8.)
                    plt.legend(loc='center left', bbox_to_anchor=(1.15, 0.48))

                    # add a scale-height-normalized flux scale on the right axis
                    Hsscaling = out['data'][planetLetter][atmosModel]['Hs']
                    # print('H scaling for this plot (%):',Hsscaling*100)
                    ax2 = ax.twinx()
                    ax2.set_ylabel('$\\Delta$ [H]')
                    axmin, axmax = ax.get_ylim()
                    rpmed = np.sqrt(np.nanmedian(1.e-2*fluxDepth_rebin))

                    # print('axmin,axmax',axmin,axmax)
                    # non-convergent spectra can be all NaN, crashing here on the sqrt
                    # should now be stopped up above in that case, but just in case add conditional
                    if np.isnan(np.nanmax(fluxDepth_rebin)):
                        log.warning('--< PROBLEM: spectrum is all NaN %s %s >--',target,planetLetter)
                        pass
                    else:
                        if axmin >= 0:
                            ax2.set_ylim((np.sqrt(1e-2*axmin) - rpmed)/Hsscaling,
                                         (np.sqrt(1e-2*axmax) - rpmed)/Hsscaling)
                        else:
                            ax2.set_ylim((-np.sqrt(-1e-2*axmin) - rpmed)/Hsscaling,
                                         (np.sqrt(1e-2*axmax) - rpmed)/Hsscaling)
                            # print('TROUBLE!! y-axis not scaled by H!!')

                    myfig.tight_layout()

                    # RID = int(os.environ.get('RUNID', None))
                    RID = os.environ.get('RUNID', None)
                    # print('RID',RID)
                    if RID:
                        RID = f'{int(RID):03}'
                    else:
                        RID = '666'
                    # print('RID',RID)

                    # plotDir = excalibur.context['data_dir'] + f"/ariel/RID{RID:03i}"
                    # plotDir = excalibur.context['data_dir'] + f"/ariel/RID{RID:03}"
                    plotDir = excalibur.context['data_dir'] + '/ariel/RID' + RID
                    if not os.path.exists(plotDir): os.mkdir(plotDir)

                    plt.savefig(plotDir +
                                '/ariel_'+atmosModel+'Atmos_' + target+'_'+planetLetter + '.png')

                    # REDUNDANT SAVE - above saves to disk; below saves as state vector
                    # plt.title('Ariel : '+target+' '+planetLetter+'; sv save for RUNID='+RID,
                    #           fontsize=16)
                    buf = io.BytesIO()
                    myfig.savefig(buf, format='png')
                    out['data'][planetLetter][atmosModel]['plot_simspectrum'] = buf.getvalue()

                    plt.close(myfig)

                    completed_at_least_one_planet = True

    # print('completed_at_least_one_planet',completed_at_least_one_planet)
    if completed_at_least_one_planet: out['STATUS'].append(True)

    return True
# ------------------------- ------------------------------------------
