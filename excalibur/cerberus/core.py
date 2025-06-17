'''cerberus core ds'''

# Heritage code shame:
# pylint: disable=too-many-arguments,too-many-branches,too-many-lines,too-many-locals,too-many-nested-blocks,too-many-positional-arguments,too-many-statements

# -- IMPORTS -- ------------------------------------------------------
import dawgie
import excalibur
import excalibur.system.core as syscore
from excalibur.target.targetlists import get_target_lists

# from excalibur.cerberus.core import savesv
from excalibur.cerberus.forward_model import (
    ctxtupdt,
    absorb,
    crbmodel,
    cloudyfmcerberus,
    clearfmcerberus,
    offcerberus,
    offcerberus1,
    offcerberus2,
    offcerberus3,
    offcerberus4,
    offcerberus5,
    offcerberus6,
    offcerberus7,
    offcerberus8,
)

from excalibur.cerberus.plotters import (
    rebin_data,
    plot_corner,
    plot_spectrumfit,
    plot_vs_prior,
    plot_walker_evolution,
    plot_fits_vs_truths,
    plot_fit_uncertainties,
    plot_mass_vs_metals,
)
from excalibur.cerberus.bounds import (
    set_prior_bound,
    add_priors,
    get_profile_limits_hstg141,
    apply_profiling,
)

import logging

import os
import pymc
from pytensor import tensor
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as img
from collections import defaultdict
from collections import namedtuple

from scipy.interpolate import interp1d as itp

log = logging.getLogger(__name__)
pymclog = logging.getLogger('pymc')
pymclog.setLevel(logging.ERROR)

CerbParams = namedtuple(
    'cerberus_params_from_runtime',
    [
        'MCMC_chain_length',
        'MCMC_sampler',
        'fitCloudParameters',
        'fitT',
        'fitCtoO',
        'fitNtoO',
    ],
)


# ----------------- --------------------------------------------------
# -- X SECTIONS LIBRARY -- -------------------------------------------
def myxsecsversion():
    '''
    Alya Al-Kibbi:111
    Changed CH4 line list to HITEMP
    Alya Al-Kibbi:112
    Built interpolator on cross sections assuming constant broadening and
    shifting effects (like exomol)
    Done to speed up processing of CH4 HITEMP line list
    '''
    return dawgie.VERSION(1, 1, 3)


# GMR: Should be in the param list
hitempdir = os.path.join(excalibur.context['data_dir'], 'CERBERUS/HITEMP')
tipsdir = os.path.join(excalibur.context['data_dir'], 'CERBERUS/TIPS')
ciadir = os.path.join(excalibur.context['data_dir'], 'CERBERUS/HITRAN/CIA')
exomoldir = os.path.join(excalibur.context['data_dir'], 'CERBERUS/EXOMOL')


def myxsecs(spc, out, verbose=False):
    '''
    G. ROUDIER: Builds Cerberus cross section library
    '''
    # logarithmicOpacitySumming = True
    logarithmic_opacity_summing = False

    # these used to be default parameters above, but are dangerous-default-values
    knownspecies = ['NO', 'OH', 'C2H2', 'N2', 'N2O', 'O3', 'O2']
    cialist = ['H2-H', 'H2-H2', 'H2-He', 'He-H']
    xmollist = ['TIO', 'H2O', 'H2CO', 'HCN', 'CO', 'CO2', 'NH3', 'CH4']

    cs = False
    planet_letters = []
    for p in spc['data'].keys():
        if (
            len(p) == 1
        ):  # filter out non-planetletter keywords, e.g. 'models','target'
            if (
                'WB' in spc['data'][p].keys()
            ):  # make sure it has a spectrum (Kepler-37e bug)
                planet_letters.append(p)
            else:
                log.warning(
                    '--< CERBERUS.XSLIB: wavelength grid is missing for %s %s >--',
                    spc['data']['target'],
                    p,
                )
    for p in planet_letters:
        out['data'][p] = {}

        # model has to be specified, if there is a list of models
        # if 'models' in spc['data'].keys():
        #    arielModel = spc['data']['models'][0]  # arbitrary model choice; all have same WB grid
        #    wgrid = np.array(spc['data'][p][arielModel]['WB'])
        # else:

        wgrid = np.array(spc['data'][p]['WB'])
        qtgrid = gettpf(knownspecies)
        library = {}

        # EDIT HERE!
        # print('cerb core  spc keys',spc['data'][p]['WB'])
        # exit()

        nugrid = (1e4 / np.copy(wgrid))[::-1]
        dwnu = np.concatenate((np.array([np.diff(nugrid)[0]]), np.diff(nugrid)))
        for myexomol in xmollist:
            # log.warning('>-- %s', str(myexomol))
            library[myexomol] = {
                'I': [],
                'nu': [],
                'T': [],
                'Itemp': [],
                'nutemp': [],
                'Ttemp': [],
                'SPL': [],
                'SPLNU': [],
            }
            thisxmdir = os.path.join(exomoldir, myexomol)
            myfiles = [f for f in os.listdir(thisxmdir) if f.endswith('K')]
            for mf in myfiles:
                xmtemp = float(mf.split('K')[0])
                with open(
                    os.path.join(thisxmdir, mf), 'r', encoding="utf-8"
                ) as fp:
                    data = fp.readlines()
                    fp.close()
                    for line in data:
                        line = np.array(line.split(' '))
                        line = line[line != '']
                        library[myexomol]['nutemp'].append(float(line[0]))
                        library[myexomol]['Itemp'].append(float(line[1]))
                        library[myexomol]['Ttemp'].append(xmtemp)
                        pass
                    pass
                pass
            for mytemp in set(library[myexomol]['Ttemp']):
                select = np.array(library[myexomol]['Ttemp']) == mytemp
                bini = []
                matnu = np.array(library[myexomol]['nutemp'])[select]
                sigma2 = np.array(library[myexomol]['Itemp'])[select]
                for nubin, mydw in zip(nugrid, dwnu):
                    select = (matnu > (nubin - mydw / 2.0)) & (
                        matnu <= nubin + mydw / 2.0
                    )
                    if logarithmic_opacity_summing:
                        # linearsum = np.sum(sigma2[select])
                        logmean = np.mean(np.log(sigma2[select]))
                        nbin = np.sum(select)
                        # print('Nbin',Nbin)
                        bini.append(nbin * np.exp(logmean))
                        # print('old,new',linearsum,Nbin * np.exp(logmean),
                        #      'Nbin,min,max',Nbin,np.min(sigma2[select]),np.max(sigma2[select]))
                    else:
                        bini.append(np.sum(sigma2[select]))
                    pass
                bini = np.array(bini) / dwnu
                library[myexomol]['nu'].extend(list(nugrid))
                library[myexomol]['I'].extend(list(bini))
                library[myexomol]['T'].extend(
                    list(np.ones(nugrid.size) * mytemp)
                )
                pass
            for iline in set(library[myexomol]['nu']):
                select = np.array(library[myexomol]['nu']) == iline
                y = np.array(library[myexomol]['I'])[select]
                x = np.array(library[myexomol]['T'])[select]
                sortme = np.argsort(x)
                x = x[sortme]
                y = y[sortme]
                myspl = itp(x, y, bounds_error=False, fill_value=0)
                library[myexomol]['SPL'].append(myspl)
                library[myexomol]['SPLNU'].append(iline)
                if verbose:
                    plt.plot(x, y, 'o')
                    xp = np.arange(101) / 100.0 * (3000.0 - np.min(x)) + np.min(
                        x
                    )
                    plt.plot(xp, myspl(xp))
                    plt.show()
                    pass
                pass
            if verbose:
                fts = 20
                plt.figure(figsize=(16, 12))
                haha = list(set(library[myexomol]['T']))
                haha = np.sort(np.array(haha))
                haha = haha[::-1]
                for temp in haha:
                    select = np.array(library[myexomol]['T']) == temp
                    plt.semilogy(
                        1e4 / (np.array(library[myexomol]['nu'])[select]),
                        np.array(library[myexomol]['I'])[select],
                        label=str(int(temp)) + 'K',
                    )
                    pass
                plt.title(myexomol)
                plt.xlabel('Wavelength $\\lambda$[$\\mu m$]', fontsize=fts + 4)
                plt.ylabel(
                    'Cross Section [$cm^{2}.molecule^{-1}$]', fontsize=fts + 4
                )
                plt.tick_params(axis='both', labelsize=fts)
                plt.legend(
                    bbox_to_anchor=(0.95, 0.0, 0.12, 1),
                    loc=5,
                    ncol=1,
                    mode='expand',
                    numpoints=1,
                    borderaxespad=0.0,
                    frameon=True,
                )
                plt.savefig(
                    excalibur.context['data_dir']
                    + '/bryden/'
                    + myexomol
                    + '_xslib.png',
                    dpi=200,
                )
                plt.show()
                pass
            pass
        for mycia in cialist:
            # log.warning('>-- %s', str(mycia))
            myfile = '_'.join((os.path.join(ciadir, mycia), '2011.cia'))
            library[mycia] = {
                'I': [],
                'nu': [],
                'T': [],
                'Itemp': [],
                'nutemp': [],
                'Ttemp': [],
                'SPL': [],
                'SPLNU': [],
            }
            with open(myfile, 'r', encoding="utf-8") as fp:
                data = fp.readlines()
                fp.close()
                # Richard et Al. 2012
                tmprtr = 666  # otherwise possibly-used-before-assignment flag
                for line in data:
                    line = np.array(line.split(' '))
                    line = line[line != '']
                    if line.size > 2:
                        tmprtr = float(line[4])
                    if line.size == 2:
                        library[mycia]['nutemp'].append(float(line[0]))
                        library[mycia]['Itemp'].append(float(line[1]))
                        library[mycia]['Ttemp'].append(tmprtr)
                        pass
                    pass
                for mytemp in set(library[mycia]['Ttemp']):
                    select = np.array(library[mycia]['Ttemp']) == mytemp
                    bini = []
                    matnu = np.array(library[mycia]['nutemp'])[select]
                    sigma2 = np.array(library[mycia]['Itemp'])[select]
                    for nubin, mydw in zip(nugrid, dwnu):
                        select = (matnu > (nubin - mydw / 2.0)) & (
                            matnu <= nubin + mydw / 2.0
                        )
                        if logarithmic_opacity_summing:
                            logmean = np.mean(np.log(sigma2[select]))
                            nbin = np.sum(select)
                            bini.append(nbin * np.exp(logmean))
                        else:
                            bini.append(np.sum(sigma2[select]))
                        pass
                    bini = np.array(bini) / dwnu
                    library[mycia]['nu'].extend(list(nugrid))
                    library[mycia]['I'].extend(list(bini))
                    library[mycia]['T'].extend(
                        list(np.ones(nugrid.size) * mytemp)
                    )
                    pass
                for iline in set(library[mycia]['nu']):
                    select = np.array(library[mycia]['nu']) == iline
                    y = np.array(library[mycia]['I'])[select]
                    x = np.array(library[mycia]['T'])[select]
                    sortme = np.argsort(x)
                    x = x[sortme]
                    y = y[sortme]
                    myspl = itp(x, y, bounds_error=False, fill_value=0)
                    library[mycia]['SPL'].append(myspl)
                    library[mycia]['SPLNU'].append(iline)
                    if verbose:
                        plt.plot(x, y, 'o')
                        xp = np.arange(101) / 100.0 * (
                            np.max(x) - np.min(x)
                        ) + np.min(x)
                        plt.plot(xp, myspl(xp))
                        plt.show()
                        pass
                    pass
                pass
            if verbose:
                for temp in set(library[mycia]['T']):
                    select = np.array(library[mycia]['T']) == temp
                    plt.semilogy(
                        1e4 / (np.array(library[mycia]['nu'])[select]),
                        np.array(library[mycia]['I'])[select],
                    )
                    pass
                plt.title(mycia)
                plt.xlabel('Wavelength $\\lambda$[$\\mu m$]')
                plt.ylabel('Line intensity $S(T)$ [$cm^{5}.molecule^{-2}$]')
                plt.show()
                pass
            pass
        for ks in knownspecies:
            # log.warning('>-- %s', str(ks))
            library[ks] = {
                'MU': [],
                'I': [],
                'nu': [],
                'S': [],
                'g_air': [],
                'g_self': [],
                'Epp': [],
                'eta': [],
                'delta': [],
            }
            myfiles = [
                f
                for f in os.listdir(os.path.join(hitempdir, ks))
                if f.endswith('.par')
            ]
            dwmin = abs(wgrid[1] - wgrid[0]) / 2.0
            dwmax = abs(wgrid[-1] - wgrid[-2]) / 2.0
            for fdata in myfiles:
                weqname = fdata.split('_')
                readit = True
                if len(weqname) > 2:
                    if float(weqname[1].split('-')[0]) != 0:
                        maxweq = 1e4 / float(weqname[1].split('-')[0])
                    else:
                        maxweq = 1e20
                    if maxweq < (np.min(wgrid) - dwmin):
                        readit = False
                    minweq = 1e4 / float(weqname[1].split('-')[1])
                    if minweq > (np.max(wgrid) + dwmax):
                        readit = False
                    pass
                if readit:
                    with open(
                        os.path.join(hitempdir, ks, fdata),
                        'r',
                        encoding="utf-8",
                    ) as fp:
                        data = fp.readlines()
                        fp.close()
                        # Rothman et Al. 2010
                        for line in data:
                            waveeq = (1e4) / float(line[3 : 3 + 12])
                            cotest = True
                            if ks == 'H2O':
                                cotest = float(line[15 : 15 + 10]) > 1e-27
                            if ks == 'CO2':
                                cotest = float(line[15 : 15 + 10]) > 1e-29
                            cmintest = waveeq < (np.max(wgrid) + dwmax)
                            cmaxtest = waveeq > (np.min(wgrid) - dwmin)
                            if cmintest and cmaxtest and cotest:
                                library[ks]['MU'].append(waveeq)
                                library[ks]['I'].append(int(line[2 : 2 + 1]))
                                library[ks]['nu'].append(
                                    float(line[3 : 3 + 12])
                                )
                                library[ks]['S'].append(
                                    float(line[15 : 15 + 10])
                                )
                                library[ks]['g_air'].append(
                                    float(line[35 : 35 + 5])
                                )
                                library[ks]['g_self'].append(
                                    float(line[40 : 40 + 5])
                                )
                                library[ks]['Epp'].append(
                                    float(line[45 : 45 + 10])
                                )
                                library[ks]['eta'].append(
                                    float(line[55 : 55 + 4])
                                )
                                library[ks]['delta'].append(
                                    float(line[59 : 59 + 8])
                                )
                                pass
                            pass
                        pass
                    pass
                pass
            if verbose:
                for i in set(library[ks]['I']):
                    select = np.array(library[ks]['I']) == i
                    plt.semilogy(
                        np.array(library[ks]['MU'])[select],
                        np.array(library[ks]['S'])[select],
                        '.',
                    )
                    pass
                plt.title(ks)
                plt.xlabel('Wavelength $\\lambda$[$\\mu m$]')
                plt.ylabel('Line intensity $S_{296K}$ [$cm.molecule^{-1}$]')
                plt.show()
                pass
            # BUILDS INTERPOLATORS SIMILAR TO EXOMOL DB DATA HANDLING
            mmr = 2.3  # Fortney 2015 for hot Jupiters
            solrad = 10.0
            hsmax = 15.0
            # increase the number of scale heights from 15 to 20, to match the Ariel forward model
            # (this is the range used for xslib; also has to be set for atmos)
            hsmax = 20.0
            nlevels = 100.0
            pgrid = np.arange(
                np.log(solrad) - hsmax,
                np.log(solrad) + hsmax / nlevels,
                hsmax / (nlevels - 1),
            )
            pgrid = np.exp(pgrid)
            pressuregrid = pgrid[::-1]
            allxsections = []
            allwavenumbers = []
            alltemperatures = []
            for tstep in np.arange(300, 2000, 100):
                # log.warning('>---- %s K', str(Tstep))
                sigma, lsig = absorb(
                    library[ks],
                    qtgrid[ks],
                    tstep,
                    pressuregrid,
                    mmr,
                    False,
                    False,
                    wgrid,
                    debug=False,
                )
                allxsections.append(sigma[0])
                allwavenumbers.append(lsig)
                alltemperatures.append(tstep)
                pass
            library[ks]['nu'] = []
            library[ks]['I'] = []
            library[ks]['T'] = []
            library[ks]['SPL'] = []
            library[ks]['SPLNU'] = []
            for indextemp, mytemp in enumerate(alltemperatures):
                bini = []
                matnu = np.array(allwavenumbers)[indextemp]
                sigma2 = np.array(allxsections)[indextemp]
                for nubin, mydw in zip(nugrid, dwnu):
                    select = (matnu > (nubin - mydw / 2.0)) & (
                        matnu <= nubin + mydw / 2.0
                    )
                    if logarithmic_opacity_summing:
                        logmean = np.mean(np.log(sigma2[select]))
                        nbin = np.sum(select)
                        bini.append(nbin * np.exp(logmean))
                    else:
                        bini.append(np.sum(sigma2[select]))
                    pass
                bini = np.array(bini) / dwnu
                library[ks]['nu'].extend(list(nugrid))
                library[ks]['I'].extend(list(bini))
                library[ks]['T'].extend(list(np.ones(nugrid.size) * mytemp))
                pass
            for iline in set(library[ks]['nu']):
                select = np.array(library[ks]['nu']) == iline
                y = np.array(library[ks]['I'])[select]
                x = np.array(library[ks]['T'])[select]
                sortme = np.argsort(x)
                x = x[sortme]
                y = y[sortme]
                myspl = itp(x, y, bounds_error=False, fill_value=0)
                library[ks]['SPL'].append(myspl)
                library[ks]['SPLNU'].append(iline)
                pass
            if verbose:
                fts = 20
                plt.figure(figsize=(16, 12))
                haha = list(set(library[ks]['T']))
                haha = np.sort(np.array(haha))
                haha = haha[::-1]
                for temp in haha:
                    select = np.array(library[ks]['T']) == temp
                    plt.semilogy(
                        1e4 / (np.array(library[ks]['nu'])[select]),
                        np.array(library[ks]['I'])[select],
                        label=str(int(temp)) + 'K',
                    )
                    pass
                plt.title(ks)
                plt.xlabel('Wavelength $\\lambda$[$\\mu m$]', fontsize=fts + 4)
                plt.ylabel(
                    'Cross Section [$cm^{2}.molecule^{-1}$]', fontsize=fts + 4
                )
                plt.tick_params(axis='both', labelsize=fts)
                plt.legend(
                    bbox_to_anchor=(0.95, 0.0, 0.12, 1),
                    loc=5,
                    ncol=1,
                    mode='expand',
                    numpoints=1,
                    borderaxespad=0.0,
                    frameon=True,
                )
                plt.show()
                pass
            pass
        out['data'][p]['XSECS'] = library
        out['data'][p]['QTGRID'] = qtgrid
        pass
    if out['data'].keys():
        cs = True
        out['STATUS'].append(True)
        pass
    return cs


# ------------------------ -------------------------------------------
# -- TOTAL PARTITION FUNCTION -- -------------------------------------
def gettpf(knownspecies, verbose=False):
    '''
    G. ROUDIER: Wrapper around HITRAN partition functions (Gamache et al. 2011)
    '''
    grid = {}
    tempgrid = list(np.arange(60.0, 3035.0, 25.0))

    for ks in knownspecies:
        grid[ks] = {'T': tempgrid, 'Q': [], 'SPL': []}
        with open(os.path.join(tipsdir, ks), 'r', encoding="utf-8") as fp:
            data = fp.readlines()
            fp.close()
            for line in data:
                grid[ks]['Q'].append([float(num) for num in line.split(',')])
                pass
            pass
        for y in grid[ks]['Q']:
            myspl = itp(tempgrid, y)
            grid[ks]['SPL'].append(myspl)
            if verbose:
                plt.plot(tempgrid, myspl(tempgrid))
                plt.plot(tempgrid, y, '+')
                pass
            pass
        if verbose:
            plt.title(ks)
            plt.xlabel('Temperature T[K]')
            plt.ylabel('Total Internal Partition Function Q')
            plt.show()
        pass
    return grid


# ------------------------------ -------------------------------------
# -- ATMOS -- --------------------------------------------------------
def atmosversion():
    '''
    Alya Al-Kibbi:121
    Changes so that retrieval is done with CH4 being from HITEMP list instead of Exomol list
    R ESTRELA:131
    Merged Spectra Capability
    '''
    return dawgie.VERSION(1, 3, 2)


def atmos(
    fin,
    xsl,
    spc,
    runtime_params,
    out,
    ext,
    hazedir=os.path.join(excalibur.context['data_dir'], 'CERBERUS/HAZE'),
    singlemod=None,
    mclen=int(1e4),
    verbose=False,
):
    '''
    G. ROUDIER: Cerberus retrieval
    '''

    okfit = False
    orbp = fin['priors'].copy()

    ssc = syscore.ssconstants(mks=True)
    crbhzlib = {'PROFILE': []}
    hazelib(crbhzlib, hazedir=hazedir, verbose=False)
    # SELECT WHICH MODELS TO RUN FOR THIS FILTER
    if ext == 'Ariel-sim':
        modfam = ['TEC']  # Ariel sims are currently only TEC equilibrium models
        # modparlbl = {'TEC':['XtoH', 'CtoO']}
        modparlbl = {'TEC': ['XtoH', 'CtoO', 'NtoO']}

        # ** select which Ariel model to fit **
        #   previously (with taurex) there were 8 options. now 4 options:
        # atmosModels = ['cerberus', 'cerberusNoclouds',
        #                'cerberuslowmmw', 'cerberuslowmmwNoclouds']
        # arielModel = 'cerberusNoclouds'
        arielmodel = 'cerberus'

        # option to fix N/O
        if not runtime_params.fitNtoO:
            modparlbl = {'TEC': ['XtoH', 'CtoO']}
        # option to fix C/O
        if not runtime_params.fitCtoO:
            modparlbl = {'TEC': ['XtoH']}

        # print('name of the forward model:',arielModel)
        # print('available models',spc['data']['models'])
        if arielmodel not in spc['data']['models']:
            log.warning('--< BIG PROB: ariel model doesnt exist!!! >--')
    else:
        modfam = ['TEC', 'PHOTOCHEM']
        modparlbl = {
            'TEC': ['XtoH', 'CtoO', 'NtoO'],
            'PHOTOCHEM': ['HCN', 'CH4', 'C2H2', 'CO2', 'H2CO'],
        }
        if not runtime_params.fitNtoO:
            modparlbl['TEC'].remove('NtoO')
        if not runtime_params.fitCtoO:
            modparlbl['TEC'].remove('CtoO')

    if (singlemod is not None) and (singlemod in modfam):
        modfam = [modfam[modfam.index(singlemod)]]

    # save the stellar params, so that analysis knows the stellar metallicity
    star_params = [
        'M*',
        'R*',
        'LOGG*',
        'RHO*',
        'FEH*',
        'T*',
        'L*',
        'Jmag',
        'Hmag',
        'Kmag',
        'spTyp',
        'AGE*',
        'dist',
    ]
    out['data']['stellar_params'] = {}
    for star_param in star_params:
        out['data']['stellar_params'][star_param] = orbp[star_param]

    # PLANET LOOP
    for p in spc['data'].keys():
        # make sure that it really is a planet letter, not another dict key
        #  (ariel has other keys, e.g. 'target', 'planets', 'models')
        # make sure it has a spectrum (Kepler-37e bug)
        if len(p) == 1 and 'WB' not in spc['data'][p].keys():
            log.warning(
                '--< CERBERUS.ATMOS: wavelength grid is missing for %s %s >--',
                spc['data']['target'],
                p,
            )
        elif len(p) == 1 and 'WB' in spc['data'][p].keys():
            if ext == 'Ariel-sim':
                if arielmodel in spc['data'][p].keys():
                    input_data = spc['data'][p][arielmodel]
                    # make sure that the wavelength is saved in usual location
                    # (the cerberus forward models expect it to come after [p])
                    # spc['data'][p]['WB'] = spc['data'][p][arielModel]['WB']
                    input_data['WB'] = spc['data'][p]['WB']
                else:
                    log.warning(
                        '--< THIS arielModel DOESNT EXIST!!! (rerun ariel task?) >--'
                    )
            else:
                input_data = spc['data'][p]

            out['data'][p] = {}
            out['data'][p]['MODELPARNAMES'] = modparlbl

            # save the planet params (mass), so that analysis can make a mass-metallicity plot
            out['data'][p]['planet_params'] = orbp[p]
            # save the tier and #-of-visits (for Ariel-sim targets), for plot labelling
            if ext == 'Ariel-sim':
                out['data'][p]['tier'] = spc['data'][p][arielmodel]['tier']
                out['data'][p]['visits'] = spc['data'][p][arielmodel]['visits']

            # eqtemp1 = orbp['T*']*np.sqrt(orbp['R*']*ssc['Rsun/AU']/(2.*orbp[p]['sma']))
            # use of L* might be better than T*,R*; it's more of a direct observable
            # eqtemp2 = 278.6 * orbp['L*']**0.25 / np.sqrt(orbp[p]['sma'])

            # These four estimates of T_equilibrium are usually about the same
            #  but sometimes they are very different
            #  system/consistency.py checks for and flags these inconsistencies
            # print('eqtemp check',eqtemp1)
            # print('eqtemp check',eqtemp2)
            # print('eqtemp check',orbp[p]['teq'])
            # print('eqtemp check',inputData['model_params']['Teq'])

            # CAREFUL:
            #  equilibrium temperatures from the archive sometimes don't match this!
            #  e.g. GJ 3053=LHS 1140 b,c (the Archive has 379K,709K vs 216,403K here)
            #  that one is wrong it seems. Lillo-Box 2020 has strange extra factor

            # print('model_params',inputData['model_params'])

            # bottom line:
            #  use the same Teq as in ariel-sim, otherwise truth/retrieved won't match
            if ext == 'Ariel-sim':
                eqtemp = input_data['model_params']['Teq']
            else:
                # (real data doesn't have any 'model_params' defined)
                # eqtemp = orbp['T*']*np.sqrt(orbp['R*']*ssc['Rsun/AU']/(2.*orbp[p]['sma']))
                eqtemp = float(orbp[p]['teq'])

            tspc = np.array(input_data['ES'])
            terr = np.array(input_data['ESerr'])
            twav = np.array(input_data['WB'])
            # twav = np.array(spc['data'][p]['WB'])

            tspecerr = abs(tspc**2 - (tspc + terr) ** 2)
            tspectrum = tspc**2
            if 'STIS-WFC3' in ext:
                filters = np.array(input_data['Fltrs'])
                cond_spec_g750 = filters == 'HST-STIS-CCD-G750L-STARE'
                # MASKING G750 WAV > 0.80
                twav_g750 = twav[cond_spec_g750]
                tspec_g750 = tspectrum[cond_spec_g750]
                tspecerr_g750 = tspecerr[cond_spec_g750]
                mask = (twav_g750 > 0.80) & (twav_g750 < 0.95)
                tspec_g750[mask] = np.nan
                tspecerr_g750[mask] = np.nan
                tspectrum[cond_spec_g750] = tspec_g750
                tspecerr[cond_spec_g750] = tspecerr_g750
            hs = input_data['Hs']
            #  Clean up
            spechs = (
                np.sqrt(tspectrum) - np.sqrt(np.nanmedian(tspectrum))
            ) / hs
            cleanup2 = abs(spechs) > 3e0  # excluding everything above +-3 Hs
            tspectrum[cleanup2] = np.nan
            tspecerr[cleanup2] = np.nan
            twav[cleanup2] = np.nan
            # cleanup = np.isfinite(tspectrum) & (tspecerr < 1e0)
            cleanup = np.isfinite(tspectrum)
            solidr = orbp[p]['rp'] * ssc['Rjup']  # MK

            for model in modfam:
                ctxtupdt(
                    cleanup=cleanup,
                    model=model,
                    p=p,
                    solidr=solidr,
                    orbp=orbp,
                    tspectrum=tspectrum,
                    xsl=xsl,
                    spc=spc,
                    modparlbl=modparlbl,
                    hzlib=crbhzlib,
                )
                out['data'][p][model] = {}

                # new method for setting priors (no change, but easier to view in bounds.py)
                prior_range_table = set_prior_bound(eqtemp)

                out['data'][p][model]['prior_ranges'] = {}
                # keep track of the bounds put on each parameter
                # this will be helpful for later plotting and analysis
                nodes = []
                with pymc.Model():

                    # set the fixed parameters (the ones that are not being fit this time)
                    fixed_params = {}
                    if not runtime_params.fitCloudParameters:
                        fixed_params['CTP'] = input_data['model_params']['CTP']
                        fixed_params['HScale'] = input_data['model_params'][
                            'HScale'
                        ]
                        fixed_params['HLoc'] = input_data['model_params'][
                            'HLoc'
                        ]
                        fixed_params['HThick'] = input_data['model_params'][
                            'HThick'
                        ]
                    if not runtime_params.fitT:
                        fixed_params['T'] = input_data['model_params']['T']
                    if not runtime_params.fitNtoO:
                        fixed_params['NtoO'] = 0.0
                    if not runtime_params.fitCtoO:
                        fixed_params['CtoO'] = input_data['model_params']['C/O']
                    # print('fixedparams',fixedParams)

                    # OFFSET BETWEEN STIS AND WFC3 filters
                    if 'STIS-WFC3' in ext:
                        cond_off0 = filters == 'HST-STIS-CCD-G430L-STARE'
                        cond_off1 = filters == 'HST-STIS-CCD-G750L-STARE'
                        cond_off2 = filters == 'HST-WFC3-IR-G102-SCAN'
                        cond_off3 = filters == 'HST-WFC3-IR-G141-SCAN'
                        valid0 = True in cond_off0
                        valid1 = True in cond_off1
                        valid2 = True in cond_off2
                        valid3 = True in cond_off3
                        if 'STIS' in filters[0]:
                            if valid0:  # G430
                                if (
                                    valid1 and valid2 and valid3
                                ):  # G430-G750-G102-G141
                                    off0_value = abs(
                                        np.nanmedian(1e2 * tspectrum[cond_off3])
                                        - np.nanmedian(
                                            1e2 * tspectrum[cond_off0]
                                        )
                                    )
                                    off1_value = abs(
                                        np.nanmedian(1e2 * tspectrum[cond_off3])
                                        - np.nanmedian(
                                            1e2 * tspectrum[cond_off1]
                                        )
                                    )
                                    off2_value = abs(
                                        np.nanmedian(1e2 * tspectrum[cond_off3])
                                        - np.nanmedian(
                                            1e2 * tspectrum[cond_off2]
                                        )
                                    )
                                    nodes.append(
                                        pymc.Uniform(
                                            'OFF0', -off0_value, off0_value
                                        )
                                    )
                                    nodes.append(
                                        pymc.Uniform(
                                            'OFF1', -off1_value, off1_value
                                        )
                                    )
                                    nodes.append(
                                        pymc.Uniform(
                                            'OFF2', -off2_value, off2_value
                                        )
                                    )
                                elif valid1 and valid2 and not valid3:
                                    off0_value = abs(
                                        np.nanmedian(1e2 * tspectrum[cond_off2])
                                        - np.nanmedian(
                                            1e2 * tspectrum[cond_off0]
                                        )
                                    )
                                    off1_value = abs(
                                        np.nanmedian(1e2 * tspectrum[cond_off2])
                                        - np.nanmedian(
                                            1e2 * tspectrum[cond_off1]
                                        )
                                    )
                                    nodes.append(
                                        pymc.Uniform(
                                            'OFF0', -off0_value, off0_value
                                        )
                                    )
                                    nodes.append(
                                        pymc.Uniform(
                                            'OFF1', -off1_value, off1_value
                                        )
                                    )
                                elif valid1 and valid3 and not valid2:
                                    off0_value = abs(
                                        np.nanmedian(1e2 * tspectrum[cond_off3])
                                        - np.nanmedian(
                                            1e2 * tspectrum[cond_off0]
                                        )
                                    )
                                    off1_value = abs(
                                        np.nanmedian(1e2 * tspectrum[cond_off3])
                                        - np.nanmedian(
                                            1e2 * tspectrum[cond_off1]
                                        )
                                    )
                                    nodes.append(
                                        pymc.Uniform(
                                            'OFF0', -off0_value, off0_value
                                        )
                                    )
                                    nodes.append(
                                        pymc.Uniform(
                                            'OFF1', -off1_value, off1_value
                                        )
                                    )
                                elif valid2 and valid3 and not valid1:
                                    off0_value = abs(
                                        np.nanmedian(1e2 * tspectrum[cond_off3])
                                        - np.nanmedian(
                                            1e2 * tspectrum[cond_off0]
                                        )
                                    )
                                    off1_value = abs(
                                        np.nanmedian(1e2 * tspectrum[cond_off3])
                                        - np.nanmedian(
                                            1e2 * tspectrum[cond_off2]
                                        )
                                    )
                                    nodes.append(
                                        pymc.Uniform(
                                            'OFF0', -off0_value, off0_value
                                        )
                                    )
                                    nodes.append(
                                        pymc.Uniform(
                                            'OFF1', -off1_value, off1_value
                                        )
                                    )
                                elif valid3 and not valid1 and not valid2:
                                    off0_value = abs(
                                        np.nanmedian(1e2 * tspectrum[cond_off3])
                                        - np.nanmedian(
                                            1e2 * tspectrum[cond_off0]
                                        )
                                    )
                                    nodes.append(
                                        pymc.Uniform(
                                            'OFF0', -off0_value, off0_value
                                        )
                                    )
                            else:
                                if valid1 and valid2 and valid3:
                                    off0_value = abs(
                                        np.nanmedian(1e2 * tspectrum[cond_off3])
                                        - np.nanmedian(
                                            1e2 * tspectrum[cond_off1]
                                        )
                                    )
                                    off1_value = abs(
                                        np.nanmedian(1e2 * tspectrum[cond_off3])
                                        - np.nanmedian(
                                            1e2 * tspectrum[cond_off2]
                                        )
                                    )
                                    nodes.append(
                                        pymc.Uniform(
                                            'OFF0', -off0_value, off0_value
                                        )
                                    )
                                    nodes.append(
                                        pymc.Uniform(
                                            'OFF1', -off1_value, off1_value
                                        )
                                    )
                                if valid1 and valid3 and not valid2:
                                    off0_value = abs(
                                        np.nanmedian(1e2 * tspectrum[cond_off3])
                                        - np.nanmedian(
                                            1e2 * tspectrum[cond_off1]
                                        )
                                    )
                                    nodes.append(
                                        pymc.Uniform(
                                            'OFF0', -off0_value, off0_value
                                        )
                                    )
                                if valid1 and valid2 and not valid3:
                                    off0_value = abs(
                                        np.nanmedian(1e2 * tspectrum[cond_off2])
                                        - np.nanmedian(
                                            1e2 * tspectrum[cond_off1]
                                        )
                                    )
                                    nodes.append(
                                        pymc.Uniform(
                                            'OFF0', -off0_value, off0_value
                                        )
                                    )
                                if valid1 and valid2 and not valid3:
                                    off0_value = abs(
                                        np.nanmedian(1e2 * tspectrum[cond_off2])
                                        - np.nanmedian(
                                            1e2 * tspectrum[cond_off1]
                                        )
                                    )
                                    nodes.append(
                                        pymc.Uniform(
                                            'OFF0', -off0_value, off0_value
                                        )
                                    )
                        if 'WFC3' in filters[0]:
                            if valid2 and valid3:
                                off0_value = abs(
                                    np.nanmedian(1e2 * tspectrum[cond_off3])
                                    - np.nanmedian(1e2 * tspectrum[cond_off2])
                                )
                                nodes.append(
                                    pymc.Uniform(
                                        'OFF0', -off0_value, off0_value
                                    )
                                )

                    # new cleaned-up version of adding on the prior bounds as pymc nodes
                    nodes, prior_ranges = add_priors(
                        prior_range_table,
                        runtime_params,
                        model,
                        modparlbl[model],
                    )

                    # before calling MCMC, save the fixed-parameter info in the context
                    ctxtupdt(
                        cleanup=cleanup,
                        model=model,
                        p=p,
                        solidr=solidr,
                        orbp=orbp,
                        tspectrum=tspectrum,
                        xsl=xsl,
                        spc=spc,
                        modparlbl=modparlbl,
                        hzlib=crbhzlib,
                        fixed_params=fixed_params,
                    )

                    # CERBERUS MCMC
                    if not runtime_params.fitCloudParameters:
                        log.warning('--< RUNNING MCMC - NO CLOUDS! >--')
                        _ = pymc.Normal(
                            'mcdata',
                            mu=clearfmcerberus(*nodes),
                            tau=1e0 / (np.nanmedian(tspecerr[cleanup]) ** 2),
                            observed=tspectrum[cleanup],
                        )
                        pass
                    else:
                        if 'STIS-WFC3' in ext:
                            if 'STIS' in filters[0]:
                                if valid0:  # G430
                                    if (
                                        valid1 and valid2 and valid3
                                    ):  # G430-G750-G102-G141
                                        _ = pymc.Normal(
                                            'mcdata',
                                            mu=offcerberus(*nodes),
                                            tau=1e0 / tspecerr[cleanup] ** 2,
                                            observed=tspectrum[cleanup],
                                        )
                                    elif valid1 and valid2 and not valid3:
                                        _ = pymc.Normal(
                                            'mcdata',
                                            mu=offcerberus1(*nodes),
                                            tau=1e0 / tspecerr[cleanup] ** 2,
                                            observed=tspectrum[cleanup],
                                        )
                                    elif valid1 and valid3 and not valid2:
                                        _ = pymc.Normal(
                                            'mcdata',
                                            mu=offcerberus2(*nodes),
                                            tau=1e0 / tspecerr[cleanup] ** 2,
                                            observed=tspectrum[cleanup],
                                        )
                                    elif valid2 and valid3 and not valid1:
                                        _ = pymc.Normal(
                                            'mcdata',
                                            mu=offcerberus3(*nodes),
                                            tau=1e0 / tspecerr[cleanup] ** 2,
                                            observed=tspectrum[cleanup],
                                        )
                                    elif valid3 and not valid1 and not valid2:
                                        _ = pymc.Normal(
                                            'mcdata',
                                            mu=offcerberus4(*nodes),
                                            tau=1e0 / tspecerr[cleanup] ** 2,
                                            observed=tspectrum[cleanup],
                                        )
                                else:
                                    if valid1 and valid2 and valid3:
                                        _ = pymc.Normal(
                                            'mcdata',
                                            mu=offcerberus5(*nodes),
                                            tau=1e0 / tspecerr[cleanup] ** 2,
                                            observed=tspectrum[cleanup],
                                        )
                                    elif valid1 and valid3 and not valid2:
                                        _ = pymc.Normal(
                                            'mcdata',
                                            mu=offcerberus6(*nodes),
                                            tau=1e0 / tspecerr[cleanup] ** 2,
                                            observed=tspectrum[cleanup],
                                        )
                                    elif valid1 and valid2 and not valid3:
                                        _ = pymc.Normal(
                                            'mcdata',
                                            mu=offcerberus7(*nodes),
                                            tau=1e0 / tspecerr[cleanup] ** 2,
                                            observed=tspectrum[cleanup],
                                        )
                            if 'WFC3' in filters[0]:
                                if valid2 and valid3:
                                    _ = pymc.Normal(
                                        'mcdata',
                                        mu=offcerberus8(*nodes),
                                        tau=1e0 / tspecerr[cleanup] ** 2,
                                        observed=tspectrum[cleanup],
                                    )
                                elif not valid2:
                                    _ = pymc.Normal(
                                        'mcdata',
                                        mu=cloudyfmcerberus(*nodes),
                                        tau=1e0
                                        / (
                                            np.nanmedian(tspecerr[cleanup]) ** 2
                                        ),
                                        observed=tspectrum[cleanup],
                                    )
                                elif not valid3:
                                    _ = pymc.Normal(
                                        'mcdata',
                                        mu=cloudyfmcerberus(*nodes),
                                        tau=1e0
                                        / (
                                            np.nanmedian(tspecerr[cleanup]) ** 2
                                        ),
                                        observed=tspectrum[cleanup],
                                    )
                                    pass
                                pass
                        if 'STIS-WFC3' not in ext:
                            log.warning('--< STANDARD MCMC (WITH CLOUDS) >--')
                            _ = pymc.Normal(
                                'mcdata',
                                mu=cloudyfmcerberus(*nodes),
                                sigma=tspecerr[cleanup],
                                observed=tspectrum[cleanup],
                            )
                            # tau=1e0/tspecerr[cleanup]**2,
                            # tau=1e0/(np.nanmedian(tspecerr[cleanup])**2),
                            # this (_mcdata) is similar to the input spectrum.
                            #   maybe it's the final forward model?  (there's only one)
                            #   oh right, this is just the definition; there's no sampling yet
                            #   so why does it bother making one call.  what param values?
                            # print('_mcdata (what is this?!)',_mcdata.eval())
                        pass

                    if runtime_params.MCMC_sampler == 'slice':
                        log.warning('>-- SLICE SAMPLER: ON  --<')
                        sampler = pymc.Slice()
                    else:
                        log.warning('>-- SLICE SAMPLER: OFF --<')
                        sampler = pymc.Metropolis()

                    log.warning(
                        '>-- MCMC nodes: %s', str([n.name for n in nodes])
                    )
                    trace = pymc.sample(
                        mclen,
                        cores=4,
                        tune=int(mclen / 4),
                        step=sampler,
                        compute_convergence_checks=False,
                        progressbar=verbose,
                    )

                    stats_summary = pymc.stats.summary(trace)

                # N_TEC = len(trace.posterior.TEC_dim_0)
                # print('# of TEC parameters',N_TEC)

                mctrace = {}
                for key in stats_summary['mean'].keys():
                    # print('key', key)
                    tracekeys = key.split('[')
                    keyname = tracekeys[0]
                    # print('tracekeys,keyname', tracekeys, keyname)
                    if len(tracekeys) > 1:
                        # print('tracekeys', tracekeys)
                        param_index = int(tracekeys[1][:-1])
                        mctrace[key] = trace.posterior[keyname][
                            :, :, param_index
                        ]
                    else:
                        mctrace[key] = trace.posterior[key]
                    print('mctrace shape', key, mctrace[key].shape)

                    # convert Nchain x Nstep 2-D posteriors to a single chain
                    # mctrace[key] = np.ravel(mctrace[key])
                    # original ravel is the correct one?  it's hard to tell with 30 steps
                    # seems better when reversed ('F' reverses the indices)
                    mctrace[key] = np.ravel(mctrace[key], order='F')
                out['data'][p][model]['MCTRACE'] = mctrace

                out['data'][p][model]['prior_ranges'] = prior_ranges
            out['data'][p]['WAVELENGTH'] = np.array(input_data['WB'])
            out['data'][p]['SPECTRUM'] = np.array(input_data['ES'])
            out['data'][p]['ERRORS'] = np.array(input_data['ESerr'])
            if ext == 'Ariel-sim':
                if 'true_spectrum' in input_data.keys():

                    out['data'][p]['TRUTH_SPECTRUM'] = np.array(
                        input_data['true_spectrum']['fluxDepth']
                    )
                    # wavelength should be the same as just above, but just in case load it here too
                    out['data'][p]['TRUTH_WAVELENGTH'] = np.array(
                        input_data['true_spectrum']['wavelength_um']
                    )
                    out['data'][p]['TRUTH_MODELPARAMS'] = input_data[
                        'model_params'
                    ]
                    # print('true modelparams in atmos:',inputData['model_params'])
            out['data'][p]['VALID'] = cleanup
            out['STATUS'].append(True)

            okfit = True
            pass
    return okfit


# ----------------------------------- --------------------------------
# -- HAZE DENSITY PROFILE LIBRARY -- ---------------------------------
def hazelib(
    sv,
    hazedir=os.path.join(excalibur.context['data_dir'], 'CERBERUS/HAZE'),
    datafile='Jup-ISS-aerosol.dat',
    verbose=False,
    fromjupiter=True,
    narrow=True,
):
    '''Haze density profiles'''
    vdensity = {
        'PRESSURE': [],
        'CONSTANT': [],
        'JMAX': [],
        'MAX': [],
        'JMEDIAN': [],
        'MEDIAN': [],
        'JAVERAGE': [],
        'AVERAGE': [],
    }
    with open(os.path.join(hazedir, datafile), 'r', encoding="utf-8") as fp:
        data = fp.readlines()
        pass
    # LATITUDE GRID
    latitude = data[0]
    latitude = np.array(latitude.split(' '))
    latitude = latitude[latitude != '']
    latitude = latitude.astype(float)  # [DEGREES]
    # DENSITY MATRIX
    pressure = []  # [mbar]
    density = []
    for line in data[1:]:
        line = np.array(line.split(' '))
        line = line[line != '']
        line = line.astype(float)
        pressure.append(line[0])
        density.append(line[1:])
        pass
    # JUPITER PROFILES
    if fromjupiter:
        pressure = 1e-3 * np.array(pressure)  # [bar]
        density = 1e-6 * np.array(density)  # [n/m^3]
        jmax = np.nanmax(density, 1)
        jmed = np.median(density, 1)
        jav = np.mean(density, 1)
        isobar = pressure * 0.0 + np.mean(jav)
        sortme = np.argsort(pressure)
        jmaxspl = itp(
            pressure[sortme],
            jmax[sortme],
            kind='linear',
            bounds_error=False,
            fill_value=0e0,
        )
        jmedspl = itp(
            pressure[sortme],
            jmed[sortme],
            kind='linear',
            bounds_error=False,
            fill_value=0e0,
        )
        javspl = itp(
            pressure[sortme],
            jav[sortme],
            kind='linear',
            bounds_error=False,
            fill_value=0e0,
        )
        pass
    else:
        if narrow:
            density = [
                0.322749941,
                0.618855379,
                0.659653289,
                0.660104488,
                0.687461411,
                0.701591071,
                0.823367371,
                0.998717644,
                1.254642603,
                1.483661838,
                1.726335787,
                1.928591783,
                2.157824745,
                2.387176443,
                2.58959867,
                2.778603657,
                2.981049632,
                3.237164569,
                3.573854191,
                3.9373308,
                4.273759202,
                4.448824507,
                4.341700309,
                4.019449062,
                3.683756827,
                3.227214438,
                2.891522204,
                2.461790549,
                1.978247447,
                1.481168369,
                1.010947518,
                0.500379957,
                0.030087865,
                -0.252101639,
            ]
            pressure = [
                1.874253147,
                1.632296367,
                1.349380195,
                1.080826407,
                0.797986227,
                0.388012349,
                -0.09324151,
                -0.461724056,
                -0.788259321,
                -1.100508193,
                -1.540042745,
                -1.922811684,
                -2.362270245,
                -2.872400855,
                -3.354110663,
                -3.849878889,
                -4.345723106,
                -4.78533365,
                -5.182996913,
                -5.524274519,
                -5.766459273,
                -5.9653289,
                -6.205005937,
                -6.40106388,
                -6.597045832,
                -6.863015911,
                -7.058997863,
                -7.282716694,
                -7.47786274,
                -7.616395156,
                -7.740945144,
                -7.851132748,
                -7.933279506,
                -7.974086915,
            ]
            pass
        else:
            density = [
                -1.222929936,
                -0.76433121,
                -0.229299363,
                0.025477707,
                0.25477707,
                0.789808917,
                1.503184713,
                2.089171975,
                2.369426752,
                2.522292994,
                2.675159236,
                2.904458599,
                3.133757962,
                3.337579618,
                3.541401274,
                3.694267516,
                3.796178344,
                3.821656051,
                3.898089172,
                3.923566879,
                4.0,
                4.025477707,
                4.050955414,
                4.127388535,
                4.152866242,
                4.178343949,
                4.127388535,
                3.974522293,
                3.796178344,
                3.541401274,
                3.261146497,
                2.904458599,
                2.624203822,
                2.267515924,
                1.961783439,
                1.630573248,
                1.299363057,
                0.968152866,
                0.636942675,
                0.280254777,
                -0.076433121,
                -0.433121019,
                -0.76433121,
                -1.070063694,
                -1.299363057,
            ]
            pressure = [
                2.812911584,
                2.809025154,
                2.63499946,
                2.327755587,
                2.156320846,
                2.083990068,
                1.976249595,
                1.835690381,
                1.528230595,
                1.086257152,
                0.678182014,
                0.337255749,
                0.030227788,
                -0.310482565,
                -0.718989528,
                -1.059268056,
                -1.534707978,
                -1.941703552,
                -2.518622477,
                -2.85782144,
                -3.33304545,
                -3.807837634,
                -4.214833207,
                -4.690057217,
                -5.097052791,
                -5.60574328,
                -6.012091115,
                -6.4175753,
                -6.788945266,
                -7.057972579,
                -7.292885674,
                -7.493252726,
                -7.626470906,
                -7.725143042,
                -7.790348699,
                -7.855338443,
                -7.954226492,
                -8.019216237,
                -8.118104286,
                -8.182878117,
                -8.281550254,
                -8.38022239,
                -8.411313829,
                -8.476519486,
                -8.508474576,
            ]
            pass
        density = 1e-6 * (10 ** (np.array(density)))
        density = np.array([density, density])
        density = density.T
        jmax = np.nanmax(density, 1)
        jmed = np.median(density, 1)
        jav = np.mean(density, 1)
        pressure = 10 ** (np.array(pressure))
        isobar = pressure * 0.0 + np.mean(jav)
        sortme = np.argsort(pressure)
        jmaxspl = itp(
            pressure[sortme],
            jmax[sortme],
            kind='linear',
            bounds_error=False,
            fill_value=0e0,
        )
        jmedspl = itp(
            pressure[sortme],
            jmed[sortme],
            kind='linear',
            bounds_error=False,
            fill_value=0e0,
        )
        javspl = itp(
            pressure[sortme],
            jav[sortme],
            kind='linear',
            bounds_error=False,
            fill_value=0e0,
        )
        pass
    if verbose:
        myfig = plt.figure(figsize=(12, 6))
        for vmr in density.T:
            plt.plot(1e6 * vmr, pressure, 'k+')
            plt.plot(1e6 * vmr, pressure)
            pass
        plt.xlabel('Aerosol Density [$n.{cm}^{-3}$]', fontsize=24)
        plt.ylabel('Pressure [bar]', fontsize=24)
        plt.semilogy()
        plt.semilogx()
        plt.gca().invert_yaxis()
        plt.tick_params(axis='both', labelsize=20)
        plt.xlim([1e-10, 1e4])
        plt.title('Latitudinal Variations', fontsize=24)
        myfig.tight_layout()

        myfig = plt.figure(figsize=(12, 6))
        plt.plot(1e6 * jmax, pressure, label='Max', color='blue')
        plt.plot(1e6 * jmed, pressure, label='Median', color='red')
        plt.plot(1e6 * jav, pressure, label='Average', color='green')
        plt.plot(1e6 * isobar, pressure, label='Constant', color='black')
        plt.legend(loc=3, frameon=False, fontsize=24)
        plt.plot(1e6 * jmaxspl(pressure), pressure, 'bo')
        plt.plot(1e6 * jmedspl(pressure), pressure, 'ro')
        plt.plot(1e6 * javspl(pressure), pressure, 'go')
        plt.xlabel('Aerosol Density [$n.{cm}^{-3}$]', fontsize=24)
        plt.ylabel('Pressure [bar]', fontsize=24)
        plt.semilogy()
        plt.semilogx()
        plt.gca().invert_yaxis()
        plt.xlim([1e-4, 1e5])
        plt.tick_params(axis='both', labelsize=20)
        plt.title('Profile Library', fontsize=24)
        myfig.tight_layout()

        plt.show()
        pass
    vdensity['PRESSURE'] = list(pressure)
    vdensity['CONSTANT'] = list(isobar)
    vdensity['JMAX'] = list(jmax)
    vdensity['MAX'].append(jmaxspl)
    vdensity['JMEDIAN'] = list(jmed)
    vdensity['MEDIAN'].append(jmedspl)
    vdensity['JAVERAGE'] = list(jav)
    vdensity['AVERAGE'].append(javspl)
    sv['PROFILE'].append(vdensity)
    return


# ---------------------------------- ---------------------------------
# -- RESULTS ----------------------------------------------------------
def resultsversion():
    '''
    V1.0.0: plot the best-fit spectrum on top of the data
    '''
    return dawgie.VERSION(1, 0, 0)


# ------------------------------ -------------------------------------
def results(trgt, filt, fin, anc, xsl, atm, out, verbose=False):
    '''
    Plot out the results from atmos()
    trgt [INPUT]: target name
    filt [INPUT]: filter
    fin [INPUT]: system.finalize.parameters
    xls [INPUT]: cerberus.xslib.data
    atm [INPUT]: cerberus.atmos.data
    out [INPUT/OUTPUT]
    verbose [OPTIONAL]: verbosity
    '''
    ssc = syscore.ssconstants(mks=True)

    if verbose:
        print('cerberus/results for target:', trgt)

    completed_at_least_one_planet = False

    # load in the table of limits used for profiling
    if filt == 'HST-WFC3-IR-G141-SCAN':
        profiling_limits = get_profile_limits_hstg141()
    else:
        profiling_limits = []

    for p in fin['priors']['planets']:
        # print('post-analysis for planet:',p)

        # TEC params - X/H, C/O, N/O
        # disEq params - HCN, CH4, C2H2, CO2, H2CO

        # check whether this planet was analyzed
        # (some planets are skipped, because they have an unbound atmosphere)
        if p not in atm.keys():
            log.warning(
                '>-- CERBERUS.RESULTS: this planet is missing cerb fit: %s %s',
                trgt,
                p,
            )

        else:
            out['data'][p] = {}

            # limit results to just the TEC model?  No, not for HST/G141
            # but do verify that TEC exists at least
            if 'TEC' not in atm[p]['MODELPARNAMES'].keys():
                log.warning('>-- %s', 'TROUBLE: theres no TEC fit!?')
                return False

            # there was a bug before where PHOTOCHEM was passed in for Ariel
            # just in case, filter out models that are missing
            models = []
            for model_name in atm[p]['MODELPARNAMES']:
                if model_name in atm[p].keys():
                    models.append(model_name)
            for model_name in models:
                all_traces = []
                all_keys = []
                for key in atm[p][model_name]['MCTRACE']:
                    # print('going through keys in MCTRACE',key)
                    all_traces.append(atm[p][model_name]['MCTRACE'][key])
                    if model_name == 'TEC':
                        if key == 'TEC[0]':
                            all_keys.append('[X/H]')
                        elif key == 'TEC[1]':
                            all_keys.append('[C/O]')
                        elif key == 'TEC[2]':
                            all_keys.append('[N/O]')
                        else:
                            all_keys.append(key)
                    elif model_name == 'PHOTOCHEM':
                        if key == 'PHOTOCHEM[0]':
                            all_keys.append('HCN')
                        elif key == 'PHOTOCHEM[1]':
                            all_keys.append('CH4')
                        elif key == 'PHOTOCHEM[2]':
                            all_keys.append('C2H2')
                        elif key == 'PHOTOCHEM[3]':
                            all_keys.append('CO2')
                        elif key == 'PHOTOCHEM[4]':
                            all_keys.append('H2CO')
                        else:
                            all_keys.append(key)
                    else:
                        all_keys.append(key)
                # print('allKeys',allKeys)

                # remove the traced phase space that is excluded by profiling
                profile_trace, applied_limits = apply_profiling(
                    trgt + ' ' + p, profiling_limits, all_traces, all_keys
                )
                keepers = np.where(profile_trace == 1)
                profiled_traces = []
                for key in atm[p][model_name]['MCTRACE']:
                    profiled_traces.append(
                        atm[p][model_name]['MCTRACE'][key][keepers]
                    )
                profiled_traces = np.array(profiled_traces)

                # make note of the bounds placed on each parameter
                if 'prior_ranges' in atm[p][model_name].keys():
                    prior_ranges = atm[p][model_name]['prior_ranges']
                else:
                    prior_ranges = {}

                fit_cloud_parameters = 'CTP' in all_keys
                fit_n_to_o = '[N/O]' in all_keys
                fit_c_to_o = '[C/O]' in all_keys
                fit_t = 'T' in all_keys

                # save the relevant info
                transitdata = {}
                transitdata['wavelength'] = atm[p]['WAVELENGTH']
                transitdata['depth'] = atm[p]['SPECTRUM'] ** 2
                transitdata['error'] = 2 * atm[p]['SPECTRUM'] * atm[p]['ERRORS']

                truth_spectrum = None
                truth_params = None
                if 'sim' in filt:
                    if 'TRUTH_SPECTRUM' in atm[p].keys():
                        # print('NOT ERROR: true spectrum found in atmos output')
                        truth_spectrum = {
                            'depth': atm[p]['TRUTH_SPECTRUM'],
                            'wavelength': atm[p]['TRUTH_WAVELENGTH'],
                        }
                        truth_params = atm[p]['TRUTH_MODELPARAMS']
                    else:
                        print(
                            'ERROR: true spectrum is missing from the atmos output'
                        )
                elif 'TRUTH_SPECTRUM' in atm[p].keys():
                    print(
                        'ERROR: true spectrum is present for non-simulated data'
                    )

                tprtrace = atm[p][model_name]['MCTRACE']['T']
                tprtrace_profiled = atm[p][model_name]['MCTRACE']['T'][keepers]
                mdplist = [
                    key
                    for key in atm[p][model_name]['MCTRACE']
                    if model_name in key
                ]
                # print('mdplist',mdplist)
                mdptrace = []
                mdptrace_profiled = []
                for key in mdplist:
                    mdptrace.append(atm[p][model_name]['MCTRACE'][key])
                for key in mdplist:
                    mdptrace_profiled.append(
                        atm[p][model_name]['MCTRACE'][key][keepers]
                    )
                if fit_cloud_parameters:
                    ctptrace = atm[p][model_name]['MCTRACE']['CTP']
                    hzatrace = atm[p][model_name]['MCTRACE']['HScale']
                    hloctrace = atm[p][model_name]['MCTRACE']['HLoc']
                    hthicktrace = atm[p][model_name]['MCTRACE']['HThick']
                    ctp = np.median(ctptrace)
                    hza = np.median(hzatrace)
                    hloc = np.median(hloctrace)
                    hthc = np.median(hthicktrace)
                    # print('fit results; CTP:',ctp)
                    # print('fit results; HScale:',hza)
                    # print('fit results; HLoc:',hloc)
                    # print('fit results; HThick:',hthc)
                    ctptrace_profiled = atm[p][model_name]['MCTRACE']['CTP'][
                        keepers
                    ]
                    hzatrace_profiled = atm[p][model_name]['MCTRACE']['HScale'][
                        keepers
                    ]
                    hloctrace_profiled = atm[p][model_name]['MCTRACE']['HLoc'][
                        keepers
                    ]
                    hthicktrace_profiled = atm[p][model_name]['MCTRACE'][
                        'HThick'
                    ][keepers]
                    ctp_profiled = np.median(ctptrace_profiled)
                    hza_profiled = np.median(hzatrace_profiled)
                    hloc_profiled = np.median(hloctrace_profiled)
                    hthc_profiled = np.median(hthicktrace_profiled)
                else:
                    ctp = atm[p]['TRUTH_MODELPARAMS']['CTP']
                    hza = atm[p]['TRUTH_MODELPARAMS']['HScale']
                    hloc = atm[p]['TRUTH_MODELPARAMS']['HLoc']
                    hthc = atm[p]['TRUTH_MODELPARAMS']['HThick']
                    ctp_profiled = ctp
                    hza_profiled = hza
                    hloc_profiled = hloc
                    hthc_profiled = hthc
                    # print(' ctp hza hloc hthc',ctp,hza,hloc,hthc)
                if fit_t:
                    tpr = np.median(tprtrace)
                    tpr_profiled = np.median(tprtrace_profiled)
                else:
                    tpr = atm[p]['TRUTH_MODELPARAMS']['T']
                    tpr_profiled = tpr
                mdp = np.median(np.array(mdptrace), axis=1)
                mdp_profiled = np.median(np.array(mdptrace_profiled), axis=1)
                # print('fit results; T:',tpr)
                # print('fit results; mdplist:',mdp)

                solidr = fin['priors'][p]['rp'] * ssc['Rjup']

                if model_name == 'TEC':
                    # if len(mdp)!=3: log.warning('--< Expecting 3 molecules for TEQ model! >--')
                    mixratio = None
                    mixratio_profiled = None
                    tceqdict = {}
                    tceqdict_profiled = {}
                    tceqdict['XtoH'] = float(mdp[0])
                    tceqdict_profiled['XtoH'] = float(mdp_profiled[0])
                    if fit_c_to_o:
                        tceqdict['CtoO'] = float(mdp[1])
                        tceqdict_profiled['CtoO'] = float(mdp_profiled[1])
                    else:
                        if ('TRUTH_MODELPARAMS' in atm[p]) and (
                            'CtoO' in atm[p]['TRUTH_MODELPARAMS']
                        ):
                            # print('truth params',atm[p]['TRUTH_MODELPARAMS'])
                            tceqdict['CtoO'] = atm[p]['TRUTH_MODELPARAMS'][
                                'CtoO'
                            ]
                        else:
                            # default is C/O=1.  Maybe the default should actually be Solar?
                            tceqdict['CtoO'] = 0.0
                        tceqdict_profiled['CtoO'] = tceqdict['CtoO']

                    if fit_n_to_o:
                        tceqdict['NtoO'] = float(mdp[2])
                        tceqdict_profiled['NtoO'] = float(mdp_profiled[2])
                    else:
                        if ('TRUTH_MODELPARAMS' in atm[p]) and (
                            'NtoO' in atm[p]['TRUTH_MODELPARAMS']
                        ):
                            # print('truth params',atm[p]['TRUTH_MODELPARAMS'])
                            tceqdict['NtoO'] = atm[p]['TRUTH_MODELPARAMS'][
                                'NtoO'
                            ]
                        else:
                            tceqdict['NtoO'] = 0.0
                        tceqdict_profiled['NtoO'] = tceqdict['NtoO']
                elif model_name == 'PHOTOCHEM':
                    if len(mdp) != 5:
                        log.warning(
                            '--< Expecting 5 molecules for PHOTOCHEM model! >--'
                        )
                    tceqdict = None
                    mixratio = {}
                    mixratio['HCN'] = float(mdp[0])
                    mixratio['CH4'] = float(mdp[1])
                    mixratio['C2H2'] = float(mdp[2])
                    mixratio['CO2'] = float(mdp[3])
                    mixratio['H2CO'] = float(mdp[4])
                    tceqdict_profiled = None
                    mixratio_profiled = {}
                    mixratio_profiled['HCN'] = float(mdp_profiled[0])
                    mixratio_profiled['CH4'] = float(mdp_profiled[1])
                    mixratio_profiled['C2H2'] = float(mdp_profiled[2])
                    mixratio_profiled['CO2'] = float(mdp_profiled[3])
                    mixratio_profiled['H2CO'] = float(mdp_profiled[4])

                else:
                    log.warning('--< Expecting TEQ or PHOTOCHEM model! >--')

                crbhzlib = {'PROFILE': []}
                hazedir = os.path.join(
                    excalibur.context['data_dir'], 'CERBERUS/HAZE'
                )
                hazelib(crbhzlib, hazedir=hazedir, verbose=False)

                param_values_median = (
                    tpr,
                    ctp,
                    hza,
                    hloc,
                    hthc,
                    tceqdict,
                    mixratio,
                )
                param_values_profiled = (
                    tpr_profiled,
                    ctp_profiled,
                    hza_profiled,
                    hloc_profiled,
                    hthc_profiled,
                    tceqdict_profiled,
                    mixratio_profiled,
                )

                # print('median fmc',np.nanmedian(fmc))
                fmc = np.zeros(transitdata['depth'].size)
                fmc = crbmodel(
                    mixratio,
                    float(hza),
                    float(ctp),
                    solidr,
                    fin['priors'],
                    xsl[p]['XSECS'],
                    xsl[p]['QTGRID'],
                    float(tpr),
                    transitdata['wavelength'],
                    hzlib=crbhzlib,
                    hzp='AVERAGE',
                    hztop=float(hloc),
                    hzwscale=float(hthc),
                    cheq=tceqdict,
                    pnet=p,
                    verbose=False,
                    debug=False,
                )
                if isinstance(fmc, tensor.variable.TensorVariable):
                    fmc = fmc.eval()  # convert tensor to numpy array
                # print('median fmc',np.nanmedian(fmc))
                # print('mean model',np.nanmean(fmc))
                # print('mean data',np.nanmean(transitdata['depth']))
                patmos_model = (
                    fmc - np.nanmean(fmc) + np.nanmean(transitdata['depth'])
                )
                # print('median pmodel',np.nanmedian(patmos_model))

                fmc_profiled = np.zeros(transitdata['depth'].size)
                fmc_profiled = crbmodel(
                    mixratio_profiled,
                    float(hza_profiled),
                    float(ctp_profiled),
                    solidr,
                    fin['priors'],
                    xsl[p]['XSECS'],
                    xsl[p]['QTGRID'],
                    float(tpr_profiled),
                    transitdata['wavelength'],
                    hzlib=crbhzlib,
                    hzp='AVERAGE',
                    hztop=float(hloc_profiled),
                    hzwscale=float(hthc_profiled),
                    cheq=tceqdict_profiled,
                    pnet=p,
                    verbose=False,
                    debug=False,
                )
                # convert tensor to numpy array
                if isinstance(fmc_profiled, tensor.variable.TensorVariable):
                    fmc_profiled = fmc_profiled.eval()
                patmos_model_profiled = (
                    fmc_profiled
                    - np.nanmean(fmc_profiled)
                    + np.nanmean(transitdata['depth'])
                )

                # calculate chi2 values to see which is the best fit
                offsets_model = (
                    patmos_model - transitdata['depth']
                ) / transitdata['error']
                chi2model = np.nansum(offsets_model**2)
                if isinstance(chi2model, tensor.variable.TensorVariable):
                    print('chi2model', chi2model.eval(), 'TENSOR YES')
                else:
                    print('chi2model', chi2model, 'TENSOR NO')

                # actually the profiled chi2 isn't used below just now, so has to be commented out
                # offsets_modelProfiled = (patmos_modelProfiled - transitdata['depth']) / transitdata['error']
                # chi2modelProfiled = np.nansum(offsets_modelProfiled**2)
                # print('chi2 after profiling',chi2modelProfiled)

                # make an array of 10 random walker results
                nrandomwalkers = 100
                nrandomwalkers = 1000
                nrandomwalkers = 20

                # fix the random seed for each target/planet, so that results are reproducable
                int_from_target = (
                    1  # arbitrary initialization for the random seed
                )
                for char in trgt + ' ' + p:
                    int_from_target = (
                        123 * int_from_target + ord(char)
                    ) % 100000
                np.random.seed(int_from_target)

                chi2best = chi2model
                patmos_best_fit = patmos_model
                param_values_best_fit = param_values_profiled
                fmcarray = []
                for _ in range(nrandomwalkers):
                    iwalker = int(len(tprtrace) * np.random.rand())
                    # iwalker = max(0, len(tprtrace) - 1 - int(1000* np.random.rand()))
                    if fit_cloud_parameters:
                        ctp = ctptrace[iwalker]
                        hza = hzatrace[iwalker]
                        hloc = hloctrace[iwalker]
                        hthc = hthicktrace[iwalker]
                    if fit_t:
                        tpr = tprtrace[iwalker]
                    mdp = np.array(mdptrace)[:, iwalker]
                    # print('shape mdp',mdp.shape)
                    # if fitCloudParameters:
                    #    print('fit results; CTP:', ctp)
                    #    print('fit results; HScale:', hza)
                    #    print('fit results; HLoc:', hloc)
                    #    print('fit results; HThick:', hthc)
                    # print('fit results; T:', tpr)
                    # print('fit results; mdplist:', mdp)

                    if model_name == 'TEC':
                        mixratio = None
                        tceqdict = {}
                        tceqdict['XtoH'] = float(mdp[0])
                        if fit_c_to_o:
                            tceqdict['CtoO'] = float(mdp[1])
                        else:
                            tceqdict['CtoO'] = atm[p]['TRUTH_MODELPARAMS'][
                                'CtoO'
                            ]
                        if fit_n_to_o:
                            tceqdict['NtoO'] = float(mdp[2])
                        else:
                            if ('TRUTH_MODELPARAMS' in atm[p]) and (
                                'NtoO' in atm[p]['TRUTH_MODELPARAMS']
                            ):
                                tceqdict['NtoO'] = atm[p]['TRUTH_MODELPARAMS'][
                                    'NtoO'
                                ]
                            else:
                                # log.warning('--< NtoO is missing from TRUTH_MODELPARAMS >--')
                                tceqdict['NtoO'] = 0.0

                    elif model_name == 'PHOTOCHEM':
                        tceqdict = None
                        mixratio = {}
                        mixratio['HCN'] = float(mdp[0])
                        mixratio['CH4'] = float(mdp[1])
                        mixratio['C2H2'] = float(mdp[2])
                        mixratio['CO2'] = float(mdp[3])
                        mixratio['H2CO'] = float(mdp[4])

                    fmcrand = np.zeros(transitdata['depth'].size)
                    fmcrand = crbmodel(
                        mixratio,
                        float(hza),
                        float(ctp),
                        solidr,
                        fin['priors'],
                        xsl[p]['XSECS'],
                        xsl[p]['QTGRID'],
                        float(tpr),
                        transitdata['wavelength'],
                        hzlib=crbhzlib,
                        hzp='AVERAGE',
                        hztop=float(hloc),
                        hzwscale=float(hthc),
                        cheq=tceqdict,
                        pnet=p,
                        verbose=False,
                        debug=False,
                    )
                    if isinstance(fmcrand, tensor.variable.TensorVariable):
                        fmcrand = fmcrand.eval()  # convert tensor to array

                    # print('len',len(fmcrand))
                    # print('median fmc', np.nanmedian(fmcrand))
                    # print('mean model', np.nanmean(fmcrand))
                    # print('stdev model', np.nanstd(fmcrand))
                    fmcarray.append(fmcrand)

                    # check to see if this model is the best one
                    patmos_modelrand = (
                        fmcrand
                        - np.nanmean(fmcrand)
                        + np.nanmean(transitdata['depth'])
                    )
                    offsets_modelrand = (
                        patmos_modelrand - transitdata['depth']
                    ) / transitdata['error']
                    chi2modelrand = np.nansum(offsets_modelrand**2)
                    # chi2modelrand = tensor.sum(offsets_modelrand**2)
                    # print('chi2 for a random walker', chi2modelrand)
                    if isinstance(
                        chi2modelrand, tensor.variable.TensorVariable
                    ):
                        print(
                            'chi2modelrand', chi2modelrand.eval(), 'TENSOR YES'
                        )
                    else:
                        print('chi2modelrand', chi2modelrand, 'TENSOR NO')
                    if isinstance(chi2best, tensor.variable.TensorVariable):
                        print('chi2best', chi2best.eval(), 'TENSOR YES')
                    else:
                        print('chi2best', chi2best, 'TENSOR NO')
                    # if chi2modelrand.eval() < chi2best.eval():
                    if chi2modelrand < chi2best:
                        # print('  using this as best', chi2modelrand)
                        chi2best = chi2modelrand
                        patmos_best_fit = patmos_modelrand
                        param_values_best_fit = (
                            tpr,
                            ctp,
                            hza,
                            hloc,
                            hthc,
                            tceqdict,
                            mixratio,
                        )

                # _______________MAKE SOME PLOTS________________
                save_dir = os.path.join(
                    excalibur.context['data_dir'], 'bryden/'
                )

                # _______________BEST-FIT SPECTRUM PLOT________________
                transitdata = rebin_data(transitdata)

                out['data'][p]['plot_spectrum_' + model_name], _ = (
                    plot_spectrumfit(
                        transitdata,
                        patmos_model,
                        patmos_model_profiled,
                        patmos_best_fit,
                        fmcarray,
                        truth_spectrum,
                        fin['priors'],
                        anc['data'][p],
                        atm[p],
                        filt,
                        model_name,
                        trgt,
                        p,
                        save_dir,
                    )
                )

                if verbose:
                    print('paramValues median  ', param_values_median)
                    print('paramValues profiled', param_values_profiled)
                    print('paramValues bestFit ', param_values_best_fit)

                # _______________CORNER PLOT________________
                out['data'][p]['plot_corner_' + model_name], _ = plot_corner(
                    all_keys,
                    all_traces,
                    profiled_traces,
                    param_values_best_fit,
                    truth_params,
                    prior_ranges,
                    filt,
                    model_name,
                    trgt,
                    p,
                    save_dir,
                )

                # _______________WALKER-EVOLUTION PLOT________________
                out['data'][p]['plot_walkerevol_' + model_name], _ = (
                    plot_walker_evolution(
                        all_keys,
                        all_traces,
                        profiled_traces,
                        truth_params,
                        prior_ranges,
                        applied_limits,
                        filt,
                        model_name,
                        trgt,
                        p,
                        save_dir,
                    )
                )

                # _______________VS-PRIOR PLOT________________
                out['data'][p]['plot_vsprior_' + model_name], _ = plot_vs_prior(
                    all_keys,
                    all_traces,
                    profiled_traces,
                    truth_params,
                    prior_ranges,
                    applied_limits,
                    filt,
                    model_name,
                    trgt,
                    p,
                    save_dir,
                )

            out['target'].append(trgt)
            out['planets'].append(p)
            completed_at_least_one_planet = True
            # print('out-data keys at end of this planet',out['data'][p].keys())

    if completed_at_least_one_planet:
        out['STATUS'].append(True)
    return out['STATUS'][-1]


# --------------------------------------------------------------------
def analysis(aspects, filt, out, verbose=False):
    '''
    Plot out the population analysis (retrieval vs truth, mass-metallicity, etc)
    aspects: cross-target information
    out [INPUT/OUTPUT]
    verbose [OPTIONAL]: verbosity
    '''
    if verbose:
        print('cerberus/analysis...')

    aspecttargets = []
    for a in aspects:
        aspecttargets.append(a)
    log.warning(
        '--< CERBERUS ANALYSIS: NUMBER OF TARGETS IN ASPECT %s >--',
        len(aspecttargets),
    )

    svname = 'cerberus.atmos'

    alltargetlists = get_target_lists()

    # set prior_ranges to avoid possible used-before-assignment problem
    # (ideally it is read in, but possibly not if there's mistake/old formatting)
    # the normal call doesn't work well here actually. and it creates nodes
    # darn.  have to just set something arbitrary
    # _, prior_ranges = addPriors(priorRangeTable, runtime_params, model, modparlbl[model])
    prior_ranges = None

    # allow for analysis of multiple target lists
    analysistargetlists = []
    # optionally specify the specific planets within multi-planet systems
    analysisplanetlist = []

    if filt == 'Ariel-sim':
        # analysistargetlists.append({
        #    'targetlistname':'Roudier+ 2022',
        #    'targets':alltargetlists['roudier62']})
        # analysistargetlists.append({
        #    'targetlistname':'MCS Nov.2023 Transit-list',
        #    'targets':alltargetlists['arielMCS_Nov2023_transit']})
        # analysistargetlists.append({
        #    'targetlistname':'MCS Nov.2023 max-visits=25',
        #    'targets':alltargetlists['arielMCS_Nov2023_maxVisits25']})
        # analysistargetlists.append({
        #    'targetlistname':'MCS Feb.2024 Transit-list',
        #    'targets':alltargetlists['arielMCS_Feb2024_transit']})
        # analysistargetlists.append({
        #    'targetlistname':'MCS Feb.2024 max-visits=25',
        #    'targets':alltargetlists['arielMCS_Feb2024_maxVisits25']})
        # analysistargetlists.append({
        #    'targetlistname':'2-year science time; Thorngren mmw (Aug.2024)',
        #    'targets':alltargetlists['ariel_Aug2024_2years']})
        #
        #  *** Tier-2 (259 planets) ***
        # analysistargetlists.append({
        #   'targetlistname':'2-year science time (Tier-2); Thorngren mmw (Nov.2024)',
        #   'targets':alltargetlists['ariel_Nov2024_2years']})
        # analysisplanetlist = {
        #    'planetlistname':'2-year science time (Tier-2); Thorngren mmw (Nov.2024)',
        #    'planets':alltargetlists['ariel_Nov2024_2years_withPlanetletters']}
        #  *** Tier-1 (626 planets) ***
        analysistargetlists.append(
            {
                'targetlistname': '2-year science time (Tier-1); Thorngren mmw (Nov.2024)',
                'targets': alltargetlists['ariel_Nov2024_2yearsTier1'],
            }
        )
        analysisplanetlist = {
            'planetlistname': '2-year science time (Tier-1); Thorngren mmw (Aug.2024)',
            'planets': alltargetlists[
                'ariel_Nov2024_2yearsTier1_withPlanetletters'
            ],
        }
    else:
        analysistargetlists.append(
            {
                'targetlistname': 'All Excalibur targets',
                'targets': alltargetlists['active'],
            }
        )
        # analysistargetlists.append({
        #    'targetlistname':'Roudier+ 2022',
        #    'targets':alltargetlists['roudier62']})
        # analysistargetlists.append({
        #    'targetlistname':'All G141 targets',
        #    'targets':alltargetlists['G141']})

    for targetlist in analysistargetlists:
        # print('  running targetlist=',targetlist['targetlistname'])
        param_names = []
        masses = []
        stellar_fehs = []
        truth_values = defaultdict(list)
        fit_values = defaultdict(list)
        fit_errors = defaultdict(list)
        fit_errors2sided = defaultdict(list)

        # FIXMEE: move to config file and fix this code
        for trgt in targetlist['targets']:
            # print('        cycling through targets',trgt)
            if trgt not in aspecttargets:
                log.warning(
                    '--< CERBERUS ANALYSIS: TARGET NOT IN ASPECT %s %s >--',
                    filt,
                    trgt,
                )
            elif svname + '.' + filt not in aspects[trgt]:
                # some targets don't have this filter; no problem
                # log.warning('--< NO CERB.ATMOS for this FILTER+TARGET %s %s >--',filt,trgt)
                pass
            elif 'STATUS' not in aspects[trgt][svname + '.' + filt]:
                log.warning(
                    '--< CERBERUS ANALYSIS: FORMAT ERROR - NO STATUS %s %s >--',
                    filt,
                    trgt,
                )
            else:
                # print('target with valid data format for this filter:',filt,trgt)
                atmos_fit = aspects[trgt][svname + '.' + filt]

                # if 'stellar_params' in atmosFit['data']:  # strange. this doesn't work
                if 'stellar_params' in atmos_fit['data'].keys():
                    stellar_feh = atmos_fit['data']['stellar_params']['FEH*']
                else:
                    stellar_feh = 0
                    log.warning(
                        '--< CERBERUS ANALYSIS: no FEH* for %s >--', trgt
                    )

                # verify SV succeeded for target
                if not atmos_fit['STATUS'][-1]:
                    log.warning(
                        '--< CERBERUS ANALYSIS: STATUS IS FALSE FOR CERB.ATMOS %s %s >--',
                        filt,
                        trgt,
                    )
                else:
                    for planet_letter in atmos_fit['data'].keys():
                        # print(trgt,atmosFit['data'][planet_letter]['MODELPARNAMES'])
                        # print(trgt,atmosFit['data'][planet_letter]['planet_params'])

                        # print('   keys:',atmosFit['data'][planet_letter].keys())
                        if (
                            planet_letter == 'stellar_params'
                        ):  # this is not a planet letter
                            pass

                        elif (
                            analysisplanetlist
                            and trgt + ' ' + planet_letter
                            not in analysisplanetlist['planets']
                        ):
                            # print(' DROP: Ariel doesnt observe this planet',trgt+' '+planet_letter)
                            pass

                        elif (
                            'TEC'
                            not in atmos_fit['data'][planet_letter][
                                'MODELPARNAMES'
                            ]
                        ):
                            log.warning(
                                '--< CERBERUS ANALYSIS: BIG PROBLEM theres no TEC model! %s %s >--',
                                filt,
                                trgt,
                            )
                        elif (
                            'prior_ranges'
                            not in atmos_fit['data'][planet_letter]['TEC']
                        ):
                            log.warning(
                                '--< CERBERUS ANALYSIS: SKIP (no prior info) - %s %s >--',
                                filt,
                                trgt,
                            )
                        else:
                            if (
                                'planet_params'
                                in atmos_fit['data'][planet_letter]
                            ):
                                masses.append(
                                    atmos_fit['data'][planet_letter][
                                        'planet_params'
                                    ]['mass']
                                )
                            else:
                                masses.append(666)

                            stellar_fehs.append(stellar_feh)

                            # (prior range should be the same for all the targets)
                            prior_ranges = atmos_fit['data'][planet_letter][
                                'TEC'
                            ]['prior_ranges']

                            all_traces = []
                            all_keys = []
                            for key in atmos_fit['data'][planet_letter]['TEC'][
                                'MCTRACE'
                            ]:
                                all_traces.append(
                                    atmos_fit['data'][planet_letter]['TEC'][
                                        'MCTRACE'
                                    ][key]
                                )

                                if key == 'TEC[0]':
                                    all_keys.append('[X/H]')
                                elif key == 'TEC[1]':
                                    all_keys.append('[C/O]')
                                elif key == 'TEC[2]':
                                    all_keys.append('[N/O]')
                                else:
                                    all_keys.append(key)

                            for key, trace in zip(all_keys, all_traces):
                                if key not in param_names:
                                    param_names.append(key)
                                med = np.median(trace)
                                fit_values[key].append(med)
                                lo = np.percentile(np.array(trace), 16)
                                hi = np.percentile(np.array(trace), 84)
                                fit_errors[key].append((hi - lo) / 2)
                                fit_errors2sided[key].append(
                                    [med - lo, hi - med]
                                )
                                if verbose:
                                    if key == '[N/O]' and (hi - lo) / 2 < 2:
                                        print(
                                            'N/O',
                                            trgt,
                                            np.median(trace),
                                            (hi - lo) / 2,
                                        )
                            if (
                                'TRUTH_MODELPARAMS'
                                in atmos_fit['data'][planet_letter].keys()
                            ) and (
                                isinstance(
                                    atmos_fit['data'][planet_letter][
                                        'TRUTH_MODELPARAMS'
                                    ],
                                    dict,
                                )
                            ):
                                truth_params = atmos_fit['data'][planet_letter][
                                    'TRUTH_MODELPARAMS'
                                ].keys()
                                # print('truth keys:',truth_params)
                            else:
                                truth_params = []

                            for trueparam, fitparam in zip(
                                ['Teq', 'metallicity', 'C/O', 'N/O', 'Mp'],
                                ['T', '[X/H]', '[C/O]', '[N/O]', 'Mp'],
                            ):
                                if trueparam in truth_params:
                                    true_value = atmos_fit['data'][
                                        planet_letter
                                    ]['TRUTH_MODELPARAMS'][trueparam]
                                    # (metallicity and C/O do not have to be converted to log-solar)
                                    # if trueparam=='metallicity':
                                    #    true_value = np.log10(true_value)
                                    # elif trueparam=='C/O':
                                    #    true_value = np.log10(true_value/0.54951)  # solar is C/O=0.55
                                    # elif trueparam=='N/O':
                                    #     true_value = true_value
                                    if (
                                        fitparam == '[N/O]'
                                        and true_value == 666
                                    ):
                                        truth_values[fitparam].append(0)
                                    else:
                                        truth_values[fitparam].append(
                                            true_value
                                        )

                                    if verbose:
                                        if (
                                            trueparam == 'Teq'
                                            and true_value > 3333
                                        ):
                                            print(
                                                'strangely high T',
                                                trgt,
                                                true_value,
                                            )
                                        if (
                                            trueparam == 'metallicity'
                                            and true_value > 66
                                        ):
                                            print(
                                                'strangely high [X/H]',
                                                trgt,
                                                true_value,
                                            )
                                            print(
                                                'atmosFit',
                                                atmos_fit['data'][
                                                    planet_letter
                                                ],
                                            )
                                        if (
                                            trueparam == 'C/O'
                                            and true_value > 0.5
                                        ):
                                            print(
                                                'strangely high [C/O]',
                                                trgt,
                                                true_value,
                                            )

                                elif trueparam == 'Mp':
                                    # if the planet mass is not in the Truth dictionary, pull it from system
                                    # print(' input keys',atmosFit['data'][planet_letter]['planet_params'])
                                    # print(' planet mass from system params:',
                                    #      atmosFit['data'][planet_letter]['planet_params']['mass'])
                                    truth_values[fitparam].append(
                                        atmos_fit['data'][planet_letter][
                                            'planet_params'
                                        ]['mass']
                                    )
                                elif trueparam == 'N/O':
                                    truth_values[fitparam].append(0)
                                else:
                                    truth_values[fitparam].append(666)

        # plot analysis of the results.  save as png and as state vector for states/view
        save_dir = os.path.join(excalibur.context['data_dir'], 'bryden/')
        fit_co_plot = False
        fit_no_plot = False
        if 'sim' in filt:
            # for simulated data, compare retrieval against the truth
            #  note that the length of plotarray depends on whether N/O and C/O are fit parameters
            # jenkins doesn't like to have a triple-packed return here because it's fussy
            plotarray = plot_fits_vs_truths(
                truth_values,
                fit_values,
                fit_errors,
                prior_ranges,
                filt,
                save_dir,
            )
            # fitTplot, fitMetalplot, fitCOplot, fitNOplot = plotarray[0],plotarray[1],plotarray[2],plotarray[3]
            fit_t_plot = plotarray[0]
            fit_metalplot = plotarray[1]
            if len(plotarray) > 2:
                fit_co_plot = plotarray[2]
            if len(plotarray) > 3:
                fit_no_plot = plotarray[3]
        else:
            # for real data, make a histogram of the retrieved uncertainties
            #  note that the length of plotarray depends on whether N/O and C/O are fit parameters
            plotarray = plot_fit_uncertainties(
                fit_values, fit_errors, prior_ranges, filt, save_dir
            )
            fit_t_plot = plotarray[0]
            fit_metalplot = plotarray[1]
            if len(plotarray) > 2:
                fit_co_plot = plotarray[2]
            if len(plotarray) > 3:
                fit_no_plot = plotarray[3]

        mass_metals_plot, _ = plot_mass_vs_metals(
            truth_values['Mp'],
            stellar_fehs,
            truth_values,
            fit_values,
            fit_errors2sided,
            prior_ranges,
            filt,
            save_dir,
        )

        # save the analysis as .csv file? (in /proj/data/spreadsheets/)
        # savesv(aspects, targetlist)

        # targetlistname = targetlist['targetlistname']

        # Add to SV
        out['data']['truths'] = dict(truth_values)
        out['data']['values'] = dict(fit_values)
        out['data']['errors'] = dict(fit_errors)
        out['data']['plot_mass_v_metals'] = mass_metals_plot
        out['data']['plot_fitT'] = fit_t_plot
        out['data']['plot_fitMetal'] = fit_metalplot
        if fit_co_plot:
            out['data']['plot_fitCO'] = fit_co_plot
        if fit_no_plot:
            out['data']['plot_fitNO'] = fit_no_plot

    out['data']['params'] = param_names
    out['data']['targetlistnames'] = [
        targetlist['targetlistname'] for targetlist in analysistargetlists
    ]

    out['STATUS'].append(True)
    return out['STATUS'][-1]


# ---------------------------------- ---------------------------------
# -- ROUDIER ET AL. 2021 RELEASE -- ----------------------------------
def rlsversion():
    '''
    GMR:110 Initial release to IPAC
    GMR:111 Removed empty keys
    '''
    return dawgie.VERSION(1, 1, 1)


def release(trgt, fin, out, verbose=False):
    '''
    GMR: Format Cerberus SV products to be released to IPAC
    trgt [INPUT]: target name
    fin [INPUT]: system.finalize.parameters
    out [INPUT/OUTPUT]
    ext [INPUT]: 'HST-WFC3-IR-G141-SCAN'
    verbose [OPTIONAL]: verbosity
    '''
    print('target name in cerb.release', trgt)
    rlsed = False
    plist = fin['priors']['planets']
    thispath = os.path.join(excalibur.context['data_dir'], 'CERBERUS')
    print('thispath', thispath)
    for p in plist:
        intxtf = os.path.join(thispath, 'P.CERBERUS.atmos', trgt + p + '.txt')
        incorrpng = os.path.join(
            thispath, 'P.CERBERUS.atmos', trgt + p + '_atmos_corr.png'
        )
        intxtpng = os.path.join(
            thispath, 'P.CERBERUS.atmos', trgt + p + '_atmos.png'
        )
        out['data'][p] = {}
        try:
            atm = np.loadtxt(intxtf)
            out['data'][p]['atmos'] = atm
            out['STATUS'].append(True)
            pass
        except FileNotFoundError:
            pass
        try:
            corrplot = img.imread(incorrpng)
            out['data'][p]['corrplot'] = corrplot
            out['STATUS'].append(True)
            pass
        except FileNotFoundError:
            pass
        try:
            modelplot = img.imread(intxtpng)
            out['data'][p]['modelplot'] = modelplot
            out['STATUS'].append(True)
            pass
        except FileNotFoundError:
            pass
        if not out['data'][p]:
            if verbose:
                log.warning('--< No data found for %s', p)
            out['data'].pop(p)
            pass
        pass
    rlsed = out['STATUS'][-1]
    if verbose:
        log.warning('--< %s', out['STATUS'])
    return rlsed
