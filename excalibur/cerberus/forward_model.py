'''cerberus forward_model ds'''

# Heritage code shame:
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments,too-many-branches,too-many-lines,too-many-locals,too-many-positional-arguments,too-many-statements

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as cst
from scipy.interpolate import interp1d as itp
import logging

import excalibur.system.core as syscore
from excalibur.util.cerberus import crbce, calcTEA, getmmw

from excalibur.cerberus.fmcontext import ctxtinit


log = logging.getLogger(__name__)

# this doesn't change results at all; just needed to avoid undefined-variable pylint
ctxt = ctxtinit()


# ----------- --------------------------------------------------------
# -- CERBERUS MODEL -- -----------------------------------------------
def crbmodel(
    temp,
    cloudtp,
    cheq=None,
    mixratio=None,
    hazescale=0.0,
    hazethick=1.0,
    hazeslope=-4.0,
    hazeloc=None,
    hazeprof='AVERAGE',
    hzlib=None,
    chemistry='TEC',
    planet=None,
    rp0=None,
    orbp=None,
    wgrid=None,
    xsecs=None,
    qtgrid=None,
    isothermal=None,
    lbroadening=None,
    lshifting=None,
    knownspecies=None,
    cialist=None,
    xmollist=None,
    nlevels=None,
    Hsmax=None,
    solrad=None,
    break_down_by_molecule=False,
    logx=False,
    verbose=False,
    debug=False,
):
    '''
    G. ROUDIER: Cerberus forward model probing up to 'Hsmax' scale heights from solid
    radius solrad evenly log divided amongst nlevels steps
    '''

    if planet is None:
        planet = ctxt.planet
    if orbp is None:
        orbp = ctxt.orbp
    if knownspecies is None:
        knownspecies = ctxt.knownspecies
    if cialist is None:
        cialist = ctxt.cialist
    else:
        cialist = ['H2-H', 'H2-H2', 'H2-He', 'He-H']
    if xmollist is None:
        xmollist = ctxt.xmollist
    else:
        xmollist = ['TIO', 'H2O', 'H2CO', 'HCN', 'CO', 'CO2', 'NH3', 'CH4']
    if nlevels is None:
        nlevels = ctxt.nlevels
    if Hsmax is None:
        Hsmax = ctxt.Hsmax
    if solrad is None:
        solrad = ctxt.solrad
    if not bool(lshifting):
        lshifting = ctxt.lshifting
    if not bool(lbroadening):
        lbroadening = ctxt.lbroadening
    if not bool(isothermal):
        isothermal = ctxt.isothermal
    if rp0 is None:
        rp0 = ctxt.rp0
    if xsecs is None:
        xsecs = ctxt.xsl['data'][ctxt.planet]['XSECS']
    if qtgrid is None:
        qtgrid = ctxt.xsl['data'][ctxt.planet]['QTGRID']
    if wgrid is None:
        wgrid = np.array(ctxt.spc['data'][ctxt.planet]['WB'])
    if hzlib is None:
        hzlib = ctxt.hzlib

    ssc = syscore.ssconstants(mks=True)
    pgrid = np.arange(
        np.log(solrad) - Hsmax,
        np.log(solrad) + Hsmax / nlevels,
        Hsmax / (nlevels - 1),
    )
    pgrid = np.exp(pgrid)
    pressure = pgrid[::-1]
    dPoverP = (pressure[1] - pressure[0]) / pressure[0]

    # print('PARAMETERS', temp, cheq['CtoO'], cheq['XtoH'])
    if not mixratio:
        fH2 = None
        fHe = None
        if cheq is None:
            log.error('neither mixratio nor cheq are defined')
        if chemistry == 'TEC':
            mixratio, fH2, fHe = crbce(
                pressure,
                temp,
                C2Or=cheq['CtoO'],
                X2Hr=cheq['XtoH'],
                N2Or=cheq['NtoO'],
            )
        elif chemistry == 'TEA':
            # 'temp' should actually be the 10 T-P profile coefficients!!
            tempCoeffs = [0, temperature, 0, 0, 0, 0, 0, 0, 0, 0]
            mixratio, fH2, fHe = calcTEA(
                tempCoeffs,
                pressure,
                metallicity=10.0**cheq['XtoH'],
                C_O=0.55 * 10.0**cheq['CtoO'],
                # N_O=?? * 10.0**cheq['NtoO'],
            )
        else:
            mixratio = {'H2O': 6}
            log.error('--< UNKNOWN CHEM MODEL: %s >--', chemistry)
        # print('mixratio',mixratio,fH2,fHe)
        mmw, fH2, fHe = getmmw(mixratio, protosolar=False, fH2=fH2, fHe=fHe)
    else:
        mmw, fH2, fHe = getmmw(mixratio)
    mmw = mmw * cst.m_p  # [kg]

    if debug:
        print('mmw', mmw * 6.022e26)

    Hs = (
        cst.Boltzmann
        * temp
        / (mmw * 1e-2 * (10.0 ** float(orbp[planet]['logg'])))
    )  # [m]

    # when the Pressure grid is log-spaced, rdz is a constant
    #  drop dz[] and dzprime[] arrays and just use this constant instead
    rdz = abs(Hs / 2.0 * np.log(1.0 + dPoverP))
    dz = 2 * rdz
    z = dz * np.linspace(0, len(pressure) - 1, len(pressure))

    rho = pressure * 1e5 / (cst.Boltzmann * temp)
    tau, tau_by_molecule, wtau = gettau(
        xsecs,
        qtgrid,
        temp,
        mixratio,
        z,
        dz,
        rho,
        rp0,
        pressure,
        wgrid,
        lbroadening,
        lshifting,
        cialist,
        fH2,
        fHe,
        xmollist,
        hazescale,
        hzlib,
        hazeprof,
        hazeslope,
        hazeloc,
        isothermal,
        hazethick,
        debug=debug,
    )

    # print('tau', tau)
    # for molecule, tau in tau_by_molecule.items:
    #    print('tau', molecule, tau)

    if not break_down_by_molecule:
        tau_by_molecule = {}
    molecules = tau_by_molecule.keys()
    # SEMI FINITE CLOUD ------------------------------------------------------------------
    reversep = np.array(pressure[::-1])
    selectcloud = pressure > 10.0**cloudtp
    blocked = False
    if np.all(selectcloud):
        tau = tau * 0
        for molecule in molecules:
            tau_by_molecule[molecule] = tau_by_molecule[molecule] * 0
        blocked = True
        pass
    if not np.all(~selectcloud) and not blocked:
        # 1) find the cloudtop location in the pressure grid
        cloudtopindex = np.max(np.arange(len(pressure))[selectcloud]) + 1

        # 2) set atmos depth to zero for all cells deeper than that
        tau[:cloudtopindex, :] = 0.0
        for molecule in molecules:
            tau_by_molecule[molecule][:cloudtopindex, :] = 0.0

        # 3) interpolate atmos depth within that cell
        for waveindex in np.arange(wtau.size):
            myspl = itp(reversep, np.asarray(tau[:, waveindex]).flatten())
            tau[cloudtopindex, waveindex] = myspl(10.0**cloudtp)
            # this line can be done for all indices (outside of loop, I mean)
            # tau[:cloudtopindex, waveindex] = 0.
            for molecule in molecules:
                myspl = itp(
                    reversep,
                    np.asarray(
                        tau_by_molecule[molecule][:, waveindex]
                    ).flatten(),
                )
                tau_by_molecule[molecule][cloudtopindex, waveindex] = myspl(
                    10.0**cloudtp
                )
                # tau_by_molecule[molecule][:cloudtopindex, waveindex] = 0.
            pass

        # adjust rp0 based on the cloudtop
        ctpdpress = 10.0**cloudtp - pressure[cloudtopindex]
        ctpdz = abs(
            Hs / 2.0 * np.log(1.0 + ctpdpress / pressure[cloudtopindex])
        )
        rp0 += z[cloudtopindex] + ctpdz
    pass

    # note that original version had rp0+z as matrix then multiply by dz after
    geometrygrid = (rp0 + z) * dz
    absorptiongrid = 1.0 - np.exp(-tau)
    atmdepth = (
        2e0 * np.asmatrix(geometrygrid) * np.asmatrix(absorptiongrid)
    ).flatten()

    model = (rp0**2 + atmdepth) / (orbp['R*'] * ssc['Rsun']) ** 2
    # model is a 1xN matrix; it needs to be a 1-d array
    #  otherwise some subsequent * or ** operations fail
    model = np.asarray(model).reshape(-1)
    model = model[::-1]

    models_by_molecule = {}
    for molecule in molecules:
        absorptiongrid = 1.0 - np.exp(-tau_by_molecule[molecule])
        atmdepth = (
            2e0 * np.asmatrix(geometrygrid) * np.asmatrix(absorptiongrid)
        ).flatten()

        models_by_molecule[molecule] = (rp0**2 + atmdepth) / (
            orbp['R*'] * ssc['Rsun']
        ) ** 2

        # convert matrix to 1-d array
        models_by_molecule[molecule] = np.asarray(
            models_by_molecule[molecule]
        ).reshape(-1)
        models_by_molecule[molecule] = models_by_molecule[molecule][::-1]

    if verbose:
        plotmodel = model.copy()
        noatm = np.nanmin(plotmodel)
        rp0hs = np.sqrt(noatm * (orbp['R*'] * ssc['Rsun']) ** 2)

        fig, ax = plt.subplots(figsize=(10, 6))
        axes = [ax, ax.twinx(), ax.twinx()]
        fig.subplots_adjust(left=0.125, right=0.775)
        axes[-1].spines['right'].set_position(('axes', 1.2))
        axes[-1].set_frame_on(True)
        axes[-1].patch.set_visible(False)
        axes[0].plot(wtau, 1e2 * plotmodel)
        axes[0].plot(wtau, plotmodel * 0 + 1e2 * noatm, '--')
        axes[0].set_xlabel('Wavelength $\\lambda$[$\\mu m$]')
        axes[0].set_ylabel('Transit Depth [%]')
        axes[0].get_yaxis().get_major_formatter().set_useOffset(False)
        yaxmin, yaxmax = axes[0].get_ylim()
        ax2min = (
            np.sqrt(1e-2 * yaxmin) * orbp['R*'] * ssc['Rsun'] - rp0hs
        ) / Hs
        ax2max = (
            np.sqrt(1e-2 * yaxmax) * orbp['R*'] * ssc['Rsun'] - rp0hs
        ) / Hs
        axes[-1].set_ylabel('Transit Depth Modulation [Hs]')
        axes[-1].set_ylim(ax2min, ax2max)
        axes[-1].get_yaxis().get_major_formatter().set_useOffset(False)
        axes[1].set_ylabel('Transit Depth Modulation [ppm]')
        axes[1].set_ylim(
            1e6 * (1e-2 * yaxmin - noatm), 1e6 * (1e-2 * yaxmax - noatm)
        )
        axes[1].get_yaxis().get_major_formatter().set_useOffset(False)
        if logx:
            plt.semilogx()
            plt.xlim([np.min(wtau), np.max(wtau)])
            pass
        plt.show()
        pass

    # print('crbmodel: model at end',model)

    if break_down_by_molecule:
        return model, models_by_molecule
    return model


# --------------------------- ----------------------------------------
# -- TAU -- ----------------------------------------------------------
def gettau(
    xsecs,
    qtgrid,
    temp,
    mixratio,
    z,
    dz,
    rho,
    rp0,
    pressure,
    wgrid,
    lbroadening,
    lshifting,
    cialist,
    fH2,
    fHe,
    xmollist,
    hazescale,
    hzlib,
    hazeprof,
    hazeslope,
    hazeloc,
    isothermal,
    hazethick,
    debug=False,
):
    '''
    G. ROUDIER: Builds optical depth matrix
    '''

    # SPHERICAL SHELL (PLANE-PARALLEL REMOVED) -------------------------------------------
    # MATRICES INIT ------------------------------------------------------------------
    Nzones = len(pressure)
    tau = np.zeros((Nzones, wgrid.size))
    # print('tau shape at the top', tau.shape)
    tau_by_molecule = {}
    # DL ARRAY, Z VERSUS ZPRIME ------------------------------------------------------
    zprime = np.broadcast_to(z, (Nzones, Nzones))
    thisz = zprime.T
    # print('zprime',zprime)
    # print('thisz',thisz)

    dl = np.sqrt(
        np.max(
            ([zprime * 0, (rp0 + zprime + dz) ** 2 - (rp0 + thisz) ** 2]),
            axis=0,
        )
    )
    dl0 = np.sqrt(
        np.max(
            [zprime * 0, (rp0 + zprime) ** 2 - (rp0 + thisz) ** 2],
            axis=0,
        )
    )
    dlarray = dl - dl0
    # print('dl', dlarray)

    # GAS ARRAY, ZPRIME VERSUS WAVELENGTH  -------------------------------------------
    for elem in mixratio:
        # tau_by_molecule[elem] = np.zeros((len(z), wgrid.size))
        mmr = 10.0 ** (mixratio[elem] - 6.0)
        # Fake use of xmollist due to changes in xslib v112
        # THIS HAS TO BE FIXED
        # if elem not in xmollist:
        if not xmollist:
            # HITEMP/HITRAN ROTHMAN ET AL. 2010 --------------------------------------
            sigma, lsig = absorb(
                xsecs[elem],
                qtgrid[elem],
                temp,
                pressure,
                mmr,
                lbroadening,
                lshifting,
                wgrid,
            )
            # sigma = np.array(sigma)  # cm^2/mol
            if True in (sigma < 0):
                sigma[sigma < 0] = 0e0
            if True in ~np.isfinite(sigma):
                sigma[~np.isfinite(sigma)] = 0e0
            sigma = sigma * 1e-4  # m^2/mol
            pass
        else:
            # EXOMOL HILL ET AL. 2013 ------------------------------------------------
            sigma, lsig = getxmolxs(temp, xsecs[elem])
            # sigma = np.array(sigma)   # cm^2/mol
            if True in (sigma < 0):
                sigma[sigma < 0] = 0e0
            if True in ~np.isfinite(sigma):
                sigma[~np.isfinite(sigma)] = 0e0
            # sigma = np.array(sigma)*1e-4  # m^2/mol
            sigma = sigma * 1e-4  # m^2/mol
            pass
        if isothermal:
            tau = tau + mmr * sigma * np.array([rho]).T
            tau_by_molecule[elem] = mmr * sigma * np.array([rho]).T
        pass
    # CIA ARRAY, ZPRIME VERSUS WAVELENGTH  -------------------------------------------
    for cia in cialist:
        if cia == 'H2-H2':
            f1 = fH2
            f2 = fH2
        elif cia == 'H2-He':
            f1 = fH2
            f2 = fHe
        elif cia == 'H2-H':
            f1 = fH2
            f2 = fH2 * 2.0
        elif cia == 'He-H':
            f1 = fHe
            f2 = fH2 * 2.0
        else:
            log.warning(
                '--< CERBERUS gettau(): UNEXPECTED MOLECULE %s >--', cia
            )
            f1 = 0
            f2 = 0
        # HITRAN RICHARD ET AL. 2012
        sigma, lsig = getciaxs(temp, xsecs[cia])  # cm^5/mol^2
        sigma = np.array(sigma) * 1e-10  # m^5/mol^2
        if True in (sigma < 0):
            sigma[sigma < 0] = 0e0
        if True in ~np.isfinite(sigma):
            sigma[~np.isfinite(sigma)] = 0e0
        tau = tau + f1 * f2 * sigma * np.array([rho**2]).T
        tau_by_molecule[cia] = f1 * f2 * sigma * np.array([rho**2]).T
    # RAYLEIGH ARRAY, ZPRIME VERSUS WAVELENGTH  --------------------------------------
    # NAUS & UBACHS 2000
    slambda0 = 750.0 * 1e-3  # microns
    sray0 = 2.52 * 1e-28 * 1e-4  # m^2/mol
    sigma = sray0 * (wgrid[::-1] / slambda0) ** (-4)
    tau = tau + fH2 * sigma * np.array([rho]).T
    tau_by_molecule['rayleigh'] = fH2 * sigma * np.array([rho]).T
    # HAZE ARRAY, ZPRIME VERSUS WAVELENGTH  ------------------------------------------
    if hzlib is None:
        slambda0 = 750.0 * 1e-3  # microns
        sray0 = 2.52 * 1e-28 * 1e-4  # m^2/mol
        sigma = sray0 * (wgrid[::-1] / slambda0) ** (hazeslope)
        hazedensity = np.ones(len(z))
        tau = tau + 10.0**hazescale * sigma * np.array([hazedensity]).T
        tau_by_molecule['haze'] = (
            10.0**hazescale * sigma * np.array([hazedensity]).T
        )
    else:
        # WEST ET AL. 2004
        sigma = (
            0.0083
            * (wgrid[::-1]) ** (hazeslope)
            * (
                1e0
                + 0.014 * (wgrid[::-1]) ** (hazeslope / 2e0)
                + 0.00027 * (wgrid[::-1]) ** (hazeslope)
            )
        )
        if hazeprof in ['MAX', 'MEDIAN', 'AVERAGE']:
            frh = hzlib['PROFILE'][0][hazeprof][0]
            rh = frh(pressure)
            rh[rh < 0] = 0.0
            haze_ref_pressure = float(pressure[rh == np.max(rh)])
            if hazeloc is None:
                hazeshift = 0e0
            else:
                hazeshift = hazeloc - np.log10(haze_ref_pressure)
            splp = np.log10(pressure[::-1])
            splrh = rh[::-1]
            thisfrh = itp(
                splp, splrh, kind='linear', bounds_error=False, fill_value=0e0
            )
            hazewdist = hazeloc - np.log10(pressure)

            if hazethick > 0:
                preval = hazeloc - hazewdist / hazethick - hazeshift
                rh = thisfrh(preval)
                rh[rh < 0] = 0e0
            else:
                rh = thisfrh(np.log10(pressure)) * 0
            if debug:
                jptprofile = 'J' + hazeprof
                jdata = np.array(hzlib['PROFILE'][0][jptprofile])
                jpres = np.array(hzlib['PROFILE'][0]['PRESSURE'])
                myfig = plt.figure(figsize=(12, 6))
                plt.plot(
                    1e6 * jdata, jpres, color='blue', label='Lavvas et al. 2017'
                )
                plt.axhline(haze_ref_pressure, linestyle='--', color='blue')
                plt.plot(
                    1e6 * rh,
                    pressure,
                    'r',
                    label='Parametrized density profile',
                )
                plt.plot(
                    1e6 * thisfrh(np.log10(pressure) - hazeshift),
                    pressure,
                    'g^',
                )
                if hazeloc is not None:
                    plt.axhline(10**hazeloc, linestyle='--', color='red')
                    pass
                plt.semilogy()
                plt.semilogx()
                plt.gca().invert_yaxis()
                plt.xlim([1e-4, np.max(1e6 * rh)])
                plt.tick_params(axis='both', labelsize=20)
                plt.xlabel('Aerosol Density [$n.{cm}^{-3}$]', fontsize=24)
                plt.ylabel('Pressure [bar]', fontsize=24)
                plt.title('Aerosol density profile', fontsize=24)
                plt.legend(
                    loc='center left',
                    frameon=False,
                    fontsize=24,
                    bbox_to_anchor=(1, 1),
                )
                myfig.tight_layout()
                plt.show()
                pass
            pass
        else:
            rh = np.array(
                [np.nanmean(hzlib['PROFILE'][0]['CONSTANT'])] * len(z)
            )
            negrh = rh < 0e0
            if True in negrh:
                rh[negrh] = 0e0
            pass
        hazecontribution = 10.0**hazescale * sigma * np.array([rh]).T
        tau = tau + hazecontribution
        tau_by_molecule['haze'] = hazecontribution
        pass

    tau = 2e0 * np.asmatrix(dlarray) * np.asmatrix(tau)

    molecules = tau_by_molecule.keys()
    for molecule in molecules:
        # print(' MOLECULE:', molecule)
        tau_by_molecule[molecule] = (
            2e0 * np.asmatrix(dlarray) * np.asmatrix(tau_by_molecule[molecule])
        )

    if debug:
        plt.figure(figsize=(12, 6))
        plt.imshow(
            np.log10(tau),
            aspect='auto',
            origin='lower',
            extent=[
                max(wgrid),
                min(wgrid),
                np.log10(max(pressure)),
                np.log10(min(pressure)),
            ],
        )
        plt.ylabel('log10(Pressure)', fontsize=24)
        plt.xlabel('Wavelength [$\\mu m$]', fontsize=24)
        plt.gca().invert_xaxis()
        plt.title('log10(Optical Depth)', fontsize=24)
        plt.tick_params(axis='both', labelsize=20)
        cbar = plt.colorbar()
        cbar.ax.tick_params(labelsize=20)
        # plt.savefig('opticalDepth1.png')  # permission denied
        # plt.savefig('/proj/sdp/bryden/opticalDepth2.png')  # no such file/dir
        # hey this works! but no way to see it.
        #  also it seems to raise bandit security error in pylint
        # plt.savefig('/tmp/opticalDepth3.png')
        # plt.savefig('/home/bryden/opticalDepth4.png') # no such file/dir
        plt.show()
        pass
    return tau, tau_by_molecule, 1e4 / lsig


# --------- ----------------------------------------------------------
# -- ATTENUATION COEFFICIENT -- --------------------------------------
def absorb(
    xsecs,
    qtgrid,
    T,
    pressure,
    mmr,
    lbroadening,
    lshifting,
    wgrid,
    iso=0,
    Tref=296.0,
    debug=False,
):
    '''
    G. ROUDIER: HITRAN HITEMP database parser
    '''

    select = np.array(xsecs['I']) == (iso + 1)
    select = np.array(xsecs['I']) == iso + 1
    S = np.array(xsecs['S'])[select]
    E = np.array(xsecs['Epp'])[select]
    gself = np.array(xsecs['g_self'])[select]
    nu = np.array(xsecs['nu'])[select]
    delta = np.array(xsecs['delta'])[select]
    eta = np.array(xsecs['eta'])[select]
    gair = np.array(xsecs['g_air'])[select]

    Qref = float(qtgrid['SPL'][iso](Tref))

    try:
        Q = float(qtgrid['SPL'][iso](T))
    except ValueError:
        Q = np.nan
    c2 = 1e2 * cst.h * cst.c / cst.Boltzmann
    tips = (Qref * np.exp(-c2 * E / T) * (1.0 - np.exp(-c2 * nu / T))) / (
        Q * np.exp(-c2 * E / Tref) * (1.0 - np.exp(-c2 * nu / Tref))
    )
    if np.all(~np.isfinite(tips)):
        tips = 0
    sigma = S * tips
    ps = mmr * pressure
    gamma = np.array(
        np.asmatrix(pressure - ps).T * np.asmatrix(gair * (Tref / T) ** eta)
        + np.asmatrix(ps).T * np.asmatrix(gself)
    )
    if lbroadening:
        if lshifting:
            matnu = np.array(
                np.asmatrix(np.ones(pressure.size)).T * np.asmatrix(nu)
                + np.asmatrix(pressure).T * np.asmatrix(delta)
            )
        else:
            matnu = np.array(nu) * np.array([np.ones(len(pressure))]).T
        pass
    else:
        matnu = np.array(nu)
    absgrid = []
    nugrid = (1e4 / wgrid)[::-1]
    dwnu = np.concatenate((np.array([np.diff(nugrid)[0]]), np.diff(nugrid)))
    if lbroadening:
        for mymatnu, mygamma in zip(matnu, gamma):
            binsigma = np.asmatrix(sigma) * np.asmatrix(
                intflor(
                    nugrid,
                    dwnu / 2.0,
                    np.array([mymatnu]).T,
                    np.array([mygamma]).T,
                )
            )
            binsigma = np.array(binsigma).flatten()
            absgrid.append(binsigma / dwnu)
            pass
        pass
    else:
        binsigma = []
        for nubin, dw in zip(nugrid, dwnu):
            select = (matnu > (nubin - dw / 2.0)) & (matnu <= nubin + dw / 2.0)
            binsigma.append(np.sum(sigma[select]))
            pass
        binsigma = np.array(binsigma) / dwnu
        absgrid.append(binsigma)
        pass
    if debug:
        plt.semilogy(1e4 / matnu.T, sigma, '.')
        plt.semilogy(wgrid[::-1], binsigma, 'o')
        plt.xlabel('Wavelength $\\lambda$[$\\mu m$]')
        plt.ylabel('Absorption Coeff [$cm^{2}.molecule^{-1}$]')
        plt.show()
        pass
    return absgrid, nugrid


# --------- ------`<----------------------------------------------------
# -- EXOMOL -- -------------------------------------------------------
def getxmolxs(temp, xsecs):
    '''
    G. ROUDIER: Wrapper around EXOMOL Cerberus library
    '''
    # sigma = np.array(list(xsecs['SPL']))
    sigma = np.array([thisspl(temp) for thisspl in xsecs['SPL']])
    nu = np.array(xsecs['SPLNU'])
    select = np.argsort(nu)
    nu = nu[select]
    sigma = sigma[select]
    return sigma, nu


# ------------------------- ------------------------------------------
# -- CIA -- ----------------------------------------------------------
def getciaxs(temp, xsecs):
    '''
    G. ROUDIER: Wrapper around CIA Cerberus library
    '''
    sigma = np.array([thisspl(temp) for thisspl in xsecs['SPL']])
    nu = np.array(xsecs['SPLNU'])
    select = np.argsort(nu)
    nu = nu[select]
    sigma = sigma[select]
    return sigma, nu


# ----------------------------- --------------------------------------
# -- PRESSURE BROADENING -- ------------------------------------------
def intflor(wave, dwave, nu, gamma):
    '''
    G. ROUDIER: Pressure Broadening
    '''
    f = (
        1e0
        / np.pi
        * (
            np.arctan((wave + dwave - nu) / gamma)
            - np.arctan((wave - dwave - nu) / gamma)
        )
    )
    return f


# -------------------------- -----------------------------------------
# -- PYMC DETERMINISTIC FUNCTIONS -- ---------------------------------
def cloudyfmcerberus(*crbinputs):
    '''
    G. ROUDIER: Wrapper around Cerberus forward model, spherical shell symmetry
    '''
    if 'T' in ctxt.fixedParams:
        tpr = ctxt.fixedParams['T']
        ctp, hazescale, hazeloc, hazethick, mdp = crbinputs
    else:
        ctp, hazescale, hazeloc, hazethick, tpr, mdp = crbinputs
    # print(
    #    ' not-fixed cloud parameters (cloudy) cloudstuff,T,mdp:',
    #    ctp,
    #    hazescale,
    #    hazeloc,
    #    hazethick,
    #    tpr,
    #    mdp
    # )

    # this extra list[] is needed for the single param case (only metallicity)
    if not isinstance(mdp, list):
        mdp = [mdp]

    fmc = np.zeros(ctxt.tspectrum.size)
    if ctxt.model in ['TEC', 'TEA']:
        tceqdict = {}
        mdpindex = 0

        if 'XtoH' in ctxt.fixedParams:
            tceqdict['XtoH'] = ctxt.fixedParams['XtoH']
        else:
            tceqdict['XtoH'] = mdp[mdpindex]
            mdpindex += 1

        if 'CtoO' in ctxt.fixedParams:
            tceqdict['CtoO'] = ctxt.fixedParams['CtoO']
        else:
            tceqdict['CtoO'] = mdp[mdpindex]
            mdpindex += 1

        if 'NtoO' in ctxt.fixedParams:
            tceqdict['NtoO'] = ctxt.fixedParams['NtoO']
        else:
            tceqdict['NtoO'] = mdp[mdpindex]
        # print(' XtoH,CtoO,NtoO =',tceqdict['XtoH'],tceqdict['CtoO'],tceqdict['NtoO'])
        fmc = crbmodel(
            tpr,
            ctp,
            hazescale=hazescale,
            hazeloc=hazeloc,
            hazethick=hazethick,
            cheq=tceqdict,
        )
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = mdp[index]

        fmc = crbmodel(
            tpr,
            ctp,
            hazescale=hazescale,
            hazeloc=hazeloc,
            hazethick=hazethick,
            mixratio=mixratio,
        )

    fmc = fmc[ctxt.cleanup]

    if len(ctxt.mcmcsig) > 0:
        fmc += np.average(ctxt.mcmcdat - fmc, weights=1 / ctxt.mcmcsig**2)

    return fmc


def clearfmcerberus(*crbinputs):
    '''
    Wrapper around Cerberus forward model - NO CLOUDS!
    (Note that this is not actually a cloud-free model; it is a fixed-cloud model!!)
    '''
    # these fixed values are probably set in ariel/core, e.g. -10 for hazescale
    ctp = ctxt.fixedParams['CTP']
    hazescale = ctxt.fixedParams['HScale']
    hazeloc = ctxt.fixedParams['HLoc']
    hazethick = ctxt.fixedParams['HThick']
    # print(' fixed cloud parameters (clear):',ctp,hazescale,hazeloc,hazethick)

    if 'T' in ctxt.fixedParams:
        tpr = ctxt.fixedParams['T']
        mdp = crbinputs[0]
        # this extra list[] is needed for the single param case (only metallicity)
        if not isinstance(mdp, list):
            mdp = [mdp]
    else:
        tpr, mdp = crbinputs
    # print(' param values inside of forward model', tpr, mdp)

    fmc = np.zeros(ctxt.tspectrum.size)
    if ctxt.model in ['TEC', 'TEA']:
        tceqdict = {}
        mdpindex = 0
        if 'XtoH' in ctxt.fixedParams:
            tceqdict['XtoH'] = ctxt.fixedParams['XtoH']
        else:
            tceqdict['XtoH'] = mdp[mdpindex]
            mdpindex += 1

        if 'CtoO' in ctxt.fixedParams:
            tceqdict['CtoO'] = ctxt.fixedParams['CtoO']
        else:
            tceqdict['CtoO'] = mdp[mdpindex]
            mdpindex += 1

        if 'NtoO' in ctxt.fixedParams:
            tceqdict['NtoO'] = ctxt.fixedParams['NtoO']
        else:
            tceqdict['NtoO'] = mdp[mdpindex]
        # print('XtoH,CtoO,NtoO =',tceqdict['XtoH'],tceqdict['CtoO'],tceqdict['NtoO'])

        # print('calculating forward model XtoH =', tceqdict['XtoH'])

        fmc = crbmodel(
            tpr,
            float(ctp),
            hazescale=float(hazescale),
            hazeloc=float(hazeloc),
            hazethick=float(hazethick),
            cheq=tceqdict,
        )
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = mdp[index]
            pass
        fmc = crbmodel(
            tpr,
            float(ctp),
            hazescale=float(hazescale),
            hazeloc=float(hazeloc),
            hazethick=float(hazethick),
            mixratio=mixratio,
        )

    fmc = fmc[ctxt.cleanup]

    # (no need for isfinite check; that's what cleanup does already)
    # if np.all(np.isfinite(ctxt.mcmcdat)):
    fmc += np.average(ctxt.mcmcdat - fmc, weights=1 / ctxt.mcmcsig**2)
    # else:
    #     fmc += np.nanmean(ctxt.mcmcdat - fmc)

    return fmc


def offcerberus(*crbinputs):
    '''
    R.ESTRELA: ADD offsets between STIS filters and STIS and WFC3 filters
    '''
    ctp, hazescale, off0, off1, off2, hazeloc, hazethick, tpr, mdp = crbinputs
    #     off0, off1, off2 = crbinputs
    #     ctp = -2.5744083
    #     hazescale = -1.425234
    #     hazeloc = -0.406851
    #     hazethick = 5.58950953
    #     tpr = 1551.41137
    #     mdp = [-1.24882918, -4.08582557, -2.4664526]
    flt = np.array(ctxt.spc['data'][ctxt.planet]['Fltrs'])
    fmc = np.zeros(ctxt.tspectrum.size)
    if ctxt.model in ['TEC', 'TEA']:
        tceqdict = {}
        tceqdict['XtoH'] = float(mdp[0])
        tceqdict['CtoO'] = float(mdp[1])
        tceqdict['NtoO'] = float(mdp[2])
        fmc = crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            cheq=tceqdict,
        )
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
        fmc = crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            mixratio=mixratio,
        )
    cond_G430 = flt[ctxt.cleanup] == 'HST-STIS-CCD-G430L-STARE'
    cond_G141 = flt[ctxt.cleanup] == 'HST-WFC3-IR-G141-SCAN'
    tspectrum_clean = ctxt.tspectrum[ctxt.cleanup]
    fmc = fmc[ctxt.cleanup] - np.nanmean(fmc[ctxt.cleanup][cond_G141])
    fmc = fmc + np.nanmean(tspectrum_clean[cond_G141])
    cond_G750 = flt[ctxt.cleanup] == 'HST-STIS-CCD-G750L-STARE'
    cond_G102 = flt[ctxt.cleanup] == 'HST-WFC3-IR-G102-SCAN'
    fmc[cond_G430] = fmc[cond_G430] - 1e-2 * float(off0)
    fmc[cond_G750] = fmc[cond_G750] - 1e-2 * float(off1)
    fmc[cond_G102] = fmc[cond_G102] - 1e-2 * float(off2)
    return fmc


def offcerberus1(*crbinputs):
    '''
    R.ESTRELA: ADD offsets between STIS filters and STIS and WFC3 filters
    '''
    ctp, hazescale, off0, off1, hazeloc, hazethick, tpr, mdp = crbinputs
    fmc = np.zeros(ctxt.tspectrum.size)
    if ctxt.model in ['TEC', 'TEA']:
        tceqdict = {}
        tceqdict['XtoH'] = float(mdp[0])
        tceqdict['CtoO'] = float(mdp[1])
        tceqdict['NtoO'] = float(mdp[2])
        fmc = crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            cheq=tceqdict,
        )
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
        fmc = crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            mixratio=mixratio,
        )
    fmc = fmc[ctxt.cleanup] - np.nanmean(fmc[ctxt.cleanup])
    fmc = fmc + np.nanmean(ctxt.tspectrum[ctxt.cleanup])
    flt = np.array(ctxt.spc['data'][ctxt.planet]['Fltrs'])
    cond_G430 = 'HST-STIS-CCD-G430L-STARE' in flt
    cond_G750 = 'HST-STIS-CCD-G750L-STARE' in flt
    fmc[cond_G430] = fmc[cond_G430] + 1e-2 * float(off0)
    fmc[cond_G750] = fmc[cond_G750] + 1e-2 * float(off1)
    return fmc


def offcerberus2(*crbinputs):
    '''
    R.ESTRELA: ADD offsets between STIS filters and STIS and WFC3 filters
    '''
    ctp, hazescale, off0, off1, hazeloc, hazethick, tpr, mdp = crbinputs
    fmc = np.zeros(ctxt.tspectrum.size)
    if ctxt.model in ['TEC', 'TEA']:
        tceqdict = {}
        tceqdict['XtoH'] = float(mdp[0])
        tceqdict['CtoO'] = float(mdp[1])
        tceqdict['NtoO'] = float(mdp[2])
        fmc = crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            cheq=tceqdict,
        )
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
        fmc = crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            mixratio=mixratio,
        )
    #    fmc = fmc[ctxt.cleanup] - np.nanmean(fmc[ctxt.cleanup])
    #    fmc = fmc + np.nanmean(ctxt.tspectrum[ctxt.cleanup])
    flt = np.array(ctxt.spc['data'][ctxt.planet]['Fltrs'])
    cond_G430 = 'HST-STIS-CCD-G430-STARE' in flt
    cond_G750 = 'HST-STIS-CCD-G750-STARE' in flt
    fmc[cond_G430] = fmc[cond_G430] + 1e-2 * float(off0)
    fmc[cond_G750] = fmc[cond_G750] + 1e-2 * float(off1)
    return fmc


def offcerberus3(*crbinputs):
    '''
    R.ESTRELA: ADD offsets between STIS filters and STIS and WFC3 filters
    '''
    ctp, hazescale, off0, off1, hazeloc, hazethick, tpr, mdp = crbinputs
    fmc = np.zeros(ctxt.tspectrum.size)
    flt = np.array(ctxt.spc['data'][ctxt.planet]['Fltrs'])
    if ctxt.model in ['TEC', 'TEA']:
        tceqdict = {}
        tceqdict['XtoH'] = float(mdp[0])
        tceqdict['CtoO'] = float(mdp[1])
        tceqdict['NtoO'] = float(mdp[2])
        fmc = crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            cheq=tceqdict,
        )
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
        fmc = crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            mixratio=mixratio,
        )
    fmc = fmc[ctxt.cleanup] - np.nanmean(fmc[ctxt.cleanup])
    fmc = fmc + np.nanmean(ctxt.tspectrum[ctxt.cleanup])
    cond_G430 = 'HST-STIS-CCD-G430-STARE' in flt
    cond_G102 = 'HST-WFC3-IR-G102-SCAN' in flt
    fmc[cond_G430] = fmc[cond_G430] + 1e-2 * float(off0)
    fmc[cond_G102] = fmc[cond_G102] + 1e-2 * float(off1)
    return fmc


def offcerberus4(*crbinputs):
    '''
    R.ESTRELA: ADD offsets between STIS filters and STIS and WFC3 filters
    '''
    ctp, hazescale, off0, hazeloc, hazethick, tpr, mdp = crbinputs
    fmc = np.zeros(ctxt.tspectrum.size)
    flt = np.array(ctxt.spc['data'][ctxt.planet]['Fltrs'])
    if ctxt.model in ['TEC', 'TEA']:
        tceqdict = {}
        tceqdict['XtoH'] = float(mdp[0])
        tceqdict['CtoO'] = float(mdp[1])
        tceqdict['NtoO'] = float(mdp[2])
        fmc = crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            cheq=tceqdict,
        )
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
        fmc = crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            mixratio=mixratio,
        )
    fmc = fmc[ctxt.cleanup] - np.nanmean(fmc[ctxt.cleanup])
    fmc = fmc + np.nanmean(ctxt.tspectrum[ctxt.cleanup])
    cond_G430 = 'HST-STIS-CCD-G430-STARE' in flt
    fmc[cond_G430] = fmc[cond_G430] + 1e-2 * float(off0)
    return fmc


def offcerberus5(*crbinputs):
    '''
    R.ESTRELA: ADD offsets between STIS filters and STIS and WFC3 filters
    '''
    ctp, hazescale, off0, off1, hazeloc, hazethick, tpr, mdp = crbinputs
    fmc = np.zeros(ctxt.tspectrum.size)
    flt = np.array(ctxt.spc['data'][ctxt.planet]['Fltrs'])
    if ctxt.model in ['TEC', 'TEA']:
        tceqdict = {}
        tceqdict['XtoH'] = float(mdp[0])
        tceqdict['CtoO'] = float(mdp[1])
        tceqdict['NtoO'] = float(mdp[2])
        fmc = crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            cheq=tceqdict,
        )
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
        fmc = crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            mixratio=mixratio,
        )
    fmc = fmc[ctxt.cleanup] - np.nanmean(fmc[ctxt.cleanup])
    fmc = fmc + np.nanmean(ctxt.tspectrum[ctxt.cleanup])
    cond_G102 = 'HST-WFC3-IR-G102-SCAN' in flt
    cond_G750 = 'HST-STIS-CCD-G750-STARE' in flt
    fmc[cond_G750] = fmc[cond_G750] + 1e-2 * float(off0)
    fmc[cond_G102] = fmc[cond_G102] + 1e-2 * float(off1)
    return fmc


def offcerberus6(*crbinputs):
    '''
    R.ESTRELA: ADD offsets between STIS filters and STIS and WFC3 filters
    '''
    ctp, hazescale, off0, hazeloc, hazethick, tpr, mdp = crbinputs
    fmc = np.zeros(ctxt.tspectrum.size)
    flt = np.array(ctxt.spc['data'][ctxt.planet]['Fltrs'])
    if ctxt.model in ['TEC', 'TEA']:
        tceqdict = {}
        tceqdict['XtoH'] = float(mdp[0])
        tceqdict['CtoO'] = float(mdp[1])
        tceqdict['NtoO'] = float(mdp[2])
        fmc = crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            cheq=tceqdict,
        )
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
        fmc = crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            mixratio=mixratio,
        )
    fmc = fmc[ctxt.cleanup] - np.nanmean(fmc[ctxt.cleanup])
    fmc = fmc + np.nanmean(ctxt.tspectrum[ctxt.cleanup])
    cond_G750 = 'HST-STIS-CCD-G750-STARE' in flt
    fmc[cond_G750] = fmc[cond_G750] + 1e-2 * float(off0)
    return fmc


def offcerberus7(*crbinputs):
    '''
    R.ESTRELA: ADD offsets between STIS filters and WFC3 filters
    '''
    ctp, hazescale, off0, hazeloc, hazethick, tpr, mdp = crbinputs
    fmc = np.zeros(ctxt.tspectrum.size)
    flt = np.array(ctxt.spc['data'][ctxt.planet]['Fltrs'])
    if ctxt.model in ['TEC', 'TEA']:
        tceqdict = {}
        tceqdict['XtoH'] = float(mdp[0])
        tceqdict['CtoO'] = float(mdp[1])
        tceqdict['NtoO'] = float(mdp[2])
        fmc = crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            cheq=tceqdict,
        )
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
        fmc = crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            mixratio=mixratio,
        )
    fmc = fmc[ctxt.cleanup] - np.nanmean(fmc[ctxt.cleanup])
    fmc = fmc + np.nanmean(ctxt.tspectrum[ctxt.cleanup])
    cond_G750 = 'HST-STIS-CCD-G750-STARE' in flt
    fmc[cond_G750] = fmc[cond_G750] + 1e-2 * float(off0)
    return fmc


def offcerberus8(*crbinputs):
    '''
    R.ESTRELA: ADD offsets between WFC3 filters
    '''
    ctp, hazescale, off0, hazeloc, hazethick, tpr, mdp = crbinputs
    fmc = np.zeros(ctxt.tspectrum.size)
    flt = np.array(ctxt.spc['data'][ctxt.planet]['Fltrs'])
    if ctxt.model in ['TEC', 'TEA']:
        tceqdict = {}
        tceqdict['XtoH'] = float(mdp[0])
        tceqdict['CtoO'] = float(mdp[1])
        tceqdict['NtoO'] = float(mdp[2])
        fmc = crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            cheq=tceqdict,
        )
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
        fmc = crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            mixratio=mixratio,
        )
    fmc = fmc[ctxt.cleanup] - np.nanmean(fmc[ctxt.cleanup])
    fmc = fmc + np.nanmean(ctxt.tspectrum[ctxt.cleanup])
    cond_G102 = 'HST-WFC3-IR-G102-SCAN' in flt
    fmc[cond_G102] = fmc[cond_G102] + 1e-2 * float(off0)
    return fmc
