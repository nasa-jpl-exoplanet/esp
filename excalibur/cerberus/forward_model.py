'''cerberus forward_model ds'''

# Heritage code shame:
# pylint: disable=invalid-name,no-member
# pylint: disable=too-many-arguments,too-many-branches,too-many-lines,too-many-locals,too-many-positional-arguments,too-many-statements

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as cst
from scipy.interpolate import interp1d as itp
import logging

import scipy.special as scipyspecial

from deprecated import deprecated

import excalibur

import excalibur.system.core as syscore

from excalibur.util.cerberus import crbce, calcTEA, getmmw

from excalibur.cerberus.fmcontext import ctxtinit

log = logging.getLogger(__name__)

# this doesn't change results at all; just needed to avoid undefined-variable pylint
ctxt = ctxtinit()


# --------------------------------------------------------------------
# -- CERBERUS FORWARD MODEL ------------------------------------------
class crbFM:
    def __init__(self):
        self.__spectrum = np.empty(0)
        self.__breakdown_by_molecule = {}
        self.__moleculeProfiles = {}
        self.__opticalDepthProfiles = {}
        self.__pressureGrid = np.empty(0)

    def crbmodel(
        self,
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
        lbroadening=None,
        lshifting=None,
        knownspecies=None,
        cialist=None,
        xmollist=None,
        atom_list=None,
        nlevels=None,
        Hsmax=None,
        solrad=None,
        break_down_by_molecule=False,
        logx=False,
        verbose=False,
        debug=False,
        atom_data=None,
        improvedBoundaryCondition=True,
        extendedBoundaryCondition=False,
    ):
        '''
        G. ROUDIER: Cerberus forward model probing up to 'Hsmax' scale heights from solid
        radius solrad evenly log divided amongst nlevels steps
        - TP profile
        - VMR profile
        - MMW profile
        '''
        if atom_list is None:
            atom_list = ['Ca', 'K', 'Na']
        if planet is None:
            planet = ctxt.planet
        if orbp is None:
            orbp = ctxt.orbp
        if knownspecies is None:
            knownspecies = ctxt.knownspecies
        if cialist is None:
            cialist = ctxt.cialist
        # else:
        #    cialist = ['H2-H', 'H2-H2', 'H2-He', 'He-H']
        if xmollist is None:
            xmollist = ctxt.xmollist
        # else:
        #    xmollist = [
        #        'TIO',
        #        'H2O',
        #        'H2CO',
        #        'HCN',
        #        'CO',
        #        'CO2',
        #        'NH3',
        #        'CH4',
        #        'C2H2',
        #    ]
        # this is passed in. why reset it here?
        #    # longer list currently used by Luke:
        #  PUT THIS INTO RUNTIME OPS!!
        #    # xmollist = ['TIO', 'H2O', 'HCN', 'CO', 'CO2', 'NH3', 'CH4', 'H2S','PH3', 'C2H2', 'OH', 'O2', 'O3', 'SO2', 'C2H6', 'C3H8', 'CH3CHO']
        # hmm some of these are actually in HITRAN, nor EXOMOL, e.g. O2 O3
        # new ones: 'H2S','PH3', 'SO2', 'C2H6', 'C3H8', 'CH3CHO'
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

        temp = np.array(temp)
        # print('  temp', temp)
        if temp.ndim:
            tpp = temp
        else:
            tpp = np.array([float(temp)] * nlevels)
            pass
        # print('  tpp', tpp)
        # option for non-isothermal T-P profile
        #  if the temperature array has just a handful of elements,
        #  then it's actually the parameters for a T-P profile
        if len(tpp) not in [int(nlevels)]:
            # print('tpp forward model with non-Isothermal T-P profile!!')
            tpp = TPprofile(temp, pressure)
        # print('  tpp', tpp)

        # verify that the temperature array has the right length (nlevels)
        if len(tpp) not in [int(nlevels)]:
            print('!!! >--< TP PROFILE != PRESSURE GRID: %s nlevels', nlevels)
            log.error(
                '!!! >--< TP PROFILE != PRESSURE GRID: %s nlevels', nlevels
            )
            pass

        mixratioprofiles = {}
        if not mixratio:
            # chemical equilibrium case
            if cheq is None:
                log.error('!!! >--< Neither mixratio nor cheq are defined')
                pass
            if chemistry.startswith('TEC'):
                mixratio, mixratioprofiles, fH2, fHe = crbce(
                    pressure,
                    tpp,
                    C2Or=cheq['CtoO'],
                    X2Hr=cheq['XtoH'],
                    N2Or=cheq['NtoO'],
                )
                mmw, fH2, fHe = getmmw(
                    mixratio,
                    protosolar=False,
                    fH2=fH2,
                    fHe=fHe,
                )
                pass

            elif chemistry.startswith('TEA'):
                #  (this one gives a div-by-0 error)
                # tempCoeffs = [0, temp, 0, 0, 0, 0, 0, 0, 0, 0]
                #  this is the correct way to pass in to Luke's _make_tp_profile
                # tempCoeffs = [0, temp, 0, 1, 0, -1, 1, 0, -1, 1]  # isothermal
                #  but now we're passing in the T array directly, not params for it
                mixratioprofiles = calcTEA(
                    tpp,
                    # tempCoeffs,
                    pressure,
                    metallicity=10.0 ** cheq['XtoH'],
                    C_O=0.55 * 10.0 ** cheq['CtoO'],
                    # N_O=?? * 10.0 ** cheq['NtoO'],
                )

                # have to take the average! (same as done in crbce)
                mixratio = {}
                for molecule in mixratioprofiles:
                    mixratio[molecule] = np.log10(
                        np.mean(10.0 ** mixratioprofiles[molecule])
                    )
                # print()
                # print('mixratio in cerb', mixratio)
                mmw, fH2, fHe = getmmw(mixratio)
                # print('TEA mmw, fH2, fHe', mmw, fH2, fHe)

                if 'ozone' in chemistry:
                    # print()
                    # print('OZONE CHECK')
                    if 'O3' not in mixratio:
                        log.error('O3 not selected for ozone model!!')

                    # add in ozone, but keep the same total metallicity
                    originalmetals = 0
                    for molecule in mixratio:
                        originalmetals += 10.0 ** mixratio[molecule]
                    # print('originalmetals', originalmetals/1.e6)

                    # print(' ozone mixratio before', mixratio['O3'])
                    # mixratio['O3'] = mixratio['O3'] * 0 + 5.0
                    # mixratio['O3'] = mixratio['O3'] * 0 + 7.0
                    # print(' NEW OZONE mixratio', mixratio['O3'])

                    # totalmetals = 0
                    # for molecule in mixratio:
                    #     totalmetals += 10.0 ** mixratio[molecule]
                    # print('totalmetals', totalmetals/1.e6)

                    newsum = originalmetals + 10.0 ** mixratio['O3']

                    for molecule in mixratio:
                        mixratio[molecule] -= np.log10(newsum / originalmetals)
                    # print(' ozone mixratio renorm', mixratio['O3'])
                    # totalmetals = 0
                    # for molecule in mixratio:
                    #    totalmetals += 10.0 ** mixratio[molecule]
                    # print('totalmetals', totalmetals/1.e6)

                    mmw, fH2, fHe = getmmw(mixratio)
                    # print('TEA mmw, fH2, fHe', mmw, fH2, fHe)

            else:
                fH2 = 0
                fHe = 0
                mixratio = {}
                mixratioprofiles = {}
                mmw = 1
                log.error('!!! >--< UNKNOWN CHEM MODEL: %s', chemistry)
                pass

            pass
        else:
            # DISEQ case
            mmw, fH2, fHe = getmmw(mixratio)

            # mixing ratio is a fixed value for all atmospheric pressures
            for molecule in mixratio:
                mixratioprofiles[molecule] = np.full(
                    (len(pressure)), mixratio[molecule]
                )

        # make sure that the mixing ratios are 1-d arrays over the pressure grid
        #  (otherwise later calls may get mis-matched broadcasting problems)
        if not fH2.ndim:
            fH2 = np.array([float(fH2)] * len(tpp))
        if not fHe.ndim:
            fHe = np.array([float(fHe)] * len(tpp))
        for molecule in mixratio:
            # mixratio[molecule] = np.array(mixratio[molecule])
            if not mixratio[molecule].ndim:
                mixratio[molecule] = np.array(
                    [float(mixratio[molecule])] * len(tpp)
                )

            # verify that the mixratio array has the right length (nlevels)
            if len(mixratio[molecule]) not in [int(nlevels)]:
                log.error(
                    '!!! >--< MIXRATIO PROFILE != PRESSURE GRID: %s nlevels',
                    nlevels,
                )

        mmw = mmw * cst.m_p  # [kg]

        if debug:
            log.info('>-- mmw: %s', mmw * 6.022e26)
            pass

        Hs = (
            cst.Boltzmann
            * tpp
            / (mmw * 1e-2 * (10.0 ** float(orbp[planet]['logg'])))
        )  # [m]

        # print('T',tpp)
        # print('mmw',mmw)
        # print('logg',10.0 ** float(orbp[planet]['logg']))
        # print('Hs!!!!!', Hs)

        # when the Pressure grid is log-spaced, rdz is a constant
        #  drop dz[] and dzprime[] arrays and just use this constant instead
        dz = 2 * abs(Hs / 2.0 * np.log(1.0 + dPoverP))

        # CB, the linspace was adapted in the case of constant dz
        z = np.concatenate(([0], np.cumsum(dz[:-1])))

        rho = pressure * 1e5 / (cst.Boltzmann * tpp)
        tau, tau_by_molecule, wtau = gettau(
            xsecs,
            qtgrid,
            tpp,
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
            atom_list,
            fH2,
            fHe,
            xmollist,
            hazescale,
            hzlib,
            hazeprof,
            hazeslope,
            hazeloc,
            hazethick,
            atom_data,
            debug=debug,
            improvedBoundaryCondition=improvedBoundaryCondition,
            extendedBoundaryCondition=extendedBoundaryCondition,
        )
        if not break_down_by_molecule:
            tau_by_molecule = {}
            pass
        molecules = tau_by_molecule.keys()
        # SEMI FINITE CLOUD ---------------------------------------------------
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
                myspl = itp(
                    reversep, np.asarray(tau[:, waveindex]).flatten()[::-1]
                )

                tau[cloudtopindex, waveindex] = myspl(10.0**cloudtp)

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
                Hs[cloudtopindex]
                / 2.0
                * np.log(1.0 + ctpdpress / pressure[cloudtopindex])
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
        plotmodel = model.copy()
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
            fig, ax1 = plt.subplots(figsize=(10, 6))
            ax2 = ax1.twiny()

            for k in mixratio.items():
                ax1.plot(pressure * 0 + k[1], pressure, label=k[0])
                pass
            ax1.legend(loc='upper left')

            ax2.plot(tpp, pressure, color='pink', label='Temperature', ls='--')
            ax2.legend(loc='center left')

            plt.semilogy()
            plt.gca().invert_yaxis()
            ax1.set_xlabel('log(VMR) [ppm]')
            ax2.set_xlabel('Temperature [K]')
            ax1.set_ylabel('Pressure [bar]')
            plt.show()

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
            ) / np.nanmean(Hs)
            ax2max = (
                np.sqrt(1e-2 * yaxmax) * orbp['R*'] * ssc['Rsun'] - rp0hs
            ) / np.nanmean(Hs)
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
        self.__spectrum = model
        self.__breakdown_by_molecule = models_by_molecule
        self.__pressureGrid = pressure
        self.__moleculeProfiles = mixratioprofiles
        self.__opticalDepthProfiles = tau  # (tau_by_molecule is also available)
        return self

    @property
    def spectrum(self):
        return self.__spectrum

    @property
    def breakdown_by_molecule(self):
        return self.__breakdown_by_molecule

    @property
    def pressureGrid(self):
        return self.__pressureGrid

    @property
    def moleculeProfiles(self):
        return self.__moleculeProfiles

    @property
    def opticalDepthProfiles(self):
        return self.__opticalDepthProfiles


@deprecated(
    'replace crbmodel() with crbFM().crbmodel() and use .spectrum method'
)
def crbmodel(temp, cloudtp, **kwargs):
    log.info('use crbFM class to get additional saved results from crbmodel')
    return crbFM().crbmodel(temp, cloudtp, **kwargs).spectrum


def TPprofile(sparseTgrid, pressures):
    '''
    interpolate a small set of temperatures over the full pressure grid
    '''

    sparsePgrid = np.linspace(
        np.log10(pressures[0]), np.log10(pressures[-1]), len(sparseTgrid)
    )

    # interp requires sparsePgrid to be increasing, so reverse its order
    temperatures = np.interp(
        np.log10(pressures), sparsePgrid[::-1], sparseTgrid[::-1]
    )
    return temperatures


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
    atom_list,
    fH2,
    fHe,
    xmollist,
    hazescale,
    hzlib,
    hazeprof,
    hazeslope,
    hazeloc,
    hazethick,
    atom_data,
    debug=False,
    improvedBoundaryCondition=True,
    extendedBoundaryCondition=False,
):
    '''
    G. ROUDIER: Builds optical depth matrix
    '''

    # SPHERICAL SHELL (PLANE-PARALLEL REMOVED) -----------------------------------
    # MATRICES INIT --------------------------------------------------------------
    Nzones = len(pressure)
    tau = np.zeros((Nzones, wgrid.size))
    tau_by_molecule = {}
    toptau_by_molecule = {}
    analytictau_by_molecule = {}
    # DL ARRAY, Z VERSUS ZPRIME --------------------------------------------------
    zprime = np.broadcast_to(z, (Nzones, Nzones))
    thisz = zprime.T

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

    # print('dlarray shape', dlarray.shape)
    # print('dlarray[0]', dlarray[0])
    # print('dlarray[-1]', dlarray[-1])

    top_rho = rho[-1]
    # bottom_rho = rho[0]
    # rp0 is in units of meters.  z,dz also
    # print('Rplanet', rp0)  #3e7 = 300km ?!  should be 4.8*Rearth = 3e4 km
    # print('top of atmosphere', z[0],z[-1])  # 4.7e6 = 4700 km
    # print('top of atmosphere', sum(dz))  # ok the last dz is not used, right?
    # Hestimate = np.exp(np.log(z[0] / z[-1]) / 20)   #
    # Rtop = rp0 + z[-1]
    NscaleHeights = np.log(pressure[0] / pressure[-1])
    Hestimate = (z[-1] - z[0]) / NscaleHeights
    # earth scale height is 8.5km.
    #  scale height goes as T/g = T Rp^2 / Mp
    #   M is 25 Mearth; R is 4.8; T is 800 so that predicts H = 22.4 km
    # print('scale height', Hestimate)  # 235 km  way off!
    # oh! ok mean molecular weight is much lower here.  it's like 3.1

    analyticIntegral = np.sqrt(2 * np.pi * Hestimate * (rp0 + z))
    # print('shape', analyticIntegral.shape)

    # GAS ARRAY, ZPRIME VERSUS WAVELENGTH  ---------------------------------------
    for elem in mixratio:
        mlp = np.array(mixratio[elem])
        if not mlp.ndim:
            mlp = np.array([float(mlp)] * len(pressure))
            pass
        if len(mlp) not in [len(pressure)]:
            log.error('!!! >--< %s VMR PROFILE NOT ON PRESSURE GRID', elem)
            pass
        mmr = 10.0 ** (mlp - 6.0)  # mmr.shape(n_pressure)
        top_mmr = mmr[-1]
        if elem not in xsecs:
            # TEA species might not have cross-sections calculated
            if elem in ['H2', 'He']:
                # ignore missing xsecs for molecules without strong features
                pass
            else:
                if elem in atom_list:
                    interp_atom = (
                        excalibur.cerberus.forward_model.ctxt.atom_xsec
                    )
                    if interp_atom is None:
                        interp_atom = atom_data
                        pass
                    if interp_atom is not None:
                        # interpolator loading
                        interpolator = interp_atom[elem]
                        T = np.repeat(temp, len(wgrid))
                        P = np.repeat(pressure, len(wgrid))
                        X_H2 = np.repeat(fH2 / (fH2 + fHe), len(wgrid))
                        wl = np.tile(wgrid, Nzones)
                        points = np.column_stack((T, P, X_H2, wl))

                        # xsec computation
                        sigma = interpolator(points)
                        sigma = np.reshape(
                            sigma, (Nzones, len(wgrid))
                        )  # cm^2/mol
                        sigma = sigma[:, ::-1].T

                        # sigma.shape(n_waves, n_pressure)
                        sigma = sigma * 1e-4  # m^2/mol
                        lsig = 1e4 / wgrid[::-1]
                        pass
                    pass
                else:
                    log.error(
                        'MISSING CROSS-SECTION: add this molecule to runtime EXOMOL  %s',
                        elem,
                    )

        else:
            # Fake use of xmollist due to changes in xslib v112
            # THIS HAS TO BE FIXED
            # if elem not in xmollist:
            if not xmollist:
                # HITEMP/HITRAN ROTHMAN ET AL. 2010 ------------------------------
                sigma, lsig = absorb(
                    xsecs[elem],
                    qtgrid[elem],
                    temp,
                    pressure,
                    mmr,
                    lbroadening,
                    lshifting,
                    wgrid,
                )  # cm^2/mol
                if True in (sigma < 0):
                    sigma[sigma < 0] = 0e0
                    pass
                if True in ~np.isfinite(sigma):
                    sigma[~np.isfinite(sigma)] = 0e0
                    pass
                sigma = sigma * 1e-4  # m^2/mol
                pass
            else:
                # EXOMOL HILL ET AL. 2013 ----------------------------------------
                sigma, lsig = getxmolxs(temp, xsecs[elem])  # cm^2/mol
                # sigma.shape(n_waves, n_pressure)
                if True in (sigma < 0):
                    sigma[sigma < 0] = 0e0
                    pass
                if True in ~np.isfinite(sigma):
                    sigma[~np.isfinite(sigma)] = 0e0
                    pass
                sigma = sigma * 1e-4  # m^2/mol
                pass

            top_sigma = sigma[:, -1]
            # this is a 1-D array, not 3-D.  just a function of wavelength
            toptau_by_molecule[elem] = top_rho * top_mmr * top_sigma
            # print('  shape check',top_rho.shape,sigma.shape,top_mmr.shape)
            # print('toptau shape', toptau_by_molecule[elem].shape) #103

            if extendedBoundaryCondition:
                # analytictau = analyticIntegral * rho * top_mmr * top_sigma
                # print(analyticIntegral.shape,
                #      rho.shape,
                #      mmr.shape,
                #      sigma.shape)

                # two options: use the full range of sigma
                #  or, to be fair, just use the 50th element to match below
                # ok also mmr should just take the top value, if it's true B.C.
                analytictau_by_molecule[elem] = (
                    analyticIntegral[:, np.newaxis]
                    * rho[:, np.newaxis]
                    * mmr[:, np.newaxis][49, :][np.newaxis, :]
                    * sigma.T[49, :][np.newaxis, :]
                )
                # print('anal shape', analytictau_by_molecule[elem].shape)

        # CB sigma (Nzones, Nzones, N_waves)
        # 1st dimension corrsponds to z
        # 2nd dimension corrsponds to z'
        # 3rd dimension corresponds to wavelength
        sigma = np.broadcast_to(
            sigma.T[None, :, :], (Nzones, Nzones, len(wgrid))
        ).copy()
        tau_by_molecule[elem] = (
            (rho * np.ones((Nzones, Nzones)))[:, :, np.newaxis]
            * (mmr * np.ones((Nzones, Nzones)))[:, :, np.newaxis]
            * sigma
        )

        tau = tau + tau_by_molecule[elem]

        pass
    # CIA ARRAY, ZPRIME VERSUS WAVELENGTH  ---------------------------------------
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
                '--< CERBERUS gettau(): UNEXPECTED CIA SPECIES %s >--', cia
            )
            f1 = 0
            f2 = 0
            pass
        # HITRAN RICHARD ET AL. 2012
        sigma, lsig = getciaxs(temp, xsecs[cia])  # cm^5/mol^2
        # sigma.shape(n_waves, n_pressure)
        sigma = np.array(sigma) * 1e-10  # m^5/mol^2
        if True in (sigma < 0):
            sigma[sigma < 0] = 0e0
            pass
        if True in ~np.isfinite(sigma):
            sigma[~np.isfinite(sigma)] = 0e0
            pass

        top_sigma = sigma[:, -1]
        top_f1 = np.array(f1)[-1]
        top_f2 = np.array(f2)[-1]

        toptau_by_molecule[cia] = top_f1 * top_f2 * top_sigma * top_rho**2

        # CB sigma (Nzones, Nzones, N_waves)
        # 1st dimension corrsponds to z
        # 2nd dimension corrsponds to z'
        # 3rd dimension corresponds to wavelength
        sigma = np.broadcast_to(
            sigma.T[None, :, :], (Nzones, Nzones, len(wgrid))
        ).copy()
        tau_by_molecule[cia] = (
            (f1 * np.ones((Nzones, Nzones)))[:, :, np.newaxis]
            * (f2 * np.ones((Nzones, Nzones)))[:, :, np.newaxis]
            * sigma
            * (rho * np.ones((Nzones, Nzones)))[:, :, np.newaxis] ** 2
        )
        tau = tau + tau_by_molecule[cia]

    # H2 RAYLEIGH ARRAY, ZPRIME VERSUS WAVELENGTH  -------------------------------
    # NAUS & UBACHS 2000
    slambda0 = 750.0 * 1e-3  # microns
    sray0 = 2.52 * 1e-28 * 1e-4  # m^2/mol
    sigma = sray0 * (wgrid[::-1] / slambda0) ** (-4e0)

    top_fH2 = np.array(fH2)[-1]
    toptau_by_molecule['rayleigh'] = top_fH2 * top_rho * sigma

    # CB sigma (Nzones, Nzones, N_waves)
    # 1st dimension corrsponds to z
    # 2nd dimension corrsponds to z'
    # 3rd dimension corresponds to wavelength
    sigma = np.broadcast_to(sigma, (Nzones, Nzones, len(wgrid))).copy()
    tau_by_molecule['rayleigh'] = (
        (fH2 * np.ones((Nzones, Nzones)))[:, :, np.newaxis]
        * (rho * np.ones((Nzones, Nzones)))[:, :, np.newaxis]
        * sigma
    )
    tau = tau + tau_by_molecule['rayleigh']

    # HAZE ARRAY, ZPRIME VERSUS WAVELENGTH  --------------------------------------
    if hzlib is None:
        slambda0 = 750.0 * 1e-3  # microns
        sray0 = 2.52 * 1e-28 * 1e-4  # m^2/mol
        sigma = sray0 * (wgrid[::-1] / slambda0) ** (hazeslope)

        toptau_by_molecule['haze'] = 10.0**hazescale * sigma

        # CB sigma (Nzones, Nzones, N_waves)
        # 1st dimension corrsponds to z
        # 2nd dimension corrsponds to z'
        # 3rd dimension corresponds to wavelength
        sigma = np.broadcast_to(sigma, (Nzones, Nzones, len(wgrid))).copy()
        tau_by_molecule['haze'] = 10.0**hazescale * sigma
        tau = tau + tau_by_molecule['haze']

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
            haze_ref_pressure = float(pressure[rh == np.max(rh)].item())
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

        top_rh = rh[-1]
        toptau_by_molecule['haze'] = 10.0**hazescale * sigma * top_rh

        # CB sigma (Nzones, Nzones, N_waves)
        # 1st dimension corrsponds to z
        # 2nd dimension corrsponds to z'
        # 3rd dimension corresponds to wavelength
        hazecontribution = (
            10.0**hazescale
            * np.broadcast_to(
                sigma * np.array([rh]).T, (Nzones, Nzones, len(wgrid))
            ).copy()
        )
        tau_by_molecule['haze'] = hazecontribution
        tau = tau + tau_by_molecule['haze']
        pass

    # CB tau (Nzones, N_waves)
    # sum over the z' for a fixed z
    tau = (2e0 * dlarray[:, :, np.newaxis] * tau).sum(axis=1)
    molecules = tau_by_molecule.keys()
    for molecule in molecules:
        tau_by_molecule[molecule] = (
            2e0 * dlarray[:, :, np.newaxis] * tau_by_molecule[molecule]
        ).sum(axis=1)
        pass

    # include the upper boundary condition on atmosphere here, after line integral
    #  use an analytic estimate for the integrated depth
    # (toptau.. is already defined above as rho*sigma at top of atmosphere)

    if improvedBoundaryCondition:
        scaleHeightsDown = np.log(pressure / pressure[-1])
        experfEquation = np.exp(scaleHeightsDown) * (
            1 - scipyspecial.erf(np.sqrt(scaleHeightsDown))
        )
        # BCintegral = np.sqrt(2 * np.pi * Hestimate * (Rtop + z)) * experfEquation
        BCintegral = np.sqrt(2 * np.pi * Hestimate * (rp0 + z)) * experfEquation

        # ok careful with array broadcasting here
        # toptau is a wavelength array (len 103)
        # BCintegral is a height array (len 100)
        # print('BCintegral', BCintegral)
        # print('BCintegral shape', BCintegral.shape)
        BCintegral = BCintegral[:, np.newaxis]
        # print('BCintegral shape', BCintegral.shape)

        # print('scaleHdown', scaleHeightsDown)
        # print('experfeq', experfEquation)
        # exit()

        for molecule in molecules:
            # print('  starting molecule=', molecule)
            # print('toptau', toptau_by_molecule[molecule])
            # print('tau shape', tau_by_molecule[molecule].shape) #100x103
            # print('tau shape', tau_by_molecule[molecule][-1].shape) #103
            # print('toptau shape', toptau_by_molecule[molecule].shape) #103

            if molecule in toptau_by_molecule:
                # fractionalChange = (
                #    toptau_by_molecule[molecule][np.newaxis, :]
                #    * BCintegral
                #    / (tau_by_molecule[molecule] + 1.0e-30)
                # )
                # hmm these are all the exact same (at same height). strange...
                # print('fractional change', fractionalChange[40,:])
                # this prints a range of heights (fixed wavelengths)
                # print('fractional change', fractionalChange[:,40])
                # check the H=10,lambda=3 case
                # print(
                #    'fractional change',
                #    molecule,
                #    fractionalChange[48, 48],
                #    np.min(fractionalChange),
                #    np.max(fractionalChange),
                # )

                tau_by_molecule[molecule] += (
                    toptau_by_molecule[molecule][np.newaxis, :] * BCintegral
                )
                tau += toptau_by_molecule[molecule][np.newaxis, :] * BCintegral
            else:
                log.warning(
                    '--< molecule missing from atmos B.C.: %s >--', molecule
                )
    if extendedBoundaryCondition:
        for molecule in molecules:
            if molecule in analytictau_by_molecule:
                # print('shape check vs analytic', molecule,
                #      tau_by_molecule[molecule].shape,
                #      analytictau_by_molecule[molecule].shape)
                # print('check vs analytic', molecule,
                #      analytictau_by_molecule[molecule] /
                #      tau_by_molecule[molecule])
                # pick a single wavelength and see how diff varies with height
                # print(
                #    'check vs analytic',
                #    molecule,
                #    analytictau_by_molecule[molecule][:, 57]
                #    / tau_by_molecule[molecule][:, 57],
                # )

                # replace the upper half with the analytic part. see how it looks
                # (and adjust the overall tau accordingly)
                tau[50:, :] -= tau_by_molecule[molecule][50:, :]
                tau[50:, :] += analytictau_by_molecule[molecule][50:, :]
                tau_by_molecule[molecule][50:, :] = analytictau_by_molecule[
                    molecule
                ][50:, :]

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
        fmc = crbFM().crbmodel(
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

        fmc = crbFM().crbmodel(
            tpr,
            ctp,
            hazescale=hazescale,
            hazeloc=hazeloc,
            hazethick=hazethick,
            mixratio=mixratio,
        )
    fmc = fmc.spectrum
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
    else:
        tpr, mdp = crbinputs
    # print('clearfmcerberus TPR = ', tpr)
    if not isinstance(mdp, list):
        mdp = [mdp]
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

        fmc = crbFM().crbmodel(
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
        fmc = crbFM().crbmodel(
            tpr,
            float(ctp),
            hazescale=float(hazescale),
            hazeloc=float(hazeloc),
            hazethick=float(hazethick),
            mixratio=mixratio,
        )
    fmc = fmc.spectrum
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
        fmc = crbFM().crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            cheq=tceqdict,
        )
        fmc = fmc.spectrum
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
        fmc = crbFM().crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            mixratio=mixratio,
        )
        fmc = fmc.spectrum
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
        fmc = crbFM().crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            cheq=tceqdict,
        )
        fmc = fmc.spectrum
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
        fmc = crbFM().crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            mixratio=mixratio,
        )
        fmc = fmc.spectrum
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
        fmc = crbFM().crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            cheq=tceqdict,
        )
        fmc = fmc.spectrum
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
        fmc = crbFM().crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            mixratio=mixratio,
        )
        fmc = fmc.spectrum
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
        fmc = crbFM().crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            cheq=tceqdict,
        )
        fmc = fmc.spectrum
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
        fmc = crbFM().crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            mixratio=mixratio,
        )
        fmc = fmc.spectrum
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
        fmc = crbFM().crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            cheq=tceqdict,
        )
        fmc = fmc.spectrum
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
        fmc = crbFM().crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            mixratio=mixratio,
        )
        fmc = fmc.spectrum
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
        fmc = crbFM().crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            cheq=tceqdict,
        )
        fmc = fmc.spectrum
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
        fmc = crbFM().crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            mixratio=mixratio,
        )
        fmc = fmc.spectrum
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
        fmc = crbFM().crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            cheq=tceqdict,
        )
        fmc = fmc.spectrum
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
        fmc = crbFM().crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            mixratio=mixratio,
        )
        fmc = fmc.spectrum
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
        fmc = crbFM().crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            cheq=tceqdict,
        )
        fmc = fmc.spectrum
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
        fmc = crbFM().crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            mixratio=mixratio,
        )
        fmc = fmc.spectrum
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
        fmc = crbFM().crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            cheq=tceqdict,
        )
        fmc = fmc.spectrum
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
        fmc = crbFM().crbmodel(
            float(tpr),
            ctp,
            hazescale=float(hazescale),
            hazeloc=hazeloc,
            hazethick=hazethick,
            mixratio=mixratio,
        )
        fmc = fmc.spectrum
    fmc = fmc[ctxt.cleanup] - np.nanmean(fmc[ctxt.cleanup])
    fmc = fmc + np.nanmean(ctxt.tspectrum[ctxt.cleanup])
    cond_G102 = 'HST-WFC3-IR-G102-SCAN' in flt
    fmc[cond_G102] = fmc[cond_G102] + 1e-2 * float(off0)
    return fmc
