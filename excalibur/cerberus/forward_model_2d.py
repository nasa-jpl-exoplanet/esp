# '''cerberus forward_model ds'''

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


def _tp_from_coeffs(
    P,
    shift,  # log-P shift
    T_base=1500,  # background scale
    alpha=0.05,  # background slope
    P_ref=1e-2,  # ref pressure
    a_rise=400,  # rise amplitude
    c_rise=-2.3,  # center of rise (log10 P)
    w_rise=0.25,  # width of rise
    a_fall=-350,  # fall amplitude
    c_fall=-1.9,  # center of fall
    w_fall=0.35,
):  # width of fall
    x = np.log10(P) + shift
    T_bg = T_base * (P / P_ref) ** alpha
    # identical to calcTEA's TP builder
    inv_rise = a_rise * (
        np.tanh((x - c_rise) / w_rise) - np.tanh(-c_rise / w_rise)
    )
    inv_fall = a_fall * (
        np.tanh((x - c_fall) / w_fall) - np.tanh(-c_fall / w_fall)
    )
    return T_bg + inv_rise + inv_fall


log = logging.getLogger(__name__)

# this doesn't change results at all; just needed to avoid undefined-variable pylint
ctxt = ctxtinit()


# ----------- --------------------------------------------------------
# -- CERBERUS MODEL -- -----------------------------------------------
def crbmodel(
    temp,  # can be: scalar T, 10-coef TP vector, or per-layer array
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

    # these used to be default parameters above, but are dangerous-default-values
    # note that these are also defined in cerberus/core/myxsecs()
    #  maybe put them inside runtime/ops.xml to ensure consistency?
    cialist = ['H2-H', 'H2-H2', 'H2-He', 'He-H']
    xmollist = [
        'TIO',
        'H2O',
        'HCN',
        'CO',
        'CO2',
        'NH3',
        'CH4',
        'H2S',
        'PH3',
        'C2H2',
        'OH',
        'O2',
        'O3',
        'SO2',
        'C2H6',
        'C3H8',
        'CH3CHO',
    ]

    ssc = syscore.ssconstants(mks=True)
    pgrid = np.arange(
        np.log(solrad) - Hsmax,
        np.log(solrad) + Hsmax / nlevels,
        Hsmax / (nlevels - 1),
    )
    pgrid = np.exp(pgrid)
    pressure = pgrid[::-1]
    dPoverP = (pressure[1] - pressure[0]) / pressure[0]
    Nz = pressure.size

    if np.isscalar(temp):
        # user passed a single T -> make a flat profile
        Tgrid = np.full(Nz, float(temp))
    else:
        tarr = np.asarray(temp, float)
        if tarr.ndim == 1 and tarr.size == 10:
            # user passed the 10 TP coefficients used for TEA -> rebuild T(P)
            Tgrid = _tp_from_coeffs(pressure, *tarr)
        elif tarr.ndim == 1 and tarr.size == Nz:
            # user passed an explicit T(P) that already matches this grid
            Tgrid = tarr
        else:
            raise ValueError(
                f"'temp' must be a scalar, a 10-element TP vector, or a {Nz}-long T(P) array"
            )
    # A scalar representative for geometry-only pieces (keeps current dz logic)
    Tgeom = float(
        np.mean(Tgrid)
    )  # not correct way to do this, going to need to adjust

    # print('PARAMETERS', temp, cheq['CtoO'], cheq['XtoH'])
    if not mixratio:
        if cheq is None:
            log.warning('neither mixratio nor cheq are defined')
        if chemistry == 'TEA':
            mixratio, fH2, fHe = calcTEA(
                pressure,
                Tgrid,
                C2Or=cheq['CtoO'],
                X2Hr=cheq['XtoH'],
                N2Or=cheq['NtoO'],
            )
        else:
            if chemistry != 'TEC':
                log.warning(
                    '--< ERROR: unknown %s chemistry model! >--', chemistry
                )
            mixratio, fH2, fHe = crbce(
                pressure,
                Tgrid,
                C2Or=cheq['CtoO'],
                X2Hr=cheq['XtoH'],
                N2Or=cheq['NtoO'],
            )
        # mmw likely in amu, may be list/array — make scalar kg for geometry
        mmw_raw, fH2, fHe = getmmw(mixratio, protosolar=False, fH2=fH2, fHe=fHe)
    else:
        mmw_raw, fH2, fHe = getmmw(mixratio)

    # ensure scalar kg for geometry (scale height)
    mmw = float(np.mean(np.asarray(mmw_raw, dtype=float))) * cst.m_p

    if debug:
        print('mmw [amu]:', float(np.mean(np.asarray(mmw_raw, dtype=float))))
        print('mmw [kg] :', mmw)
    # gravity (log10 g in cgs) for scale height
    if isinstance(orbp, dict) and planet in orbp:
        if 'logg' in orbp[planet]:
            logg = float(orbp[planet]['logg'])  # expected: log10(g[cm s^-2])
        elif 'g' in orbp[planet]:
            # fallback if only g in m/s^2 is available
            logg = float(
                np.log10(orbp[planet]['g'] * 100.0)
            )  # m/s^2 -> cm/s^2, then log10
        else:
            raise KeyError(
                f"orbp['{planet}'] must have 'logg' (cgs) or 'g' (m/s^2)."
            )
    else:
        raise KeyError(
            "Parameter 'orbp' must be a dict with a key for the selected 'planet'."
        )

    Hs = (cst.Boltzmann * Tgeom) / (mmw * 1e-2 * (10.0**logg))

    # when the Pressure grid is log-spaced, rdz is a constant
    #  drop dz[] and dzprime[] arrays and just use this constant instead
    rdz = abs(Hs / 2.0 * np.log(1.0 + dPoverP))
    dz = 2 * rdz
    z = dz * np.linspace(0, len(pressure) - 1, len(pressure))

    rho = pressure * 1e5 / (cst.Boltzmann * Tgrid)

    tau, tau_by_molecule, wtau = gettau(
        xsecs,
        qtgrid,
        Tgrid,  # <-- per-layer temperatures
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
        axes[-1].set_ylim(np.nanmin(ax2min), np.nanmax(ax2max))
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
    tau_by_molecule = {}
    # DL ARRAY, Z VERSUS ZPRIME ------------------------------------------------------
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

    # Handle temperature: scalar or per-layer
    tarr = np.asarray(temp, dtype=float)
    layered_T = tarr.ndim == 1 and tarr.size == Nzones
    if layered_T:
        Tlayers = tarr
    else:
        Tlayers = np.full(Nzones, float(tarr))

    # GAS ARRAY, ZPRIME VERSUS WAVELENGTH  -------------------------------------------
    for elem in mixratio:
        val = mixratio[elem]  # expected log10(ppm)

        if np.ndim(val) == 0:
            # single number everywhere
            logppm_profile = float(val) * np.ones(Nzones, dtype=float)
        else:
            # 1-D profile from chemistry (may NOT match Nzones)
            val_arr = np.asarray(val, dtype=float).ravel()
            if val_arr.size == Nzones:
                logppm_profile = val_arr
            else:
                # interpolate the chem profile to the FM grid in log-pressure space
                logP_dst = np.log10(np.asarray(pressure, dtype=float))
                order = np.argsort(logP_dst)
                logP_sorted = logP_dst[order]
                logP_src = np.linspace(
                    logP_sorted.min(), logP_sorted.max(), val_arr.size
                )
                prof_sorted = np.interp(logP_sorted, logP_src, val_arr)
                logppm_profile = np.empty_like(logP_dst)
                logppm_profile[order] = prof_sorted

        # convert log10(ppm) → mole fraction
        mmr = 10.0 ** (logppm_profile - 6.0)  # ppm → fraction

        # --- Cross sections ---
        # NOTE: historical hack in original code:
        #   if not xmollist: HITRAN/HITEMP path
        #   else: EXOMOL path (always taken because xmollist is a non-empty list)
        # We preserve that behavior.
        if not xmollist:
            # HITRAN/HITEMP: evaluate per layer when needed
            if layered_T:
                sigma_layers = []
                for i in range(Nzones):
                    # per-layer pressure and mmr to get correct broadening
                    sigma_i, lsig = absorb(
                        xsecs[elem],
                        qtgrid[elem],
                        float(Tlayers[i]),
                        np.array([pressure[i]]),
                        float(mmr[i]),
                        lbroadening,
                        lshifting,
                        wgrid,
                    )
                    sigma_layers.append(np.asarray(sigma_i).ravel())
                sigma_mat = np.vstack(sigma_layers) * 1e-4  # m^2/mol
            else:
                sigma, lsig = absorb(
                    xsecs[elem],
                    qtgrid[elem],
                    float(Tlayers[0]),
                    pressure,
                    mmr,
                    lbroadening,
                    lshifting,
                    wgrid,
                )
                sigma = np.array(sigma) * 1e-4  # m^2/mol
                sigma_mat = (
                    sigma if sigma.ndim == 2 else np.tile(sigma, (Nzones, 1))
                )
        else:
            # EXOMOL: evaluate σ(T) per layer (temperature dependence only)
            sigma, lsig = getxmolxs(
                Tlayers, xsecs[elem]
            )  # can return (layers×λ)
            sigma = np.array(sigma)
            if sigma.ndim == 1:
                sigma_mat = np.tile(sigma, (Nzones, 1))
            else:
                sigma_mat = sigma
            sigma_mat = sigma_mat * 1e-4  # m^2/mol

        # Clean up any negatives/nans/infs
        sigma_mat = np.where(np.isfinite(sigma_mat), sigma_mat, 0.0)
        sigma_mat[sigma_mat < 0.0] = 0.0

        # Build contribution (layers × λ)
        rho_mat = rho[:, None]
        if np.isscalar(mmr):
            mmr_col = np.full(Nzones, mmr)[:, None]
        else:
            mmr_col = np.asarray(mmr, float).reshape(Nzones, 1)

        contrib = mmr_col * sigma_mat * rho_mat
        tau += contrib
        tau_by_molecule[elem] = contrib

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

        # Evaluate CIA σ(T) per layer
        sigma, lsig = getciaxs(Tlayers, xsecs[cia])  # shape: (layers×λ) or (λ,)
        sigma = np.array(sigma) * 1e-10  # m^5/mol^2
        if sigma.ndim == 1:
            sigma_mat = np.tile(sigma, (Nzones, 1))
        else:
            sigma_mat = sigma

        sigma_mat = np.where(np.isfinite(sigma_mat), sigma_mat, 0.0)
        sigma_mat[sigma_mat < 0.0] = 0.0

        rho2_col = (rho**2)[:, None]

        if np.isscalar(f1):
            f1_arr = np.full(Nzones, f1)
        else:
            f1_arr = np.asarray(f1)
        if np.isscalar(f2):
            f2_arr = np.full(Nzones, f2)
        else:
            f2_arr = np.asarray(f2)
        fprod_col = (f1_arr * f2_arr)[:, None]

        contrib = fprod_col * sigma_mat * rho2_col
        tau += contrib
        tau_by_molecule[cia] = contrib

    # RAYLEIGH ARRAY, ZPRIME VERSUS WAVELENGTH  --------------------------------------
    # NAUS & UBACHS 2000 (temperature-independent form used here)
    slambda0 = 750.0 * 1e-3  # microns
    sray0 = 2.52 * 1e-28 * 1e-4  # m^2/mol
    sigma = sray0 * (wgrid[::-1] / slambda0) ** (-4)

    if sigma.ndim == 1:
        sigma_mat = np.tile(sigma, (Nzones, 1))
    else:
        sigma_mat = sigma

    rho_mat = rho[:, None]

    if np.isscalar(fH2):
        fH2_arr = np.full(Nzones, fH2)
    else:
        fH2_arr = np.asarray(fH2)
        if fH2_arr.shape[0] != Nzones:
            raise ValueError(f"fH2 must have length {Nzones}")
    fH2_mat = fH2_arr[:, None]

    contrib = fH2_mat * sigma_mat * rho_mat
    tau += contrib
    tau_by_molecule['rayleigh'] = contrib

    # HAZE ARRAY, ZPRIME VERSUS WAVELENGTH  ------------------------------------------
    if hzlib is None:
        slambda0 = 750.0 * 1e-3  # microns
        sray0 = 2.52 * 1e-28 * 1e-4  # m^2/mol
        sigma = sray0 * (wgrid[::-1]) ** (hazeslope)
        hazedensity = np.ones(len(z))
        hazecontribution = 10.0**hazescale * sigma * np.array([hazedensity]).T
        tau = tau + hazecontribution
        tau_by_molecule['haze'] = hazecontribution
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

    # Line-of-sight integration
    tau = 2e0 * np.asmatrix(dlarray) * np.asmatrix(tau)

    molecules = list(tau_by_molecule.keys())
    for molecule in molecules:
        tau_by_molecule[molecule] = (
            2e0 * np.asmatrix(dlarray) * np.asmatrix(tau_by_molecule[molecule])
        )

    # Choose a wavelength grid to report back
    wtau = (
        1e4 / lsig
    )  # μm   (lsig is wavenumber grid from last species evaluated)

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
        plt.show()
        pass
    return tau, tau_by_molecule, wtau


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
    (T is expected to be scalar here; per-layer handling is done in gettau)
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
    Accepts scalar T or per-layer T array; returns σ(T) sorted by ν.
    '''
    temps = np.asarray(temp, dtype=float)
    nu = np.array(xsecs['SPLNU'])
    select = np.argsort(nu)
    nu = nu[select]

    if temps.ndim == 0:
        sigma = np.array([thisspl(float(temps)) for thisspl in xsecs['SPL']])
        sigma = sigma[select]
        return sigma, nu
    else:
        sig_layers = []
        for Ti in temps.ravel():
            sigma_i = np.array([thisspl(float(Ti)) for thisspl in xsecs['SPL']])
            sig_layers.append(sigma_i[select])
        sigma_mat = np.vstack(sig_layers)  # layers × ν
        return sigma_mat, nu


# ------------------------- ------------------------------------------
# -- CIA -- ----------------------------------------------------------
def getciaxs(temp, xsecs):
    '''
    G. ROUDIER: Wrapper around CIA Cerberus library
    Accepts scalar T or per-layer T array; returns σ(T) sorted by ν.
    '''
    temps = np.asarray(temp, dtype=float)
    nu = np.array(xsecs['SPLNU'])
    select = np.argsort(nu)
    nu = nu[select]

    if temps.ndim == 0:
        sigma = np.array([thisspl(float(temps)) for thisspl in xsecs['SPL']])
        sigma = sigma[select]
        return sigma, nu
    else:
        sig_layers = []
        for Ti in temps.ravel():
            sigma_i = np.array([thisspl(float(Ti)) for thisspl in xsecs['SPL']])
            sig_layers.append(sigma_i[select])
        sigma_mat = np.vstack(sig_layers)  # layers × ν
        return sigma_mat, nu


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

    if 'T' in ctxt.fixedParams:
        tpr = ctxt.fixedParams['T']
        mdp = crbinputs[0]
        # this extra list[] is needed for the single param case (only metallicity)
        if not isinstance(mdp, list):
            mdp = [mdp]
    else:
        tpr, mdp = crbinputs

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
    fmc += np.average(ctxt.mcmcdat - fmc, weights=1 / ctxt.mcmcsig**2)
    return fmc


def offcerberus(*crbinputs):
    '''
    R.ESTRELA: ADD offsets between STIS filters and STIS and WFC3 filters
    '''
    ctp, hazescale, off0, off1, off2, hazeloc, hazethick, tpr, mdp = crbinputs
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
