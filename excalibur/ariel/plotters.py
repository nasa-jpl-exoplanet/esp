'''ariel core ds'''

# Heritage code shame:
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments,too-many-locals,too-many-positional-arguments,too-many-branches,too-many-statements

# -- IMPORTS -- ------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from excalibur.util.plotters import add_scale_height_labels, save_plot_tosv

# ------------------------- ------------------------------------------


def plot_spectrum(
    target,
    planet_letter,
    tier,
    visits,
    wavelength_um,
    fluxDepth,
    fluxDepth_observed,
    uncertainties_percent,
    molecules,
    fluxDepth_by_molecule,
    Hsscaling,
    plottype='Ariel',
    verbose=False,
):
    # PLOT THE SPECTRA
    myfig, ax = plt.subplots(figsize=(8, 4))
    if plottype == 'Ariel':
        tierlabel = 'Tier-' + str(tier) + ' '
    else:
        tierlabel = ''
    plt.title(
        plottype
        + ' simulation : '
        + target
        + ' '
        + planet_letter
        + ' : '
        + tierlabel
        + str(visits)
        + ' visits',
        fontsize=16,
    )
    plt.xlabel(str('Wavelength [$\\mu m$]'), fontsize=14)
    plt.ylabel(str('$(R_p/R_*)^2$ [%]'), fontsize=14)

    # plot the true (model) spectrum - raw
    #  (this is for testing higher resolution cross-sections
    #   but currently raw and rebin are the same thing)
    # plt.plot(
    #    wavelength_um_raw,
    #    fluxDepth_raw,
    #    color='palegoldenrod',
    #    ls='-',
    #    lw=1,
    #    zorder=1,
    #    label='truth raw',
    # )
    # plot the true (model) spectrum - binned
    plt.plot(
        wavelength_um,
        fluxDepth,
        color='k',
        ls='-',
        lw=1,
        zorder=3,
        # label='truth binned',
        label='truth',
    )
    yrange = plt.ylim()
    # plot the simulated data points
    plt.scatter(
        wavelength_um,
        fluxDepth_observed,
        marker='o',
        s=20,
        color='None',
        edgecolor='k',
        zorder=4,
        label='sim data',
    )
    plt.errorbar(
        wavelength_um,
        fluxDepth_observed,
        yerr=uncertainties_percent,
        linestyle='None',
        lw=0.2,
        color='grey',
        zorder=2,
    )
    plt.ylim(yrange)
    if plottype == 'Ariel':
        plt.xlim(0.0, 8.0)
    else:
        xlims = plt.xlim()
        plt.xlim(0.0, xlims[1])
    xlims = plt.xlim()

    # calculate the baseline and total range
    #  use this to omit molecules that are not contributing
    baseline = 666
    maxdepth = -666
    for imole, molecule in enumerate(molecules):
        baseline = np.min(np.append(baseline, fluxDepth_by_molecule[molecule]))
        maxdepth = np.max(np.append(maxdepth, fluxDepth_by_molecule[molecule]))
    negligible_molecules = ''
    negligible_molecules_more = ''
    Nnegligible = 0
    if plottype == 'Ariel':
        negligibletextloc = 6.8
    else:
        negligibletextloc = xlims[1] * 0.85
    for imole, molecule in enumerate(molecules):
        colorlist = [
            'red',
            'orange',
            'palegreen',
            'lightseagreen',
            'blueviolet',
            'fuchsia',
        ]
        stylelist = ['--'] * 6 + [':'] * 6 + ['-.'] * 6
        feature_strength = (
            np.max(fluxDepth_by_molecule[molecule]) - baseline
        ) / (maxdepth - baseline)
        # cut off anything less than 1% of maximum contribution
        if feature_strength < 0.01:
            # print('  dropped:', molecule, feature_strength)
            Nnegligible += 1
            if Nnegligible < 10:
                negligible_molecules += ' ' + molecule
            else:
                negligible_molecules_more += ' ' + molecule
        else:
            # print(' OK:', molecule, feature_strength)
            plt.plot(
                wavelength_um,
                fluxDepth_by_molecule[molecule],
                color=colorlist[imole % len(colorlist)],
                ls=stylelist[imole % len(stylelist)],
                lw=1,
                zorder=2,
                label=molecule,
            )
        extra = (maxdepth - baseline) / 13
        plt.ylim((baseline - extra, maxdepth + extra))
        yrange = plt.ylim()
        plt.text(
            negligibletextloc,
            yrange[0] + (yrange[1] - yrange[0]) * (-0.13),
            'negligible contribution:',
            fontsize=8,
        )
        # there's formating problems when too many negligible molecules
        # better to split it up over two lines
        plt.text(
            negligibletextloc,
            yrange[0] + (yrange[1] - yrange[0]) * (-0.18),
            negligible_molecules,
            fontsize=8,
        )
        plt.text(
            negligibletextloc,
            yrange[0] + (yrange[1] - yrange[0]) * (-0.23),
            negligible_molecules_more,
            fontsize=8,
        )
    if plottype == 'Ariel':
        plt.xlim(0.0, 8.0)
    else:
        # plt.xlim(0.3, 1.1)
        xlims = plt.xlim()
        plt.xlim(0.0, xlims[1])

    plt.legend(loc='center left', bbox_to_anchor=(1.16, 0.48))

    # add a scale-height-normalized flux scale on the right axis
    # print('H scaling for this plot (%):',Hsscaling*100)
    add_scale_height_labels({'Hs': [Hsscaling]}, 1.0e-2 * fluxDepth, ax, myfig)

    # option to save plot to disk
    # plot_dir = (
    #     excalibur.context['data_dir'] + '/ariel/savedplots'
    # )
    # if not os.path.exists(plot_dir):
    #    os.mkdir(plot_dir)
    # plt.savefig(plot_dir + '/ariel_' +
    #             target + '_' + planet_letter + '.png')

    savedFigure = save_plot_tosv(myfig)
    if verbose:
        plt.show()
    plt.close(myfig)
    return savedFigure


# ------------------------- ------------------------------------------


def plot_spectrum_topmolecules(
    target,
    planet_letter,
    tier,
    visits,
    wavelength_um,
    wavelengthedge_low,
    wavelengthedge_high,
    fluxDepth,
    fluxDepth_observed,
    uncertainties_percent,
    molecules,
    fluxDepth_by_molecule,
    Hsscaling,
    plottype='Ariel',
    verbose=False,
):
    # PLOT THE SPECTRA
    myfig, ax = plt.subplots(figsize=(8, 4))
    if plottype == 'Ariel':
        tierlabel = 'Tier-' + str(tier) + ' '
    else:
        tierlabel = ''
    plt.title(
        plottype
        + ' simulation : '
        + target
        + ' '
        + planet_letter
        + ' : '
        + tierlabel
        + str(visits)
        + ' visits',
        fontsize=16,
    )
    plt.xlabel(str('Wavelength [$\\mu m$]'), fontsize=14)
    plt.ylabel(str('$(R_p/R_*)^2$ [%]'), fontsize=14)

    # plot the true (model) spectrum - binned
    plt.plot(
        wavelength_um,
        fluxDepth,
        color='k',
        ls='-',
        lw=1,
        zorder=3,
        label='truth',
    )
    yrange = plt.ylim()
    # plot the simulated data points
    plt.scatter(
        wavelength_um,
        fluxDepth_observed,
        marker='o',
        s=20,
        color='None',
        edgecolor='k',
        zorder=4,
        label='sim data',
    )
    plt.errorbar(
        wavelength_um,
        fluxDepth_observed,
        yerr=uncertainties_percent,
        linestyle='None',
        lw=0.2,
        color='grey',
        zorder=2,
    )
    plt.ylim(yrange)

    colorlist = [
        'red',
        'palegreen',
        'lightseagreen',
        'blueviolet',
        'fuchsia',
    ]
    moleculeColorMatch = {
        'H2O': 'steelblue',
        'CO': 'goldenrod',
        'CO2': 'gold',
        'CH4': 'darkorange',
        'haze': 'deepskyblue',
    }

    dominantMolecule_byWavelength = []
    # print('molecules', molecules)
    for iwave in range(len(wavelength_um)):
        # print('wave', iwave, wavelength_um[iwave])
        # print('  len check', len(fluxDepth_by_molecule['H2O']))
        fluxesThisWavelength = np.array(
            [fluxDepth_by_molecule[molecule][iwave] for molecule in molecules]
        )
        # print(' fluxes', fluxesThisWavelength)
        wheremax = np.where(
            fluxesThisWavelength == np.max(fluxesThisWavelength)
        )[0][0]
        # print(' where max', wheremax)
        dominantMolecule_byWavelength.append(molecules[wheremax])
    # print('dominantMolecule_byWavelength', dominantMolecule_byWavelength)

    moleculeSpacing = (yrange[1] - yrange[0]) / 12
    nextMoleculeYpos = yrange[1] + moleculeSpacing
    moleculesFound = []
    moleculeYpos = []
    for iwave in range(len(wavelength_um)):
        thisMolecule = dominantMolecule_byWavelength[iwave]
        if thisMolecule not in moleculesFound:
            # print(' found a new molecule',thisMolecule)
            moleculesFound.append(thisMolecule)
            moleculeYpos.append(nextMoleculeYpos)
            nextMoleculeYpos += moleculeSpacing

        imole = moleculesFound.index(thisMolecule)
        moleculeColor = moleculeColorMatch.get(
            thisMolecule, colorlist[imole % len(colorlist)]
        )
        # print('thismolecule,color:', thisMolecule, moleculeColor)
        plt.plot(
            [wavelengthedge_low[iwave], wavelengthedge_high[iwave]],
            [moleculeYpos[imole], moleculeYpos[imole]],
            color=moleculeColor,
            ls='-',
            lw=2,
            zorder=2,
        )
        moleculetextloc = wavelengthedge_high[-1] * 1.03
        plt.text(
            moleculetextloc,
            moleculeYpos[imole],
            thisMolecule,
            color=moleculeColor,
            fontsize=12,
            verticalalignment='center',
        )
    # if plottype == 'Ariel':
    #     plt.xlim(0, 8.5)
    xlims = plt.xlim()
    plt.xlim(0, xlims[1]*1.08)
    plt.ylim(yrange[0], nextMoleculeYpos)
    # plt.legend(loc='center left', bbox_to_anchor=(1.16, 0.48))

    # add a scale-height-normalized flux scale on the right axis
    # print('H scaling for this plot (%):',Hsscaling*100)
    add_scale_height_labels({'Hs': [Hsscaling]}, 1.0e-2 * fluxDepth, ax, myfig)

    savedFigure = save_plot_tosv(myfig)
    if verbose:
        plt.show()
    plt.close(myfig)
    return savedFigure


# ------------------------- ------------------------------------------


def plot_depthprobed(
    target,
    planet_letter,
    model_params,
    wavelength_um,
    pressure,
    opacityProfiles,
    # molecules,
    verbose=False,
):

    opacityProfiles = np.array(opacityProfiles)
    npressure, nwave = opacityProfiles.shape
    # print('# wavelengths, pressures', nwave, npressure)

    # convert the local opacity to a transmission map (as func of wave)
    # transmission is 1 at top of atmos, 0 at the bottom
    throughput = np.zeros((npressure + 1, nwave))
    throughput[-1, :] = 1
    for ipressure in range(npressure)[::-1]:  # start at high index = atmos top
        throughput[ipressure, :] = throughput[ipressure + 1, :] * np.exp(
            -opacityProfiles[ipressure, :]
        )
    # cut off that top buffer edge
    throughput = throughput[:-1, :]

    myfig, _ = plt.subplots(figsize=(8, 4))
    myfig.subplots_adjust(top=0.92, bottom=0.13, left=0.09, right=0.98)
    plt.title(
        'Atmospheric depth probed : ' + target + ' ' + planet_letter,
        fontsize=16,
    )
    plt.xlabel('Wavelength [$\\mu m$]', fontsize=14)
    plt.ylabel('log(Pressure) [bar]', fontsize=14)

    plt.contourf(
        wavelength_um[::-1],  # cerberus fm calculates on a backward wave grid
        np.log10(pressure),
        throughput,
        zorder=2,
        colors='grey',
        levels=np.array([0.159, 0.841]),
    )
    plt.contour(
        wavelength_um[::-1],  # cerberus fm calculates on a backward wave grid
        np.log10(pressure),
        throughput,
        zorder=3,
        colors='k',
        levels=np.array([0.5]),
    )

    # reverse the y-axis so 10 bar is at the bottom
    ylims = plt.ylim()
    plt.ylim(ylims[1], ylims[0])

    if 'CTP' in model_params:
        # dashed line showing the cloud deck
        # print('CTP:', model_params['CTP'])
        plt.plot(
            wavelength_um,
            [model_params['CTP']] * len(wavelength_um),
            'k:',
            zorder=3,
        )
        plt.text(
            6.0,
            model_params['CTP'] + (ylims[0] - ylims[1]) / 25.0,
            'top of cloud deck',
            fontsize=8,
        )

    savedFigure = save_plot_tosv(myfig)
    if verbose:
        plt.show()
    plt.close(myfig)
    return savedFigure


# ----------------- --------------------------------------------


def plot_vertical_profiles(
    target,
    planet_letter,
    molecule_profiles,
    pressure,
    verbose=False,
):
    colorlist = [
        'red',
        'orange',
        'palegreen',
        'lightseagreen',
        'blueviolet',
        'fuchsia',
    ]
    stylelist = ['-'] * 6 + ['--'] * 6 + [':'] * 6 + ['-.'] * 6
    # set mixing ratio boundaries for main molecules, trace molecules, and absolute cutoff
    floor_ppm = 1e-12
    split_ppm = 1e-6
    xmax_ppm = 1e6

    myfig, (ax_main, ax_trace) = plt.subplots(
        2, 1, sharey=True, gridspec_kw={"height_ratios": [3, 1]}, figsize=(8, 4)
    )
    for imole, molecule in enumerate(molecule_profiles.keys()):
        x = 10 ** molecule_profiles[molecule]
        x = x.astype(float)
        # main molecules, shown in main plot
        x_main = x.copy()
        x_main[(x_main < split_ppm)] = np.nan
        ax_main.plot(
            x_main,
            pressure,
            label=molecule,
            color=colorlist[imole % len(colorlist)],
            ls=stylelist[imole % len(stylelist)],
        )
        # trace molecules shown in smaller snippet
        x_trace = x.copy()
        x_trace[(x_trace >= split_ppm) | (x_trace < floor_ppm)] = np.nan
        ax_trace.plot(
            x_trace,
            pressure,
            color=colorlist[imole % len(colorlist)],
            ls=stylelist[imole % len(stylelist)],
        )
    ax_main.set_yscale('log')
    ax_trace.set_yscale('log')
    ax_main.set_xscale('log')
    ax_trace.set_xscale('log')
    ax_main.set_xlim(split_ppm, xmax_ppm)
    ax_trace.set_xlim(floor_ppm, split_ppm)
    pmin = float(np.nanmin(pressure))
    pmax = float(np.nanmax(pressure))
    ax_main.set_ylim(pmax, pmin)
    # 2 column legend if more than 10 molecules
    if len(molecule_profiles.keys()) > 10:
        ncols = 2
    else:
        ncols = 1
    # ax_main.legend(bbox_to_anchor=(1.0, 1.0), ncol=ncols)
    ax_main.legend(loc='center left', bbox_to_anchor=(1.16, 0.48), ncol=ncols)
    ax_trace.set_xlabel('Volume Mixing Ratio [ppm]', fontsize=14)
    ax_main.set_ylabel('Pressure [bar]', fontsize=14)
    ax_main.set_title(
        'Vertical Chemical Profile : ' + target + ' ' + planet_letter,
        fontsize=16,
    )
    plt.tight_layout()
    savedFigure = save_plot_tosv(myfig)
    if verbose:
        plt.show()
    plt.close(myfig)
    return savedFigure


# --------------------------------------------------------------------
