'''ariel core ds'''

# Heritage code shame:
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments,too-many-locals

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
    wavelength_um_rebin,
    fluxDepth_rebin,
    fluxDepth_observed,
    uncertainties_percent,
    molecules,
    fluxDepth_by_molecule_rebin,
    Hsscaling,
    atmosModel,
    verbose=False,
):
    # PLOT THE SPECTRA
    myfig, ax = plt.subplots(figsize=(8, 4))
    plt.title(
        'Ariel simulation : '
        + target
        + ' '
        + planet_letter
        + ' : Tier-'
        + str(tier)
        + ' '
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
    #    wavelength_um,
    #    fluxDepth,
    #    color='palegoldenrod',
    #    ls='-',
    #    lw=1,
    #    zorder=1,
    #    label='truth raw',
    # )
    # plot the true (model) spectrum - binned
    plt.plot(
        wavelength_um_rebin,
        fluxDepth_rebin,
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
        wavelength_um_rebin,
        fluxDepth_observed,
        marker='o',
        s=20,
        color='None',
        edgecolor='k',
        zorder=4,
        label='sim data',
    )
    plt.errorbar(
        wavelength_um_rebin,
        fluxDepth_observed,
        yerr=uncertainties_percent,
        linestyle='None',
        lw=0.2,
        color='grey',
        zorder=2,
    )
    plt.ylim(yrange)

    # calculate the baseline and total range
    #  use this to omit molecules that are not contributing
    baseline = 666
    maxdepth = -666
    for imole, molecule in enumerate(molecules):
        baseline = np.min(
            np.append(baseline, fluxDepth_by_molecule_rebin[molecule])
        )
        maxdepth = np.max(
            np.append(maxdepth, fluxDepth_by_molecule_rebin[molecule])
        )
    negligible_molecules = ''
    negligible_molecules_more = ''
    Nnegligible = 0
    for imole, molecule in enumerate(molecules):
        colorlist = [
            'red',
            'orange',
            'palegreen',
            'lightseagreen',
            'blueviolet',
            'fuchsia',
        ]
        stylelist = [
            '--',
            '--',
            '--',
            '--',
            '--',
            '--',
            ':',
            ':',
            ':',
            ':',
            ':',
            ':',
            '-.',
            '-.',
            '-.',
            '-.',
            '-.',
            '-.',
        ]
        feature_strength = (
            np.max(fluxDepth_by_molecule_rebin[molecule]) - baseline
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
                wavelength_um_rebin,
                fluxDepth_by_molecule_rebin[molecule],
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
            6.9,
            yrange[0] + (yrange[1] - yrange[0]) * (-0.13),
            'negligible contribution:',
            fontsize=8,
        )
        # there's formating problems when too many negligible molecules
        # better to split it up over two lines
        plt.text(
            6.8,
            yrange[0] + (yrange[1] - yrange[0]) * (-0.18),
            negligible_molecules,
            fontsize=8,
        )
        plt.text(
            6.8,
            yrange[0] + (yrange[1] - yrange[0]) * (-0.23),
            negligible_molecules_more,
            fontsize=8,
        )
    plt.xlim(0.0, 8.0)
    plt.legend(loc='center left', bbox_to_anchor=(1.16, 0.48))

    # add a scale-height-normalized flux scale on the right axis
    # print('H scaling for this plot (%):',Hsscaling*100)
    add_scale_height_labels(
        {'Hs': [Hsscaling]}, 1.0e-2 * fluxDepth_rebin, ax, myfig
    )

    # option to save plot to disk
    # plot_dir = (
    #     excalibur.context['data_dir'] + '/ariel/savedplots'
    # )
    # if not os.path.exists(plot_dir):
    #    os.mkdir(plot_dir)
    # plt.savefig(plot_dir +
    #             '/ariel_' + atmosModel + 'Atmos_' +
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
    wavelength_um_rebin,
    fluxDepth_rebin,
    fluxDepth_observed,
    uncertainties_percent,
    molecules,
    fluxDepth_by_molecule_rebin,
    Hsscaling,
    atmosModel,
    verbose=False,
):
    # PLOT THE SPECTRA
    myfig, ax = plt.subplots(figsize=(8, 4))
    plt.title(
        'Ariel simulation : '
        + target
        + ' '
        + planet_letter
        + ' : Tier-'
        + str(tier)
        + ' '
        + str(visits)
        + ' visits',
        fontsize=16,
    )
    plt.xlabel(str('Wavelength [$\\mu m$]'), fontsize=14)
    plt.ylabel(str('$(R_p/R_*)^2$ [%]'), fontsize=14)

    # plot the true (model) spectrum - binned
    plt.plot(
        wavelength_um_rebin,
        fluxDepth_rebin,
        color='k',
        ls='-',
        lw=1,
        zorder=3,
        label='truth',
    )
    yrange = plt.ylim()
    # plot the simulated data points
    plt.scatter(
        wavelength_um_rebin,
        fluxDepth_observed,
        marker='o',
        s=20,
        color='None',
        edgecolor='k',
        zorder=4,
        label='sim data',
    )
    plt.errorbar(
        wavelength_um_rebin,
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
        'orange',
        'palegreen',
        'lightseagreen',
        'blueviolet',
        'fuchsia',
    ]

    dominantMolecule_byWavelength = []
    print('molecules', molecules)
    for iwave in range(len(wavelength_um_rebin)):
        print('wave', iwave, wavelength_um_rebin[iwave])
        print('  len check', len(fluxDepth_by_molecule_rebin['H2O']))
        fluxesThisWavelength = np.array(
            [
                fluxDepth_by_molecule_rebin[molecule][iwave]
                for molecule in molecules
            ]
        )
        print(' fluxes', fluxesThisWavelength)
        wheremax = np.where(
            fluxesThisWavelength == np.max(fluxesThisWavelength)
        )[0]
        print(' where max', wheremax)
        # print(' index',fluxesThisWavelength.index(wheremax))
        dominantMolecule_byWavelength.append(np.array(molecules)[wheremax])
    print('dominantMolecule_byWavelength', dominantMolecule_byWavelength)

    # feature_strength = fluxDepth_by_molecule_rebin[molecule]
    # plt.plot(
    #    wavelength_um_rebin,
    #    fluxDepth_by_molecule_rebin[molecule],
    #    color=colorlist[imole % len(colorlist)],
    #    ls=stylelist[imole % len(stylelist)],
    #    lw=1,
    #    zorder=2,
    #    label=molecule,
    # )

    # extra = (maxdepth - baseline) / 13
    # plt.ylim((baseline - extra, maxdepth + extra))
    # yrange = plt.ylim()
    # plt.text(
    #    6.8,
    #    yrange[0] + (yrange[1] - yrange[0]) * (-0.18),
    #    negligible_molecules,
    #    fontsize=8,
    # )
    # plt.xlim(0.0, 8.0)
    plt.legend(loc='center left', bbox_to_anchor=(1.16, 0.48))

    # add a scale-height-normalized flux scale on the right axis
    # print('H scaling for this plot (%):',Hsscaling*100)
    add_scale_height_labels(
        {'Hs': [Hsscaling]}, 1.0e-2 * fluxDepth_rebin, ax, myfig
    )

    # option to save plot to disk
    # plot_dir = (
    #     excalibur.context['data_dir'] + '/ariel/savedplots'
    # )
    # if not os.path.exists(plot_dir):
    #    os.mkdir(plot_dir)
    # plt.savefig(plot_dir +
    #             '/ariel_' + atmosModel + 'Atmos_' +
    #             target + '_' + planet_letter + '.png')

    savedFigure = save_plot_tosv(myfig)
    if verbose:
        plt.show()
    plt.close(myfig)
    return savedFigure


# ------------------------- ------------------------------------------


def plot_depthprobed(
    target,
    planet_letter,
    tier,
    visits,
    wavelength_um,
    pressure,
    opticalDepthProfiles,
    # molecules,
    verbose=False,
):
    # PLOT THE SPECTRA
    myfig, ax = plt.subplots(figsize=(8, 4))
    plt.title(
        'Ariel simulation : '
        + target
        + ' '
        + planet_letter
        + ' : Tier-'
        + str(tier)
        + ' '
        + str(visits)
        + ' visits',
        fontsize=16,
    )
    plt.xlabel('Wavelength [$\\mu m$]', fontsize=14)
    plt.ylabel('Pressure [bar]', fontsize=14)

    plt.contour(
        wavelength_um,
        pressure,
        opticalDepthProfiles,
    )

    savedFigure = save_plot_tosv(myfig)
    if verbose:
        plt.show()
    plt.close(myfig)
    return savedFigure


# ------------------------- ------------------------------------------
