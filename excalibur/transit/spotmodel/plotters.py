import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pyplot
from matplotlib.cm import get_cmap, ScalarMappable
from matplotlib.colors import Normalize

from excalibur.util.plotters import save_plot_tosv


def plot_transit_depths(f_spot_array, tempSpot_array, lambdaEff_nm, D_lambda, D_lambda_juststar):

    # Configura os parâmetros visuais padrão
    plt.rcParams['axes.linewidth'] = (
        2  # Aumenta a espessura das bordas dos eixos
    )

    # Caminho para o arquivo de resultados
    # file_path = f"simulation_results_{target}.txt"
    # Carregar os dados (assumindo valores separados por vírgula e com cabeçalho)
    # data = np.loadtxt(file_path, delimiter=',', skiprows=1)
    # f_spot_array, tempSpot_array, lambdaEff_nm, D_lambda = data.T

    # Cria a figura com dois subplots (empilhados verticalmente)
    fig, (ax1, ax2) = plt.subplots(
        nrows=2, ncols=1, figsize=(10, 6), sharex=True
    )

    ###############################################################################
    # Primeiro subplot: Linhas coloridas conforme o filling factor (f_spot)
    ###############################################################################
    # Agrupa as simulações únicas (combinações de f_spot e tempSpot)

    # print('f_spot_array', f_spot_array)
    # print('tempSpot_array', tempSpot_array)

    # Escolhe um colormap e cria o normalizador para f_spot
    cmap_f = get_cmap("winter_r")
    norm_f = Normalize(f_spot_array.min(), f_spot_array.max())

    # Loop para plotar cada simulação
    for i_ff, ff in enumerate(f_spot_array):
        for i_spottemp, spottemp in enumerate(tempSpot_array):
            # Define a cor com base no filling factor
            color = cmap_f(norm_f(ff))
            ax1.plot(
                lambdaEff_nm,
                D_lambda[i_ff, i_spottemp],
                marker='o',
                linestyle='-',
                color=color,
                alpha=1,
            )
    ax1.plot(
        lambdaEff_nm,
        D_lambda_juststar[0,0],
        marker='x',
        linestyle='-',
        color='k',
        alpha=1,
    )

    # Adiciona a colorbar para o filling factor
    sm_f = ScalarMappable(cmap=cmap_f, norm=norm_f)
    sm_f.set_array(f_spot_array)
    cbar1 = fig.colorbar(sm_f, ax=ax1, orientation='vertical', pad=0.02)
    cbar1.set_label('Filling Factor', fontsize=16, labelpad=15)
    cbar1.ax.tick_params(labelsize=16)
    cbar1.ax.yaxis.label.set_rotation(270)
    cbar1.ax.yaxis.label.set_verticalalignment('bottom')

    # Personaliza o primeiro subplot
    ax1.set_title('', fontsize=16)
    ax1.set_xlabel('', fontsize=16)
    ax1.set_ylabel('Transit Depth [ppm]', fontsize=16)
    ax1.tick_params(
        axis="x",
        direction="in",
        labelsize=16,
        width=2,
        length=7,
        pad=3,
        top=True,
    )
    ax1.tick_params(
        axis="y",
        direction="in",
        labelsize=16,
        width=2,
        length=7,
        pad=3,
        right=True,
    )

    ###############################################################################
    # Segundo subplot: Linhas coloridas conforme a spot temperature (T_spot)
    ###############################################################################
    # Filtra os valores válidos de temperatura (não NaN) para a normalização
    valid_tempSpot = tempSpot_array[~np.isnan(tempSpot_array)]
    if len(valid_tempSpot) == 0:
        raise ValueError("All temperature values are NaN. Check the data.")

    # Escolhe um colormap e cria o normalizador para T_spot
    cmap_t = get_cmap("cool_r")
    norm_t = Normalize(valid_tempSpot.min(), valid_tempSpot.max())

    # Loop para plotar cada simulação (desconsiderando as que possuem temperatura NaN)
    for i_ff, ff in enumerate(f_spot_array):
        for i_spottemp, spottemp in enumerate(tempSpot_array):
            if np.isnan(spottemp):
                continue

            # Define a cor com base na spot temperature
            color = cmap_t(norm_t(spottemp))
            ax2.plot(
                lambdaEff_nm,
                D_lambda[i_ff, i_spottemp],
                marker='o',
                linestyle='-',
                color=color,
                alpha=1,
            )
    ax2.plot(
        lambdaEff_nm,
        D_lambda_juststar[0,0],
        marker='x',
        linestyle='-',
        color='k',
        alpha=1,
    )

    # Adiciona a colorbar para a spot temperature
    sm_t = ScalarMappable(cmap=cmap_t, norm=norm_t)
    sm_t.set_array(valid_tempSpot)
    cbar2 = fig.colorbar(sm_t, ax=ax2, orientation='vertical', pad=0.02)
    cbar2.set_label('Spot Temperature [K]', fontsize=16, labelpad=15)
    cbar2.ax.tick_params(labelsize=16)
    cbar2.ax.yaxis.label.set_rotation(270)
    cbar2.ax.yaxis.label.set_verticalalignment('bottom')

    # Personaliza o segundo subplot
    ax2.set_title('', fontsize=16)
    ax2.set_xlabel('Wavelength [nm]', fontsize=16)
    ax2.set_ylabel('Transit Depth [ppm]', fontsize=16)
    ax2.tick_params(
        axis="x",
        direction="in",
        labelsize=16,
        width=2,
        length=7,
        pad=3,
        top=True,
    )
    ax2.tick_params(
        axis="y",
        direction="in",
        labelsize=16,
        width=2,
        length=7,
        pad=3,
        right=True,
    )

    # Ajusta o layout e exibe a figura
    plt.tight_layout()
    # plt.show()
    savedplot = save_plot_tosv(fig)
    plt.close(fig)
    return savedplot


def plot_lightcurves(
    num_wavelengths,
    stack_tempoHoras,
    stack_curvaLuz,
    lambdaEff_nm,
    tempoTransito,
):

    fig = plt.figure(figsize=(6, 5))

    count4 = 0
    palette = pyplot.cm.cool_r(np.linspace(0, 1, num_wavelengths))
    count_palette = num_wavelengths - 1

    while count4 < num_wavelengths:
        pyplot.plot(
            stack_tempoHoras[count4],
            stack_curvaLuz[count4],
            label=int(lambdaEff_nm[count4]),
            color=palette[count_palette],
            linewidth=1.5,
        )
        count4 += 1
        count_palette -= 1

    min_D = np.min(stack_curvaLuz)

    pyplot.axis([-tempoTransito / 5, tempoTransito / 5, min_D, 1.00005])
    pyplot.xlabel('Time from transit center (hr)', fontsize=16)
    pyplot.ylabel('Relative flux', fontsize=16)
    pyplot.tick_params(
        axis="x", direction="in", labelsize=16, length=7, width=2, top=True
    )
    pyplot.tick_params(
        axis="y",
        direction="in",
        labelsize=16,
        length=7,
        width=2,
        right=True,
    )
    pyplot.rcParams['axes.linewidth'] = (
        2  # Aumenta a espessura das bordas dos eixos
    )
    pyplot.tight_layout()
    # pyplot.show()
    savedplot = save_plot_tosv(fig)
    plt.close(fig)
    return savedplot
