'''transit Spotmodel ds'''

# Viktor Y. D. Sumida e-mail: viktorsumida@icloud.com
# PhD student in Astrophysics at Mackenzie University
# Contact via e-mail: viktorsumida@icloud.com

import numpy as np

from excalibur.transit.spotmodel.star import Star
from excalibur.transit.spotmodel.eclipse_nv1 import Eclipse
from excalibur.transit.spotmodel.Planeta import Planeta
from excalibur.transit.spotmodel.plotters import plot_lightcurves


def planck(wavelength, Temp):

    h_Planck = 6.626e-34
    c_Light = 2.998e8
    k_Boltzmann = 1.38e-23

    aux = h_Planck * c_Light / (wavelength * k_Boltzmann * Temp)
    intensity = (2 * h_Planck * (c_Light**2)) / (
        (wavelength**5) * (np.exp(aux) - 1.0)
    )
    return intensity


# Variáveis globais para armazenar o baseline em memória
BASELINE_WAVE = None
BASELINE_DLAMDA = None


class SpotModel:

    def __init__(self, params):
        self.target = params['target']
        self.num_wavelengths = params['num_wavelengths']
        self.c1 = params['c1']
        self.c2 = params['c2']
        self.c3 = params['c3']
        self.c4 = params['c4']
        self.lambdaEff = params['lambdaEff']
        self.raio = params['raio']
        self.intensidadeMaxima = params['intensidadeMaxima']
        self.tamanhoMatriz = params['tamanhoMatriz']
        self.raioStar = params['raioStar']
        self.ecc = params['ecc']
        self.anom = params['anom']
        self.tempStar = params['tempStar']
        self.starspots = params['starspots']
        self.quantidade = params['quantidade']
        # Mantenha lat e longt como arrays para suportar múltiplas manchas
        self.lat = params['lat']
        self.longt = params['longt']
        self.r = params['r']
        self.semiEixoUA = params['semiEixoUA']
        self.massStar = params['massStar']
        self.plot_anim = params['plot_anim']
        self.periodo = params['periodo']
        self.anguloInclinacao = params['anguloInclinacao']
        self.raioPlanetaRj = params['raioPlanetaRj']
        self.plot_graph = params['plot_graph']
        self.plot_star = params['plot_star']
        self.tempSpot = params['tempSpot']

        raioSun = self.raioStar
        raioStar = self.raioStar * 696340  # converting from solar radii to km

        # Cria a estrela
        lambdaEfetivo = 0.0028976 / self.tempStar
        intensidadeEstrelaLambdaMaximo = planck(lambdaEfetivo, self.tempStar)

        raioPlanJup = self.raioPlanetaRj

        # print('stellar effective temperature:', self.tempStar)
        # print('active region effective temperature:', self.tempSpot)

        intensidadeMancha = np.zeros(self.num_wavelengths)
        intensidadeManchaNormalizada = np.zeros(self.num_wavelengths)
        intensidadeManchaRazao = np.zeros(self.num_wavelengths)

        # Remova a conversão de lat e longt para float, mantendo-os como arrays:
        # self.r já é um escalar e não precisa ser modificado.

        massPlaneta = 0  # EXCALIBUR does not handle planet mass
        planeta_ = Planeta(
            self.semiEixoUA,
            raioPlanJup,
            self.periodo,
            self.anguloInclinacao,
            self.ecc,
            self.anom,
            raioStar,
            massPlaneta,
        )

        stack_curvaLuz = None
        stack_tempoHoras = None

        count3 = 0
        while count3 < self.num_wavelengths:
            # print(
            #    'Starting the simulation '
            #    + str(count3 + 1)
            #    + ' of '
            #    + str(self.num_wavelengths)
            # )

            intensidadeEstrelaLambda = planck(
                self.lambdaEff[count3] * 1.0e-6, self.tempStar
            )
            intensidadeEstrelaLambdaNormalizada = (
                float(self.intensidadeMaxima)
                * intensidadeEstrelaLambda
                / intensidadeEstrelaLambdaMaximo
            )
            Nx = self.tamanhoMatriz
            Ny = self.tamanhoMatriz

            estrela_ = Star(
                self.raio,
                raioSun,
                intensidadeEstrelaLambdaNormalizada,
                self.c1[count3],
                self.c2[count3],
                self.c3[count3],
                self.c4[count3],
                self.tamanhoMatriz,
            )

            estrela_.resetManchas()

            if not self.starspots:
                # spotless star
                print('Spotless star created')
                estrela_matriz = estrela_.getMatrizEstrela()

            elif self.starspots:

                mancha = [None] * self.quantidade

                intensidadeMancha[count3] = planck(
                    self.lambdaEff[count3] * 1.0e-6, self.tempSpot
                )
                intensidadeManchaNormalizada[count3] = (
                    intensidadeMancha[count3]
                    * intensidadeEstrelaLambdaNormalizada
                    / intensidadeEstrelaLambda
                )
                intensidadeManchaRazao[count3] = (
                    intensidadeManchaNormalizada[count3]
                    / intensidadeEstrelaLambdaNormalizada
                )
                count2 = 0
                while count2 < self.quantidade:
                    # Usa os valores individuais de lat e longt para cada mancha
                    mancha[count2] = Star.Mancha(
                        intensidadeManchaRazao[count3],
                        self.r,
                        self.lat[count2],
                        self.longt[count2],
                    )
                    #  Geoff: THIS IS NOT USED
                    # raioMancha = self.r * raioStar

                    #  Geoff: THIS IS NOT USED
                    # area = np.pi * (raioMancha**2)

                    # print('Number of spots:', self.quantidade)
                    # print('each spot radius:', self.r)
                    # print('spot radius (sum):', self.r * quantidade)
                    # print('each spot coverage area ff (sum):', self.r ** 2)
                    # print('spot coverage area ff (sum):', (self.r ** 2) * quantidade)

                    estrela_.addMancha(mancha[count2])
                    count2 += 1

                estrela_.criaEstrelaManchada()
                estrela_matriz = estrela_.getMatrizEstrela()

            if self.plot_star:
                estrela_.Plotar(self.tamanhoMatriz, estrela_.getMatrizEstrela())

            # Eclipse (sem CME ou luas)
            eclipse_ = Eclipse(Nx, Ny, self.raio, estrela_, planeta_)
            eclipse_.setEstrela(estrela_matriz)
            eclipse_.criarEclipse(self.plot_anim, self.plot_graph)

            tempoHoras = 1
            eclipse_.geraTempoHoras(tempoHoras)

            tempoTransito = eclipse_.getTempoTransito()
            curvaLuz = eclipse_.getCurvaLuz()
            tempoHoras = eclipse_.getTempoHoras()

            if count3 == 0:
                tamanho_curva = len(curvaLuz)
                tamanho_tempo = len(tempoHoras)
                stack_curvaLuz = np.zeros((self.num_wavelengths, tamanho_curva))
                stack_tempoHoras = np.zeros(
                    (self.num_wavelengths, tamanho_tempo)
                )

            stack_curvaLuz[count3, :] = eclipse_.getCurvaLuz()
            stack_tempoHoras[count3, :] = eclipse_.getTempoHoras()

            count3 += 1
            # print(f"Finishing iteration {count3}")

        #######################################################################
        # Plotting the light curves
        #######################################################################

        lambdaEff_nm = [0.0] * self.num_wavelengths
        D_lambda = [0.0] * self.num_wavelengths
        count4 = 0
        while count4 < self.num_wavelengths:
            lambdaEff_nm[count4] = self.lambdaEff[count4] * 1000
            index_midTrans = int(np.floor(len(stack_curvaLuz[count4]) / 2))
            D_lambda_mid = stack_curvaLuz[count4]
            # Transit depth in ppm
            D_lambda[count4] = (1.0 - D_lambda_mid[index_midTrans]) * 1.0e6
            count4 += 1
        D_lambda = np.array(D_lambda)
        lambdaEff_nm = np.array(lambdaEff_nm)

        # asdf
        self.plot_lightcurves = plot_lightcurves(
            self.num_wavelengths,
            stack_tempoHoras,
            stack_curvaLuz,
            lambdaEff_nm,
            tempoTransito,
        )

        def salvar_dados_simulacao(f_spot, tempSpot, lambdaEff_nm, D_lambda):
            f_spot_array = np.full(len(lambdaEff_nm), f_spot)
            tempSpot_array = np.full(len(lambdaEff_nm), tempSpot)

            # print('f_spot,tempSpot,wavelength,D_lambda')
            # print(np.transpose(
            #    [f_spot_array, tempSpot_array, lambdaEff_nm, D_lambda]))

        salvar_dados_simulacao(
            self.r**2 * self.quantidade, self.tempSpot, lambdaEff_nm, D_lambda
        )

        # save the results (only the depth is used; the rest is redundant)
        self.ff = self.r**2 * self.quantidade
        self.T = self.tempSpot
        self.wavearray = lambdaEff_nm
        self.depth = D_lambda
