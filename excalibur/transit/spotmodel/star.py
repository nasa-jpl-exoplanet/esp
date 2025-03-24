import math
import matplotlib.pyplot as plt
import numpy as np
from typing import List
# from numba import njit, prange


class Star:
    '''
    A classe Estrela recebe como objeto o raio, intensidade máxima, coeficientes de escurecimento de limbo.
    A estrela é construída em uma matriz de tamanho default (por exemplo, 856).

    Parâmetros:
      - raio: O raio da estrela em pixels.
      - raioSun: O raio da estrela em unidades de Rsun.
      - intensidadeMaxima: Intensidade no centro da estrela.
      - coeficienteHum, coeficienteDois, coeficienteTres, coeficienteQuatro: coeficientes do modelo de limbo.
      - tamanhoMatriz: Tamanho da matriz em que a estrela será construída.
      - A estrela construída com os coeficientes de limbo é armazenada em 'estrelaMatriz'.
    '''

    def __init__(
        self,
        raio,
        raioSun,
        intensidadeMaxima,
        coeficienteHum,
        coeficienteDois,
        coeficienteTres,
        coeficienteQuatro,
        tamanhoMatriz,
    ):
        self.raio = raio  # em pixel
        self.raioSun = raioSun * 696340  # em relação ao raio do Sol
        self.intensidadeMaxima = intensidadeMaxima
        self.coeficienteHum = coeficienteHum
        self.coeficienteDois = coeficienteDois
        self.coeficienteTres = coeficienteTres
        self.coeficienteQuatro = coeficienteQuatro
        self.tamanhoMatriz = tamanhoMatriz

        # Inicializa atributos necessários
        self.Nx = self.tamanhoMatriz
        self.Ny = self.tamanhoMatriz
        self.color = "hot"
        self.manchas: List[Star.Mancha] = []  # Inicializa sem manchas
        self.faculas: List[Star.Facula] = []
        self.estrelaMatriz = Star.criaEstrela(
            self.Nx,
            self.Ny,
            self.tamanhoMatriz,
            self.raio,
            self.intensidadeMaxima,
            self.coeficienteHum,
            self.coeficienteDois,
            self.coeficienteTres,
            self.coeficienteQuatro,
        )

    class Mancha:
        def __init__(self, intensidade, raio, latitude, longitude):
            self.intensidade = intensidade  # intensidade relativa (menor que 1)
            self.raio = raio  # relativo ao raio da estrela
            self.latitude = latitude
            self.longitude = longitude
            # print(f"Mancha criada: intensidade={self.intensidade}, r={self.raio}, lat={self.latitude}, longt={self.longitude}")

    class Facula:
        def __init__(self, raio, intensidade, latitude, longitude):
            self.intensidade = intensidade  # intensidade relativa (maior que 1)
            self.raio = raio  # relativo ao raio da estrela
            self.latitude = latitude
            self.longitude = longitude

    # @njit(parallel=True)
    def criaEstrela(
        lin,
        col,
        tamanhoMatriz,
        raio,
        intensidadeMaxima,
        coeficienteHum,
        coeficienteDois,
        coeficienteTres,
        coeficienteQuatro,
    ):
        """
        Cria a matriz da estrela utilizando o modelo de escurecimento de limbo de 4 parâmetros (Claret)
        com aceleração via JIT.

        Parâmetros:
          - lin, col: dimensões da matriz;
          - tamanhoMatriz: tamanho total da matriz (assume-se quadrada);
          - raio: raio da estrela em pixels;
          - intensidadeMaxima: intensidade máxima no centro da estrela;
          - coeficienteHum, coeficienteDois, coeficienteTres, coeficienteQuatro: coeficientes do modelo.

        Retorna:
          - numpy.ndarray com a matriz da estrela.
        """
        estrela = np.zeros((lin, col), dtype=np.int32)
        center = tamanhoMatriz / 2.0

        # for i in prange(lin):
        for i in range(lin):
            for j in range(col):
                dx = i - center
                dy = j - center
                dist = math.sqrt(dx * dx + dy * dy)
                if dist <= raio:
                    cosTheta = math.sqrt(1.0 - (dist / raio) ** 2)
                    part1 = (
                        1
                        - coeficienteHum * (1 - math.sqrt(cosTheta))
                        - coeficienteDois * (1 - cosTheta)
                    )
                    part2 = -coeficienteTres * (
                        1 - cosTheta**1.5
                    ) - coeficienteQuatro * (1 - cosTheta**2.0)
                    intensidade = intensidadeMaxima * (part1 + part2)
                    if intensidade < 0:
                        intensidade = 0
                    estrela[i, j] = int(intensidade)

        return estrela

    '''
    Ruidos podem ser Manchas ou Fáculas
    '''

    def criaRuidos(self, ruidos):
        # print("Iniciando aplicação de manchas...")

        for ruido in ruidos:
            raio_mancha_pixel = (
                self.raio * ruido.raio
            )  # raio em função do raio da estrela em pixels

            # Converte coordenadas de posicionamento da mancha de graus para radianos
            degreeToRadian = np.pi / 180.0
            latitudeMancha = ruido.latitude * degreeToRadian
            longitudeMancha = ruido.longitude * degreeToRadian

            # Calcula a posição da mancha em pixels, relativa ao centro da estrela
            ys = self.raio * np.sin(latitudeMancha)
            xs = self.raio * np.cos(latitudeMancha) * np.sin(longitudeMancha)
            anguloHelio = np.arccos(
                np.cos(latitudeMancha) * np.cos(longitudeMancha)
            )

            # Efeito de projeção (elipticidade)
            yy = (
                ys + self.Ny / 2
            )  # posição em pixel relativa à origem da matriz
            xx = xs + self.Nx / 2

            kk = np.arange(self.Ny * self.Nx)
            vx = kk - self.Nx * np.int64(1.0 * kk / self.Nx) - xx
            vy = kk / self.Ny - yy

            # Ângulo de rotação da mancha
            anguloRot = np.abs(np.arctan(ys / xs))
            if latitudeMancha * longitudeMancha > 0:
                anguloRot = -anguloRot
            elif latitudeMancha * longitudeMancha == 0:
                anguloRot = 0

            (ii,) = np.where(
                (
                    (
                        (vx * np.cos(anguloRot) - vy * np.sin(anguloRot))
                        / np.cos(anguloHelio)
                    )
                    ** 2
                    + (vx * np.sin(anguloRot) + vy * np.cos(anguloRot)) ** 2
                )
                < raio_mancha_pixel**2
            )

            spot = np.ones(self.Ny * self.Nx)
            spot[ii] = ruido.intensidade
            spot = spot.reshape([self.Ny, self.Nx])

            self.estrelaMatriz = self.estrelaMatriz * spot

        return self.estrelaMatriz

    # Inserção de manchas
    def addMancha(self, mancha: Mancha):
        '''
        Cria e adiciona uma Mancha à estrela.
        '''
        self.manchas.append(mancha)

    def criaEstrelaManchada(self):
        self.estrelaMatriz = self.criaRuidos(self.manchas)
        return self.estrelaMatriz

    # Inserção de Fáculas
    def addFacula(self, facula: Facula):
        '''
        Cria e adiciona uma Fácula à estrela.
        '''
        if facula.intensidade < 1:
            print("O valor da intensidade da Fácula deve ser maior que 1.")
            return
        self.faculas.append(facula)

    def criaEstrelaComFaculas(self):
        self.criaRuidos(self.faculas)

    # Getters e Setters
    def getNx(self):
        return self.Nx

    def getNy(self):
        return self.Ny

    def getRaioStar(self):
        return self.raio

    def getMatrizEstrela(self):
        return self.estrelaMatriz

    def getu1(self):
        return self.coeficienteHum

    def getu2(self):
        return self.coeficienteDois

    def getTamanhoMatriz(self):
        return self.tamanhoMatriz

    def getRaioSun(self):
        return self.raioSun

    def getIntensidadeMaxima(self):
        return self.intensidadeMaxima

    def getError(self):
        return self.error

    def setStarName(self, starName):
        self.starName = starName

    def getStarName(self):
        return self.starName

    def getCadence(self):
        return self.cadence

    def setCadence(self, cadence):
        self.cadence = cadence

    def Plotar(self, tamanhoMatriz, estrela):
        plt.axis([0, tamanhoMatriz, 0, tamanhoMatriz])
        plt.imshow(estrela, self.color)
        plt.gca().invert_yaxis()  # Corrige o eixo Y invertido
        plt.show()

    def resetManchas(self):
        self.manchas = []  # Reseta as manchas para uma lista vazia
