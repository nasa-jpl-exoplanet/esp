import numpy as np

np.seterr(divide='ignore', invalid='ignore')


def Validar(msg):
    '''função criada para validar entradas, por exemplo numeros nao float/int ou negativos'''
    valor = 0
    while True:
        try:
            n = float(input(msg))
            if n >= 0:
                valor = n
                return valor
            else:
                print("\033[0;31mErro! Digite uma entrada válida\033[m")
        except Exception as erro:
            print(
                f'\033[0;31mO valor digitado é inválido. Por favor, digite novamente. O tipo de problema encontrado foi{erro.__class__}\n\n\033[m'
            )


def ValidarEscolha(msg):
    '''função criada para validar escolhas (1 ou 2)'''
    valor = 0
    while True:
        try:
            n = int(input(msg))
            if (n == 1) or (n == 2):
                valor = n
                return valor
            else:
                print("\033[0;31mErro! Digite uma entrada válida\033[m")
        except Exception as erro:
            print(
                f'\033[0;31mO valor digitado é inválido. Por favor, digite novamente. O tipo de problema encontrado foi{erro.__class__}\n\n\033[m'
            )


def calSemiEixo(periodo, mass):
    '''
    funcao que calcula o semieixo do planeta de acordo com o peridodo atraves da 3a lei de Kepler
    parametros:
    periodo :: periodo do planeta em dias
    G :: constante gravitacional universal
    Pi:: numero de pi
    periodos:: periodo convertido convertido para segundos
    mass:: massa da estrela em relacao a massa do sol
    massestrela:: conversao da massa da estrela
    a :: semi eixo orbital retornado
    '''
    # print(
    # '''
    #                              3a LEI DE KEPLER
    # \033[1;35m------------------------------------------------------------------------------
    # períodos**2= ((4*(pi))**2/G*(massaestrela+massaplaneta))*(semieixoorbital***3)
    # G=9,806 65 m/s²,
    # Pi=3.14159265359
    # -------------------------------------------------------------------------------
    # A seguir, digite a massa da estrela em Kg para que a 3a Lei de Kepler seja apli-
    # cada e entao, o Semi Eixo orbital seja calculado.
    # \033[m''')
    G = 6.674184 * (10 ** (-11))  # constante gravitacao universal
    periodoSeg = (
        periodo * 86400
    )  # transformando o periodo que é dado em dias em segundos
    massEstrela = mass * (1.989 * (10**30))
    a = (((periodoSeg**2) * G * massEstrela) / (4 * (np.pi**2))) ** (1 / 3)
    print('a: ', a)
    return a


def calculaLat(semiEixoRaioStar, anguloInclinacao):
    '''Funcao que calcula latitude para que a mancha seja influente na curva de luz'''
    dtor = np.pi / 180
    lat = -(
        np.arcsin(semiEixoRaioStar * np.cos(anguloInclinacao * dtor)) / dtor
    )
    return lat
