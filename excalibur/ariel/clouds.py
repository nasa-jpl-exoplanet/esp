'''ariel planet_clouds ds'''

# Heritage code shame:
# pylint: disable=invalid-name


# -- IMPORTS -- ------------------------------------------------------
import numpy as np
import csv
import logging

# from excalibur.cerberus.bounds import setPriorBound
import excalibur

log = logging.getLogger(__name__)


# ______________________________________________________


def readCloudTable(
    dataDir=excalibur.context['data_dir'] + '/ariel/',
    filename='estrela22_cloudFits.csv',
):
    '''
    Read in the best-fit cloud parameters from Estrela 2022
    '''

    with open(dataDir + filename, 'r', encoding='ascii') as file:
        csvFile = csv.DictReader(file)

        cloudParamList = []
        # print('starting to read',filename)
        for line in csvFile:
            cloudParamList.append(line)

    for cloudParams in cloudParamList:
        # take the log of the parameters, to match the expectation in cerberus
        cloudParams['CTP'] = np.log10(float(cloudParams['CTP']))
        cloudParams['HScale'] = np.log10(float(cloudParams['HScale']))
        cloudParams['HLoc'] = np.log10(float(cloudParams['HLoc']))
        cloudParams['HThick'] = np.log10(float(cloudParams['HThick']))

    return cloudParamList


# ______________________________________________________


def randomCloudParameters():
    '''
    Randomly pull a set of 4 cloud parameters from the best-fits to HST/G141 data (Estrela 2022)
    '''

    cloudParamTable = readCloudTable()
    # print(len(cloudParamTable))
    # print(cloudParamTable[10])

    randomPlanetIndex = int(np.random.random() * len(cloudParamTable))
    # print('randomPlanetIndex',randomPlanetIndex)
    # if randomPlanetIndex >= len(cloudParamTable):
    #    exit('ERROR: out of array len')

    randomCloudParams = cloudParamTable[randomPlanetIndex]

    # verify that the cloud parameters lie within the prior range used for cerberus fitting
    # priorbounds = setPriorBound()
    # for cloudParams in cloudParamTable:
    #    if cloudParams['CTP'] < priorbounds['CTP'][0] or \
    #       cloudParams['CTP'] > priorbounds['CTP'][1]:
    #        print('CTP is outside prior range',cloudParams['CTP'],priorbounds['CTP'])
    #    if cloudParams['HScale'] < priorbounds['HScale'][0] or \
    #       cloudParams['HScale'] > priorbounds['HScale'][1]:
    #        print('HScale is outside prior range',cloudParams['HScale'],priorbounds['HScale'])
    #    if cloudParams['HLoc'] < priorbounds['HLoc'][0] or \
    #       cloudParams['HLoc'] > priorbounds['HLoc'][1]:
    #        print('HLoc is outside prior range',cloudParams['HLoc'],priorbounds['HLoc'])
    #    if cloudParams['HThick'] < priorbounds['HThick'][0] or \
    #       cloudParams['HThick'] > priorbounds['HThick'][1]:
    #        print('HThick is outside prior range',cloudParams['HThick'],priorbounds['HThick'])

    return randomCloudParams


# ______________________________________________________


def fixedCloudParameters():
    '''
    All planets are assigned the same median cloud parameters (from Estrela 2022)
    '''
    medianCloudParams = {}
    medianCloudParams['CTP'] = -1.52  # 0.03 bar
    medianCloudParams['HScale'] = -2.10  # 0.008
    medianCloudParams['HLoc'] = -2.30  # 0.005 bar location (meaningless)
    medianCloudParams['HThick'] = 9.76  # very wide vertical thickness

    return medianCloudParams


# ______________________________________________________
