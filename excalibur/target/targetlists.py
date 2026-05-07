'''target targetlists ds'''

# Heritage code shame:
# pylint: disable=invalid-name
# pylint: disable=too-many-branches

import excalibur
import os
import logging

log = logging.getLogger(__name__)


# --------------------------------------------------------------------
def get_target_lists():
    '''
    Load in various target lists
    '''

    targetlistdir = excalibur.context['data_dir'] + '/targetlists/'

    targetlist_filenames = {
        'active': 'all',
        'roudier62': 'roudier62',
        'G141': 'HST_G141',
        'JWST': 'JWST',
        'Spitzer': 'Spitzer',
        'ariel_Aug2024_2years': 'stars_tier2_thorngren_aug2024',
        'ariel_Aug2024_2years_withPlanetletters': 'planets_tier2_thorngren_aug2024',
        'ariel_Nov2024_2years': 'stars_tier2_thorngren_nov2024',
        'ariel_Nov2024_2years_withPlanetletters': 'planets_tier2_thorngren_nov2024',
        'ariel_Nov2024_2yearsTier1': 'stars_tier1_thorngren_nov2024',
        'ariel_Nov2024_2yearsTier1_withPlanetletters': 'planets_tier1_thorngren_nov2024',
        'ariel_massesNeeded': 'planets_arielROSESmassneeded',
    }

    targetlists = {}
    for listname, filename in targetlist_filenames.items():
        targetlists[listname] = []

        if os.path.isfile(targetlistdir + filename + '.txt'):
            with open(
                targetlistdir + filename + '.txt', 'r', encoding='utf-8'
            ) as data:
                lines = data.readlines()
                for line in lines:
                    if line.startswith('#'):
                        # print('skipping a comment')
                        pass
                    else:
                        target = line.replace("'", "")
                        target = target.split('#')[0].split(',')[0].strip()
                        targetlists[listname].append(target)
        else:
            print('target list missing:', targetlistdir + filename)
        # print(targetlists[listname])

    # for targetlist, targets in targetlists.items():
    #    print('# of targets:', targetlist, len(targets))

    return targetlists


# --------------------------------------------------------------------


def read_ArielMCS_info(filename='Ariel_MCS_Known_2024-02-14.csv'):
    '''
    Load in the MCS table with all the Ariel target info
    The most recent file is 2/14/24, but also consider the 11/20/23 file
    '''

    arielDir = excalibur.context['data_dir'] + '/ariel/'

    listofDictionaries = []

    if not os.path.isfile(arielDir + filename):
        log.warning('--< PROBLEM: Ariel MCS table not found : %s >--', filename)
    else:
        with open(arielDir + filename, 'r', encoding='ascii') as file:
            csvFile = csv.DictReader(file)

            # print('starting to read',filename)
            for line in csvFile:
                listofDictionaries.append(line)

    return listofDictionaries


# --------------------------------------------------------------------


def targetlist_ArielMCSknown(
    filedate='Nov2023', maxVisits=666, transitCategoryOnly=False
):
    '''
    Select a batch of targets from the Ariel MCS list of known planets

    If transitCategoryOnly, include all planets in the 'transit' or 'either' category
    Otherwise make the selection based on the number of transit visits (<= maxVisits)
    '''

    if filedate == 'Nov2023':
        filename = 'Ariel_MCS_Known_2023-11-20.csv'
    elif filedate == 'Feb2024':
        filename = 'Ariel_MCS_Known_2024-02-14.csv'
    else:
        log.warning('--< PROBLEM: Unknown date for the Ariel MCS table >--')
        filename = 'Ariel_MCS_Known_2023-11-20.csv'

    targetinfo = read_ArielMCS_info(filename=filename)

    aliases = arielAliases()

    targetList = []
    for target in targetinfo:
        selected = False
        if transitCategoryOnly:
            if (
                target['Preferred Method'] == 'Transit'
                or target['Preferred Method'] == 'Either'
            ):
                selected = True
        else:
            if float(target['Tier 2 Transits']) <= maxVisits:
                selected = True

        if selected:
            # special case for TOI-216.01 and TOI-216.02
            if target['Planet Name'].endswith('.01'):
                target['Planet Name'] = target['Planet Name'][:-3] + ' b'
                # print('fixing planet name from .01 to',target['Planet Name'])
            elif target['Planet Name'].endswith('.02'):
                target['Planet Name'] = target['Planet Name'][:-3] + ' c'
                # print('fixing planet name from .02 to',target['Planet Name'])

            # translate a handful of MCS target names to Excalibur aliases
            if target['Star Name'] in aliases:
                alias = aliases[target['Star Name']]
                target['Star Name'] = alias
                target['Planet Name'] = alias + target['Planet Name'][-2:]

            # save the star names (not the planet names).  that's what's needed for pipeline call
            # targetList.append(target['Planet Name'])
            if target['Star Name'] in targetList:
                # print('skipping a multi-planet entry for',target['Star Name'])
                pass
            else:
                targetList.append(target['Star Name'])

    return targetList


# --------------------------------------------------------------------


def arielAliases():
    '''
    Some of the names in Ariel MCS are different from in Excalibur.
    Translate them using these aliases.
    (first column is MCS name, second column is Excalibur name)
    '''

    aliases = {
        'HD 3167': 'K2-96',
        'L 168-9': 'GJ 4332',
        'LHS 1140': 'GJ 3053',
        'LHS 475': 'GJ 4102',
        'TOI-836': 'LTT 5972',
        'WASP-94 A': 'WASP-94',
        'GJ 143': 'HD 21749',
        'HIP 41378': 'K2-93',
    }

    return aliases
