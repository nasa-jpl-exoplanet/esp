'''util logg ds'''

# -- IMPORTS -- ------------------------------------------------------
import numpy


# -------------- -----------------------------------------------------
# -- SV VALIDITY -- --------------------------------------------------
def calculate_logg(mass, radius, sscmks, units='solar'):
    '''calculate log(g).  units should be solar for a star, Jupiter for a planet'''

    if units == 'solar':
        mass_mks = float(mass) * sscmks['Msun']
        radius_mks = float(radius) * sscmks['Rsun']
    else:
        mass_mks = float(mass) * sscmks['Mjup']
        radius_mks = float(radius) * sscmks['Rjup']

    g_mks = sscmks['G'] * mass_mks / radius_mks**2
    logg = numpy.log10(g_mks)

    return logg
