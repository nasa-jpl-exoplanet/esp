'''system constants'''


# ------------- ------------------------------------------------------
# -- SOLAR SYSTEM CONSTANTS -- ---------------------------------------
def ssconstants(mks=False, cgs=False):
    '''
    G. ROUDIER: IAU 2012
    '''
    if mks and cgs:
        ssc = {'Idiot': True}
        pass
    elif mks:
        ssc = {
            'Rsun': 6.957e8,
            'Msun': 1.9884158605722263e30,
            'Lsun': 3.828e26,
            'Rjup': 7.1492e7,
            'Mjup': 1.8985233541508517e27,
            'Rearth': 6.371e6,
            'Mearth': 5.972168e24,
            'AU': 1.495978707e11,
            'G': 6.67428e-11,
            'c': 2.99792e8,
            'Rgas': 8.314462e3,
            'Rjup/Rsun': 1.0276268506540176e-1,
            'Rsun/AU': 4.650467260962158e-3,
            'Tsun': 5772,
        }
        pass
    elif cgs:
        ssc = {
            'Rsun': 6.957e10,
            'Msun': 1.9884158605722263e33,
            'Lsun': 3.828e33,
            'Rjup': 7.1492e9,
            'Mjup': 1.8985233541508517e30,
            'Rearth': 6.371e8,
            'Mearth': 5.972168e27,
            'AU': 1.495978707e13,
            'G': 6.67428e-8,
            'c': 2.99792e10,
            'Rgas': 8.314462e7,
            'Rjup/Rsun': 1.0276268506540176e-1,
            'Rsun/AU': 4.650467260962158e-3,
            'Tsun': 5772,
        }
        pass
    else:
        # hopefully never gets here (mks or cgs should be specified)
        ssc = {
            'Rjup/Rsun': 1.0276268506540176e-1,
            'Rsun/AU': 4.650467260962158e-3,
        }
        pass

    ssc['day'] = 24.0 * 60.0 * 60.0
    return ssc
