'''system overwriter ds'''

# Heritage code shame:
# pylint: disable=too-many-lines,too-many-statements

import numpy
import copy
from importlib import import_module as fetch
from excalibur.system.autofill import (
    derive_LOGGplanet_from_R_and_M,
    derive_Teqplanet_from_Lstar_and_sma,
)


def fix_default_reference(target):
    '''
    Set the default reference by hand
    (as opposed to the below code, which sets individual parameters by hand)
    '''

    set_refs_by_hand = {
        # 'GJ 9827':'Bonomo et al. 2023',    # not Rice et al. 2019   # problem is actually scrape/collect
        # not really much improvement on these:
        # 'GJ 3053':'Lillo-Box et al. 2020',   # residuals are better, mainly from 1 point itk
        # 'HD 191939':'Orell-Miquel et al. 2022',  # no change in residuals (887 was a bit better)
        # 'K2-3':'Diamond-Lowe et al. 2022',  # why are residuals a bit better? should be same
        # 'LTT 1445 A':'Oddo et al. 2023',  # 887 is a lot lower residual, very diff looking. why?
    }

    # default_reference = False
    # if target in set_refs_by_hand: default_reference = set_refs_by_hand[target]
    default_reference = set_refs_by_hand.get(target, False)

    return default_reference


# -- PRIORITY PARAMETERS -- ------------------------------------------
def ppar():
    '''
    'starID': stellar ID from the targetlist() i.e: 'WASP-12'
    'planet': planet letter i.e: 'b'
    '[units]': value in units indicated inside brackets
    ref: reference
    overwrite[starID] =
    {
        'R*':[Rsun], 'R*_ref':ref,
        'T*':[K], 'T*_lowerr':[K], 'T*_uperr':[K], 'T*_ref':ref,
        'FEH*':[dex], 'FEH*_lowerr':[dex], 'FEH*_uperr':[dex], 'FEH*_ref':ref,
        'LOGG*':[dex CGS], 'LOGG*_lowerr':[dex CGS], 'LOGG*_uperr':[dex CGS], 'LOGG*_ref':ref,
        planet:
        {
        'inc':[degrees], 'inc_lowerr':[degrees], 'inc_uperr':[degrees], 'inc_ref':ref,
        't0':[JD], 't0_lowerr':[JD], 't0_uperr':[JD], 't0_ref':ref,
        'sma':[AU], 'sma_lowerr':[AU], 'sma_uperr':[AU], 'sma_ref':ref,
        'period':[days], 'period_ref':ref,
        'ecc':[], 'ecc_ref':ref,
        'rp':[Rjup], 'rp_ref':ref
        }
    }
    '''
    sscmks = fetch('excalibur.system.core').ssconstants(cgs=True)

    overwrite = {}
    # this one is somewhat higher than bonomo (0.6).  drop here ig
    # overwrite['GJ 9827'] = {
    #    'R*':0.637, 'R*_uperr':0.063, 'R*_lowerr':-0.063,
    #    'R*_ref':'Prieto-Arranz et al. 2018',
    # }
    overwrite['HAT-P-11'] = {
        #  archive gives Yee 2018 but it's rounded to 0.68 for some reason
        'R*': 0.683,
        'R*_uperr': 0.009,
        'R*_lowerr': -0.009,
        'R*_ref': 'Yee et al. 2018',
        # Yee 2018 table 2 doesn't seem to have any error bars.  archive has 7e-5; we have 3e-7
        #  there's no error bar on t0 either. strange.  oh it says that both come from Huber 2017
        'b': {
            'period': 4.887802443,
            'period_uperr': 3e-7,
            'period_lowerr': -3e-7,
            'period_ref': 'Yee et al. 2018',
        },
    }
    # overwrite['HAT-P-17'] = {
    # stassun 2017 default is similar (0.87 vs 0.838).  let's drop this one
    # 'R*':0.838, 'R*_uperr':0.021, 'R*_lowerr':-0.021,
    # 'R*_ref':'Howard et al. 2012',
    # 'b':{
    # stassun 2017 inc is the same, but with shifted errors (+-0.15 instead of +0.2-0.1)
    # 'inc':89.2,
    # 'inc_uperr':0.2, 'inc_lowerr':-0.1,
    # 'inc_ref':'Howard et al. 2012',
    # NOT SURE ABOUT THIS ONE.  currently we use Kokuri 2023.  should be a lot better
    # 't0':2454801.16945,
    # 't0_uperr':0.0002, 't0_lowerr':-0.0002,
    # 't0_ref':'Howard et al. 2012',
    # this sma doesn't at all match the one derived from P,M* (0.0926).  drop this one
    # 'sma':0.06,
    # 'sma_uperr':0.0014, 'sma_lowerr':-0.0014,
    # 'sma_ref':'Howard et al. 2012 + high ecc',
    # this is the same as stassun 2017, more or less (10.33852)
    # 'period':10.338523,
    # 'period_ref':'Howard et al. 2012'}
    # }
    # drop the sma fill in; better to use a consistently derived number
    # also drop t0; the default is now a newer kokori (2023)
    # overwrite['HAT-P-3'] = {
    #    'b':{'sma':0.03878,  # this is 0.0406 when derived from P,M* (the default)
    #         'sma_uperr':0.00065, 'sma_lowerr':-0.00065,
    #         'sma_ref':'Kokori et al. 2021',
    #         't0':2455694.72623,
    #         't0_uperr':8e-05, 't0_lowerr':-8e-05,
    #         't0_ref':'Kokori et al. 2021'}
    # }
    # this R* is very different from stassun (2.05); keep the self-consistent one
    #  same for sma
    # overwrite['HAT-P-41'] = {
    #     'R*':1.683, 'R*_uperr':0.058, 'R*_lowerr':-0.036,
    #     'R*_ref':'Hartman et al. 2012',
    #     'b':{'sma':0.04258,
    #          'sma_uperr':0.00047, 'sma_lowerr':-0.00048,
    #          'sma_ref':'Kokori et al. 2021'}
    # }
    # default (barros 2017) has similar values. no real reason to set here by hand
    # overwrite['HD 106315'] = {
    #     'R*':1.18, 'R*_uperr':0.11, 'R*_lowerr':-0.11,
    #     'R*_ref':'Crossfield et al. 2017',
    #     'T*':6290, 'T*_uperr':60, 'T*_lowerr':-60,
    #     'T*_ref':'Crossfield et al. 2017',
    #     'FEH*':-0.24, 'FEH*_uperr':0.04, 'FEH*_lowerr':-0.04,
    #     'FEH*_ref':'Crossfield et al. 2017',
    #     'LOGG*':4.29, 'LOGG*_uperr':0.07, 'LOGG*_lowerr':-0.07,
    #     'LOGG*_ref':'Crossfield et al. 2017',
    # }
    # let's keep the derived value of 0.0831
    # overwrite['HD 97658'] = {
    #    'b':{'sma':0.0805,
    #         'sma_uperr':0.001, 'sma_lowerr':-0.001,
    #         'sma_ref':'ExoFOP-TESS TOI'}}

    overwrite['K2-33'] = {
        'FEH*': 0.0,
        'FEH*_uperr': 0.13,
        'FEH*_lowerr': -0.14,
        'FEH*_units': '[dex]',
        'FEH*_ref': 'Mann et al. 2016',
    }

    # default (becker 2019) has similar values. no real reason to set here by hand
    #  hmm, the new one is missing Rp though. strange
    # overwrite['K2-93'] = {
    #     'R*':1.4, 'R*_uperr':0.19, 'R*_lowerr':-0.19,
    #     'R*_ref':'Vanderburg et al. 2016',
    #     'T*':6199, 'T*_uperr':50, 'T*_lowerr':-50,
    #     'T*_ref':'Vanderburg et al. 2016',
    #     'FEH*':-0.11, 'FEH*_uperr':0.08, 'FEH*_lowerr':-0.08,
    #     'FEH*_ref':'Vanderburg et al. 2016',
    #     'LOGG*':4.18, 'LOGG*_uperr':0.1, 'LOGG*_lowerr':-0.1,
    #     'LOGG*_ref':'Vanderburg et al. 2016',
    #     'b':{'inc':88.4,
    #          'inc_uperr':1.6, 'inc_lowerr':-1.6,
    #          'inc_ref':'Vanderburg et al. 2016',
    #          't0':2457152.2844,
    #          't0_uperr':0.0021, 't0_lowerr':-0.0021,
    #          't0_ref':'Vanderburg et al. 2016',
    #          'sma':0.12695775622426692,
    #          'sma_uperr':0.0029, 'sma_lowerr':-0.0029,
    #          'sma_ref':'Vanderburg et al. 2016',
    #          'ecc':0,
    #          'ecc_ref':'Vanderburg et al. 2016',
    #          'rp':0.2561240978011526,
    #          'rp_ref':'Vanderburg et al. 2016',
    #          'period':15.5712,
    #          'period_ref':'Vanderburg et al. 2016',
    # our assumed M-R relation gives 0.03103
    # 'mass':0.0258,
    # 'mass_uperr':0.0110,
    # 'mass_lowerr':-0.0077,
    # 'mass_ref':'Pearson 2019',
    # 'mass_units':'Jupiter mass',
    # 'logg':2.9791,
    # 'logg_lowerr':-0.1539, 'logg_uperr':0.1539,
    # 'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]'
    #          }}
    # default (bonomo 2023) has similar values. mostly identical in fact
    # overwrite['K2-96'] = {
    #     'R*':0.872, 'R*_uperr':0.057, 'R*_lowerr':-0.057,
    #     'R*_ref':'Christiansen et al. 2017',
    #     'T*':5261, 'T*_uperr':60, 'T*_lowerr':-60,
    #     'T*_ref':'Christiansen et al. 2017',
    #     'FEH*':0.04, 'FEH*_uperr':0.05, 'FEH*_lowerr':-0.05,
    #     'FEH*_ref':'Christiansen et al. 2017',
    #     'LOGG*':4.47, 'LOGG*_uperr':0.05, 'LOGG*_lowerr':-0.05,
    #     'LOGG*_ref':'Christiansen et al. 2017',
    #     'b':{'inc':83.4,
    #          'inc_uperr':4.6, 'inc_lowerr':-7.7,
    #          'inc_ref':'Christiansen et al. 2017',
    #          't0':2457394.37454,
    #          't0_uperr':0.00043, 't0_lowerr':-0.00043,
    #          't0_ref':'Christiansen et al. 2017',
    #          'sma':0.01815,
    #          'sma_uperr':0.00023, 'sma_lowerr':-0.00023,
    #          'sma_ref':'Christiansen et al. 2017',
    #          'ecc':0,
    #          'ecc_ref':'Christiansen et al. 2017',
    #          'rp':0.152,
    #          'rp_ref':'Christiansen et al. 2017',
    #          'period':0.959641,
    #          'period_ref':'Christiansen et al. 2017',
    #          'logg':3.22900042686,
    #          'logg_lowerr':0.3, 'logg_uperr':-0.3,
    #          'logg_ref':'System Prior Auto Fill', 'logg_units':'log10[cm.s-2]'},
    #     'c':{'inc':89.30,
    #          'inc_uperr':0.5, 'inc_lowerr':-0.96,
    #          'inc_ref':'Christiansen et al. 2017',
    #          't0':2457394.9788,
    #          't0_uperr':0.0012, 't0_lowerr':-0.0012,
    #          't0_ref':'Christiansen et al. 2017',
    #          'sma':0.1795,
    #          'sma_uperr':0.0023, 'sma_lowerr':-0.0023,
    #          'sma_ref':'Christiansen et al. 2017',
    #          'ecc':0,
    #          'ecc_ref':'Christiansen et al. 2017',
    #          'rp':0.269,
    #          'rp_ref':'Christiansen et al. 2017',
    #          'period':29.8454,
    #          'period_ref':'Christiansen et al. 2017',
    #          'logg':3.02377443746,
    #          'logg_lowerr':0.3, 'logg_uperr':-0.3,
    #          'logg_ref':'System Prior Auto Fill', 'logg_units':'log10[cm.s-2]'}}
    overwrite['KELT-1'] = {
        # default FEH is now 0.131 (from TICv8)
        # 'FEH*':0.009, 'FEH*_uperr':0.073, 'FEH*_lowerr':-0.073,
        # 'FEH*_units':'[dex]', 'FEH*_ref':'Siverd et al. 2012',
        'b': {
            # hmm, these values are off a bit, even though its the same reference
            # might as well stick with the originals
            # 'inc':87.6,
            # 'inc_uperr':1.4, 'inc_lowerr':-1.9,
            # 'inc_ref':'Siverd et al. 2012',
            't0': 2455933.61,
            't0_uperr': 0.00041,
            't0_lowerr': -0.00039,
            't0_ref': 'Siverd et al. 2012 + GMR',
            # 'sma':0.02470,
            # 'sma_uperr':0.00039, 'sma_lowerr':-0.00039,
            # 'sma_ref':'Siverd et al. 2012',
            # 'ecc':0, 'ecc_ref':'Siverd et al. 2012',
            "Spitzer_IRAC1_subarray": [
                0.3499475318779155,
                -0.13450119362315333,
                0.07098128685193948,
                -0.019248332190717504,
            ],
            "Spitzer_IRAC2_subarray": [
                0.34079591311025204,
                -0.21763621595372798,
                0.1569303075862828,
                -0.048363772020055255,
            ],
        }
    }

    overwrite['Kepler-16'] = {
        'R*': 0.665924608009903,
        'R*_uperr': 0.0013,
        'R*_lowerr': -0.0013,
        'R*_ref': 'Oroz + GMR',
        # the Triaud 2022 period (226 days, with no accompaning transit midtime) is no good
        #  none of the HST/G141 falls within the transit; data.timing is empty
        #  (in it's defense, that publication gives an errorbar of 1.7 days!)
        # the only other reference is the discovery paper with 228.776+-0.03
        #  no idea how they got such a small error bar; off by >100-sigma!
        # it's really hard to get the HST just right.  how did they schedule it?!
        #  HST is 12 orbits since the published T_0 in Jan.2020
        #   so a 0.01 error in period translates to a 3 hour shift in Jun.2017 (HST)
        # seems like t0=225.165 but it's still off a bit
        # let's just use the original params from a year ago:
        'b': {
            'inc': 89.7511397641686,
            'inc_uperr': 0.0323,
            'inc_lowerr': -0.04,
            'inc_ref': 'Oroz + GMR',
            't0': 2457914.235774330795,
            't0_uperr': 0.004,
            't0_lowerr': -0.004,
            't0_ref': 'Oroz',
            # this period is the default. can be dropped here
            # 'period':228.776,
            # 'period_uperr':0.03,'period_lowerr':-0.03,
            # 'period_ref':'Doyle et al. 2011'
        },
    }
    # overwrite['Kepler-1625'] = {
    # logg is very similar to our derived result
    # 'b':{'logg':3.0132,
    #     'logg_uperr':0.15, 'logg_lowerr':-0.15,
    #     'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]',
    # this sma would be a bit higher than our derived value (0.8407)
    #     'sma':0.898,
    #     'sma_uperr':0.1,
    #     'sma_lowerr':-0.1,
    #     'sma_ref':'Morton et al. 2016'}}
    # might as well use default Chakrabarty & Sengupta 2019 rather than two refs here
    overwrite['WASP-12'] = {
        #     'R*':1.59, 'R*_uperr':0.18, 'R*_lowerr':-0.18,
        #     'R*_ref':'Stassun et al. 2017',
        #     'T*':6300, 'T*_uperr':150, 'T*_lowerr':-150,
        #     'T*_ref':'Stassun et al. 2017',
        # same FEH as ozul emden 2017 (but larger uncertainty here)
        # 'FEH*':0.3, 'FEH*_uperr':0.14, 'FEH*_lowerr':-0.17,
        # 'FEH*_ref':'Stassun et al. 2017',
        #     'LOGG*':4.38, 'LOGG*_uperr':0.1, 'LOGG*_lowerr':-0.1,
        #     'LOGG*_ref':'Stassun et al. 2017',
        'b': {
            #         'inc':82.5,
            #         'inc_uperr':0.75, 'inc_lowerr':-0.75,
            #         'inc_ref':'Stassun et al. 2017',
            #         't0':2456176.6683,
            #         't0_uperr':0.000078, 't0_lowerr':-0.000078,
            #         't0_ref':'Collins et al. 2017',
            #         'sma':0.0234,
            #         'sma_uperr':0.00056, 'sma_lowerr':-0.0005,
            #         'sma_ref':'Collins et al. 2017',
            #         'period':1.09142245,
            #         'period_ref':'Stassun et al. 2017',
            # limb darkening
            "Spitzer_IRAC1_subarray": [
                0.36885966190119,
                -0.148367404490232,
                0.07112997446947285,
                -0.014533906130047942,
            ],
            "Spitzer_IRAC2_subarray": [
                0.33948631752691805,
                -0.19254408234857706,
                0.1277084571166541,
                -0.037068426815200436,
            ],
        }
    }
    # this is significantly different from Mancini 2017 (0.01) but no reason to think it's better
    # overwrite['WASP-39'] = {
    #    'FEH*':-0.10, 'FEH*_uperr':0.1, 'FEH*_lowerr':-0.1,
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Faedi et al. 2011',
    #    }
    overwrite['WASP-43'] = {
        # FEH is exactly same as bonomo 2017 default
        # 'FEH*':-0.05, 'FEH*_uperr':0.17, 'FEH*_lowerr':-0.17,
        # 'FEH*_units':'[dex]', 'FEH*_ref':'Hellier et al. 2011',
        "b": {
            "Spitzer_IRAC1_subarray": [
                0.5214015151262713,
                -0.116913722511716,
                -0.0025615252155260474,
                0.008679785618454554,
            ],
            "Spitzer_IRAC2_subarray": [
                0.43762215323543396,
                -0.17305029863164503,
                0.09760807455104326,
                -0.029028877897651247,
            ],
        }
    }
    # overwrite['WASP-6'] = {
    # drop t0. there's a newer Kokori (2023) as default
    #     'b':{'t0':2455591.28967,
    #          't0_uperr':7e-05, 't0_lowerr':-7e-05,
    #          't0_ref':'Kokori et al. 2021',
    # use derived sma as default
    # 'sma':0.04217,
    # 'sma_uperr':0.00079,
    # 'sma_lowerr':-0.0012,
    # 'sma_ref':'Kokori et al. 2021'}
    # }
    # drop this one; it's the same as the default
    # overwrite['XO-2'] = {
    #     'b':{'t0':2454508.73829,
    #          't0_uperr':0.00014, 't0_lowerr':-0.00016,
    #          't0_ref':'Crouzet et al. 2012'}
    # }
    # overwrite['XO-3'] = {
    #     'b':{
    # inclination is quite different from stassun 2017 (79.32)
    #  seems best to keep a value that is consistent with the other system params
    # 'inc': 84.20,
    # 'inc_uperr':0.54,
    # 'inc_lowerr':-0.54,
    # 'inc_ref':'Bonomo et al. 2017',
    # omega is similar to bonomo 2017 (349.35)
    # 'omega':347,
    # 'omega_lowerr':-3,
    # 'omega_uperr':3,
    # 'omega_ref':'Wong et al. 2014',
    # ecc is the same as the current value, but errors here much smaller, and ref missing
    # 'ecc':0.29,
    # 'ecc_lowerr':-0.01,
    # 'ecc_uperr':0.01,
    #    }
    # }
    overwrite['HAT-P-23'] = {
        'b': {
            "Spitzer_IRAC1_subarray": [
                0.4028945813566236,
                -0.1618193396025557,
                0.08312362942354319,
                -0.019766348298489313,
            ],
            "Spitzer_IRAC2_subarray": [
                0.3712209668752866,
                -0.1996422788905644,
                0.12409504521199885,
                -0.033786702881953186,
            ],
        }
    }
    overwrite['WASP-14'] = {
        'b': {
            "Spitzer_IRAC1_subarray": [
                0.3556193331718539,
                -0.13491841927882636,
                0.06201863236774508,
                -0.012634699997427995,
            ],
            "Spitzer_IRAC2_subarray": [
                0.3352914789225599,
                -0.1977755003834447,
                0.13543229121842332,
                -0.040489856045654665,
            ],
        }
    }
    overwrite['WASP-34'] = {
        'b': {
            "Spitzer_IRAC1_subarray": [
                0.428983331524869,
                -0.18290950217251944,
                0.09885346596732751,
                -0.025116667946425204,
            ],
            "Spitzer_IRAC2_subarray": [
                0.3809141367993901,
                -0.19189122283729515,
                0.11270592554391648,
                -0.02937059129121932,
            ],
        }
    }
    # barsato 2019 is about the same, but improved precision
    # overwrite['Kepler-9'] = {
    #     'b':{
    #        'mass':0.137,
    #        'mass_uperr':0.008,
    #        'mass_lowerr':-0.01,
    #        'mass_ref':'Hadden & Lithwick et al. 2017',
    #        'mass_units':'Jupiter mass',
    #    }}
    # overwrite['LHS 3844'] = {
    #    'FEH*':0, 'FEH*_uperr':0.25, 'FEH*_lowerr':-0.25,
    #    'FEH*_units':'[dex]', 'FEH*_ref':"Kyle's best guess",
    #    # this mass is somewhat higher than our assumed M-R relation (0.0063)
    #    # 'b':{
    #    #    'mass':0.0118,
    #    #    'mass_uperr':0.0051,
    #    #    'mass_lowerr':-0.0036,
    #    #    'mass_ref':'Pearson 2019',
    #    #    'mass_units':'Jupiter mass',
    #    #    'logg':3.323,
    #    #    'logg_lowerr':-0.15, 'logg_uperr':0.16,
    #    #    'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]'
    #    # }
    # }
    # real masses now from Mayo 2023 (0.01353, 0.05695, 0.00944)
    # omega is quite different (124 vs 24.6) but it's best to stay consistent with Mayo
    #  and both have large uncertainty anyway.  more or less junk
    # overwrite['K2-136'] = {
    #    'b':{
    #        'mass':0.0105,
    #        'mass_uperr':0.0045,
    #        'mass_lowerr':-0.0032,
    #        'mass_ref':'Pearson 2019',
    #        'mass_units':'Jupiter mass',
    #        'logg':3.513,
    #        'logg_lowerr':-0.15, 'logg_uperr':0.15,
    #        'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]'
    #    },
    #    'c':{
    #        'mass':0.0272,
    #        'mass_uperr':0.0116,
    #        'mass_lowerr':-0.0081,
    #        'mass_ref':'Pearson 2019',
    #        'mass_units':'Jupiter mass',
    #        'logg':2.9791,
    #        'logg_lowerr':-0.15, 'logg_uperr':0.15,
    #        'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]',
    #        'omega':24.6,
    #        'omega_lowerr':-74.0,
    #        'omega_uperr':141.0,
    #        'omega_ref':"Mann et al. 2017"
    #    },
    #    'd':{
    #        'mass':0.0130,
    #        'mass_uperr':0.0056,
    #        'mass_lowerr':-0.0039,
    #        'mass_ref':'Pearson 2019',
    #        'mass_units':'Jupiter mass',
    #        'logg':3.2749,
    #        'logg_lowerr':-0.15, 'logg_uperr':0.15,
    #        'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]'
    #    }}
    # overwrite['K2-25'] = {
    #    'b':{
    # there is a real mass measurement now; it's a lot lower (0.07709)  stefansson 2020
    # 'mass':0.0335,
    # 'mass_uperr':0.0143,
    # 'mass_lowerr':-0.0100,
    # 'mass_ref':'Pearson 2019',
    # 'mass_units':'Jupiter mass',
    # 'logg':2.948,
    # 'logg_lowerr':-0.15, 'logg_uperr':0.15,
    # 'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]',
    # sma is similar to stefansson
    # 'sma':0.029535117574370662,
    # 'sma_uperr':0.0021173577439160705,
    # 'sma_lowerr':-0.0021173577439160705,
    # 'sma_ref':'Mann et al. 2016',
    # omega is different now (120+12-14) (from stefansson 2020)
    # 'omega':62,
    # 'omega_lowerr':-39,
    # 'omega_uperr':44,
    # 'omega_ref':"Mann et al. 2016"
    #    }}
    # our assumed M-R relation gives 0.0310
    # overwrite['K2-124'] = {
    #     'b':{
    #         'mass':0.0259,
    #         'mass_uperr':0.0110,
    #         'mass_lowerr':-0.0077,
    #         'mass_ref':'Pearson 2019',
    #         'mass_units':'Jupiter mass',
    #         'logg':2.9791,
    #         'logg_lowerr':-0.15, 'logg_uperr':0.15,
    #         'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]'
    #     }}
    # our assumed M-R relation gives 0.0205
    # overwrite['K2-167'] = {
    #     'b':{
    #         'mass':0.0246,
    #         'mass_uperr':0.0105,
    #         'mass_lowerr':-0.0074,
    #         'mass_ref':'Pearson 2019',
    #         'mass_units':'Jupiter mass',
    #         'logg':2.9871,
    #         'logg_lowerr':-0.15, 'logg_uperr':0.15,
    #         'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]'
    #     }}
    # overwrite['K2-19'] = {
    # there are real mass measurements now by petigura 2020 (0.102,0.034,0.031) very different now
    # 'b':{
    # 'mass':0.2554,
    # 'mass_uperr':0.109,
    # 'mass_lowerr':-0.0764,
    # 'mass_ref':'Pearson 2019',
    # 'mass_units':'Jupiter mass',
    # 'logg':3.1228,
    # 'logg_lowerr':-0.15, 'logg_uperr':0.15,
    # 'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]',
    # there's a newer one now Nespral 2017. 271-12.  quite different but not much effect
    # 'omega':179,
    # 'omega_lowerr':-52,
    # 'omega_uperr':52,
    # 'omega_ref':"Barros et al. 2015"
    # },
    # 'c':{
    # 'mass':0.0683,
    # 'mass_uperr':0.0292,
    # 'mass_lowerr':-0.0204,
    # 'mass_ref':'Pearson 2019',
    # 'mass_units':'Jupiter mass',
    # 'logg':2.9542,
    # 'logg_lowerr':-0.15, 'logg_uperr':0.15,
    # 'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]'
    # },
    # 'd':{
    # 'mass':0.011,
    # 'mass_uperr':0.0047,
    # 'mass_lowerr':-0.0033,
    # 'mass_ref':'Pearson 2019',
    # 'mass_units':'Jupiter mass',
    # 'logg':3.4207,
    # 'logg_lowerr':-0.15, 'logg_uperr':0.15,
    # 'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]'
    # }}
    # our assumed M-R relation gives 0.0133,0.0236
    # overwrite['K2-21'] = {
    #     'b':{
    #         'mass':0.0153,
    #         'mass_uperr':0.0065,
    #         'mass_lowerr':-0.0046,
    #         'mass_ref':'Pearson 2019',
    #         'mass_units':'Jupiter mass',
    #         'logg':3.1488,
    #         'logg_lowerr':-0.15, 'logg_uperr':0.15,
    #         'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]',
    #     },
    #     'c':{
    #         'mass':0.0210,
    #         'mass_uperr':0.0089,
    #         'mass_lowerr':-0.0063,
    #         'mass_ref':'Pearson 2019',
    #         'mass_units':'Jupiter mass',
    #         'logg':3.0235,
    #         'logg_lowerr':-0.15, 'logg_uperr':0.15,
    #         'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]',
    #     }}
    # our assumed M-R relation gives 0.0263
    # overwrite['K2-212'] = {
    #     'b':{
    #         'mass':0.0213,
    #         'mass_uperr':0.0091,
    #         'mass_lowerr':-0.0064,
    #         'mass_ref':'Pearson 2019',
    #         'mass_units':'Jupiter mass',
    #         'logg':3.0195,
    #         'logg_lowerr':-0.15, 'logg_uperr':0.15,
    #         'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]'
    #     }}
    # our assumed M-R relation gives 0.0268
    # overwrite['K2-26'] = {
    #     'b':{
    #         'mass':0.0230,
    #         'mass_uperr':0.0098,
    #         'mass_lowerr':-0.0069,
    #         'mass_ref':'Pearson 2019',
    #         'mass_units':'Jupiter mass',
    #         'logg':3.0014,
    #         'logg_lowerr':-0.15, 'logg_uperr':0.15,
    #         'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]'
    #     }}
    # our assumed M-R relation gives 0.0208
    # overwrite['K2-28'] = {
    #     'b':{
    #         'mass':0.0193,
    #         'mass_uperr':0.0082,
    #         'mass_lowerr':-0.0058,
    #         'mass_ref':'Pearson 2019',
    #         'mass_units':'Jupiter mass',
    #         'logg':3.0487,
    #         'logg_lowerr':-0.15, 'logg_uperr':0.15,
    #         'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]'
    #     }}
    # our assumed M-R relation gives 0.2985 (quite a bit lower)
    # overwrite['K2-289'] = {
    #     'b':{
    #         'mass':0.4283,
    #         'mass_uperr':0.1824,
    #         'mass_lowerr':-0.1279,
    #         'mass_ref':'Pearson 2019',
    #         'mass_units':'Jupiter mass',
    #         'logg':3.2067,
    #         'logg_lowerr':-0.15, 'logg_uperr':0.15,
    #         'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]'
    #     }}
    # our assumed M-R relation gives 1
    # overwrite['K2-52'] = {
    #     'b':{
    #         'mass':1.0821,
    #         'mass_uperr':0.4619,
    #         'mass_lowerr':-0.3237,
    #         'mass_ref':'Pearson 2019',
    #         'mass_units':'Jupiter mass',
    #         'logg':3.0168,
    #         'logg_lowerr':-0.15, 'logg_uperr':0.15,
    #         'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]'
    #     }}
    overwrite['K2-55'] = {
        'FEH*': 0,
        'FEH*_uperr': 0.25,
        'FEH*_lowerr': -0.25,
        'FEH*_units': '[dex]',
        'FEH*_ref': "Kyle's best guess",
        # our assumed M-R relation gives 0.0497
        # 'b':{
        #     'mass':0.0539,
        #     'mass_uperr':0.0230,
        #     'mass_lowerr':-0.0161,
        #     'mass_ref':'Pearson 2019',
        #     'mass_units':'Jupiter mass',
        #     'logg':2.9414,
        #     'logg_lowerr':-0.15, 'logg_uperr':0.15,
        #     'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]'}
    }
    # our assumed M-R relation gives 0.0270,0.0104,0.0116
    # overwrite['K2-58'] = {
    #     'b':{
    #         'mass':0.0221,
    #         'mass_uperr':0.0094,
    #         'mass_lowerr':-0.0066,
    #         'mass_ref':'Pearson 2019',
    #         'mass_units':'Jupiter mass',
    #         'logg':3.0107,
    #         'logg_lowerr':-0.15, 'logg_uperr':0.15,
    #         'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]'
    #     },
    #     'c':{
    #         'mass':0.0138,
    #         'mass_uperr':0.0059,
    #         'mass_lowerr':-0.0041,
    #         'mass_ref':'Pearson 2019',
    #         'mass_units':'Jupiter mass',
    #         'logg':3.2137,
    #         'logg_lowerr':-0.15, 'logg_uperr':0.15,
    #         'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]'
    #     },
    #     'd':{
    #         'mass':0.0144,
    #         'mass_uperr':0.0061,
    #         'mass_lowerr':-0.0043,
    #         'mass_ref':'Pearson 2019',
    #         'mass_units':'Jupiter mass',
    #         'logg':3.1854,
    #         'logg_lowerr':-0.1540, 'logg_uperr':0.15,
    #         'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]'
    #     }}
    # our assumed M-R relation gives 0.02895
    # overwrite['K2-79'] = {
    #     'b':{
    #         'mass':0.0382,
    #         'mass_uperr':0.0163,
    #         'mass_lowerr':-0.0114,
    #         'mass_ref':'Pearson 2019',
    #         'mass_units':'Jupiter mass',
    #         'logg':2.9409,
    #         'logg_lowerr':-0.1541, 'logg_uperr':0.1541,
    #         'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]'
    #     }}
    # our assumed M-R relation gives 0.180 (quite a bit lower)
    # overwrite['K2-87'] = {
    #    'b':{
    #        'mass':0.2410,
    #        'mass_uperr':0.1029,
    #        'mass_lowerr':-0.0721,
    #        'mass_ref':'Pearson 2019',
    #        'mass_units':'Jupiter mass',
    #        'logg':3.1135,
    #        'logg_lowerr':-0.1544, 'logg_uperr':0.1544,
    #        'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]'
    #    }}
    # our assumed M-R relation gives 0.0197
    # overwrite['K2-9'] = {
    #    'b':{
    #        'mass':0.0187,
    #        'mass_uperr':0.008,
    #        'mass_lowerr':-0.0056,
    #        'mass_ref':'Pearson 2019',
    #        'mass_units':'Jupiter mass',
    #        'logg':3.0604,
    #        'logg_lowerr':-0.1540, 'logg_uperr':0.1540,
    #        'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]'
    #    }}
    # our assumed M-R relation gives 0.0256,0.0063 (lower, especially c)
    # overwrite['K2-90'] = {
    #     'b':{
    #         'mass':0.0222,
    #         'mass_uperr':0.0094,
    #         'mass_lowerr':-0.0066,
    #         'mass_ref':'Pearson 2019',
    #         'mass_units':'Jupiter mass',
    #         'logg':3.0095,
    #         'logg_lowerr':-0.1539, 'logg_uperr':0.1539,
    #         'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]'
    #     },
    #     'c':{
    #         'mass':0.0119,
    #         'mass_uperr':0.0051,
    #         'mass_lowerr':-0.0036,
    #         'mass_ref':'Pearson 2019',
    #         'mass_units':'Jupiter mass',
    #         'logg':3.3391,
    #         'logg_lowerr':-0.1548, 'logg_uperr':0.1548,
    #         'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]',
    #     }}
    # our assumed M-R relation gives 0.0389
    # overwrite['K2-95'] = {
    #    'b':{
    #        'mass':0.0420,
    #        'mass_uperr':0.0179,
    #        'mass_lowerr':-0.0125,
    #        'mass_ref':'Pearson 2019',
    #        'mass_units':'Jupiter mass',
    #        'logg':2.9385,
    #        'logg_lowerr':-0.1541, 'logg_uperr':0.1541,
    #        'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]'
    #    }}
    # bonomo 2023 masses are somewhat larger (0.01564 and 0.03502)
    # overwrite['K2-96'] = {
    #    'b':{
    #        'mass':0.0143,
    #        'mass_uperr':0.0061,
    #        'mass_lowerr':-0.0043,
    #        'mass_ref':'Pearson 2019',
    #        'mass_units':'Jupiter mass',
    #        'logg':3.1884,
    #        'logg_lowerr':-0.1543, 'logg_uperr':0.1543,
    #        'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]'
    #    },
    #    'c':{
    #        'mass':0.0252,
    #        'mass_uperr':0.0107,
    #        'mass_lowerr':-0.0075,
    #        'mass_ref':'Pearson 2019',
    #        'mass_units':'Jupiter mass',
    #        'logg':2.9825,
    #        'logg_lowerr':-0.1539, 'logg_uperr':0.1539,
    #        'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]'
    #    }}
    # Grunblatt 2018 default is somewhat lower than this (0.48)
    # overwrite['K2-97'] = {
    #    'b':{
    #        'mass':0.6808,
    #        'mass_uperr':0.2892,
    #        'mass_lowerr':-0.2030,
    #        'mass_ref':'Pearson 2019',
    #        'mass_units':'Jupiter mass',
    #        'logg':3.2746,
    #        'logg_lowerr':-0.1538, 'logg_uperr':0.1538,
    #        'logg_ref':'Pearson 2019', 'logg_units':'log10[cm.s-2]'
    #    }}
    overwrite["CoRoT-2"] = {
        "b": {
            "Spitzer-IRAC-IR-45-SUB": {
                "rprs": 0.15417,
                "ars": 6.60677,
                "inc": 88.08,
                "ref": "KAP",
            },
            "Spitzer_IRAC1_subarray": [
                0.4423526389772671,
                -0.20200004648957037,
                0.11665312313321362,
                -0.03145249632862833,
            ],
            "Spitzer_IRAC2_subarray": [
                0.38011095282800866,
                -0.18511089957904475,
                0.10671540314411156,
                -0.027341272754041506,
            ],
        }
    }
    overwrite['GJ 1132'] = {
        'b': {
            "Spitzer_IRAC1_subarray": [
                0.8808056407530139,
                -0.7457051451918199,
                0.4435989599088468,
                -0.11533694981224148,
            ],
            "Spitzer_IRAC2_subarray": [
                0.8164503264173831,
                -0.7466319094521022,
                0.448251686664617,
                -0.11545411611119284,
            ],
        }
    }

    # Limb darkening coefficients computed with
    # https://github.com/ucl-exoplanets/ExoTETHyS
    overwrite["HD 189733"] = {
        "b": {
            "Spitzer_IRAC1_subarray": [
                0.44354479528696034,
                -0.08947631097404696,
                -0.00435085810513761,
                0.00910594090013591,
            ],
            "Spitzer_IRAC2_subarray": [
                0.38839253774457744,
                -0.17732423236450323,
                0.11810733772498544,
                -0.03657823168637474,
            ],
        }
    }

    overwrite["HD 209458"] = {
        "b": {
            "Spitzer_IRAC1_subarray": [
                0.3801802812964466,
                -0.14959456473437277,
                0.08226460839494475,
                -0.02131689251855459,
            ],
            "Spitzer_IRAC2_subarray": [
                0.36336971411871394,
                -0.21637839776809858,
                0.14792839620718698,
                -0.04286208270953803,
            ],
        }
    }

    overwrite["KELT-9"] = {
        "b": {
            "Spitzer_IRAC1_subarray": [
                0.3313765755137107,
                -0.3324189051186633,
                0.2481428693159906,
                -0.0730038845221279,
            ],
            "Spitzer_IRAC2_subarray": [
                0.3081800911324395,
                -0.32930034680921816,
                0.24236433024537915,
                -0.07019145797258527,
            ],
        }
    }

    overwrite["WASP-33"] = {
        "b": {
            "Spitzer_IRAC1_subarray": [
                0.3360838105569875,
                -0.20369556446757797,
                0.14180512020307806,
                -0.04279692505871632,
            ],
            "Spitzer_IRAC2_subarray": [
                0.3250093115902963,
                -0.2634497438671309,
                0.19583740005736275,
                -0.05877816796715111,
            ],
        }
    }

    overwrite["WASP-103"] = {
        "b": {
            "Spitzer_IRAC1_subarray": [
                0.3758384436095625,
                -0.1395975318171088,
                0.0693688736769638,
                -0.0162794345748232,
            ],
            "Spitzer_IRAC2_subarray": [
                0.35938501700892445,
                -0.20947695473210462,
                0.14144809404384948,
                -0.040990804709028064,
            ],
        }
    }

    overwrite["KELT-16"] = {
        "b": {
            "Spitzer_IRAC1_subarray": [
                0.36916688544219783,
                -0.13957103936844534,
                0.07119044558535764,
                -0.018675120861031937,
            ],
            "Spitzer_IRAC2_subarray": [
                0.3549546430076809,
                -0.2155295179664244,
                0.15043075368738368,
                -0.04514276133343067,
            ],
        }
    }

    overwrite["WASP-121"] = {
        "b": {
            "Spitzer_IRAC1_subarray": [
                0.35343787303428653,
                -0.13444332321181401,
                0.06955169678670275,
                -0.018419427272667512,
            ],
            "Spitzer_IRAC2_subarray": [
                0.34304917676671737,
                -0.21584682478514353,
                0.15457937580646092,
                -0.04742589029069567,
            ],
        }
    }

    overwrite["KELT-20"] = {
        "b": {
            "Spitzer_IRAC1_subarray": [
                0.3344300791273318,
                -0.28913882534855895,
                0.21010188872209157,
                -0.061790086627358104,
            ],
            "Spitzer_IRAC2_subarray": [
                0.31807958384684176,
                -0.31082384355713244,
                0.22719653693837855,
                -0.06610531710958574,
            ],
        }
    }

    overwrite["HAT-P-7"] = {
        "b": {
            "Spitzer_IRAC1_subarray": [
                0.3637280943230625,
                -0.1424523111830543,
                0.06824894539731155,
                -0.014578311756816686,
            ],
            "Spitzer_IRAC2_subarray": [
                0.3397626817587848,
                -0.1976647471694525,
                0.13403605366799531,
                -0.039618202725551235,
            ],
        }
    }

    overwrite["WASP-76"] = {
        "b": {
            "Spitzer_IRAC1_subarray": [
                0.37508785151968,
                -0.15123541065822635,
                0.07733376118565834,
                -0.018551329687575616,
            ],
            "Spitzer_IRAC2_subarray": [
                0.3503514529020057,
                -0.20455423732189046,
                0.13804368117965113,
                -0.040335740501177636,
            ],
        }
    }

    overwrite["WASP-19"] = {
        "b": {
            "Spitzer_IRAC1_subarray": [
                0.4485159019167513,
                -0.20853342549558768,
                0.12340401852800424,
                -0.034379873907955626,
            ],
            "Spitzer_IRAC2_subarray": [
                0.3825821042080216,
                -0.18816485899811294,
                0.11078803979402557,
                -0.029013926574001703,
            ],
        }
    }

    overwrite["KELT-7"] = {
        "b": {
            "Spitzer_IRAC1_subarray": [
                0.34264654838082775,
                -0.15263897274083513,
                0.09529544911293918,
                -0.028445262603159684,
            ],
            "Spitzer_IRAC2_subarray": [
                0.3362025591248816,
                -0.23293403482921576,
                0.17083787539542783,
                -0.05207951610486718,
            ],
        }
    }

    overwrite["KELT-14"] = {
        "b": {
            "Spitzer_IRAC1_subarray": [
                0.4493742805668009,
                -0.24639842685467592,
                0.1640017017159755,
                -0.04777987900661082,
            ],
            "Spitzer_IRAC2_subarray": [
                0.3941297227246355,
                -0.22099935114922328,
                0.13393719173432234,
                -0.03441830202500127,
            ],
        }
    }

    overwrite["WASP-74"] = {
        "b": {
            "Spitzer_IRAC1_subarray": [
                0.3953736210736681,
                -0.16565075764571777,
                0.09318866035061182,
                -0.024247454399023875,
            ],
            "Spitzer_IRAC2_subarray": [
                0.3696898156305795,
                -0.2125028155143125,
                0.13938377401131619,
                -0.03911691640241164,
            ],
        }
    }

    overwrite["HD 149026"] = {
        "b": {
            "Spitzer_IRAC1_subarray": [
                0.367160938263667,
                -0.1305588325879479,
                0.06612034580484898,
                -0.018337084470844645,
            ],
            "Spitzer_IRAC2_subarray": [
                0.3601039845595216,
                -0.22202760617949738,
                0.15729906194707344,
                -0.047668070962046706,
            ],
        }
    }

    overwrite["TrES-3"] = {
        "b": {
            "Spitzer_IRAC1_subarray": [
                0.43597266162376314,
                -0.19001215158398352,
                0.1056322815109545,
                -0.027744065630210032,
            ],
            "Spitzer_IRAC2_subarray": [
                0.3810582548901189,
                -0.18972122146795323,
                0.1119886599627006,
                -0.029375587180729256,
            ],
        }
    }

    overwrite["WASP-77 A"] = {
        "b": {
            "Spitzer_IRAC1_subarray": [
                0.44975696488393374,
                -0.17728194779592824,
                0.09158922569805375,
                -0.02479615071127561,
            ],
            "Spitzer_IRAC2_subarray": [
                0.3906922322805967,
                -0.20039168705178179,
                0.13035218295191758,
                -0.03796619908753919,
            ],
        }
    }

    overwrite["WASP-95"] = {
        "b": {
            "Spitzer_IRAC1_subarray": [
                0.413096705663071,
                -0.17163366883006448,
                0.09079314140373773,
                -0.022304168790776884,
            ],
            "Spitzer_IRAC2_subarray": [
                0.3772265443367232,
                -0.2006526678014234,
                0.12252287674677063,
                -0.03283398178266608,
            ],
        }
    }

    overwrite["WASP-140"] = {
        "b": {
            "Spitzer_IRAC1_subarray": [
                0.4294348108523915,
                -0.10313120787437623,
                0.016430633217998755,
                0.0008036783543789007,
            ],
            "Spitzer_IRAC2_subarray": [
                0.39023842724510155,
                -0.19846841760266795,
                0.13665600099047498,
                -0.04254450064955062,
            ],
        }
    }

    overwrite["WASP-52"] = {
        "b": {
            "Spitzer_IRAC1_subarray": [
                0.4542826787797213,
                -0.10475364767168102,
                0.01183804531437866,
                0.0029171050937958822,
            ],
            "Spitzer_IRAC2_subarray": [
                0.3906750566059761,
                -0.17705502335732198,
                0.11703224188529365,
                -0.03571381970099965,
            ],
        }
    }

    overwrite["GJ 1214"] = {
        "b": {
            # switch from the most recent ephemeris (Kokori 2022 = exoClock) back to the default
            # 'period':1.58040433,
            # 'period_uperr':1.3e-7, 'period_lowerr':-1.3e-7,
            # 'period_ref':'Cloutier et al. 2021',
            # 't0':2455701.413328,
            # 't0_uperr':0.000066, 't0_lowerr':-0.000059,
            # 't0_ref':'Cloutier et al. 2021',
            "Spitzer_IRAC1_subarray": [
                0.9083242210542111,
                -0.7976911808204602,
                0.4698074336560188,
                -0.12001861589169728,
            ],
            "Spitzer_IRAC2_subarray": [
                0.8239880988090422,
                -0.760781868877928,
                0.4513165756893245,
                -0.11497950826716168,
            ],
        }
    }

    # overwrite['WASP-87'] = {
    #    'FEH*':0, 'FEH*_uperr':0.25, 'FEH*_lowerr':-0.25,
    #    'FEH*_units':'[dex]', 'FEH*_ref':"Default to solar metallicity"}

    # only one system (this one) is missing an H_mag
    #  make a guess at it based on V=16.56,I=15.30
    overwrite['OGLE-TR-056'] = {
        'Jmag': 14.5,
        'Jmag_uperr': 1,
        'Jmag_lowerr': -1,
        'Jmag_units': '[mag]',
        'Jmag_ref': 'Geoff guess',
        'Hmag': 14,
        'Hmag_uperr': 1,
        'Hmag_lowerr': -1,
        'Hmag_units': '[mag]',
        'Hmag_ref': 'Geoff guess',
        'Kmag': 14,
        'Kmag_uperr': 1,
        'Kmag_lowerr': -1,
        'Kmag_units': '[mag]',
        'Kmag_ref': 'Geoff guess',
    }

    # there's a bug in the Archive where this planet's radius is
    #  only given in Earth units, not our standard Jupiter units
    #  0.37+-0.18 REarth = 0.033+-0.16
    # mass is normally filled in via an assumed MRrelation; needs to be done here instead
    # logg is normally calculated from M+R; needs to be done here instead
    # Weiss 2024 has this as radius=blank and flagA=candidate planet that might be noise
    # 8/28/24 drop this overwrite info
    #   1) the radius is only from an older KOI table and is only SNR=2
    #   2) the mass is a recent measurement (but probably it's just the TTV measure?)
    #   3) 0.4 REarth and 8 MEarth is pretty crazy
    # overwrite['Kepler-37'] = {
    #    'e':{'rp':0.033,
    #         'rp_uperr':0.016, 'rp_lowerr':-0.016,
    #         'rp_units':'[Jupiter radius]',
    #         'rp_ref':'Q1-Q8 KOI Table',
    #         # there's a real mass measurement now from Weiss 2024 (0.0255)
    #         #   hold on this is 100x larger!
    #         # 'mass':0.0002,
    #         # 'mass_uperr':0.0002, 'mass_lowerr':-0.0001,
    #         # 'mass_units':'[Jupiter mass]',
    #         # 'mass_ref':'Assumed mass/radius relation',
    #         # 'logg':2.7,
    #         # 'logg_uperr':0.3, 'logg_lowerr':-0.2,
    #         # 'logg_units':'log10[cm.s-2]',
    #         # 'logg_ref':'Assumed mass/radius relation'
    #         # we still need to set logg, because there is originally no radius, so it's blank
    #         'logg':4.76,  # gees it's 100x larger gravity now.  this is junk
    #         'logg_uperr':0.3, 'logg_lowerr':-0.2,
    #         'logg_units':'log10[cm.s-2]',
    #         'logg_ref':'from Mp and Rp'
    #         }}

    # for the newly added comfirmed-planet Ariel targets, some metallicities are missing
    #  oh that's funny. the Chen 2021 compilation has zero for these (with no error bar)
    #   but the source is listed as flag 5, which is the Exoplanet Archive.  fake news!
    # overwrite['HATS-50'] = {
    #    'FEH*':0.0, 'FEH*_uperr':0.25, 'FEH*_lowerr':-0.25,
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Default to solar metallicity'}
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Chen et al. 2021'}
    # overwrite['HATS-51'] = {
    #    'FEH*':0.0, 'FEH*_uperr':0.25, 'FEH*_lowerr':-0.25,
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Default to solar metallicity'}
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Chen et al. 2021'}
    overwrite['HATS-52'] = {
        'FEH*': -0.09,
        'FEH*_uperr': 0.17,
        'FEH*_lowerr': -0.17,
        'FEH*_units': '[dex]',
        'FEH*_ref': 'Magrini et al. 2022',
    }
    # overwrite['HATS-53'] = {
    #    'FEH*':0.0, 'FEH*_uperr':0.25, 'FEH*_lowerr':-0.25,
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Default to solar metallicity'}
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Chen et al. 2021'}
    # overwrite['HATS-58 A'] = {
    #    'FEH*':0.0, 'FEH*_uperr':0.25, 'FEH*_lowerr':-0.25,
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Default to solar metallicity'}
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Chen et al. 2021'}
    overwrite['K2-129'] = {
        'FEH*': 0.105,
        'FEH*_uperr': 0.235,
        'FEH*_lowerr': -0.235,
        'FEH*_units': '[dex]',
        'FEH*_ref': 'Hardagree-Ullman et al. 2020',
    }
    # overwrite['LHS 1678'] = {
    #    'FEH*':0.0, 'FEH*_uperr':0.25, 'FEH*_lowerr':-0.25,
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Default to solar metallicity'}
    # Ammons 2006 has 0.76+-2.26, which is absurd
    # overwrite['TIC 257060897'] = {
    #    'FEH*':0.0, 'FEH*_uperr':0.25, 'FEH*_lowerr':-0.25,
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Default to solar metallicity'}
    # overwrite['TOI-122'] = {
    #    'FEH*':0.0, 'FEH*_uperr':0.25, 'FEH*_lowerr':-0.25,
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Default to solar metallicity'}
    # overwrite['TOI-1227'] = {
    #    'FEH*':0.0, 'FEH*_uperr':0.25, 'FEH*_lowerr':-0.25,
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Default to solar metallicity'}
    # overwrite['TOI-1442'] = {
    #    'FEH*':0.0, 'FEH*_uperr':0.25, 'FEH*_lowerr':-0.25,
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Default to solar metallicity'}
    # overwrite['TOI-1693'] = {
    #    'FEH*':0.0, 'FEH*_uperr':0.25, 'FEH*_lowerr':-0.25,
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Default to solar metallicity'}
    # overwrite['TOI-237'] = {
    #    'FEH*':0.0, 'FEH*_uperr':0.25, 'FEH*_lowerr':-0.25,
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Default to solar metallicity'}
    # overwrite['TOI-2411'] = {
    #    'FEH*':0.0, 'FEH*_uperr':0.25, 'FEH*_lowerr':-0.25,
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Default to solar metallicity'}
    # overwrite['TOI-2427'] = {
    #    'FEH*':0.0, 'FEH*_uperr':0.25, 'FEH*_lowerr':-0.25,
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Default to solar metallicity'}
    # overwrite['TOI-451'] = {
    #    'FEH*':0.0, 'FEH*_uperr':0.25, 'FEH*_lowerr':-0.25,
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Default to solar metallicity'}
    # overwrite['TOI-540'] = {
    #    'FEH*':0.0, 'FEH*_uperr':0.25, 'FEH*_lowerr':-0.25,
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Default to solar metallicity'}
    # there's a value for this one now (osborne 2023)
    # overwrite['TOI-544'] = {
    #     'FEH*':0.0, 'FEH*_uperr':0.25, 'FEH*_lowerr':-0.25,
    #     'FEH*_units':'[dex]', 'FEH*_ref':'Default to solar metallicity'}
    # overwrite['TOI-833'] = {
    #    'FEH*':0.0, 'FEH*_uperr':0.25, 'FEH*_lowerr':-0.25,
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Default to solar metallicity'}
    # this one is missing JHK photometry.  not sure why; it's in 2MASS/Simbad
    overwrite['K2-295'] = {
        'Jmag': 11.807,
        'Jmag_uperr': 0.027,
        'Jmag_lowerr': -0.027,
        'Jmag_units': '[mag]',
        'Jmag_ref': '2MASS',
        'Hmag': 11.259,
        'Hmag_uperr': 0.023,
        'Hmag_lowerr': -0.023,
        'Hmag_units': '[mag]',
        'Hmag_ref': '2MASS',
        'Kmag': 11.135,
        'Kmag_uperr': 0.025,
        'Kmag_lowerr': -0.025,
        'Kmag_units': '[mag]',
        'Kmag_ref': '2MASS',
    }
    # this one is missing R* and M*.  That's unusual!
    # ah wait it does actually have a log-g measure of 4.1 (lower than Solar)
    # arg this one is really tricky.  planet semi-major axis is undefined without M*
    overwrite['WASP-110'] = {
        # R* has a value now (0.86 from 'ExoFOP-TESS TOI')
        # 'R*':1.0, 'R*_uperr':0.25, 'R*_lowerr':-0.25,
        # 'R*_ref':'Default to solar radius',
        'M*': 1.0,
        'M*_uperr': 0.25,
        'M*_lowerr': -0.25,
        'M*_ref': 'Default to solar mass',
        # 1/22/25 There is a LOGG* value now, so no need to fill in M* or RHO*
        # RHO* derivation (from R* and M*) comes before this, so we have to set it here
        # 'RHO*':1.4, 'RHO*_uperr':0.25, 'RHO*_lowerr':-0.25,
        # 'RHO*_ref':'Default to solar density',
        #  older density above gives an inconsistency flag (doesn't match current radius)
        # redo density for R* = 0.864 RSun  (solar density is 1.41)
        'RHO*': 2.2,
        'RHO*_uperr': 0.25,
        'RHO*_lowerr': -0.25,
        'RHO*_ref': 'Assumes solar mass',
        # L* needed for teq (actually it's set below)
        # L* is derived from R*,T*, now that R* has a default value
        # 'L*':1.0, 'L*_uperr':0.25, 'L*_lowerr':-0.25,
        # 'L*_ref':'Default to solar luminosity',
        # 'LOGG*':4.3, 'LOGG*_uperr':0.1, 'LOGG*_lowerr':-0.1,
        # 'LOGG*_ref':'Default to solar log(g)',
        # make sure LOGG* matches M* above (and R*=0.86 now)
        'LOGG*': 4.57,
        'LOGG*_uperr': 0.1,
        'LOGG*_lowerr': -0.1,
        'LOGG*_ref': 'derived from M*,R*',
        # Period is 3.87 days
        # teq derivation (from L* and sma) comes before this, so we have to set it here
        'b': {
            'sma': 0.05,
            'sma_uperr': 0.01,
            'sma_lowerr': -0.01,
            'sma_ref': 'Assume solar mass',
            # teq is derived from R*,T* now that R* has a default value
            # 'teq':1245, 'teq_uperr':100, 'teq_lowerr:':-100,
            # 'teq_units':'[K]', 'teq_ref':'derived from L*,sma',
        },
    }

    # this one is weird. there's a metallicity value in the Archive for 'c' but not for 'b'
    # looks like this is fixed by 2024 archive update; Capistrant et al. 2024 is similar (0.03)
    # overwrite['HD 63433'] = {
    #    'FEH*':0.05, 'FEH*_uperr':0.05, 'FEH*_lowerr':-0.05,
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Dai et al. 2020'}

    # overwrite['TOI-1411'] = {
    #    'FEH*':0.0, 'FEH*_uperr':0.25, 'FEH*_lowerr':-0.25,
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Default to solar metallicity'}

    # why isn't this in the archive?  non-hipparcos, but still..
    overwrite['TRAPPIST-1'] = {
        'dist': (1000.0 / 80.2123),
        'dist_uperr': 0.01,
        'dist_lowerr': -0.01,
        'dist_units': '[pc]',
        'dist_ref': 'Gaia EDR3',
    }
    # had to use vizier for this one; not in simbad for some reason
    overwrite['NGTS-10'] = {
        'dist': (1000.0 / 3.8714),
        'dist_uperr': 12.0,
        'dist_lowerr': -12.0,
        'dist_units': '[pc]',
        'dist_ref': 'Gaia EDR3',
    }
    overwrite['Kepler-1314'] = {
        'dist': (1000.0 / 7.0083),
        'dist_uperr': 2.0,
        'dist_lowerr': -2.0,
        'dist_units': '[pc]',
        'dist_ref': 'Gaia EDR3',
    }
    # also had to use vizier for this one
    # it's in Gaia, but there's no parallax
    # wikipedia has it at 980pc from 2011 schneider site from buchhave 2011 discovery paper
    # the paper says it is from Girardi isochrone fitting
    overwrite['Kepler-14'] = {
        'dist': 980.0,
        'dist_uperr': 100.0,
        'dist_lowerr': -100.0,
        'dist_units': '[pc]',
        'dist_ref': 'Buchhave et al. 2011',
    }

    # 11/10/23 period update to match G141 phase
    overwrite['HAT-P-26'] = {
        # 'b':{'period':4.234520,  # this is the default. decreasing it a bit
        'b': {
            'period': 4.2345002,
            'period_uperr': 7e-7,
            'period_lowerr': -7e-7,
            'period_ref': 'Kokori et al. 2022',
        }
    }

    # 11/12/23 period updates to match G141 phase
    # overwrite['HAT-P-18'] = {
    # 'b':{'period':5.508023,  # this is the default. increasing it a bit
    # hmm, these are about the same.  what about t0?  308P+1.4min diff
    # yeah it's about the same itk. drop it
    # 'b':{'period':5.5080287,
    # 'period_uperr':1.4e-6, 'period_lowerr':-1.4e-6,
    # 'period_ref':'Ivshina & Winn 2022',
    # this is similar to Hartman 2011 120+-56
    # 'omega':104.0,
    # 'omega_lowerr':-50,
    # 'omega_uperr':50,
    # 'omega_ref':"Esposito et al. 2014"
    # }
    # }

    # some of the new JWST targets are missing mandatory parameters
    #  (without these system.finalize will crash)

    # not much in Vizier.  there's two StarHorse metallicities 0.0698 and -0.101483
    # overwrite['GJ 4102'] = {
    #    'FEH*':0.0, 'FEH*_uperr':0.25, 'FEH*_lowerr':-0.25,
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Default to solar metallicity'}
    # even less in Vizier for this white dwarf.  e.g. C/He and Ca/He both blank
    # overwrite['WD 1856'] = {
    #    'FEH*':0.0, 'FEH*_uperr':0.25, 'FEH*_lowerr':-0.25,
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Default to solar metallicity'}

    # for the 75 new Ariel targets in the Feb.14,2024 Edwards target list
    #  3 are missing the mandatory stellar metallicity (2 aren't even in SIMBAD!)
    #
    overwrite['TOI-2445'] = {
        'FEH*': -0.140,
        'FEH*_uperr': 0.25,
        'FEH*_lowerr': -0.25,
        'FEH*_units': '[dex]',
        'FEH*_ref': 'Sprague et al. 2022',
    }
    overwrite['TOI-2459'] = {  # CD-39 1993
        'FEH*': 0.01,
        'FEH*_uperr': 0.25,
        'FEH*_lowerr': -0.25,
        'FEH*_units': '[dex]',
        'FEH*_ref': 'Bochanski et al. 2018',
    }
    overwrite['TOI-5803'] = {  # TYC 556-982-1
        'FEH*': 0.02,
        'FEH*_uperr': 0.25,
        'FEH*_lowerr': -0.25,
        'FEH*_units': '[dex]',
        'FEH*_ref': 'Ammons et al. 2006',
    }
    # and another 100+ new Ariel targets considered (eclipse targets plus Nov.2023 targets)
    #  8 more are missing the mandatory stellar metallicity (toi-4308 is not in simbad even)
    #
    # overwrite['Gaia-1'] = {
    #    'FEH*':0.0, 'FEH*_uperr':0.25, 'FEH*_lowerr':-0.25,
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Default to solar metallicity'}
    overwrite['Gaia-2'] = {
        'FEH*': -0.49,
        'FEH*_uperr': 0.25,
        'FEH*_lowerr': -0.25,
        'FEH*_units': '[dex]',
        'FEH*_ref': 'Ammons et al. 2006',
    }
    overwrite['HIP 9618'] = {
        'FEH*': -0.07,
        'FEH*_uperr': 0.25,
        'FEH*_lowerr': -0.25,
        'FEH*_units': '[dex]',
        'FEH*_ref': 'Xiang et al. 2019',
    }
    overwrite['K2-321'] = {
        'FEH*': -0.05,
        'FEH*_uperr': 0.25,
        'FEH*_lowerr': -0.25,
        'FEH*_units': '[dex]',
        'FEH*_ref': 'Ding et al. 2022',
    }
    # overwrite['K2-417'] = {
    #    'FEH*':0.0, 'FEH*_uperr':0.25, 'FEH*_lowerr':-0.25,
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Default to solar metallicity'}
    overwrite['TOI-206'] = {
        'FEH*': 0.057,
        'FEH*_uperr': 0.25,
        'FEH*_lowerr': -0.25,
        'FEH*_units': '[dex]',
        'FEH*_ref': 'Sprague et al. 2022',
    }
    # overwrite['TOI-4308'] = {
    #    'FEH*':0.0, 'FEH*_uperr':0.25, 'FEH*_lowerr':-0.25,
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Default to solar metallicity'}
    overwrite['TOI-4342'] = {
        'FEH*': -0.090,
        'FEH*_uperr': 0.25,
        'FEH*_lowerr': -0.25,
        'FEH*_units': '[dex]',
        'FEH*_ref': 'Yu et al. 2023',
    }

    # stellar distance isn't an excalibur-mandatory parameter, but it's used by ArielRad
    #  so try to fill it in if it's blank
    #
    # this one is in simbad.   parallax = 3.5860 [0.0397]
    overwrite['TOI-3540 A'] = {
        'dist': 278.9,
        'dist_uperr': 3.1,
        'dist_lowerr': -3.1,
        'dist_units': '[pc]',
        'dist_ref': 'Gaia EDR3',
    }
    # this one is in simbad.   parallax = 2.9341 [0.0939]
    overwrite['TOI-2977'] = {
        'dist': 340.82,
        'dist_uperr': 10.9,
        'dist_lowerr': -10.9,
        'dist_units': '[pc]',
        'dist_ref': 'Gaia EDR3',
    }
    # this one is NOT in simbad.  Gaia Plx missing in Vizier.   parallax =
    #  oh wait, it's already set above. nevermind
    # overwrite['Kepler-14'] = {
    #    'dist':, 'dist_uperr':, 'dist_lowerr':,
    #    'dist_units':'[pc]', 'dist_ref':'Gaia EDR3'}
    # this one is in simbad.   parallax = 7.0083 [0.1077]
    #  drop it - it's already done above here
    # overwrite['Kepler-1314'] = {
    #    'dist':142.69, 'dist_uperr':2.2, 'dist_lowerr':-2.2,
    #    'dist_units':'[pc]', 'dist_ref':'Gaia EDR3'}
    # this one is NOT in simbad, but Vizier has it.   parallax = 3.8714+-0.1759
    #  drop it - it's already done above here
    # overwrite['NGTS-10'] = {
    #    'dist':258.30, 'dist_uperr':11.7, 'dist_lowerr':-11.7,
    #    'dist_units':'[pc]', 'dist_ref':'Gaia EDR3'}
    # this one is in simbad.   parallax = 80.2123 [0.0716]
    #  drop it - it's already done above here
    # overwrite['TRAPPIST-1'] = {
    #    'dist':12.467, 'dist_uperr':0.01, 'dist_lowerr':-0.01,
    #    'dist_units':'[pc]', 'dist_ref':'Gaia EDR3'}

    # for debugging:
    # overwrite["WASP-69"] = {
    #    "b": {
    #        't0':2455748.8342,
    #        't0_uperr':0.00018, 't0_lowerr':-0.00018,
    #        't0_ref':'Bonomo et al. 2017',
    #    }}

    # (for Raissa's paper on L 98-59)
    # make planets c and d consistent with planet b
    # 12/23/2024  fixed now with the self-consistent parameter selection
    # overwrite['L 98-59'] = {
    #    'c':{
    #        'rp':0.1236, 'rp_uperr':0.0085, 'rp_lowerr':-0.0067,
    #        'rp_ref':'Demangeon et al. 2021',
    #        'mass':0.00698, 'mass_uperr':0.00082, 'mass_lowerr':-0.00079,
    #        'mass_ref':'Demangeon et al. 2021',
    #        'logg_ref':'Demangeon et al. 2021',
    #        'teq_ref':'Demangeon et al. 2021',
    #        'sma':0.0304, 'sma_uperr':0.0011, 'sma_lowerr':-0.0012,
    #        'sma_ref':'Demangeon et al. 2021',
    #        'period':3.6906777, 'period_uperr':1.6e-6, 'period_lowerr':-2.6e-6,
    #        'period_ref':'Demangeon et al. 2021',
    #        't0':2458367.27375, 't0_uperr':0.00013, 't0_lowerr':-0.00022,
    #        't0_ref':'Demangeon et al. 2021',
    #        'inc':88.11, 'inc_uperr':0.36, 'inc_lowerr':-0.16,
    #        'inc_ref':'Demangeon et al. 2021',
    #        'ecc':0.103, 'ecc_uperr':0.045, 'ecc_lowerr':-0.058,
    #        'ecc_ref':'Demangeon et al. 2021',
    #        'omega':261., 'omega_uperr':20., 'omega_lowerr':-10.,
    #        'omega_ref':'Demangeon et al. 2021',
    #        'impact':0.601, 'impact_uperr':0.081, 'impact_lowerr':-0.066,
    #        'impact_ref':'Demangeon et al. 2021',
    #        'trandur':1.346, 'trandur_uperr':0.122, 'trandur_lowerr':-0.069,
    #        'trandur_ref':'Demangeon et al. 2021',
    #        'ars':19.00, 'ars_uperr':1.20, 'ars_lowerr':-0.80,
    #        'ars_ref':'Demangeon et al. 2021',
    #        'rprs':0.04088, 'rprs_uperr':0.00068, 'rprs_lowerr':-0.00056,
    #        'rprs_ref':'Demangeon et al. 2021',
    #        },
    #    'd':{
    #        'rp':0.1357, 'rp_uperr':0.0106, 'rp_lowerr':-0.0087,
    #        'rp_ref':'Demangeon et al. 2021',
    #        'mass':0.00610, 'mass_uperr':0.00088, 'mass_lowerr':-0.00088,
    #        'mass_ref':'Demangeon et al. 2021',
    #        'logg_ref':'Demangeon et al. 2021',
    #        'teq_ref':'Demangeon et al. 2021',
    #        'sma':0.0486, 'sma_uperr':0.0018, 'sma_lowerr':-0.0019,
    #        'sma_ref':'Demangeon et al. 2021',
    #        'period':7.4507245, 'period_uperr':8.1e-6, 'period_lowerr':-4.6e-6,
    #        'period_ref':'Demangeon et al. 2021',
    #        't0':2458362.73974, 't0_uperr':0.00031, 't0_lowerr':-0.00040,
    #        't0_ref':'Demangeon et al. 2021',
    #        'inc':88.449, 'inc_uperr':0.058, 'inc_lowerr':-0.111,
    #        'inc_ref':'Demangeon et al. 2021',
    #        'ecc':0.0740, 'ecc_uperr':0.0570, 'ecc_lowerr':-0.0460,
    #        'ecc_ref':'Demangeon et al. 2021',
    #        'omega':180., 'omega_uperr':27., 'omega_lowerr':-50.,
    #        'omega_ref':'Demangeon et al. 2021',
    #        'impact':0.922, 'impact_uperr':0.059, 'impact_lowerr':-0.059,
    #        'impact_ref':'Demangeon et al. 2021',
    #        'trandur':0.840, 'trandur_uperr':0.150, 'trandur_lowerr':-0.200,
    #        'trandur_ref':'Demangeon et al. 2021',
    #        'ars':33.7, 'ars_uperr':1.9, 'ars_lowerr':-1.7,
    #        'ars_ref':'Demangeon et al. 2021',
    #        'rprs':0.04480, 'rprs_uperr':0.00106, 'rprs_lowerr':-0.00100,
    #        'rprs_ref':'Demangeon et al. 2021',
    #         }}

    # g = sscmks['G'] * float(overwrite['L 98-59']['c']['mass'])*sscmks['Mjup'] / \
    #     (float(overwrite['L 98-59']['c']['rp'])*sscmks['Rjup'])**2
    # overwrite['L 98-59']['c']['logg'] = numpy.log10(g)
    # g = sscmks['G'] * float(overwrite['L 98-59']['d']['mass'])*sscmks['Mjup'] / \
    #    (float(overwrite['L 98-59']['d']['rp'])*sscmks['Rjup'])**2
    # overwrite['L 98-59']['d']['logg'] = numpy.log10(g)
    # print('logg before',overwrite['L 98-59']['c']['logg'])
    # print('logg before',overwrite['L 98-59']['d']['logg'])

    # system_info = copy.deepcopy(overwrite['L 98-59'])
    # system_info['L*'] = [0.0113]
    # system_info['L*_uperr'] = [0.0004]
    # system_info['L*_lowerr'] = [-0.0004]
    # for p in ['c','d']:
    #    system_info[p]['teq'] = ['']
    #    system_info[p]['teq_uperr'] = ['']
    #    system_info[p]['teq_lowerr'] = ['']
    #    system_info[p]['teq_ref'] = ['']
    #
    #    logg_derived, logg_lowerr_derived, logg_uperr_derived, logg_ref_derived = \
    #        derive_LOGGplanet_from_R_and_M(system_info, p)
    #    overwrite['L 98-59'][p]['logg'] = logg_derived[0]
    #    overwrite['L 98-59'][p]['logg_lowerr'] = logg_lowerr_derived[0]
    #    overwrite['L 98-59'][p]['logg_uperr'] = logg_uperr_derived[0]
    #    overwrite['L 98-59'][p]['logg_ref'] = logg_ref_derived[0]
    #
    #    teq_derived, teq_lowerr_derived, teq_uperr_derived, teq_ref_derived = \
    #        derive_Teqplanet_from_Lstar_and_sma(system_info, p)
    #    overwrite['L 98-59'][p]['teq'] = teq_derived[0]
    #    overwrite['L 98-59'][p]['teq_lowerr'] = teq_lowerr_derived[0]
    #    overwrite['L 98-59'][p]['teq_uperr'] = teq_uperr_derived[0]
    #    overwrite['L 98-59'][p]['teq_ref'] = teq_ref_derived[0]

    # print('overwrite final',overwrite['L 98-59'])
    # print(derive_LOGGplanet_from_R_and_M(overwrite['L 98-59'], 'c'))
    # print(derive_LOGGplanet_from_R_and_M(overwrite['L 98-59'], 'd'))
    # print('logg after',overwrite['L 98-59']['c']['logg'])
    # print('logg after',overwrite['L 98-59']['d']['logg'])
    # print('')

    # these had FEH* upper limits (e.g. >-0.5 for TOI-2081) that were removed
    # overwrite['TOI-2081'] = {
    #    'FEH*':0.0, 'FEH*_uperr':0.25, 'FEH*_lowerr':-0.25,
    #    'FEH*_units':'[dex]', 'FEH*_ref':'Default to solar metallicity'}
    # (not sure what the second one was; modified output to include target name)

    # Crossfield 2016 has zero for the transit duration
    #  (van Eylen 2016 is the default publication, but it has blank transit duration)
    # TICv8 has 8.184+-0.2.  (made up error bar below)
    overwrite['K2-39'] = {
        'b': {
            'trandur': 8.79,
            'trandur_uperr': 1.0,
            'trandur_lowerr': -1.0,
            'trandur_units': '[hour]',
            'trandur_ref': 'Vandenburg et al. 2016',
        }
    }

    # Galazutdinov 2023 has stellar FEH* = 7.79+-0.12; results in mmw = 4 million for planet
    overwrite['TOI-1408'] = {
        'FEH*': 0.25,
        'FEH*_uperr': 0.06,
        'FEH*_lowerr': -0.06,
        'FEH*_units': '[dex]',
        'FEH*_ref': 'Korth et al. 2024',
    }

    # the new batch of targets (mainly from Burt et al.) has some missing parallaxi
    overwrite['K2-65'] = {
        'dist': (1000.0 / 15.8468),
        'dist_uperr': 0.1,
        'dist_lowerr': -0.1,
        'dist_units': '[pc]',
        'dist_ref': 'Gaia EDR3',
    }
    overwrite['KOI-1257'] = {
        'dist': (1000.0 / 0.4387),
        'dist_uperr': 1000.0,
        'dist_lowerr': 1000.0,
        'dist_units': '[pc]',
        'dist_ref': 'Gaia DR3',
    }
    overwrite['Kepler-565'] = {
        'dist': (1000.0 / 0.7233),
        'dist_uperr': 100.0,
        'dist_lowerr': 100.0,
        'dist_units': '[pc]',
        'dist_ref': 'Gaia DR3',
    }
    overwrite['Kepler-621'] = {
        'dist': (1000.0 / 1.2779),
        'dist_uperr': 100.0,
        'dist_lowerr': -100.0,
        'dist_units': '[pc]',
        'dist_ref': 'Gaia DR3',
    }
    overwrite['Kepler-799'] = {
        'dist': (1000.0 / 0.6499),
        'dist_uperr': 50.0,
        'dist_lowerr': -50.0,
        'dist_units': '[pc]',
        'dist_ref': 'Gaia DR3',
    }
    overwrite['Kepler-808'] = {
        'dist': (1000.0 / 2.8936),
        'dist_uperr': 10.0,
        'dist_lowerr': -10.0,
        'dist_units': '[pc]',
        'dist_ref': 'Gaia DR3',
    }
    overwrite['WTS-2'] = {
        'dist': (1000.0 / 1.4183),
        'dist_uperr': 40.0,
        'dist_lowerr': -40.0,
        'dist_units': '[pc]',
        'dist_ref': 'Gaia DR3',
    }
    overwrite['SPECULOOS-3'] = {
        'dist': (1000.0 / 59.7005),
        'dist_uperr': 0.01,
        'dist_lowerr': -0.01,
        'dist_units': '[pc]',
        'dist_ref': 'Gaia EDR3',
    }

    overwrite['testJup'] = {
        'R*': 1,
        'M*': 1,
        'RHO*': 1.41,
        'LOGG*': 4.438,
        'L*': 1,
        'T*': 5888,
        'FEH*': 0,
        'Jmag': 5,
        'Hmag': 5,
        'Kmag': 5,
        'x': {  # must match added_planet_letter in core.py
            'teq': 880.0,
            'inc': 90,
            't0': 0,
            'sma': 0.1,
            'period': 10,
            'ecc': 0,
            'rp': 1,
            'mass': 1,
            'logg': 3.394,
        },
    }

    # (for Ines's paper on L 98-59)
    # 08/13/2025 : update planet b's parameters to the Cadieux et al. 2025 parameters
    overwrite['L 98-59'] = {
        'R*': 0.3155,
        'R*_uperr': 0.0062,
        'R*_lowerr': -0.0062,
        'R*_ref': 'Cadieux et al. 2025',
        'M*': 0.2923,
        'M*_uperr': 0.0067,
        'M*_lowerr': -0.0067,
        'M*_ref': 'Cadieux et al. 2025',
        'RHO*': 13.0,
        'RHO*_uperr': 1.1,
        'RHO*_lowerr': -0.9,
        'RHO*_ref': 'Cadieux et al. 2025',
        'LOGG*': 4.91,
        'LOGG*_uperr': 0.02,
        'LOGG*_lowerr': -0.02,
        'LOGG*_ref': 'Cadieux et al. 2025',
        'L*': 0.0123,
        'L*_uperr': 0.0009,
        'L*_lowerr': -0.0011,
        'L*_ref': 'Cadieux et al. 2025',
        'T*': 3415,
        'T*_uperr': 60,
        'T*_lowerr': -60,
        'T*_ref': 'Cadieux et al. 2025',
        'AGE*': 4.94,
        'AGE*_uperr': 0.28,
        'AGE*_lowerr': -0.28,
        'AGE*_ref': 'Cadieux et al. 2025',
        'b': {
            'rp': 0.837,
            'rp_uperr': 0.0019,
            'rp_lowerr': -0.0019,
            'rp_ref': 'Cadieux et al. 2025',
            'mass': 0.0014,
            'mass_uperr': 0.0003,
            'mass_lowerr': -0.0003,
            'mass_ref': 'Cadieux et al. 2025',
            'logg_ref': 'Cadieux et al. 2025',
            'teq_ref': 'Cadieux et al. 2025',
            'sma': 0.0223,
            'sma_uperr': 0.0007,
            'sma_lowerr': -0.0007,
            'sma_ref': 'Cadieux et al. 2025',
            'period': 2.2531140,
            'period_uperr': 0.0000004,
            'period_lowerr': -0.0000004,
            'period_ref': 'Cadieux et al. 2025',
            't0': 2458366.17056,
            't0_uperr': 0.00013,
            't0_lowerr': -0.00022,
            't0_ref': 'Cadieux et al. 2025',
            'inc': 88.08,
            'inc_uperr': 0.23,
            'inc_lowerr': -0.20,
            'inc_ref': 'Cadieux et al. 2025',
            'ecc': 0.031,
            'ecc_uperr': 0.017,
            'ecc_lowerr': -0.016,
            'ecc_ref': 'Cadieux et al. 2025',
            'impact': 0.51,
            'impact_uperr': 0.04,
            'impact_lowerr': -0.05,
            'impact_ref': 'Cadieux et al. 2025',
            'trandur': 1.01,
            'trandur_uperr': 0.03,
            'trandur_lowerr': -0.03,
            'trandur_ref': 'Cadieux et al. 2025',
        },
    }

    g = (
        sscmks['G']
        * float(overwrite['L 98-59']['b']['mass'])
        * sscmks['Mjup']
        / (float(overwrite['L 98-59']['b']['rp']) * sscmks['Rjup']) ** 2
    )
    overwrite['L 98-59']['b']['logg'] = numpy.log10(g)

    system_info = copy.deepcopy(overwrite['L 98-59'])
    system_info['L*'] = [0.0123]
    system_info['L*_uperr'] = [0.0009]
    system_info['L*_lowerr'] = [-0.0011]

    system_info['b']['teq'] = ['']
    system_info['b']['teq_uperr'] = ['']
    system_info['b']['teq_lowerr'] = ['']
    system_info['b']['teq_ref'] = ['']

    logg_derived, logg_lowerr_derived, logg_uperr_derived, logg_ref_derived = (
        derive_LOGGplanet_from_R_and_M(system_info, 'b')
    )
    overwrite['L 98-59']['b']['logg'] = logg_derived[0]
    overwrite['L 98-59']['b']['logg_lowerr'] = logg_lowerr_derived[0]
    overwrite['L 98-59']['b']['logg_uperr'] = logg_uperr_derived[0]
    overwrite['L 98-59']['b']['logg_ref'] = logg_ref_derived[0]

    teq_derived, teq_lowerr_derived, teq_uperr_derived, teq_ref_derived = (
        derive_Teqplanet_from_Lstar_and_sma(system_info, 'b')
    )
    overwrite['L 98-59']['b']['teq'] = teq_derived[0]
    overwrite['L 98-59']['b']['teq_lowerr'] = teq_lowerr_derived[0]
    overwrite['L 98-59']['b']['teq_uperr'] = teq_uperr_derived[0]
    overwrite['L 98-59']['b']['teq_ref'] = teq_ref_derived[0]

    return overwrite


# -------------------------------------------------------------------
