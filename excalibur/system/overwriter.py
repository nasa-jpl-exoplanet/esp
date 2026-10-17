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

    overwrite['K2-33'] = {
        'FEH*': 0.0,
        'FEH*_uperr': 0.13,
        'FEH*_lowerr': -0.14,
        'FEH*_units': '[dex]',
        'FEH*_ref': 'Mann et al. 2016',
    }    
    overwrite['KELT-1'] = {
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
    overwrite['WASP-12'] = {
        'b': {
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
    overwrite['WASP-43'] = {
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
    overwrite['K2-55'] = {
        'FEH*': 0,
        'FEH*_uperr': 0.25,
        'FEH*_lowerr': -0.25,
        'FEH*_units': '[dex]',
        'FEH*_ref': "Kyle's best guess",
        # our assumed M-R relation gives 0.0497
    }
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
        # --< GMR: JWST NIRSPEC
        'R*': 0.215,
        'R*_uperr': 0.008,
        'R*_lowerr': -0.008,
        'R*_ref': 'https://doi.org/10.3847/1538-3881/ac1584',
        'T*': 3250,
        'T*_uperr': 100,
        'T*_lowerr': -100,
        'T*_ref': 'https://doi.org/10.3847/1538-3881/ac1584',
        'LOGG*': 5.026,
        'LOGG*_uperr': 0.04,
        'LOGG*_lowerr': -0.04,
        'LOGG*_ref': 'https://doi.org/10.3847/1538-3881/ac1584',
        'FEH*': 0.29,
        'FEH*_uperr': 0.12,
        'FEH*_lowerr': -0.12,
        'FEH*_ref': 'https://doi.org/10.3847/1538-3881/ac1584',
        "b": {
            'period': 1.580404341,
            'period_uperr': 0.000000079,
            'period_lowerr': -0.000000079,
            'period_ref': 'https://doi.org/10.3847/2041-8213/ad7fef',
            't0': 2460143.93029503,
            't0_uperr': 0.00000358,
            't0_lowerr': -0.00000358,
            't0_ref': 'https://doi.org/10.3847/2041-8213/ad7fef',
            'inc': 89.32,
            'inc_uperr': 0.03,
            'inc_lowerr': -0.03,
            'inc_ref': 'https://doi.org/10.3847/2041-8213/ad7fef',
            'sma': 0.015267716541101813,
            'sma_uperr': 2e-5,
            'sma_lowerr': -2e-5,
            'sma_ref': 'https://doi.org/10.3847/2041-8213/ad7fef',
            # >--
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
        },
    }

    overwrite["GJ 357"] = {
        # --< CB: JWST NIRSPEC
        "b": {
            't0': 60282.3413,
            't0_uperr': 3e-5,
            't0_lowerr': -3e-5,
            't0_ref': 'https://iopscience.iop.org/article/10.3847/1538-3881/adee92/pdf',
        },
    }

    overwrite["K2-18"] = {
        # --< CB: JWST NIRSPEC
        "b": {
            't0': 59964.969453,
            't0_uperr': 0.0001,
            't0_lowerr': -0.0001,
            't0_ref': 'https://iopscience.iop.org/article/10.3847/2041-8213/acf577/pdf',
        },
    }

    overwrite["HAT-P-11"] = {
        # --< CB: JWST NIRSPEC
        'R*': 0.872,
        'b': {
            'period': 4.88781501,
            'period_uperr': 6.8e-7,
            'period_lowerr': -6.8e-7,
            'period_ref': 'Winn et al 2010',
            'inc': 89.17,
            'inc_uperr': 0.46,
            'inc_lowerr': -0.60,
            'inc_ref': 'Winn et al 2010',
            'sma': 0.05196,
            'sma_ref': 'Winn et al 2010',
            't0': 60504.965,
            't0_uperr': 0.0001,
            't0_lowerr': -0.0001,
        },
    }

    overwrite["HAT-P-14"] = {
        # --< CB: JWST NIRSPEC
        "b": {
            't0': 59729.203001524,
            't0_uperr': 0.000010,
            't0_lowerr': -0.000010,
            't0_ref': 'https://iopscience.iop.org/article/10.1088/1538-3873/aca3d3/pdf',
        },
    }

    overwrite["HAT-P-26"] = {
        # --< CB: JWST NIRSPEC
        "b": {
            'period': 4.2344923,
            'period_ref': 'https://iopscience.iop.org/article/10.3847/1538-3881/ae0929/pdf',
            't0': 60110.30657380267,
            't0_uperr': 0.000036,
            't0_lowerr': -0.000036,
            't0_ref': 'https://iopscience.iop.org/article/10.3847/1538-3881/ae0929/pdf',
            'sma': 0.0459,
            'sma_ref': 'https://iopscience.iop.org/article/10.3847/1538-3881/ae0929/pdf',
            'inc': 87.8,
            'inc_ref': 'https://iopscience.iop.org/article/10.3847/1538-3881/ae0929/pdf',
            'ecc': 0.0,
            'ecc_ref': 'https://iopscience.iop.org/article/10.3847/1538-3881/ae0929/pdf',
        },
    }

    overwrite["HD 15337"] = {
        # --< CB: JWST NIRSPEC
        "b": {
            't0': 60142.14151,
            't0_uperr': 0.00011,
            't0_lowerr': -0.00011,
            't0_ref': 'https://ntrs.nasa.gov/api/citations/20260005584/downloads/20260005584-TOI_402_01.pdf',
        },
        "c": {
            't0': 60149.25631,
            't0_uperr': 0.00012,
            't0_lowerr': -0.00012,
            't0_ref': 'https://arxiv.org/pdf/2602.22327',
            'inc': 88.0,
            'inc_uperr': 0.1,
            'inc_lowerr': -0.1,
            'inc_ref': 'https://arxiv.org/pdf/2602.22327',
            'sma': 0.1113,
            'sma_ref': 'https://arxiv.org/pdf/2602.22327',
        },
    }

    overwrite["LHS 3844"] = {
        # --< CB: JWST NIRSPEC
        "b": {
            't0': 58325.22129184415,
            't0_uperr': 0.0001,
            't0_lowerr': -0.0001,
            't0_ref': 'file:///Users/cbernard/Downloads/s41550-026-02860-3-3.pdf',
        },
    }

    overwrite["LTT 1445 A"] = {
        # --< CB: JWST NIRSPEC
        "c": {
            't0': 58412.076569440804,
            't0_uperr': 0.00078,
            't0_lowerr': -0.00074,
            't0_ref': 'CB',
        },
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

    # for the newly added comfirmed-planet Ariel targets, some metallicities are missing
    #  oh that's funny. the Chen 2021 compilation has zero for these (with no error bar)
    #   but the source is listed as flag 5, which is the Exoplanet Archive.  fake news!
    
    overwrite['HATS-52'] = {
        'FEH*': -0.09,
        'FEH*_uperr': 0.17,
        'FEH*_lowerr': -0.17,
        'FEH*_units': '[dex]',
        'FEH*_ref': 'Magrini et al. 2022',
    }

    overwrite['K2-129'] = {
        'FEH*': 0.105,
        'FEH*_uperr': 0.235,
        'FEH*_lowerr': -0.235,
        'FEH*_units': '[dex]',
        'FEH*_ref': 'Hardagree-Ullman et al. 2020',
    }

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
    # wikipedia has it at 980pc from 2011 schneider site
    # from buchhave 2011 discovery paper
    # the paper says it is from Girardi isochrone fitting
    overwrite['Kepler-14'] = {
        'dist': 980.0,
        'dist_uperr': 100.0,
        'dist_lowerr': -100.0,
        'dist_units': '[pc]',
        'dist_ref': 'Buchhave et al. 2011',
    }

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

    overwrite['TOI-206'] = {
        'FEH*': 0.057,
        'FEH*_uperr': 0.25,
        'FEH*_lowerr': -0.25,
        'FEH*_units': '[dex]',
        'FEH*_ref': 'Sprague et al. 2022',
    }

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
            'rp': 0.0747,
            'rp_uperr': 0.0017,
            'rp_lowerr': -0.0017,
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
    overwrite['L 98-59']['b']['logg'] = float(logg_derived[0])
    overwrite['L 98-59']['b']['logg_lowerr'] = float(logg_lowerr_derived[0])
    overwrite['L 98-59']['b']['logg_uperr'] = float(logg_uperr_derived[0])
    overwrite['L 98-59']['b']['logg_ref'] = logg_ref_derived[0]

    teq_derived, teq_lowerr_derived, teq_uperr_derived, teq_ref_derived = (
        derive_Teqplanet_from_Lstar_and_sma(system_info, 'b')
    )
    overwrite['L 98-59']['b']['teq'] = float(teq_derived[0])
    overwrite['L 98-59']['b']['teq_lowerr'] = float(teq_lowerr_derived[0])
    overwrite['L 98-59']['b']['teq_uperr'] = float(teq_uperr_derived[0])
    overwrite['L 98-59']['b']['teq_ref'] = teq_ref_derived[0]

    # nov.2025 addition of 450 new systems from the Archive
    #  3 are missing distances
    overwrite['Kepler-1676'] = {
        'dist': (1000.0 / 1.1204),
        'dist_uperr': 50.0,
        'dist_lowerr': -50.0,
        'dist_units': '[pc]',
        'dist_ref': 'Gaia EDR3',
    }
    overwrite['Kepler-477'] = {
        'dist': (1000.0 / 2.1511),
        'dist_uperr': 5.0,
        'dist_lowerr': -5.0,
        'dist_units': '[pc]',
        'dist_ref': 'Gaia EDR3',
    }
    overwrite['Kepler-478'] = {
        'dist': (1000.0 / 1.3568),
        'dist_uperr': 200.0,
        'dist_lowerr': -200.0,
        'dist_units': '[pc]',
        'dist_ref': 'Gaia EDR3',
    }
    #  1 is missing star mass
    #   actually TESS input catalog (Stassun 2019 and Paegert 2021)
    #    coincidentally has it as 1.000, same as our default value
    #   other vizier refs have it as 0.939 and 0.948
    overwrite['EPIC 212624936'] = {
        'M*': 1.000,
        'M*_uperr': 0.25,
        'M*_lowerr': -0.25,
        'M*_ref': 'Staussun et al. 2019',
        'RHO*': 1.13,  # from R*=0.93
        'RHO*_uperr': 0.25,
        'RHO*_lowerr': -0.25,
        'RHO*_ref': 'derived from M*,R*',
    }

    # missing t0 for Kepler-453 b (so it comes up as no planets)
    # strange, it actually is in the default publication. binary oddness
    overwrite['Kepler-453'] = {
        'b': {
            't0': 2455069.020,
            't0_uperr': 0.054,
            't0_lowerr': -0.054,
            't0_ref': 'Welsh et al. 2015',  # Table 3
        }
    }

    # GMR: BJT_TBD_to_MJD_UTC
    overwrite['WASP-47'] = {
        'b': {
            't0': 56982.4774228593,
        },
        'd': {
            't0': 56987.87473479712,
        },
    }

    return overwrite


# -------------------------------------------------------------------
