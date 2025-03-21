'''cerberus forward_model ds'''

# Heritage code shame:
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments,too-many-branches,too-many-lines,too-many-locals,too-many-positional-arguments,too-many-statements

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as cst
from scipy.interpolate import interp1d as itp
import logging

# import sys
# import pytensor
# from pytensor.ifelse import ifelse
from pytensor import tensor

import excalibur

# 2/13/25 Geoff: this line has been here for a while. maybe not required anymore? We'll see...
# import excalibur.cerberus.forward_model  # is this needed for the context updater?
import excalibur.system.core as syscore
from excalibur.util.cerberus import crbce, getmmw

# -- GLOBAL CONTEXT FOR PYMC DETERMINISTICS ---------------------------------------------
from collections import namedtuple

log = logging.getLogger(__name__)

CONTEXT = namedtuple(
    'CONTEXT',
    [
        'cleanup',
        'model',
        'p',
        'solidr',
        'orbp',
        'tspectrum',
        'xsl',
        'spc',
        'modparlbl',
        'hzlib',
        'fixedParams',
    ],
)
ctxt = CONTEXT(
    cleanup=None,
    model=None,
    p=None,
    solidr=None,
    orbp=None,
    tspectrum=None,
    xsl=None,
    spc=None,
    modparlbl=None,
    hzlib=None,
    fixedParams=None,
)


def ctxtupdt(
    cleanup=None,
    model=None,
    p=None,
    solidr=None,
    orbp=None,
    tspectrum=None,
    xsl=None,
    spc=None,
    modparlbl=None,
    hzlib=None,
    fixed_params=None,
):
    '''
    G. ROUDIER: Update context
    '''
    excalibur.cerberus.forward_model.ctxt = CONTEXT(
        cleanup=cleanup,
        model=model,
        p=p,
        solidr=solidr,
        orbp=orbp,
        tspectrum=tspectrum,
        xsl=xsl,
        spc=spc,
        fixedParams=fixed_params,
        modparlbl=modparlbl,
        hzlib=hzlib,
    )
    return


# ----------- --------------------------------------------------------
# -- CERBERUS MODEL -- -----------------------------------------------
def crbmodel(
    mixratio,
    rayleigh,
    cloudtp,
    rp0,
    orbp,
    xsecs,
    qtgrid,
    temp,
    wgrid,
    lbroadening=False,
    lshifting=False,
    # nlevels=100,
    # TEMPORARY REDUCTION IN ATMOS RESOLUTION WHILE DEBUGGING PYMC
    nlevels=7,
    # increase the number of scale heights from 15 to 20, to match the Ariel forward model
    Hsmax=20.0,
    solrad=10.0,
    hzlib=None,
    hzp=None,
    hzslope=-4.0,
    hztop=None,
    hzwscale=1e0,
    cheq=None,
    logx=False,
    pnet='b',
    break_down_by_molecule=False,
    verbose=False,
    debug=False,
):
    '''
    G. ROUDIER: Cerberus forward model probing up to 'Hsmax' scale heights from solid
    radius solrad evenly log divided amongst nlevels steps
    '''

    # these used to be default parameters above, but are dangerous-default-values
    # note that these are also defined in cerberus/core/myxsecs()
    #  maybe put them inside runtime/ops.xml to ensure consistency?
    cialist = ['H2-H', 'H2-H2', 'H2-He', 'He-H']
    xmollist = ['TIO', 'H2O', 'H2CO', 'HCN', 'CO', 'CO2', 'NH3', 'CH4']

    ssc = syscore.ssconstants(mks=True)
    pgrid = np.arange(
        np.log(solrad) - Hsmax,
        np.log(solrad) + Hsmax / nlevels,
        Hsmax / (nlevels - 1),
    )
    # print('pgrid before exponential',pgrid)
    pgrid = np.exp(pgrid)
    # dp = np.diff(pgrid[::-1])
    p = pgrid[::-1]
    print('pressure', len(p), p)
    # print('delta-pressure',len(dp),dp)
    dPoverP = (p[1] - p[0]) / p[0]

    print()
    if not mixratio:
        # print('SHOULD BE HERE')
        if cheq is None:
            log.warning('neither mixratio nor cheq are defined')
        mixratio, fH2, fHe = crbce(
            p, temp, C2Or=cheq['CtoO'], X2Hr=cheq['XtoH'], N2Or=cheq['NtoO']
        )
        mmw, fH2, fHe = getmmw(mixratio, protosolar=False, fH2=fH2, fHe=fHe)
    else:
        # print('SHOULD NOT BE HERE')
        mmw, fH2, fHe = getmmw(mixratio)
    mmw = mmw * cst.m_p  # [kg]
    print('mmw', mmw.eval() * 6.022e26)

    Hs = (
        cst.Boltzmann
        * temp
        / (mmw * 1e-2 * (10.0 ** float(orbp[pnet]['logg'])))
    )  # [m]

    # when the Pressure grid is log-spaced, rdz is a constant
    #  drop dz[] and dzprime[] arrays and just use this constant instead
    rdz = abs(Hs / 2.0 * np.log(1.0 + dPoverP))

    z = [0]
    # dz = []
    # addz = []
    # for press, dpress in zip(p[:-1], dp):
    # print('dp/p',dpress/press,np.log(1. + dpress/press),dPoverP)
    # for press in p[:-1]:
    # rdz = abs(Hs/2.*np.log(1. + dpress/press))
    # rdz = abs(Hs/2.*np.log(1. + dPoverP))
    # print('press,rdz',press,rdz.eval())
    # if addz:
    # dz.append(addz[-1]/2. + rdz)
    # else:
    # tem = tensor.dscalar()
    # print('tem',tem)
    # tem = 2.*rdz
    # print('tem',tem)
    # dz.append(tem)
    # dz.append(2.*rdz)
    # print('temp',dz[-1])
    # addz.append(2.*rdz)
    # z.append(z[-1]+addz[-1])
    for _ in p[:-1]:
        z.append(z[-1] + 2.0 * rdz)
    # dz.append(addz[-1])
    # print()
    # print('len check on z',len(z))
    # print('len check on dz',len(dz))
    # print()
    # print('z going into gettau',z)
    # print('z going into gettau',[zlist.eval() for zlist in z[1:]])
    # print('dz going into gettau',dz)
    # print('dz going into gettau',[zlist.eval() for zlist in dz])
    # print()

    z = np.linspace(0, len(p) - 1, len(p))
    z = 2 * rdz * z
    # print('z', z.eval())
    # print('z redo',[zlist.eval() for zlist in z[1:]])

    # simplify dz.  will help tremendously below, where tensor keeps crashing
    # dz = np.array([dz[0].eval()]*len(dz))
    # print('new dz',dz)
    # print(' rdz',rdz.eval()*2)
    dz = 2 * rdz
    # print('new dz', dz.eval())

    rho = p * 1e5 / (cst.Boltzmann * temp)
    tau, tau_by_molecule, wtau = gettau(
        xsecs,
        qtgrid,
        temp,
        mixratio,
        z,
        dz,
        rho,
        rp0,
        p,
        wgrid,
        lbroadening,
        lshifting,
        cialist,
        fH2,
        fHe,
        xmollist,
        rayleigh,
        hzlib,
        hzp,
        hzslope,
        hztop,
        hzwscale=hzwscale,
        debug=debug,
    )

    print('tau', tau.eval())
    for molecule in tau_by_molecule:
        print('tau', molecule, tau_by_molecule[molecule].eval())

    if not break_down_by_molecule:
        tau_by_molecule = {}
    molecules = tau_by_molecule.keys()
    # SEMI FINITE CLOUD ------------------------------------------------------------------
    reversep = np.array(p[::-1])
    selectcloud = p > 10.0**cloudtp
    blocked = False
    if all(selectcloud):
        tau = tau * 0
        for molecule in molecules:
            tau_by_molecule[molecule] = tau_by_molecule[molecule] * 0
        blocked = True
        pass
    if not all(~selectcloud) and not blocked:
        cloudindex = np.max(np.arange(len(p))[selectcloud]) + 1
        for index in np.arange(wtau.size):
            # asdf
            print('reversep', reversep)
            print('reversep', reversep.shape)
            print('tau check', index, tau[:, index].eval())
            print('tau check', index, tau[:, index].eval().shape)
            myspl = itp(reversep, tau[:, index].eval())
            # print('myspl',myspl)  # interp1d object
            print(' cloudtp', cloudtp)
            powcheck = 10.0**cloudtp
            print(' powcheck', powcheck)
            print('  indices', cloudindex, index)
            print(' DOES INTERP WORK?', myspl(10.0**cloudtp))
            # tau[cloudindex, index] = myspl(10.0**cloudtp)  # fails. says to use .set or .inc
            print('  tau[]', tau[cloudindex, index].eval())
            tau = tau[cloudindex, index].set(myspl(10.0**cloudtp))
            print('  tau[]', tau[cloudindex, index].eval())
            # tau[:cloudindex, index] = 0.0
            tau = tau[:cloudindex, index].set(0.0)
            for molecule in molecules:
                myspl = itp(
                    reversep, tau_by_molecule[molecule][:, index].eval()
                )
                tau_by_molecule[molecule] = tau_by_molecule[molecule][
                    cloudindex, index
                ].set(myspl(10.0**cloudtp))
                tau_by_molecule[molecule] = tau_by_molecule[molecule][
                    :cloudindex, index
                ].set(0.0)
            pass
        ctpdpress = 10.0**cloudtp - p[cloudindex]
        ctpdz = abs(Hs / 2.0 * np.log(1.0 + ctpdpress / p[cloudindex]))
        rp0 += z[cloudindex] + ctpdz
        pass
    # atmdepth = (2e0 * np.array(
    #    np.asmatrix((rp0 + np.array(z)) * np.array(dz))
    #    * np.asmatrix(1.0 - np.exp(-tau))).flatten())
    matrix1 = (rp0 + z) * dz
    matrix2 = 1.0 - tensor.exp(-tau)
    print('matrix1 shape', matrix1.eval().shape)
    print('matrix2 shape', matrix2.eval().shape)
    atmdepth = 2e0 * tensor.nlinalg.matrix_dot(matrix1, matrix2)
    print('result shape', atmdepth.eval().shape)
    #  flatten is not needed I think.  it's already 1-D, no?
    # atmdepth = tensor.flatten(atmdepth)
    # print('result shape',atmdepth.eval().shape)
    # atmdepth = (2e0 * np.array(
    #    np.asmatrix((rp0 + np.array(z)) * np.array(dz))
    #    * np.asmatrix(1.0 - np.exp(-tau))).flatten())
    model = (rp0**2 + atmdepth) / (orbp['R*'] * ssc['Rsun']) ** 2
    print('final model result!', model.eval())
    print('final model result!', model.eval().shape)
    models_by_molecule = {}
    for molecule in molecules:
        # atmdepth = (
        #    2e0
        #    * np.array(
        #        np.asmatrix((rp0 + np.array(z)) * np.array(dz))
        #        * np.asmatrix(1.0 - np.exp(-tau_by_molecule[molecule]))
        # ).flatten()
        # )
        atmdepth = 2e0 * tensor.nlinalg.matrix_dot(
            (rp0 + z) * dz, 1.0 - tensor.exp(-tau_by_molecule[molecule])
        )
        models_by_molecule[molecule] = (rp0**2 + atmdepth) / (
            orbp['R*'] * ssc['Rsun']
        ) ** 2
        models_by_molecule[molecule] = models_by_molecule[molecule][::-1]
    if verbose:
        noatm = rp0**2 / (orbp['R*'] * ssc['Rsun']) ** 2
        noatm = tensor.nanmin(model)
        rp0hs = tensor.sqrt(noatm * (orbp['R*'] * ssc['Rsun']) ** 2)

        fig, ax = plt.subplots(figsize=(10, 6))
        axes = [ax, ax.twinx(), ax.twinx()]
        fig.subplots_adjust(left=0.125, right=0.775)
        axes[-1].spines['right'].set_position(('axes', 1.2))
        axes[-1].set_frame_on(True)
        axes[-1].patch.set_visible(False)
        axes[0].plot(wtau, 1e2 * model)
        axes[0].plot(wtau, model * 0 + 1e2 * noatm, '--')
        axes[0].set_xlabel('Wavelength $\\lambda$[$\\mu m$]')
        axes[0].set_ylabel('Transit Depth [%]')
        axes[0].get_yaxis().get_major_formatter().set_useOffset(False)
        yaxmin, yaxmax = axes[0].get_ylim()
        ax2min = (
            np.sqrt(1e-2 * yaxmin) * orbp['R*'] * ssc['Rsun'] - rp0hs
        ) / Hs
        ax2max = (
            np.sqrt(1e-2 * yaxmax) * orbp['R*'] * ssc['Rsun'] - rp0hs
        ) / Hs
        axes[-1].set_ylabel('Transit Depth Modulation [Hs]')
        axes[-1].set_ylim(ax2min, ax2max)
        axes[-1].get_yaxis().get_major_formatter().set_useOffset(False)
        axes[1].set_ylabel('Transit Depth Modulation [ppm]')
        axes[1].set_ylim(
            1e6 * (1e-2 * yaxmin - noatm), 1e6 * (1e-2 * yaxmax - noatm)
        )
        axes[1].get_yaxis().get_major_formatter().set_useOffset(False)
        if logx:
            plt.semilogx()
            plt.xlim([np.min(wtau), np.max(wtau)])
            pass
        plt.show()
        pass
    if break_down_by_molecule:
        return model[::-1], models_by_molecule
    return model[::-1]


# --------------------------- ----------------------------------------
# -- TAU -- ----------------------------------------------------------
def gettau(
    xsecs,
    qtgrid,
    temp,
    mixratio,
    z,
    dz,
    rho,
    rp0,
    p,
    wgrid,
    lbroadening,
    lshifting,
    cialist,
    fH2,
    fHe,
    xmollist,
    rayleigh,
    hzlib,
    hzp,
    hzslope,
    hztop,
    isothermal=True,
    hzwscale=1e0,
    debug=False,
):
    '''
    G. ROUDIER: Builds optical depth matrix
    '''
    # SPHERICAL SHELL (PLANE-PARALLEL REMOVED) -------------------------------------------
    # MATRICES INIT ------------------------------------------------------------------
    Nzones = len(p)
    # print('z (at start of gettau)',[zz.eval() for zz in z])
    print('z (at start of gettau)', z.eval())
    # tau = np.zeros((len(z), wgrid.size))
    tau = np.zeros((Nzones, wgrid.size))
    print('tau shape at the top', tau.shape)
    tau_by_molecule = {}
    # DL ARRAY, Z VERSUS ZPRIME ------------------------------------------------------
    dlarray = []
    # these are coming in as lists created with append loops.  want arrays
    zprime = z
    # zprime = np.array(z)
    # dzprime = np.array(dz)
    # zprime = z
    # dzprime = dz
    for iz, thisz in enumerate(z):
        print()
        print('LOOP iz,thisz', iz, thisz.eval())
        # dltemp = (rp0 + zprime + dzprime)**2 - (rp0 + thisz)**2
        # print('rp0',rp0)
        # print('zprime',zprime)  # messed up.  first element (zero) is diff
        # print('dzprime',dzprime) # messed up. first and last are diff (mul vs add)
        # print('dzprime',[dzp.eval() for dzp in dzprime]) # ok
        # print('dltemp',dltemp)
        # print('dltemp',[dlt.eval()for dlt in dltemp]) # ok
        # print('dltemp',dltemp.eval())  # fails.  no attribute eval for a numpy array. huh? oh i see it's an array of tensors ig
        # dl = tensor.sqrt((rp0))
        # print('dl sqrt works1?',dl.eval()) # ok
        # pytensor.dprint(zprime)   # same as a regular print; doesn't show chart
        # pytensor.dprint(dzprime)  # same as a regular print; doesn't show chart
        # dl = tensor.sqrt((dzprime))    # <--- *** THIS FAILS ***
        # dl = [tensor.sqrt(dzp) for dzp in dzprime]
        # print('dl sqrt works1b?',[ddd.eval() for ddd in dl])  # ok
        # dl = [tensor.sqrt(dzp) for dzp in zprime]
        # print('dl sqrt works1c?',dl) # ok

        # dl11 = np.sqrt((rp0 + zprime + dzprime)**2)
        # print('dl sqrt works first half?',dl11)
        # dl22 = np.sqrt((rp0 + thisz)**2)
        # print('dl sqrt works second part?',dl22)
        # print('   failing?  loop it')
        # testt = []
        # for dl11element in dl11:
        #    testt.append(dl11element - dl22)
        # print('   schould work',testt)
        # print('   worked',[t.eval() for t in testt])
        # print('   failing?')
        # testt = dl11 - dl22  # fails
        # print('   schould fail')

        # dl = tensor.sqrt((rp0 + zprime + dzprime)**2 - (rp0 + thisz)**2) # fails

        #  hmm weird.  this works on first pass of the loop, but not second
        # dl = np.sqrt((rp0 + zprime + dzprime)**2 - (rp0 + thisz)**2)  # second iteration trouble
        # dl = np.sqrt((rp0 + zprime + dzprime)**2) - np.sqrt((rp0 + thisz)**2)  # fails
        #  try making the list (first half) into an array first
        # dl = np.sqrt(np.array((rp0 + zprime + dzprime))**2) - np.sqrt((rp0 + thisz)**2) # fails
        # it's as if it's treating first part as a list.  so loop through it
        # dl = []
        # for zpr,dzpr in zip(zprime,dzprime):
        # dl = np.zeros(len(zprime))
        # dl = np.zeros(Nzones)
        print(' zprime', zprime.eval())
        # print(' dz', dz.eval())
        # print('dl',dl)
        # WHAT ABOUT THIS NOW? Yes!
        # dl = np.sqrt((rp0 + zprime + dz) ** 2)  # works!
        dl = np.sqrt((rp0 + zprime + dz) ** 2 - (rp0 + thisz) ** 2)
        # for izz in range(Nzones):
        #    print(' subloop',izz,dl)
        #    dl[izz] = rp0
        #    print('  zprime!!!!',zprime)
        #    print('  zprime!!!!',zprime[0])
        #    print('  zprime!!!!',zprime[0])
        #    print('  zprime!!!!',zprime[1])
        #    dl[izz] = zprime[izz]
        #    # print('  dzprime!!!!',dz)
        #    # print('  dzprime!!!!',dz[0])
        #    # print('  dzprime!!!!',dzprime[0].eval())
        #    # print('  dzprime!!!!',dzprime[1])
        #    # print('  dzprime!!!!',dzprime[1].eval())
        #    print('  dz',dz.eval())
        #    print('   izz',izz)
        #    dl[izz] = dz
        #    print('   izz',izz)
        #    dl[izz] = rp0 + zprime[izz] + dz
        #    dl[izz] = np.sqrt((rp0 + zprime[izz] + dz)**2)
        #    dl[izz] = np.sqrt((rp0 + zprime[izz] + dz)**2 - (rp0 + thisz)**2)
        # print('dl sqrt works subtract both',dl)
        # print('dl sqrt works subtract both',[d.eval() for d in dl])
        print('dl sqrt works subtract both', dl.eval()/dz.eval())

        # for d in dl: print('loop check1',d.eval())
        # for id,d in enumerate(dl): print('loop check2',id,d.eval())

        # dl =(rp0 + zprime + dzprime)**2 - (rp0 + thisz)**2
        # print('dl sqrt works subtract sqrs?',dl)

        # print('dl works py?',[d.eval() for d in dl])
        # dl[:iz] = tensor.dscalar(0e0)  # fails.  must be a string
        #  oh! it expects a string because it's the name.  how to put in a value?
        # print('iz,len(dl)',iz,len(dl))
        # print('dl1',dl)
        # dl[:iz] = tensor.dscalar()  # I think you can leave the name blank
        # print('dl2',dl)
        # dl[:iz] = 0.
        # dl[:iz+1] = 0.  # FAILS NOW (with dz/z update)
        # print('dl after zeros at start',dl)
        # print('dl just the zeros check',dl[:iz+1])
        # dl = tensor.sqrt((rp0 + zprime + dzprime)**2 - (rp0 + thisz)**2)
        # dl[:iz] *= 0.
        # print('dl4',dl)
        # oof nasty bug here!
        #  sometimes equal terms are off by the instrument precision
        #  so sqrt(1e15 - 1e15) = sqrt(-1) = NaN
        # take absolute value, just to be sure there's no problem
        # dl[iz:] = dl[iz:] - np.sqrt((rp0 + zprime[iz:])**2 - (rp0 + thisz)**2)
        #  convert this previous line to a tensor form
        #   uh hold on, what's a tensor here?
        # print('dl',dl)  # array of tensors
        # print('zprime',zprime)  # array of tensors, except the first one is zero!
        # print('thisz',thisz)  # float
        # print('rp0',rp0)  # float
        # print()
        # dl = np.sqrt((rp0 + zprime + dzprime)**2 - (rp0 + thisz)**2)
        # print('dl works?',[d.eval() for d in dl])
        # test = (rp0 + thisz)**2
        # print('test works1?',test)
        # print('test works1?',[d.eval() for d in test]) # fails. it's a float not a list

        # test = (rp0 + zprime[iz:])**2
        # print('test works2?',test)
        # print('test works2?',[d.eval() for d in test[1:]])  # no work for first cell
        # test = (rp0 + zprime[iz:])**2 - (rp0 + thisz)**2
        # print('test works3?',test)
        # print('test works3?',[d.eval() for d in test])

        # print(' *** thisz vs zprime check! ***',zprime[iz],thisz)
        # test = (np.abs((rp0 + zprime[iz+1:])**2 - (rp0 + thisz)**2))
        # print('test works4?',test)
        # test = (np.abs((rp0 + zprime[iz:])**2 - (rp0 + thisz)**2))
        # print('test works5?',test)
        # test = np.sqrt(np.abs((rp0 + zprime[iz+1:])**2 - (rp0 + thisz)**2))
        # print('test works6?',test)
        # print('test works4?',[d.eval() for d in test])

        # print()
        # dl[iz:] = dl[iz:] - np.sqrt(np.abs((rp0 + zprime[iz:])**2 - (rp0 + thisz)**2))
        # if 0:
        #    for did, d in enumerate(dl):
        #        print('loop check2', id, d.eval())
        #        if did < iz + 1:
        #            dl[did] = dl[did] * 0.0
        #            print('hup', dl[did].eval())
        #        else:
        #            dl[did] = dl[did] - np.sqrt(
        #                np.abs((rp0 + zprime[did]) ** 2 - (rp0 + thisz) ** 2)
        #            )
        #            print('dop', dl[did].eval())

        # (above) dl = np.sqrt((rp0 + zprime + dz) ** 2 - (rp0 + thisz) ** 2)

        dl = dl - np.sqrt(np.abs((rp0 + zprime) ** 2 - (rp0 + thisz) ** 2))
        # print('dl with negative still',[dd.eval() for dd in dl])
        print('dl with negative still', dl.eval()/dz.eval())

        # dl0 = np.sqrt(np.abs((rp0 + zprime) ** 2 - (rp0 + thisz) ** 2))
        # print(' dl0 old', dl0.eval()/dz.eval())
        # dl0 = np.sqrt(np.max([zprime*0,(rp0 + zprime)**2 - (rp0 + thisz)**2])) # fails
        dl = np.sqrt(
            tensor.max(
                [zprime * 0, (rp0 + zprime + dz) ** 2 - (rp0 + thisz) ** 2], axis=0
            )
        )
        dl0 = np.sqrt(
            tensor.max(
                [zprime * 0, (rp0 + zprime) ** 2 - (rp0 + thisz) ** 2], axis=0
            )
        )
        print(' dl1 ', dl.eval()/dz.eval())
        print(' dl0 ', dl0.eval()/dz.eval())
        dl = dl - dl0
        print(' dl new', dl.eval()/dz.eval())
        # asdfasdf

        # dl = ifelse(zprime > thisz, dl - dl0, dl * 0)  # fails
        # print('dl with nan fixed?',dl.eval())

        # if zprime > thisz:  # fails
        #    dl = dl - dl0
        # else:
        #    dl = dl * 0
        # print(' dl old part',iz,dl[:iz+1])
        # print(' dl new part',iz,[d.eval() for d in dl[iz+1:]])
        dlarray.append(dl)
    print('MADE IT!!!')
    # print('dlarray',[d.eval() for d in dlarray])
    # print('dlarray',dlarray.eval())  # fails.  it's a list
    # print()
    dlarray = np.array(dlarray)
    # print('dlarray', [d.eval() for d in dlarray])
    # print('dlarray',dlarray.eval()) still fails even after change list to array

    # dlarray is a list of tensors which are 1-d arrays (yuch!)
    # convert it to a 2-d tensor
    dlarraymod = []
    for dla in dlarray:
        dlarraymod.append(dla.eval())
    dlarray = np.array(dlarraymod)
    print('dlarray up top', dlarray)

    # GAS ARRAY, ZPRIME VERSUS WAVELENGTH  -------------------------------------------
    for elem in mixratio:
        # tau_by_molecule[elem] = np.zeros((len(z), wgrid.size))
        mmr = 10.0 ** (mixratio[elem] - 6.0)
        # Fake use of xmollist due to changes in xslib v112
        # THIS HAS TO BE FIXED
        # if elem not in xmollist:
        if not xmollist:
            # HITEMP/HITRAN ROTHMAN ET AL. 2010 --------------------------------------
            sigma, lsig = absorb(
                xsecs[elem],
                qtgrid[elem],
                temp,
                p,
                mmr,
                lbroadening,
                lshifting,
                wgrid,
                debug=False,
            )
            # sigma = np.array(sigma)  # cm^2/mol
            if True in (sigma < 0):
                sigma[sigma < 0] = 0e0
            if True in ~np.isfinite(sigma):
                sigma[~np.isfinite(sigma)] = 0e0
            sigma = sigma * 1e-4  # m^2/mol
            pass
        else:
            # EXOMOL HILL ET AL. 2013 ------------------------------------------------
            sigma, lsig = getxmolxs(temp, xsecs[elem])
            # sigma = np.array(sigma)   # cm^2/mol
            if True in (sigma < 0):
                sigma[sigma < 0] = 0e0
            if True in ~np.isfinite(sigma):
                sigma[~np.isfinite(sigma)] = 0e0
            # sigma = np.array(sigma)*1e-4  # m^2/mol
            sigma = sigma * 1e-4  # m^2/mol
            pass
        if isothermal:
            tau = tau + mmr * sigma * np.array([rho]).T
            tau_by_molecule[elem] = mmr * sigma * np.array([rho]).T
            # print('    mmr', mmr)  # tensor
            # print('    sigma', sigma) # 100 floats
            # print('    rho', rho)  # 7 floats
            # print('tau for this molecule', elem, tau_by_molecule[elem].eval())
        pass
    # CIA ARRAY, ZPRIME VERSUS WAVELENGTH  -------------------------------------------
    for cia in cialist:
        if cia == 'H2-H2':
            f1 = fH2
            f2 = fH2
        elif cia == 'H2-He':
            f1 = fH2
            f2 = fHe
        elif cia == 'H2-H':
            f1 = fH2
            f2 = fH2 * 2.0
        elif cia == 'He-H':
            f1 = fHe
            f2 = fH2 * 2.0
        else:
            log.warning(
                '--< CERBERUS gettau(): UNEXPECTED MOLECULE %s >--', cia
            )
            f1 = 0
            f2 = 0
        # HITRAN RICHARD ET AL. 2012
        sigma, lsig = getciaxs(temp, xsecs[cia])  # cm^5/mol^2
        sigma = np.array(sigma) * 1e-10  # m^5/mol^2
        if True in (sigma < 0):
            sigma[sigma < 0] = 0e0
        if True in ~np.isfinite(sigma):
            sigma[~np.isfinite(sigma)] = 0e0
        tau = tau + f1 * f2 * sigma * np.array([rho**2]).T
        tau_by_molecule[cia] = f1 * f2 * sigma * np.array([rho**2]).T
        # print('tau for this molecule', cia, tau_by_molecule[cia].eval())
    # RAYLEIGH ARRAY, ZPRIME VERSUS WAVELENGTH  --------------------------------------
    # NAUS & UBACHS 2000
    slambda0 = 750.0 * 1e-3  # microns
    sray0 = 2.52 * 1e-28 * 1e-4  # m^2/mol
    sigma = sray0 * (wgrid[::-1] / slambda0) ** (-4)
    tau = tau + fH2 * sigma * np.array([rho]).T
    tau_by_molecule['rayleigh'] = fH2 * sigma * np.array([rho]).T
    # print('tau for this molecule', 'rayleigh', tau_by_molecule['rayleigh'].eval())
    # HAZE ARRAY, ZPRIME VERSUS WAVELENGTH  ------------------------------------------
    if hzlib is None:
        slambda0 = 750.0 * 1e-3  # microns
        sray0 = 2.52 * 1e-28 * 1e-4  # m^2/mol
        sigma = sray0 * (wgrid[::-1] / slambda0) ** (hzslope)
        hazedensity = np.ones(len(z))
        tau = tau + (10.0**rayleigh) * sigma * np.array([hazedensity]).T
        tau_by_molecule['haze'] = (
            (10.0**rayleigh) * sigma * np.array([hazedensity]).T
        )
    else:
        # WEST ET AL. 2004
        sigma = (
            0.0083
            * (wgrid[::-1]) ** (hzslope)
            * (
                1e0
                + 0.014 * (wgrid[::-1]) ** (hzslope / 2e0)
                + 0.00027 * (wgrid[::-1]) ** (hzslope)
            )
        )
        if hzp in ['MAX', 'MEDIAN', 'AVERAGE']:
            frh = hzlib['PROFILE'][0][hzp][0]
            rh = frh(p)
            rh[rh < 0] = 0.0
            refhzp = float(p[rh == np.max(rh)])
            if hztop is None:
                hzshift = 0e0
            else:
                hzshift = hztop - np.log10(refhzp)
            splp = np.log10(p[::-1])
            splrh = rh[::-1]
            thisfrh = itp(
                splp, splrh, kind='linear', bounds_error=False, fill_value=0e0
            )
            hzwdist = hztop - np.log10(p)
            if hzwscale > 0:
                preval = hztop - hzwdist / hzwscale - hzshift
                rh = thisfrh(preval)
                rh[rh < 0] = 0e0
            else:
                rh = thisfrh(np.log10(p)) * 0
            if debug:
                jptprofile = 'J' + hzp
                jdata = np.array(hzlib['PROFILE'][0][jptprofile])
                jpres = np.array(hzlib['PROFILE'][0]['PRESSURE'])
                myfig = plt.figure(figsize=(12, 6))
                plt.plot(
                    1e6 * jdata, jpres, color='blue', label='Lavvas et al. 2017'
                )
                plt.axhline(refhzp, linestyle='--', color='blue')
                plt.plot(1e6 * rh, p, 'r', label='Parametrized density profile')
                plt.plot(1e6 * thisfrh(np.log10(p) - hzshift), p, 'g^')
                if hztop is not None:
                    plt.axhline(10**hztop, linestyle='--', color='red')
                    pass
                plt.semilogy()
                plt.semilogx()
                plt.gca().invert_yaxis()
                plt.xlim([1e-4, np.max(1e6 * rh)])
                plt.tick_params(axis='both', labelsize=20)
                plt.xlabel('Aerosol Density [$n.{cm}^{-3}$]', fontsize=24)
                plt.ylabel('Pressure [bar]', fontsize=24)
                plt.title('Aerosol density profile', fontsize=24)
                plt.legend(
                    loc='center left',
                    frameon=False,
                    fontsize=24,
                    bbox_to_anchor=(1, 1),
                )
                myfig.tight_layout()
                plt.show()
                pass
            pass
        else:
            rh = np.array(
                [np.nanmean(hzlib['PROFILE'][0]['CONSTANT'])] * len(z)
            )
            negrh = rh < 0e0
            if True in negrh:
                rh[negrh] = 0e0
            pass
        # print('lower haze',rayleigh)
        # print('lower haze',sigma)
        print('lower haze',rh)
        hazecontribution = 10.0**rayleigh * sigma * np.array([rh]).T
        # convert the haze contribution to a tensor
        #  otherwise haze will be different type than all other contributions
        # hazecontribution = tensor.as_tensor(hazecontribution)
        tau = tau + hazecontribution
        tau_by_molecule['haze'] = hazecontribution
        pass
    # careful, haze is the weird one. it's float, not tensor
    # print('tau for this molecule', 'haze', tau_by_molecule['haze'])
    # print('OVERALL TAU!!!', tau.eval())

    # any trouble with this matrix multiplication?
    print('dlarray', dlarray)
    print('dlarray', dlarray.shape)  # 7x7 elements
    print('dlarray[0][0]', dlarray[0][0])
    print('tau', tau.eval().shape)  # 7x103 elements

    # module 'pytensor.tensor' has no attribute 'as_matrix'. Did you mean: 'bmatrix'?
    #    tau = 2e0 * np.array(tensor.as_matrix(dlarray) * tensor.as_matrix(tau))
    # tau = 2e0 * dlarray * tau  # fails
    # tau = 2e0 * tensor.nlinalg.matrix_dot(dlarray, tau)  # keyerror: 'object'

    print('dlarray', dlarray)
    # TypeError: Unsupported dtype for TensorType: object
    # testtest = tensor.as_tensor(dlarray)
    # print('dlarray',testtest)

    # dlarraymod = []
    # for dla in dlarray: dlarraymod.append(dla.eval())
    # dlarraymod = np.array(dlarraymod)
    # print('dlarraymod',dlarraymod)
    # print('dlarraymod shape',dlarraymod.shape)
    # testtest = tensor.as_tensor(dlarraymod)
    # print('dlarray as tensor',testtest)
    # print('OK!')
    # print()

    # tau = 2e0 * np.array(np.asmatrix(dlarray) * np.asmatrix(tau))
    # print('tau after matrixing',tau.shape)
    tau = 2e0 * tensor.nlinalg.matrix_dot(tensor.as_tensor(dlarray), tau)
    print('tau shape after matrixing', tau.eval().shape)

    molecules = tau_by_molecule.keys()
    for molecule in molecules:
        print(' MOLECULE:', molecule)
        # careful here.  most of the time this is a tensor
        #  but for haze it comes out as a numpy array
        #  (haze is based on a float-based density profile)
        # print('tau[0][0] for molecule',
        #      molecule, tau_by_molecule[molecule][0][0])
        # this is only an issue for printing.  this handles either case fine
        tau_by_molecule[molecule] = 2e0 * tensor.nlinalg.matrix_dot(
            tensor.as_tensor(dlarray), tau_by_molecule[molecule]
        )
        print(
            'tau[0][0] for molecule',
            molecule,
            tau_by_molecule[molecule][0][0].eval(),
        )
    if debug:
        plt.figure(figsize=(12, 6))
        plt.imshow(
            tensor.log10(tau).eval(),
            aspect='auto',
            origin='lower',
            extent=[max(wgrid), min(wgrid), np.log10(max(p)), np.log10(min(p))],
        )
        plt.ylabel('log10(Pressure)', fontsize=24)
        plt.xlabel('Wavelength [$\\mu m$]', fontsize=24)
        plt.gca().invert_xaxis()
        plt.title('log10(Optical Depth)', fontsize=24)
        plt.tick_params(axis='both', labelsize=20)
        cbar = plt.colorbar()
        cbar.ax.tick_params(labelsize=20)
        # asdf
        # plt.savefig('opticalDepth1.png')  # permission denied
        # plt.savefig('/proj/sdp/bryden/opticalDepth2.png')  # no such file/dir
        plt.savefig(
            '/tmp/opticalDepth3.png'
        )  # hey this works! but where is it?!
        # plt.savefig('/home/bryden/opticalDepth4.png') # no such file/dir
        plt.show()
        pass
    return tau, tau_by_molecule, 1e4 / lsig


# --------- ----------------------------------------------------------
# -- ATTENUATION COEFFICIENT -- --------------------------------------
def absorb(
    xsecs,
    qtgrid,
    T,
    p,
    mmr,
    lbroadening,
    lshifting,
    wgrid,
    iso=0,
    Tref=296.0,
    debug=False,
):
    '''
    G. ROUDIER: HITRAN HITEMP database parser
    '''
    select = np.array(xsecs['I']) == iso + 1
    S = np.array(xsecs['S'])[select]
    E = np.array(xsecs['Epp'])[select]
    gself = np.array(xsecs['g_self'])[select]
    nu = np.array(xsecs['nu'])[select]
    delta = np.array(xsecs['delta'])[select]
    eta = np.array(xsecs['eta'])[select]
    gair = np.array(xsecs['g_air'])[select]
    Qref = float(qtgrid['SPL'][iso](Tref))
    try:
        Q = float(qtgrid['SPL'][iso](T))
    except ValueError:
        Q = np.nan
    c2 = 1e2 * cst.h * cst.c / cst.Boltzmann
    tips = (Qref * np.exp(-c2 * E / T) * (1.0 - np.exp(-c2 * nu / T))) / (
        Q * np.exp(-c2 * E / Tref) * (1.0 - np.exp(-c2 * nu / Tref))
    )
    if np.all(~np.isfinite(tips)):
        tips = 0
    sigma = S * tips
    ps = mmr * p
    gamma = np.array(
        np.asmatrix(p - ps).T * np.asmatrix(gair * (Tref / T) ** eta)
        + np.asmatrix(ps).T * np.asmatrix(gself)
    )
    if lbroadening:
        if lshifting:
            matnu = np.array(
                np.asmatrix(np.ones(p.size)).T * np.asmatrix(nu)
                + np.asmatrix(p).T * np.asmatrix(delta)
            )
        else:
            matnu = np.array(nu) * np.array([np.ones(len(p))]).T
        pass
    else:
        matnu = np.array(nu)
    absgrid = []
    nugrid = (1e4 / wgrid)[::-1]
    dwnu = np.concatenate((np.array([np.diff(nugrid)[0]]), np.diff(nugrid)))
    if lbroadening:
        for mymatnu, mygamma in zip(matnu, gamma):
            binsigma = np.asmatrix(sigma) * np.asmatrix(
                intflor(
                    nugrid,
                    dwnu / 2.0,
                    np.array([mymatnu]).T,
                    np.array([mygamma]).T,
                )
            )
            binsigma = np.array(binsigma).flatten()
            absgrid.append(binsigma / dwnu)
            pass
        pass
    else:
        binsigma = []
        for nubin, dw in zip(nugrid, dwnu):
            select = (matnu > (nubin - dw / 2.0)) & (matnu <= nubin + dw / 2.0)
            binsigma.append(np.sum(sigma[select]))
            pass
        binsigma = np.array(binsigma) / dwnu
        absgrid.append(binsigma)
        pass
    if debug:
        plt.semilogy(1e4 / matnu.T, sigma, '.')
        plt.semilogy(wgrid[::-1], binsigma, 'o')
        plt.xlabel('Wavelength $\\lambda$[$\\mu m$]')
        plt.ylabel('Absorption Coeff [$cm^{2}.molecule^{-1}$]')
        plt.show()
        pass
    return absgrid, nugrid


# --------- ----------------------------------------------------------
# -- EXOMOL -- -------------------------------------------------------
def getxmolxs(temp, xsecs):
    '''
    G. ROUDIER: Wrapper around EXOMOL Cerberus library
    '''
    print('--getxmolxs start--')
    # print('xsecs',xsecs['SPL'])
    # sigma = np.array([thisspl for thisspl in xsecs['SPL']])
    # unneccessary-comprehension error here.  but itk maybe needed for tensor version?
    sigma = np.array(list(xsecs['SPL']))
    # print('sigma3', sigma[3])  # scipy.interp object
    # print('sigma3', xsecs['SPL'][3])  # scipy.interp object
    #    print('nu', xsecs['SPLNU'])      # a list of floats
    # print('nu3', xsecs['SPLNU'][3])  # a float
    # temp is not a tensor, for ariel-sim call. for cerberus?  asdf
    # print('temperature', temp)
    #    for thisspl in xsecs['SPL']:
    # print('does this single call work?', thisspl)
    #        print('does this single call work?', thisspl(temp))  # <-- was failing before
    sigma = np.array([thisspl(temp) for thisspl in xsecs['SPL']])
    nu = np.array(xsecs['SPLNU'])
    select = np.argsort(nu)
    nu = nu[select]
    sigma = sigma[select]
    print('--getxmolxs end--')
    return sigma, nu


# ------------------------- ------------------------------------------
# -- CIA -- ----------------------------------------------------------
def getciaxs(temp, xsecs):
    '''
    G. ROUDIER: Wrapper around CIA Cerberus library
    '''
    sigma = np.array([thisspl(temp) for thisspl in xsecs['SPL']])
    nu = np.array(xsecs['SPLNU'])
    select = np.argsort(nu)
    nu = nu[select]
    sigma = sigma[select]
    return sigma, nu


# ----------------------------- --------------------------------------
# -- PRESSURE BROADENING -- ------------------------------------------
def intflor(wave, dwave, nu, gamma):
    '''
    G. ROUDIER: Pressure Broadening
    '''
    f = (
        1e0
        / np.pi
        * (
            np.arctan((wave + dwave - nu) / gamma)
            - np.arctan((wave - dwave - nu) / gamma)
        )
    )
    return f


# -------------------------- -----------------------------------------
# -- PYMC DETERMINISTIC FUNCTIONS -- ---------------------------------
# @tco.as_op(itypes=[tt.dscalar, tt.dscalar, tt.dscalar, tt.dscalar, tt.dscalar,
#                   tt.dvector],
#           otypes=[tt.dvector])
def cloudyfmcerberus(*crbinputs):
    '''
    G. ROUDIER: Wrapper around Cerberus forward model, spherical shell symmetry
    '''
    ctp, hza, hzloc, hzthick, tpr, mdp = crbinputs
    print(
        ' not-fixed cloud parameters (cloudy):',
        tpr.eval(),
        ctp.eval(),
        hza.eval(),
        hzloc.eval(),
        hzthick.eval(),
    )
    fmc = np.zeros(ctxt.tspectrum.size)
    if ctxt.model == 'TEC':
        tceqdict = {}
        mdpindex = 0
        if 'XtoH' in ctxt.fixedParams:
            tceqdict['XtoH'] = ctxt.fixedParams['XtoH']
        else:
            tceqdict['XtoH'] = mdp[mdpindex].eval()
            # tceqdict['XtoH'] = mdp[mdpindex]     make sure to fix/test this next!!!
            # asdf
            mdpindex += 1

        if 'CtoO' in ctxt.fixedParams:
            tceqdict['CtoO'] = ctxt.fixedParams['CtoO']
        else:
            tceqdict['CtoO'] = mdp[mdpindex].eval()
            # tceqdict['CtoO'] = mdp[mdpindex]
            mdpindex += 1

        if 'NtoO' in ctxt.fixedParams:
            tceqdict['NtoO'] = ctxt.fixedParams['NtoO']
        else:
            tceqdict['NtoO'] = mdp[mdpindex].eval()
            # tceqdict['NtoO'] = mdp[mdpindex]
        # print('XtoH,CtoO,NtoO =',tceqdict['XtoH'],tceqdict['CtoO'],tceqdict['NtoO'])

        #        fmc = crbmodel(None, hza.eval(), ctp.eval(), ctxt.solidr, ctxt.orbp,
        #                       ctxt.xsl['data'][ctxt.p]['XSECS'],
        #                       ctxt.xsl['data'][ctxt.p]['QTGRID'],
        #                       tpr.eval(), np.array(ctxt.spc['data'][ctxt.p]['WB']),
        #                       hzlib=ctxt.hzlib,  hzp='AVERAGE', hztop=hzloc.eval(),
        #                       hzwscale=hzthick.eval(), cheq=tceqdict, pnet=ctxt.p,
        #                       verbose=False, debug=False)
        fmc = crbmodel(
            None,
            hza,
            ctp,
            ctxt.solidr,
            ctxt.orbp,
            ctxt.xsl['data'][ctxt.p]['XSECS'],
            ctxt.xsl['data'][ctxt.p]['QTGRID'],
            tpr,
            np.array(ctxt.spc['data'][ctxt.p]['WB']),
            hzlib=ctxt.hzlib,
            hzp='AVERAGE',
            hztop=hzloc,
            hzwscale=hzthick,
            cheq=tceqdict,
            pnet=ctxt.p,
            verbose=False,
            debug=False,
        )
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = mdp[index].eval()

        fmc = crbmodel(
            mixratio,
            hza.eval(),
            ctp.eval(),
            ctxt.solidr,
            ctxt.orbp,
            ctxt.xsl['data'][ctxt.p]['XSECS'],
            ctxt.xsl['data'][ctxt.p]['QTGRID'],
            tpr.eval(),
            np.array(ctxt.spc['data'][ctxt.p]['WB']),
            hzlib=ctxt.hzlib,
            hzp='AVERAGE',
            hztop=hzloc.eval(),
            hzwscale=hzthick.eval(),
            cheq=None,
            pnet=ctxt.p,
            verbose=False,
            debug=False,
        )

    fmc = fmc[ctxt.cleanup] - np.nanmean(fmc[ctxt.cleanup])
    fmc = fmc + np.nanmean(ctxt.tspectrum[ctxt.cleanup])
    return fmc


# @tco.as_op(itypes=[tt.dscalar, tt.dvector],
#           otypes=[tt.dvector])
def clearfmcerberus(*crbinputs):
    '''
    Wrapper around Cerberus forward model - NO CLOUDS!
    '''
    # ctp = 3.    # cloud deck is very deep - 1000 bars
    # hza = -10.  # small number means essentially no haze
    # hzloc = 0.
    # hzthick = 0.
    ctp = ctxt.fixedParams['CTP']
    hza = ctxt.fixedParams['HScale']
    hzloc = ctxt.fixedParams['HLoc']
    hzthick = ctxt.fixedParams['HThick']
    # print(' fixed cloud parameters (clear):',ctp,hza,hzloc,hzthick)

    tpr, mdp = crbinputs

    # if you don't want to fit Teq, fix it here. otherwise modify decorators somehow
    # tpr = 593.5
    # if 'T' in ctxt.fixedParams: exit('fixing T is tricky; need to adjust the decorators')
    # if 'T' in ctxt.fixedParams:
    #     tpr = ctxt.fixedParams['T']
    #     mdp = crbinputs
    # else:
    #    tpr, mdp = crbinputs

    fmc = np.zeros(ctxt.tspectrum.size)
    if ctxt.model == 'TEC':
        tceqdict = {}
        mdpindex = 0
        if 'XtoH' in ctxt.fixedParams:
            tceqdict['XtoH'] = ctxt.fixedParams['XtoH']
        else:
            tceqdict['XtoH'] = float(mdp[mdpindex])
            mdpindex += 1

        if 'CtoO' in ctxt.fixedParams:
            tceqdict['CtoO'] = ctxt.fixedParams['CtoO']
        else:
            tceqdict['CtoO'] = float(mdp[mdpindex])
            mdpindex += 1

        if 'NtoO' in ctxt.fixedParams:
            tceqdict['NtoO'] = ctxt.fixedParams['NtoO']
        else:
            tceqdict['NtoO'] = float(mdp[mdpindex])
        # print('XtoH,CtoO,NtoO =',tceqdict['XtoH'],tceqdict['CtoO'],tceqdict['NtoO'])

        fmc = crbmodel(
            None,
            float(hza),
            float(ctp),
            ctxt.solidr,
            ctxt.orbp,
            ctxt.xsl['data'][ctxt.p]['XSECS'],
            ctxt.xsl['data'][ctxt.p]['QTGRID'],
            float(tpr),
            np.array(ctxt.spc['data'][ctxt.p]['WB']),
            hzlib=ctxt.hzlib,
            hzp='AVERAGE',
            hztop=float(hzloc),
            hzwscale=float(hzthick),
            cheq=tceqdict,
            pnet=ctxt.p,
            verbose=False,
            debug=False,
        )
        pass
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
            pass
        fmc = crbmodel(
            mixratio,
            float(hza),
            float(ctp),
            ctxt.solidr,
            ctxt.orbp,
            ctxt.xsl['data'][ctxt.p]['XSECS'],
            ctxt.xsl['data'][ctxt.p]['QTGRID'],
            float(tpr),
            np.array(ctxt.spc['data'][ctxt.p]['WB']),
            hzlib=ctxt.hzlib,
            hzp='AVERAGE',
            hztop=float(hzloc),
            hzwscale=float(hzthick),
            cheq=None,
            pnet=ctxt.p,
            verbose=False,
            debug=False,
        )
        pass
    fmc = fmc[ctxt.cleanup] - np.nanmean(fmc[ctxt.cleanup])
    fmc = fmc + np.nanmean(ctxt.tspectrum[ctxt.cleanup])
    return fmc


# @tco.as_op(itypes=[tt.dscalar, tt.dscalar, tt.dscalar, tt.dscalar, tt.dscalar,
#                   tt.dscalar, tt.dscalar, tt.dscalar, tt.dvector],
#           otypes=[tt.dvector])
def offcerberus(*crbinputs):
    '''
    R.ESTRELA: ADD offsets between STIS filters and STIS and WFC3 filters
    '''
    ctp, hza, off0, off1, off2, hzloc, hzthick, tpr, mdp = crbinputs
    #     off0, off1, off2 = crbinputs
    #     ctp = -2.5744083
    #     hza = -1.425234
    #     hzloc = -0.406851
    #     hzthick = 5.58950953
    #     tpr = 1551.41137
    #     mdp = [-1.24882918, -4.08582557, -2.4664526]
    wbb = np.array(ctxt.spc['data'][ctxt.p]['WB'])
    flt = np.array(ctxt.spc['data'][ctxt.p]['Fltrs'])
    #  cond_wav = (wbb < 0.56) | (wbb > 1.02)
    fmc = np.zeros(ctxt.tspectrum.size)
    if ctxt.model == 'TEC':
        tceqdict = {}
        tceqdict['XtoH'] = float(mdp[0])
        tceqdict['CtoO'] = float(mdp[1])
        tceqdict['NtoO'] = float(mdp[2])
        fmc = crbmodel(
            None,
            float(hza),
            float(ctp),
            ctxt.solidr,
            ctxt.orbp,
            ctxt.xsl['data'][ctxt.p]['XSECS'],
            ctxt.xsl['data'][ctxt.p]['QTGRID'],
            float(tpr),
            wbb,
            hzlib=ctxt.hzlib,
            hzp='AVERAGE',
            hztop=float(hzloc),
            hzwscale=float(hzthick),
            cheq=tceqdict,
            pnet=ctxt.p,
            verbose=False,
            debug=False,
        )
        pass
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
            pass
        fmc = crbmodel(
            mixratio,
            float(hza),
            float(ctp),
            ctxt.solidr,
            ctxt.orbp,
            ctxt.xsl['data'][ctxt.p]['XSECS'],
            ctxt.xsl['data'][ctxt.p]['QTGRID'],
            float(tpr),
            np.array(ctxt.spc['data'][ctxt.p]['WB']),
            hzlib=ctxt.hzlib,
            hzp='AVERAGE',
            hztop=float(hzloc),
            hzwscale=float(hzthick),
            cheq=None,
            pnet=ctxt.p,
            verbose=False,
            debug=False,
        )
        pass
    cond_G430 = flt[ctxt.cleanup] == 'HST-STIS-CCD-G430L-STARE'
    cond_G141 = flt[ctxt.cleanup] == 'HST-WFC3-IR-G141-SCAN'
    tspectrum_clean = ctxt.tspectrum[ctxt.cleanup]
    fmc = fmc[ctxt.cleanup] - np.nanmean(fmc[ctxt.cleanup][cond_G141])
    fmc = fmc + np.nanmean(tspectrum_clean[cond_G141])
    #     fmc = fmc[ctxt.cleanup] - np.nanmean(fmc[ctxt.cleanup])
    #     fmc = fmc + np.nanmean(ctxt.tspectrum[ctxt.cleanup])
    ww = wbb
    ww = ww[ctxt.cleanup]
    cond_G750 = flt[ctxt.cleanup] == 'HST-STIS-CCD-G750L-STARE'
    cond_G102 = flt[ctxt.cleanup] == 'HST-WFC3-IR-G102-SCAN'
    fmc[cond_G430] = fmc[cond_G430] - 1e-2 * float(off0)
    fmc[cond_G750] = fmc[cond_G750] - 1e-2 * float(off1)
    fmc[cond_G102] = fmc[cond_G102] - 1e-2 * float(off2)
    return fmc


# @tco.as_op(itypes=[tt.dscalar, tt.dscalar, tt.dscalar, tt.dscalar, tt.dscalar,
#                   tt.dscalar, tt.dscalar, tt.dvector],
#           otypes=[tt.dvector])
def offcerberus1(*crbinputs):
    '''
    R.ESTRELA: ADD offsets between STIS filters and STIS and WFC3 filters
    '''
    ctp, hza, off0, off1, hzloc, hzthick, tpr, mdp = crbinputs
    wbb = np.array(ctxt.spc['data'][ctxt.p]['WB'])
    fmc = np.zeros(ctxt.tspectrum.size)
    if ctxt.model == 'TEC':
        tceqdict = {}
        tceqdict['XtoH'] = float(mdp[0])
        tceqdict['CtoO'] = float(mdp[1])
        tceqdict['NtoO'] = float(mdp[2])
        fmc = crbmodel(
            None,
            float(hza),
            float(ctp),
            ctxt.solidr,
            ctxt.orbp,
            ctxt.xsl['data'][ctxt.p]['XSECS'],
            ctxt.xsl['data'][ctxt.p]['QTGRID'],
            float(tpr),
            wbb,
            hzlib=ctxt.hzlib,
            hzp='AVERAGE',
            hztop=float(hzloc),
            hzwscale=float(hzthick),
            cheq=tceqdict,
            pnet=ctxt.p,
            verbose=False,
            debug=False,
        )
        pass
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
            pass
        fmc = crbmodel(
            mixratio,
            float(hza),
            float(ctp),
            ctxt.solidr,
            ctxt.orbp,
            ctxt.xsl['data'][ctxt.p]['XSECS'],
            ctxt.xsl['data'][ctxt.p]['QTGRID'],
            float(tpr),
            np.array(ctxt.spc['data'][ctxt.p]['WB']),
            hzlib=ctxt.hzlib,
            hzp='AVERAGE',
            hztop=float(hzloc),
            hzwscale=float(hzthick),
            cheq=None,
            pnet=ctxt.p,
            verbose=False,
            debug=False,
        )
        pass
    fmc = fmc[ctxt.cleanup] - np.nanmean(fmc[ctxt.cleanup])
    fmc = fmc + np.nanmean(ctxt.tspectrum[ctxt.cleanup])
    ww = wbb
    ww = ww[ctxt.cleanup]
    flt = np.array(ctxt.spc['data'][ctxt.p]['Fltrs'])
    cond_G430 = 'HST-STIS-CCD-G430L-STARE' in flt
    cond_G750 = 'HST-STIS-CCD-G750L-STARE' in flt
    fmc[cond_G430] = fmc[cond_G430] + 1e-2 * float(off0)
    fmc[cond_G750] = fmc[cond_G750] + 1e-2 * float(off1)
    return fmc


# @tco.as_op(itypes=[tt.dscalar, tt.dscalar, tt.dscalar, tt.dscalar, tt.dscalar,
#                   tt.dscalar, tt.dscalar, tt.dvector],
#           otypes=[tt.dvector])
def offcerberus2(*crbinputs):
    '''
    R.ESTRELA: ADD offsets between STIS filters and STIS and WFC3 filters
    '''
    ctp, hza, off0, off1, hzloc, hzthick, tpr, mdp = crbinputs
    wbb = np.array(ctxt.spc['data'][ctxt.p]['WB'])
    fmc = np.zeros(ctxt.tspectrum.size)
    if ctxt.model == 'TEC':
        tceqdict = {}
        tceqdict['XtoH'] = float(mdp[0])
        tceqdict['CtoO'] = float(mdp[1])
        tceqdict['NtoO'] = float(mdp[2])
        fmc = crbmodel(
            None,
            float(hza),
            float(ctp),
            ctxt.solidr,
            ctxt.orbp,
            ctxt.xsl['data'][ctxt.p]['XSECS'],
            ctxt.xsl['data'][ctxt.p]['QTGRID'],
            float(tpr),
            wbb,
            hzlib=ctxt.hzlib,
            hzp='AVERAGE',
            hztop=float(hzloc),
            hzwscale=float(hzthick),
            cheq=tceqdict,
            pnet=ctxt.p,
            verbose=False,
            debug=False,
        )
        pass
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
            pass
        fmc = crbmodel(
            mixratio,
            float(hza),
            float(ctp),
            ctxt.solidr,
            ctxt.orbp,
            ctxt.xsl['data'][ctxt.p]['XSECS'],
            ctxt.xsl['data'][ctxt.p]['QTGRID'],
            float(tpr),
            np.array(ctxt.spc['data'][ctxt.p]['WB']),
            hzlib=ctxt.hzlib,
            hzp='AVERAGE',
            hztop=float(hzloc),
            hzwscale=float(hzthick),
            cheq=None,
            pnet=ctxt.p,
            verbose=False,
            debug=False,
        )
        pass
    fmc = fmc[ctxt.cleanup] - np.nanmean(fmc[ctxt.cleanup])
    fmc = fmc + np.nanmean(ctxt.tspectrum[ctxt.cleanup])
    ww = wbb
    ww = ww[ctxt.cleanup]
    flt = np.array(ctxt.spc['data'][ctxt.p]['Fltrs'])
    cond_G430 = 'HST-STIS-CCD-G430-STARE' in flt
    cond_G750 = 'HST-STIS-CCD-G750-STARE' in flt
    fmc[cond_G430] = fmc[cond_G430] + 1e-2 * float(off0)
    fmc[cond_G750] = fmc[cond_G750] + 1e-2 * float(off1)
    return fmc


# @tco.as_op(itypes=[tt.dscalar, tt.dscalar, tt.dscalar, tt.dscalar, tt.dscalar,
#                   tt.dscalar, tt.dscalar, tt.dvector],
#           otypes=[tt.dvector])
def offcerberus3(*crbinputs):
    '''
    R.ESTRELA: ADD offsets between STIS filters and STIS and WFC3 filters
    '''
    ctp, hza, off0, off1, hzloc, hzthick, tpr, mdp = crbinputs
    wbb = np.array(ctxt.spc['data'][ctxt.p]['WB'])
    fmc = np.zeros(ctxt.tspectrum.size)
    flt = np.array(ctxt.spc['data'][ctxt.p]['Fltrs'])
    if ctxt.model == 'TEC':
        tceqdict = {}
        tceqdict['XtoH'] = float(mdp[0])
        tceqdict['CtoO'] = float(mdp[1])
        tceqdict['NtoO'] = float(mdp[2])
        fmc = crbmodel(
            None,
            float(hza),
            float(ctp),
            ctxt.solidr,
            ctxt.orbp,
            ctxt.xsl['data'][ctxt.p]['XSECS'],
            ctxt.xsl['data'][ctxt.p]['QTGRID'],
            float(tpr),
            wbb,
            hzlib=ctxt.hzlib,
            hzp='AVERAGE',
            hztop=float(hzloc),
            hzwscale=float(hzthick),
            cheq=tceqdict,
            pnet=ctxt.p,
            verbose=False,
            debug=False,
        )
        pass
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
            pass
        fmc = crbmodel(
            mixratio,
            float(hza),
            float(ctp),
            ctxt.solidr,
            ctxt.orbp,
            ctxt.xsl['data'][ctxt.p]['XSECS'],
            ctxt.xsl['data'][ctxt.p]['QTGRID'],
            float(tpr),
            np.array(ctxt.spc['data'][ctxt.p]['WB']),
            hzlib=ctxt.hzlib,
            hzp='AVERAGE',
            hztop=float(hzloc),
            hzwscale=float(hzthick),
            cheq=None,
            pnet=ctxt.p,
            verbose=False,
            debug=False,
        )
        pass
    fmc = fmc[ctxt.cleanup] - np.nanmean(fmc[ctxt.cleanup])
    fmc = fmc + np.nanmean(ctxt.tspectrum[ctxt.cleanup])
    ww = wbb
    ww = ww[ctxt.cleanup]
    cond_G430 = 'HST-STIS-CCD-G430-STARE' in flt
    cond_G102 = 'HST-WFC3-IR-G102-SCAN' in flt
    fmc[cond_G430] = fmc[cond_G430] + 1e-2 * float(off0)
    fmc[cond_G102] = fmc[cond_G102] + 1e-2 * float(off1)
    return fmc


# @tco.as_op(itypes=[tt.dscalar, tt.dscalar, tt.dscalar, tt.dscalar, tt.dscalar, tt.dscalar,
#                   tt.dvector],
#           otypes=[tt.dvector])
def offcerberus4(*crbinputs):
    '''
    R.ESTRELA: ADD offsets between STIS filters and STIS and WFC3 filters
    '''
    ctp, hza, off0, hzloc, hzthick, tpr, mdp = crbinputs
    wbb = np.array(ctxt.spc['data'][ctxt.p]['WB'])
    fmc = np.zeros(ctxt.tspectrum.size)
    flt = np.array(ctxt.spc['data'][ctxt.p]['Fltrs'])
    if ctxt.model == 'TEC':
        tceqdict = {}
        tceqdict['XtoH'] = float(mdp[0])
        tceqdict['CtoO'] = float(mdp[1])
        tceqdict['NtoO'] = float(mdp[2])
        fmc = crbmodel(
            None,
            float(hza),
            float(ctp),
            ctxt.solidr,
            ctxt.orbp,
            ctxt.xsl['data'][ctxt.p]['XSECS'],
            ctxt.xsl['data'][ctxt.p]['QTGRID'],
            float(tpr),
            wbb,
            hzlib=ctxt.hzlib,
            hzp='AVERAGE',
            hztop=float(hzloc),
            hzwscale=float(hzthick),
            cheq=tceqdict,
            pnet=ctxt.p,
            verbose=False,
            debug=False,
        )
        pass
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
            pass
        fmc = crbmodel(
            mixratio,
            float(hza),
            float(ctp),
            ctxt.solidr,
            ctxt.orbp,
            ctxt.xsl['data'][ctxt.p]['XSECS'],
            ctxt.xsl['data'][ctxt.p]['QTGRID'],
            float(tpr),
            np.array(ctxt.spc['data'][ctxt.p]['WB']),
            hzlib=ctxt.hzlib,
            hzp='AVERAGE',
            hztop=float(hzloc),
            hzwscale=float(hzthick),
            cheq=None,
            pnet=ctxt.p,
            verbose=False,
            debug=False,
        )
        pass
    fmc = fmc[ctxt.cleanup] - np.nanmean(fmc[ctxt.cleanup])
    fmc = fmc + np.nanmean(ctxt.tspectrum[ctxt.cleanup])
    ww = wbb
    ww = ww[ctxt.cleanup]
    cond_G430 = 'HST-STIS-CCD-G430-STARE' in flt
    fmc[cond_G430] = fmc[cond_G430] + 1e-2 * float(off0)
    return fmc


# @tco.as_op(itypes=[tt.dscalar, tt.dscalar, tt.dscalar, tt.dscalar, tt.dscalar,
#                   tt.dscalar, tt.dscalar, tt.dvector],
#           otypes=[tt.dvector])
def offcerberus5(*crbinputs):
    '''
    R.ESTRELA: ADD offsets between STIS filters and STIS and WFC3 filters
    '''
    ctp, hza, off0, off1, hzloc, hzthick, tpr, mdp = crbinputs
    wbb = np.array(ctxt.spc['data'][ctxt.p]['WB'])
    fmc = np.zeros(ctxt.tspectrum.size)
    flt = np.array(ctxt.spc['data'][ctxt.p]['Fltrs'])
    if ctxt.model == 'TEC':
        tceqdict = {}
        tceqdict['XtoH'] = float(mdp[0])
        tceqdict['CtoO'] = float(mdp[1])
        tceqdict['NtoO'] = float(mdp[2])
        fmc = crbmodel(
            None,
            float(hza),
            float(ctp),
            ctxt.solidr,
            ctxt.orbp,
            ctxt.xsl['data'][ctxt.p]['XSECS'],
            ctxt.xsl['data'][ctxt.p]['QTGRID'],
            float(tpr),
            wbb,
            hzlib=ctxt.hzlib,
            hzp='AVERAGE',
            hztop=float(hzloc),
            hzwscale=float(hzthick),
            cheq=tceqdict,
            pnet=ctxt.p,
            verbose=False,
            debug=False,
        )
        pass
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
            pass
        fmc = crbmodel(
            mixratio,
            float(hza),
            float(ctp),
            ctxt.solidr,
            ctxt.orbp,
            ctxt.xsl['data'][ctxt.p]['XSECS'],
            ctxt.xsl['data'][ctxt.p]['QTGRID'],
            float(tpr),
            np.array(ctxt.spc['data'][ctxt.p]['WB']),
            hzlib=ctxt.hzlib,
            hzp='AVERAGE',
            hztop=float(hzloc),
            hzwscale=float(hzthick),
            cheq=None,
            pnet=ctxt.p,
            verbose=False,
            debug=False,
        )
        pass
    fmc = fmc[ctxt.cleanup] - np.nanmean(fmc[ctxt.cleanup])
    fmc = fmc + np.nanmean(ctxt.tspectrum[ctxt.cleanup])
    ww = wbb
    ww = ww[ctxt.cleanup]
    cond_G102 = 'HST-WFC3-IR-G102-SCAN' in flt
    cond_G750 = 'HST-STIS-CCD-G750-STARE' in flt
    fmc[cond_G750] = fmc[cond_G750] + 1e-2 * float(off0)
    fmc[cond_G102] = fmc[cond_G102] + 1e-2 * float(off1)
    return fmc


# @tco.as_op(itypes=[tt.dscalar, tt.dscalar, tt.dscalar, tt.dscalar, tt.dscalar, tt.dscalar,
#                   tt.dvector],
#           otypes=[tt.dvector])
def offcerberus6(*crbinputs):
    '''
    R.ESTRELA: ADD offsets between STIS filters and STIS and WFC3 filters
    '''
    ctp, hza, off0, hzloc, hzthick, tpr, mdp = crbinputs
    wbb = np.array(ctxt.spc['data'][ctxt.p]['WB'])
    fmc = np.zeros(ctxt.tspectrum.size)
    flt = np.array(ctxt.spc['data'][ctxt.p]['Fltrs'])
    if ctxt.model == 'TEC':
        tceqdict = {}
        tceqdict['XtoH'] = float(mdp[0])
        tceqdict['CtoO'] = float(mdp[1])
        tceqdict['NtoO'] = float(mdp[2])
        fmc = crbmodel(
            None,
            float(hza),
            float(ctp),
            ctxt.solidr,
            ctxt.orbp,
            ctxt.xsl['data'][ctxt.p]['XSECS'],
            ctxt.xsl['data'][ctxt.p]['QTGRID'],
            float(tpr),
            wbb,
            hzlib=ctxt.hzlib,
            hzp='AVERAGE',
            hztop=float(hzloc),
            hzwscale=float(hzthick),
            cheq=tceqdict,
            pnet=ctxt.p,
            verbose=False,
            debug=False,
        )
        pass
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
            pass
        fmc = crbmodel(
            mixratio,
            float(hza),
            float(ctp),
            ctxt.solidr,
            ctxt.orbp,
            ctxt.xsl['data'][ctxt.p]['XSECS'],
            ctxt.xsl['data'][ctxt.p]['QTGRID'],
            float(tpr),
            np.array(ctxt.spc['data'][ctxt.p]['WB']),
            hzlib=ctxt.hzlib,
            hzp='AVERAGE',
            hztop=float(hzloc),
            hzwscale=float(hzthick),
            cheq=None,
            pnet=ctxt.p,
            verbose=False,
            debug=False,
        )
        pass
    fmc = fmc[ctxt.cleanup] - np.nanmean(fmc[ctxt.cleanup])
    fmc = fmc + np.nanmean(ctxt.tspectrum[ctxt.cleanup])
    ww = wbb
    ww = ww[ctxt.cleanup]
    cond_G750 = 'HST-STIS-CCD-G750-STARE' in flt
    fmc[cond_G750] = fmc[cond_G750] + 1e-2 * float(off0)
    return fmc


# @tco.as_op(itypes=[tt.dscalar, tt.dscalar, tt.dscalar, tt.dscalar, tt.dscalar, tt.dscalar,
#                   tt.dvector],
#           otypes=[tt.dvector])
def offcerberus7(*crbinputs):
    '''
    R.ESTRELA: ADD offsets between STIS filters and WFC3 filters
    '''
    ctp, hza, off0, hzloc, hzthick, tpr, mdp = crbinputs
    wbb = np.array(ctxt.spc['data'][ctxt.p]['WB'])
    fmc = np.zeros(ctxt.tspectrum.size)
    flt = np.array(ctxt.spc['data'][ctxt.p]['Fltrs'])
    if ctxt.model == 'TEC':
        tceqdict = {}
        tceqdict['XtoH'] = float(mdp[0])
        tceqdict['CtoO'] = float(mdp[1])
        tceqdict['NtoO'] = float(mdp[2])
        fmc = crbmodel(
            None,
            float(hza),
            float(ctp),
            ctxt.solidr,
            ctxt.orbp,
            ctxt.xsl['data'][ctxt.p]['XSECS'],
            ctxt.xsl['data'][ctxt.p]['QTGRID'],
            float(tpr),
            wbb,
            hzlib=ctxt.hzlib,
            hzp='AVERAGE',
            hztop=float(hzloc),
            hzwscale=float(hzthick),
            cheq=tceqdict,
            pnet=ctxt.p,
            verbose=False,
            debug=False,
        )
        pass
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
            pass
        fmc = crbmodel(
            mixratio,
            float(hza),
            float(ctp),
            ctxt.solidr,
            ctxt.orbp,
            ctxt.xsl['data'][ctxt.p]['XSECS'],
            ctxt.xsl['data'][ctxt.p]['QTGRID'],
            float(tpr),
            np.array(ctxt.spc['data'][ctxt.p]['WB']),
            hzlib=ctxt.hzlib,
            hzp='AVERAGE',
            hztop=float(hzloc),
            hzwscale=float(hzthick),
            cheq=None,
            pnet=ctxt.p,
            verbose=False,
            debug=False,
        )
        pass
    fmc = fmc[ctxt.cleanup] - np.nanmean(fmc[ctxt.cleanup])
    fmc = fmc + np.nanmean(ctxt.tspectrum[ctxt.cleanup])
    ww = wbb
    ww = ww[ctxt.cleanup]
    cond_G750 = 'HST-STIS-CCD-G750-STARE' in flt
    fmc[cond_G750] = fmc[cond_G750] + 1e-2 * float(off0)
    return fmc


# @tco.as_op(itypes=[tt.dscalar, tt.dscalar, tt.dscalar, tt.dscalar, tt.dscalar,
#                   tt.dscalar, tt.dscalar, tt.dvector],
#           otypes=[tt.dvector])
def offcerberus8(*crbinputs):
    '''
    R.ESTRELA: ADD offsets between WFC3 filters
    '''
    ctp, hza, off0, hzloc, hzthick, tpr, mdp = crbinputs
    wbb = np.array(ctxt.spc['data'][ctxt.p]['WB'])
    fmc = np.zeros(ctxt.tspectrum.size)
    flt = np.array(ctxt.spc['data'][ctxt.p]['Fltrs'])
    if ctxt.model == 'TEC':
        tceqdict = {}
        tceqdict['XtoH'] = float(mdp[0])
        tceqdict['CtoO'] = float(mdp[1])
        tceqdict['NtoO'] = float(mdp[2])
        fmc = crbmodel(
            None,
            float(hza),
            float(ctp),
            ctxt.solidr,
            ctxt.orbp,
            ctxt.xsl['data'][ctxt.p]['XSECS'],
            ctxt.xsl['data'][ctxt.p]['QTGRID'],
            float(tpr),
            wbb,
            hzlib=ctxt.hzlib,
            hzp='AVERAGE',
            hztop=float(hzloc),
            hzwscale=float(hzthick),
            cheq=tceqdict,
            pnet=ctxt.p,
            verbose=False,
            debug=False,
        )
        pass
    else:
        mixratio = {}
        for index, key in enumerate(ctxt.modparlbl[ctxt.model]):
            mixratio[key] = float(mdp[index])
            pass
        fmc = crbmodel(
            mixratio,
            float(hza),
            float(ctp),
            ctxt.solidr,
            ctxt.orbp,
            ctxt.xsl['data'][ctxt.p]['XSECS'],
            ctxt.xsl['data'][ctxt.p]['QTGRID'],
            float(tpr),
            np.array(ctxt.spc['data'][ctxt.p]['WB']),
            hzlib=ctxt.hzlib,
            hzp='AVERAGE',
            hztop=float(hzloc),
            hzwscale=float(hzthick),
            cheq=None,
            pnet=ctxt.p,
            verbose=False,
            debug=False,
        )
        pass
    fmc = fmc[ctxt.cleanup] - np.nanmean(fmc[ctxt.cleanup])
    fmc = fmc + np.nanmean(ctxt.tspectrum[ctxt.cleanup])
    ww = wbb
    ww = ww[ctxt.cleanup]
    cond_G102 = 'HST-WFC3-IR-G102-SCAN' in flt
    fmc[cond_G102] = fmc[cond_G102] + 1e-2 * float(off0)
    return fmc
