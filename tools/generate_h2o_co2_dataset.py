import os
import pickle
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))

import excalibur.system.core as syscore
import excalibur.cerberus.core as crbcore
import excalibur.cerberus.states as crbstt
import excalibur.cerberus.forward_model as crbfm


CORNICHONPATH = Path('/Users/siddharthmanikandan/Documents/Excalibur-Pickles/')
TARGET = 'TOI-260'
PLANET = 'b'
FLT = 'JWST-NIRSPEC-NRS-F290LP-G395H'
OUT_DIR = ROOT / 'datasets'
OUT_DIR.mkdir(exist_ok=True)


def load_setup():
    name = ['system', 'Finalize', 'parameters']
    rid = 1519
    fullname = [TARGET, str(rid)]
    fullname.extend(name)
    strname = '.'.join(fullname)
    sysfin_path = CORNICHONPATH / f'{strname}.pkl'
    if sysfin_path.exists():
        with open(sysfin_path, 'rb') as f:
            sysfin = pickle.load(f)
    else:
        raise FileNotFoundError(f'Missing system pickle: {sysfin_path}')

    name = ['transit', 'Spectrum', FLT]
    rid = 1565
    fullname = [TARGET, str(rid)]
    fullname.extend(name)
    strname = '.'.join(fullname)
    trnspc_path = CORNICHONPATH / f'{strname}.pkl'
    if trnspc_path.exists():
        with open(trnspc_path, 'rb') as f:
            trnspc = pickle.load(f)
    else:
        raise FileNotFoundError(f'Missing transit pickle: {trnspc_path}')

    spc = {'data': {PLANET: {'WB': np.concatenate([trnspc['data'][PLANET][det]['1']['WB'] for det in ['NRS1', 'NRS2']])}}}

    rtp = crbcore.CerbXSlibParams(
        knownspecies=['NO', 'OH', 'C2H2', 'N2', 'N2O', 'O3', 'O2'],
        cialist=['H2-H', 'H2-H2', 'H2-He', 'He-H'],
        xmollist=['TIO', 'H2O', 'H2CO', 'HCN', 'CO', 'CO2', 'NH3', 'CH4', 'C2H2', 'C2H6', 'C3H8', 'CH3CHO', 'SO2', 'H2S'],
        nlevels=100,
        solrad=10,
        Hsmax=25,
        lbroadening=False,
        lshifting=False,
    )

    name = ['cerberus', 'XSLib', FLT]
    rid = 6666
    fullname = [TARGET, str(rid)]
    fullname.extend(name)
    strname = '.'.join(fullname)
    xsl_path = CORNICHONPATH / f'{strname}.pkl'
    if xsl_path.exists():
        with open(xsl_path, 'rb') as f:
            crbxsl = pickle.load(f)
    else:
        raise FileNotFoundError(f'Missing xsec pickle: {xsl_path}')

    return sysfin, spc, rtp, crbxsl


def make_mixratio(value_h2o=-10.0, value_co2=-10.0):
    species = ['H2O', 'CO2', 'CH4', 'CO', 'HCN', 'C2H2', 'NH3', 'SO2', 'H2S']
    noabs = -10.0
    mix = {mol: noabs for mol in species}
    mix['H2O'] = value_h2o
    mix['CO2'] = value_co2
    return mix


def make_spectrum(sysfin, spc, rtp, crbxsl, h2o_value, co2_value):
    ssc = syscore.ssconstants(mks=True)
    Tfact = np.sqrt(sysfin['priors']['R*'] * ssc['Rsun/AU'] / (2.0 * sysfin['priors'][PLANET]['sma']))
    Tav = float(sysfin['priors']['T*'] * Tfact)
    Tav = np.asarray([Tav] * 100, dtype=float)

    pgrid = np.arange(np.log(1) - 20, np.log(1) + 20 / 100, 20 / (100 - 1))
    pressure = np.exp(pgrid[::-1])
    Tav[pressure < 10 ** -3] = 1000

    crbhzlib = {'PROFILE': []}
    crbcore.hazelib(crbhzlib, hazedir=str(CORNICHONPATH), verbose=False)

    haze = -10
    hzwscale = 3.1
    hztop = -3.43
    hzp = 'AVERAGE'
    ctp = 1
    solidr = sysfin['priors'][PLANET]['rp'] * ssc['Rjup']
    tceqdict = None

    mixratio = make_mixratio(value_h2o=h2o_value, value_co2=co2_value)

    result = crbfm.crbFM().crbmodel(
        Tav,
        ctp,
        cheq=tceqdict,
        mixratio=mixratio,
        hazescale=haze,
        hazethick=hzwscale,
        hazeloc=hztop,
        hazeprof=hzp,
        hzlib=crbhzlib,
        planet=PLANET,
        rp0=solidr,
        orbp=sysfin['priors'],
        wgrid=spc['data'][PLANET]['WB'],
        xsecs=crbxsl['data'][PLANET]['XSECS'],
        qtgrid=crbxsl['data'][PLANET]['QTGRID'],
        cialist=rtp.cialist,
        xmollist=rtp.xmollist,
        nlevels=100,
        Hsmax=25,
        solrad=1,
        break_down_by_molecule=True,
        logx=False,
        verbose=False,
        debug=False,
    )
    return np.asarray(result['spectrum'], dtype=float), spc['data'][PLANET]['WB']


def generate_dataset():
    sysfin, spc, rtp, crbxsl = load_setup()
    abundance_levels = np.array([0.0, 1.0, 2.0, 3.0, 4.0], dtype=float)

    wavelengths = spc['data'][PLANET]['WB']
    spectra = []
    labels = []

    for h2o_value in abundance_levels:
        for co2_value in abundance_levels:
            spectrum, _ = make_spectrum(sysfin, spc, rtp, crbxsl, h2o_value, co2_value)
            spectra.append(spectrum)
            labels.append(np.array([h2o_value, co2_value], dtype=float))

    spectra = np.stack(spectra, axis=0)
    labels = np.stack(labels, axis=0)

    out_path = OUT_DIR / 'h2o_co2_synthetic_spectra.npz'
    np.savez_compressed(
        out_path,
        wavelengths=wavelengths,
        spectra=spectra,
        labels=labels,
        h2o_levels=abundance_levels,
        co2_levels=abundance_levels,
    )

    print(f'Created dataset: {out_path}')
    print(f'Array shapes: spectra={spectra.shape}, labels={labels.shape}, wavelengths={wavelengths.shape}')
    print('First label example:', labels[0])
    print('First spectrum sample:', spectra[0][:5])


if __name__ == '__main__':
    generate_dataset()
