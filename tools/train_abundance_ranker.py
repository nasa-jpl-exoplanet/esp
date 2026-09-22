import json
import pickle
import sys
from pathlib import Path

import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
import joblib

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import excalibur.system.core as syscore
import excalibur.cerberus.core as crbcore
import excalibur.cerberus.forward_model as crbfm

CORNICHONPATH = Path('/Users/siddharthmanikandan/Documents/Excalibur-Pickles/')
TARGET = 'TOI-260'
PLANET = 'b'
FLT = 'JWST-NIRSPEC-NRS-F290LP-G395H'
OUT_DIR = ROOT / 'datasets'
OUT_DIR.mkdir(exist_ok=True)

SPECIES = ['H2O', 'CO2', 'CH4', 'CO', 'NH3', 'HCN', 'C2H2']
LEVELS = np.array([-10.0, -8.0, -6.0, -4.0, -2.0, 0.0, 2.0], dtype=float)
DEFAULT_NUM_SAMPLES = 1000


def load_setup():
    name = ['system', 'Finalize', 'parameters']
    rid = 1519
    fullname = [TARGET, str(rid)]
    fullname.extend(name)
    strname = '.'.join(fullname)
    sysfin_path = CORNICHONPATH / f'{strname}.pkl'
    if not sysfin_path.exists():
        raise FileNotFoundError(f'Missing system pickle: {sysfin_path}')
    with open(sysfin_path, 'rb') as f:
        sysfin = pickle.load(f)

    name = ['transit', 'Spectrum', FLT]
    rid = 1565
    fullname = [TARGET, str(rid)]
    fullname.extend(name)
    strname = '.'.join(fullname)
    trnspc_path = CORNICHONPATH / f'{strname}.pkl'
    if not trnspc_path.exists():
        raise FileNotFoundError(f'Missing transit pickle: {trnspc_path}')
    with open(trnspc_path, 'rb') as f:
        trnspc = pickle.load(f)

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
    if not xsl_path.exists():
        raise FileNotFoundError(f'Missing xsec pickle: {xsl_path}')
    with open(xsl_path, 'rb') as f:
        crbxsl = pickle.load(f)

    return sysfin, spc, rtp, crbxsl


def make_mixratio(values):
    noabs = -10.0
    mix = {mol: noabs for mol in SPECIES}
    for mol, value in zip(SPECIES, values):
        mix[mol] = float(value)
    return mix


def make_spectrum(sysfin, spc, rtp, crbxsl, mixratio_values):
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

    result = crbfm.crbFM().crbmodel(
        Tav,
        ctp,
        cheq=tceqdict,
        mixratio=make_mixratio(mixratio_values),
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
    spectrum = np.asarray(result['spectrum'], dtype=float)
    return spectrum


def generate_dataset(num_samples=DEFAULT_NUM_SAMPLES, seed=42):
    rng = np.random.default_rng(seed)
    sysfin, spc, rtp, crbxsl = load_setup()
    wavelengths = spc['data'][PLANET]['WB']

    spectra = []
    labels = []

    for _ in range(num_samples):
        values = rng.uniform(-9.0, 2.5, size=len(SPECIES))
        strong_species = rng.choice(len(SPECIES), size=2, replace=False)
        values[strong_species] += rng.uniform(1.0, 3.0, size=2)
        values = np.clip(values, -10.0, 3.0)

        spectrum = make_spectrum(sysfin, spc, rtp, crbxsl, values)
        noise = 1.0 + rng.normal(0.0, 0.02, size=spectrum.shape)
        spectrum = spectrum * noise
        spectra.append(spectrum)
        labels.append(values.copy())

    spectra = np.stack(spectra, axis=0)
    labels = np.stack(labels, axis=0)

    dataset_path = OUT_DIR / 'abundance_ranker_dataset.npz'
    np.savez_compressed(
        dataset_path,
        wavelengths=wavelengths,
        spectra=spectra,
        labels=labels,
        species=np.array(SPECIES),
    )

    info = {
        'num_samples': int(num_samples),
        'num_features': int(spectra.shape[1]),
        'num_species': len(SPECIES),
        'species': SPECIES,
        'dataset_path': str(dataset_path),
    }
    with open(OUT_DIR / 'abundance_ranker_dataset.json', 'w', encoding='utf-8') as f:
        json.dump(info, f, indent=2)

    return dataset_path, spectra, labels, wavelengths


def train_model(dataset_path):
    data = np.load(dataset_path)
    X = data['spectra']
    y = data['labels']

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )

    model = RandomForestRegressor(
        n_estimators=500,
        random_state=42,
        n_jobs=-1,
        min_samples_leaf=1,
        max_features='sqrt',
    )
    model.fit(X_train, y_train)
    pred = model.predict(X_test)

    metrics = {
        'species': SPECIES,
        'r2_by_species': {
            species: float(r2_score(y_test[:, i], pred[:, i]))
            for i, species in enumerate(SPECIES)
        },
        'rmse_by_species': {
            species: float(np.sqrt(mean_squared_error(y_test[:, i], pred[:, i])))
            for i, species in enumerate(SPECIES)
        },
    }

    with open(OUT_DIR / 'abundance_ranker_metrics.json', 'w', encoding='utf-8') as f:
        json.dump(metrics, f, indent=2)

    model_path = OUT_DIR / 'abundance_ranker_model.joblib'
    joblib.dump(model, model_path)

    return model, pred, metrics


def topk_from_prediction(predicted_abundances, species, k=5):
    order = np.argsort(predicted_abundances)[::-1]
    return [species[i] for i in order[:k]], order[:k]


def evaluate_topk(model, dataset_path, k=5):
    data = np.load(dataset_path)
    X = data['spectra']
    y = data['labels']
    pred = model.predict(X)

    topk_true = []
    topk_pred = []
    for true_row, pred_row in zip(y, pred):
        true_order = np.argsort(true_row)[::-1]
        pred_order = np.argsort(pred_row)[::-1]
        topk_true.append([SPECIES[i] for i in true_order[:k]])
        topk_pred.append([SPECIES[i] for i in pred_order[:k]])

    overlap = []
    for truth, forecast in zip(topk_true, topk_pred):
        overlap.append(len(set(truth) & set(forecast)) / k)

    results = {
        'k': k,
        'mean_topk_overlap': float(np.mean(overlap)),
        'example_true_top5': topk_true[:3],
        'example_pred_top5': topk_pred[:3],
    }
    with open(OUT_DIR / 'abundance_ranker_topk.json', 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2)

    return results


def main():
    dataset_path, spectra, labels, wavelengths = generate_dataset(num_samples=DEFAULT_NUM_SAMPLES, seed=21)
    print(f'Created dataset: {dataset_path}')
    print(f'Array shapes: spectra={spectra.shape}, labels={labels.shape}, wavelengths={wavelengths.shape}')

    model, pred, metrics = train_model(dataset_path)
    print('Model fit complete; per-species R2 values:')
    for species, score in metrics['r2_by_species'].items():
        print(f'  {species}: R2={score:.3f}')

    ranking_metrics = evaluate_topk(model, dataset_path, k=5)
    print('Top-5 overlap summary:')
    print(ranking_metrics)

    print('Prepared outputs:')
    print('  - datasets/abundance_ranker_dataset.npz')
    print('  - datasets/abundance_ranker_model.joblib')
    print('  - datasets/abundance_ranker_metrics.json')
    print('  - datasets/abundance_ranker_topk.json')


if __name__ == '__main__':
    main()
