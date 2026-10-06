'''
GMR: Tools for excalibur dev on mentor
'''
import os
import sys
import pickle
import dawgie
import logging
import importlib

log = logging.getLogger(__name__)

class DuckDS:
    '''
    GMR:If it quacks like a duck it is a dawgie ds
    '''
    def __init__(self, alg:dawgie.Algorithm):
        self._alg = alg
        pass
    def update(self):
        for sv in self._alg.state_vectors():
            for k, v in sv.items():
                if not isinstance(v, dawgie.Value):
                    log.critical(f'--< {self._alg.name()}["{sv.name()}"]["{k}"] does not extend dawgie.Value >--')
                    pass
                pass
            pass
        pass
    pass

def envvar(filename:str):
    '''
    GMR:Based on Al's sv_loading.ipynb
    '''
    with open(filename, 'rt') as file:
        for line in file.readlines():
            key,value = line.replace ('export ','').replace ('\\\n','').strip().split('=')
            os.environ[key] = os.path.expandvars(os.path.expanduser(value))
            pass
        pass
    return

def thisenv(repository_root, myenv, mainpipeline=True):
    '''
    GMR:Based on Al's sv_loading.ipynb
    '''
    envvar(os.path.join(repository_root, '.docker/.env'))
    envvar(os.path.join(repository_root, myenv))
    if mainpipeline: envvar('/proj/sdp/ops/db-read-access') 
    os.environ["LDTK_ROOT"] = '/proj/sdp/data/ldtk'
    # GMR: This code now sits inside excalibur/utils/
    # It means that the following line should have been done from outside already
    # sys.path.append(repository_root)
    # GMR: on the fly imports
    # ConnectionRefusedError: [Errno 111] Connection refused
    import dawgie.db
    import dawgie.context
    import dawgie.security
    dawgie.security.initialize(
        path=os.path.expandvars(os.path.expanduser(dawgie.context.guest_public_keys)),
        myname=dawgie.context.ssl_pem_myname,
        myself=os.path.expandvars(os.path.expanduser(dawgie.context.ssl_pem_myself)),
        system=dawgie.context.ssl_pem_file
    )
    _ = dawgie.db.reopen()
    return

def loadSV(nms, trg, rid, xcd, SVroot=False):
    '''
    GMR:Returns a database product (SV)
    [I]:nms:[LIST]:SV name (['transit', 'Spectrum', 'JWST-NIRSPEC-NRS-F290LP-G395H'])
    [I]:trg:[STR]:Host star name ('HAT-P-26')
    [I]:rid:[INT]:Excalibur database runID number (666)
    [I]:xcd:[STR]:Path to excalibur code ('${HOME}/esp/excalibur')
    [OPT]:SVroot:[BOOL]:Returns SV root instead of SV[filter]
    '''
    strtask, stralgo, strsv = nms
    # GMR: on the fly imports
    # ConnectionRefusedError: [Errno 111] Connection refused
    import dawgie.pl
    import dawgie.pl.scan
    dawgie.context.ae_base_path = os.path.expandvars(xcd)
    dawgie.context.ae_base_package = 'excalibur'
    dawgie.pl.scan.for_factories(dawgie.context.ae_base_path,
                                 dawgie.context.ae_base_package)
    taskmod = importlib.import_module('excalibur.' + strtask)        
    algomod = importlib.import_module('excalibur.' + strtask + '.algorithms')
    algo = getattr(algomod, stralgo)
    alg = algo()
    for s in alg.state_vectors():
        dawgie.db.connect(alg,
                          taskmod.task(strtask, 0, rid, trg),
                          trg).load(dawgie.SV_REF(taskmod.task, alg, s))
        pass
    if SVroot: out = alg
    else: out = alg.sv_as_dict()[strsv]
    return out

def pickleSV(trg, nms, rid, odr, esp, clncrn=False, saveme=None):
    '''
    GMR:Returns either a read pickle or loadSV output
    [I]:nms:[LIST]:SV name (['transit', 'Spectrum', 'JWST-NIRSPEC-NRS-F290LP-G395H'])
    [I]:odr:[STR]:Output directory for saving pickles
    [OPT]:clncrn:[BOOL]:Cleans up pickles before reading/saving them
    '''
    fullname = [trg, str(rid)]
    fullname.extend(nms)
    strname = '.'.join(fullname)
    mycornichon = odr + strname + '.pkl'
    if clncrn:
        if os.path.isfile(mycornichon):
            os.remove(mycornichon)
            pass
        pass
    if saveme is not None:
        with open(mycornichon, 'wb') as f:
            pickle.dump(saveme, f, pickle.HIGHEST_PROTOCOL)
            pass
        return saveme
    else:
        if os.path.isfile(mycornichon):
            logging.info('>-- FROM CORNICHON %s', strname)
            with open(mycornichon, 'rb') as f:
                out = pickle.load(f)
                pass
            pass
        else:
            logging.info('>-- FROM DATABASE %s', strname)
            out = loadSV(nms, trg, rid, os.path.join(esp, 'excalibur'), SVroot=True)
            with open(mycornichon, 'wb') as f:
                logging.info('>-- PICKLING %s', strname)
                pickle.dump(out, f, pickle.HIGHEST_PROTOCOL)
                pass
            pass
        return out
    pass

def nospam():
    '''
    GMR:For sanity
    '''
    logging.basicConfig(level=logging.INFO)
    logging.getLogger("arviz").setLevel(logging.WARNING)
    logging.getLogger("dawgie").setLevel(logging.CRITICAL)
    logging.getLogger("numexpr").setLevel(logging.WARNING)
    logging.getLogger("pytensor").setLevel(logging.CRITICAL)
    return
