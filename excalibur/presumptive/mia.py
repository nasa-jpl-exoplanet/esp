'''checks that the pipeline worker capacity'''

import os
import pickle
import requests

from . import email
from . import perform

from urllib.parse import urljoin
from urllib.parse import urlparse


def farm(args):
    '''check that all farms have presented the correct workers

    Go to each node and determine what the ID of the current worker. Compare
    the remote IDs to the ops worker (the on mentor0 where it is built). If
    they are not the versions, then perform retrenchment on that node.
    '''
    wid = perform.worker_id()
    fallow = list(
        filter(lambda node, wid=wid: wid != perform.worker_id(node), args.nodes)
    )
    if fallow:
        email.send(
            args,
            f'''
Found {len(fallow)} fallow farms on mentor(s) {' '.join(f'mentor{n}' for n in fallow)}. Going to retrench the workers.
            ''',
        )
        for node in fallow:
            perform.retrenchment(node)
        return 1
    return 0


def worker(args):
    '''check that individual workers have not gone AWOL

    Check the number of workers is the correct total. It should be straight
    forward math given the number of nodes and replicas. Since it is one
    docker compose file, it should be the same replicas for each node. Easy
    product.

    The trick comes in that a couple of workers might be in transit at any
    time it is checked. Because of this, we really only care if there are two
    consecutive reads with the same idle and crew that does not add up to the
    easy product.

    It will miss situations where the worker dies and reincarnation does not
    take place. How often is this?
    '''
    fn = f'/tmp/{urlparse(args.url).netloc}.mia.pkl'
    total = len(args.nodes) * args.replicas
    response = requests.get(
        urljoin(args.url, 'app/schedule/crew'), timeout=90, verify=args.ca_path
    )
    response.raise_for_status()
    current = response.json()
    current['idle'] = int(current['idle'])
    current['busy'] = sorted(
        task.split('duration')[0].strip() for task in current['busy']
    )
    alive = current['idle'] + len(current['busy'])
    if alive != total:
        previous = {}
        if os.path.isfile(fn):
            with open(fn, 'br') as file:
                previous = pickle.load(file)
        if current == previous:
            if alive < total:
                # have workers on strike that have not returned to the workforce
                email.send(
                    args,
                    f'''
 Some {total - alive} workers are refusing reincarnation. Please check the docker  compose file because this should not happen. Will attempt to repatriate it.
                ''',
                )
                perform.repatriation(perform.fallen(args.nodes)[1])
            else:
                # have undead workers that are both idle and supposedly working
                if current['idle'] == total:
                    email.send(
                        args,
                        f'''
 The pipeline is reporting {alive - total} undead processes and no other tasks. Resetting the pipeline to purge the undead workers and restore the pipeline availability.
                    ''',
                    )
                    perform.reset()
                else:
                    fallen, where = perform.fallen(args.nodes)
                    if fallen:
                        email.send(
                            args,
                            f'''
 The pipeline is reporting {alive - total - fallen} undead workers and {fallen} workers. Attempting to repatriate the fallen.
                            ''',
                        )
                        perform.repatriation(where)
            pass
        else:
            with open(fn, 'bw') as file:
                pickle.dump(current, file)
        return 1
    if os.path.isfile(fn):
        os.unlink(fn)
    return 0
