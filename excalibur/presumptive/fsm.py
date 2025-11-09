'''checks that the pipeline is available

1. if the pipeline is running, it does nothing
2. if the pipeline is not running, it checks the previous state to see how
   long it has not been running
3. if the pipeline has not been in the running state for some period of time,
   then restart it
'''

import os
import pickle
import requests

from . import email
from . import perform

from datetime import datetime
from datetime import timedelta

from urllib.parse import urljoin
from urllib.parse import urlparse


def check(args):
    fn = f'/tmp/{urlparse(args.url).netloc}.fsm.pkl'
    response = requests.get(urljoin(args.url, 'app/pl/status'), timeout=90)
    response.raise_for_status()
    current = response.json()
    reset = False
    if current['status'] != 'active' or current['name'] != 'running':
        # set previous to current if there is not a previous on disk
        if os.path.isfile(fn):
            with open(fn, 'br') as file:
                previous = pickle.load(file)
        else:
            previous = current.copy()
            previous['when'] = datetime.now()
        # if there is a change in state from the previous, reset the timer
        if (
            current['status'] != previous['status']
            or current['name'] != previous['name']
        ):
            previous = current.copy()
            previous['when'] = datetime.now()
        # have we surpased the desired wait time
        duration = datetime.now() - previous['when']
        if duration.total_seconds() > args.threshold:
            previous['when'] = datetime.now() + timedelta(days=1000)
            if current['status'] == 'active' and current['name'] == 'loading':
                msg = f'''
 The pipeline is stuck in "loading" for {duration.total_seconds()} seconds. Restarting the pipeline does not make sense because there are probably messages in /proj/sdp/data/logs/ops.log that will indicate why it has not finished loading. A pipeline restart should result in the same condition. Hence, you need to read the logs, fix the bug, and then restart the pipeline either through a github.com merge or manually.
                '''
            else:
                msg = f'''
 It has been {duration.total_seconds()} seconds since first detected the change to status "{current['status']}" and state "{current['name']}". The duration {duration.total_seconds()} seconds is greater than the desired threshold {args.threshold} seconds given to me. Restarting the pipeline now.
            '''
                reset = True
            email.send(args, msg)
            if reset:
                perform.reboot()
        # save new previous
        with open(fn, 'bw') as file:
            pickle.dump(previous, file)
        return 1
    if os.path.isfile(fn):
        os.unlink(fn)
    return 0
