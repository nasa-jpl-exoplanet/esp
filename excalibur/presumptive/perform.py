'''utility module to perform actions as sdppiped on mentor0'''

import subprocess
import time


def reboot():
    '''shutdown then restart the pipeline

    First, shutdown all of the workers. Then, shutdown the operational pipeline.
    Wait some amount of time for the system to stabilize (paranoia more than
    need). Restart just the operational pipeline. Once the pipeline is running,
    the presumptive despot should notice the missing workers and retrench them.
    '''
    subprocess.run('${HOME}/farm.sh down', check=True, shell=True)
    subprocess.run('${HOME}/stop_ops.sh', check=True, shell=True)
    time.sleep(5)  # paranoia
    subprocess.run('${HOME}/run_ops.sh', check=True, shell=True)


def repatriation():
    # find workers that have stopped and restart them individually
    pass

def reset():
    '''do a pp_reset.sh ops'''
    subprocess.run(
        '/proj/sdp/ops/ae/tools/pp_reset.sh ops', check=True, shell=True
    )


def retrenchment(deploy=False):
    '''try to bring back all of the workers'''
    subprocess.run('${HOME}/farm.sh down', check=True, shell=True)
    time.sleep(5)  # paranoia
    if deploy:
        subprocess.run('${HOME}/deploy.sh', check=True, shell=True)
        time.sleep(5)  # paranoia
    subprocess.run('${HOME}/farm.sh up', check=True, shell=True)
