'''utility module to perform actions as sdppiped on mentor0'''

import subprocess
import time


def fallen():
    '''report the number of workers that are not running'''
    return 0


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
    '''restart the stopped workers without restarting the whole farm'''
    subprocess.run('${HOME}/repatriate.sh', check=True, shell=True)
    pass


def reset():
    '''do a pp_reset.sh ops'''
    subprocess.run(
        '/proj/sdp/ops/ae/tools/pp_reset.sh ops', check=True, shell=True
    )


def retrenchment(node):
    '''try to bring back all of the workers'''
    subprocess.run(
        f'ssh -o ConnectTimeout=60 mentor{node} "${{HOME}}/deploy.sh 1"',
        check=True,
        shell=True,
    )
    subprocess.run(
        f'ssh -o ConnectTimeout=60 mentor{node} "${{HOME}}/run_workers.sh"',
        check=True,
        shell=True,
    )


def worker_id(node=None):
    cmd = 'docker images -q esp/worker'
    if node:
        cmd = f'ssh -o ConnectTimeout=60 mentor{node} "{cmd}"'
    result = subprocess.run(cmd, capture_output=True, check=True, shell=True)
    return result.stdout.decode('utf-8').strip()
