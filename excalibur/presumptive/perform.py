'''utility module to perform actions as sdppiped on mentor0'''

import subprocess
import time


def fallen(nodes):
    '''report the number of workers that are not running'''
    cmd = 'docker ps -aq --filter "name=ops-workers*" --filter "status=exited"'
    count = 0
    downed = []
    for node in nodes:
        result = subprocess.run(
            f'ssh -o ConnectTimeout=60 mentor{node} {cmd}',
            capture_output=True,
            check=True,
            shell=True,
        )
        ids = result.stdout.decode('utf-8').strip().split()
        count += len(ids)
        if ids:
            downed.append(node)
    return count, downed


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


def repatriation(nodes):
    '''restart the stopped workers without restarting the whole farm

    Had to move the command to a script because this was not working as an
    inline set:
    cmd = 'docker ps -aq --filter "name=ops-workers*" --filter "status=exited"'
    cmd = f'docker rm $({cmd}) && ${{HOME}}/run_workers.sh'
    '''
    cmd = './repatriate.sh'
    for node in nodes:
        subprocess.run(
            f'ssh -o ConnectTimeout=60 mentor{node} {cmd}',
            check=True,
            shell=True,
        )
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
