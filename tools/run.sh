#! /usr/bin/env bash
#
# Script runs the work flow locally. This is not a general script that runs any
# workflow like act. It is clearning and build to process these workflows since
# they are all all Ubuntu based and require a virtual Python 3.12 environment.

do_job()
{
    python -m venv $3/venv
    . $3/venv/bin/activate
    python -m pip install pyyaml
    python <<EOF > $3/do.sh
import os
import sys
import yaml

with open ("$1", 'rt') as file: y = yaml.safe_load(file)
lc = 0
steps = []
if "$2" == 'PyTesting':
    print ('docker pull postgres:latest')
    print ('docker run --detach --env POSTGRES_PASSWORD=password '
           '--env POSTGRES_USER=tester --name ga_postgres --network host '
           '--rm  postgres:latest')
    print ('sleep 5')
    print ('docker exec -i ga_postgres createdb -U tester testspace')
for step in y['jobs']['$2']['steps']:
    if 'run' in step:
        for line in filter (len, step['run'].split('\n')):
            if line.startswith ('sudo'): continue
            if line.startswith ('createdb'): continue
            if '&&' in line: line = line[:line.find('&&')]
            if '||' in line: line = line[:line.find('||')]
            cmd = line.split()[0]
            lc += 1
            if cmd == 'black' and 'KEEP_CHANGES' in os.environ:
                line = line.replace ('--check --diff ','')
            print (line,f'&& echo "result of github action step: success,{cmd}" '
                   f'|| echo "result of github action step: failure,{cmd}"')
if lc: print (f'echo "github actions expected steps: {lc}"')
if "$2" == 'PyTesting':
    print ('docker stop ga_postgres')

EOF
    chmod 755 $3/do.sh
    $3/do.sh
}

within()
{
    [[ $# -eq 1 ]] && return 0
    base=$1
    shift
    while [[ $# -gt 0 ]]
    do
        [[ "$base" = "$1" ]] && return 0
        shift
    done
    return 1
}

if [ -z "$(which python)" ]
then
    echo "'python' is not defined in your path. It should execute Python 3.12 or later."
    exit -1
fi
python <<EOF
import sys
if sys.version_info.major >= 3 and sys.version_info.minor >= 12:
    sys.exit(0)
else: sys.exit(-2)
EOF
if [ $? -ne 0 ]
then
    echo "'python' is an older version. Need to use python 3.12 or later"
    python --version
    exit -2
fi
expected_jobs=""
root=$(realpath $(dirname $0)/..)
this=$(realpath $(dirname $0)/../.github/workflows)
wdir=$this/local
export LDTK_ROOT=$wdir/.ldtk
cd $root
rm -f $root/*.rpt.txt
trap "rm -rf $wdir $root/Python/build $root/Python/dawgie.egg-info $root/Python/dawgie/fe/requirements.txt" EXIT
for yaml in $(ls $this/*.yaml)
do
    echo "yaml: $(basename $yaml)"
    for job in $(python <<EOF
import sys
import yaml

with open ("$yaml", 'rt') as file: y = yaml.safe_load(file)
print(' '.join(y['jobs'].keys()))
EOF
                 )
    do
        if $(within $job $*)
        then
            echo "   queuing job: $job"
            expected_jobs="${expected_jobs} $job"
            mkdir -p $wdir/$job
            do_job $yaml $job $wdir/$job > $root/$job.rpt.txt 2>&1 &
        else
            echo "   skipping $job"
        fi
    done
done
echo "waiting for jobs to complete"
wait
rm -rf $wdir $root/build $root/excalibur.egg-info ${root}/.coverage ${root}/.pytest_cache
declare -i summary=0
for job in $expected_jobs
do
    if [ -f $root/$job.rpt.txt ]
    then
        python <<EOF
import os,sys
def extract (frm:str, stmnt:str)->[str]:
    result = []
    idx = frm.find(stmnt)
    while idx >= 0:
       result.append (frm[idx+len(stmnt):frm.find('\n',idx)].strip())
       idx = frm.find(stmnt, idx+len(stmnt))
    return result
with open ("$root/$job.rpt.txt", 'rt') as file: content = file.read()
count = int(extract (frm=content,
                     stmnt='github actions expected steps: ')[0])
results = extract (frm=content,
                   stmnt='result of github action step: ')
summary = False
if len(results) != count:
    print ('Failure: $job did not produce all of the exptected messages. See $root/$job.rpt.txt for details.')
else:
    summary = all(r.startswith('success') for r in results) 
    if summary: print ('Success: $job')
    else: print ('Failure: $job - see $root/$job.rpt.txt for details')
for result in results:
    state,cmd = result.split(',')
    print (' ', state+':', cmd)
if summary:
    if 'KEEP_REPORTS' not in os.environ: os.unlink("$root/$job.rpt.txt")
    sys.exit(0)
else: sys.exit(1)
EOF
        summary=$summary+$?
    else
        echo "Failure: $job did not create expected report."
        summary=1
    fi
done
[ -z "${KEEP_CHANGES+x}" ] && rm -f $root/.github/workflows/bandit_full.json
[[ $summary -eq 0 ]] && echo "Summary: All test verifications were successful" || echo "Summary: Some or all test verifications failed." 
trap - EXIT
