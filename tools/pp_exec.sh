#! /usr/bin/env bash

usage()
{    
    echo "usage: $0 <task name> <target name> <environment_profile>"
    echo "   <task name> the python package name after excalbur as in cerberus"
    echo "               given the full package nameexcalibur.cerberus"
    echo ""
    echo "   <target name> is the exactly that"
    echo ""
    echo "   <environment profile> : the name of the environment variable set to"
    echo "                           used when starting the private pipeline."
    echo ""
    echo "   All of <environment profile>,  <task name>, and <target name> are"
    echo "   required."
    echo ""
    echo "   Environment variable RUNID can be used to set the runid for this"
    echo "   job. It too is optional and defaults to 17 if not set."
    echo ""
    echo "example: $0 system 'GJ 1214' alsMT"
    echo "example: $0 system 'GJ 1214' # defaults to username for <environment profile>"
    echo ""
}

[[ ${1:-""} == "-?" ]] && usage && exit 0
[[ ${1:-""} == "-h" ]] && usage && exit 0
[[ ${1:-""} == "--help" ]] && usage && exit 0
[[ $# -lt 2 ]] && usage && exit -1
[[ $# -gt 3 ]] && usage && exit -1

ep=${3:-${USER}}
root=$(realpath $(dirname $0)/..)

if [ -f $ep ]
then
    . $ep
else
    if [ -f $root/envs/$ep ]
    then
        . $root/envs/$ep
    else
        echo "Could not resolve $ep"
        exit -1
    fi
fi

docker compose \
       --file $root/.docker/compose.yaml \
       --project-name $(basename ${ep,,}) \
       exec \
          -e DISPLAY=${DISPLAY} \
          -e RUNID=${RUNID:-17} \
          -e TARGET_NAME="${2}" \
          pipeline python3 -m excalibur.${1}
