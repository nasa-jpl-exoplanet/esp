#! /usr/bin/env bash

usage()
{
    echo "usage: $0 <environment profile>"
    echo "   <environment profile> : the name of the environment variable set to"
    echo "                           used when starting the private pipeline."
    echo ""
    echo "example: $0 alsMT"
    echo "example: $0 # defaults to username for <environment profile>" 
    echo ""
}

[[ ${1:-""} == "-?" ]] && usage && exit 0
[[ ${1:-""} == "-h" ]] && usage && exit 0
[[ ${1:-""} == "--help" ]] && usage && exit 0
[[ $# -gt 1 ]] && usage && exit -1

ep=${1:-${USER}}
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
       up --detach pipeline
