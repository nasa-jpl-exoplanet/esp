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

myself=${DAWGIE_SSL_PEM_MYSELF:-/proj/sdp/${EXCALIBUR_USER:-${USER:-nouser}}/certs/excalibur_identity.pem}

curl -XPOST --cert ${myself} "https://${EXCALIBUR_HOST:-localhost}:${DAWGIE_SFE_PORT:-9995}/app/reset?archive=${DAWGIE_ARCHIVE:-false}"
