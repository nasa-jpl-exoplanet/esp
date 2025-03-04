#! /usr/bin/env bash

root=$(realpath $(dirname $0)/..)
export DAWGIE_DB_HOST=excalibur.jpl.nasa.gov \
       DAWGIE_DB_IMPL=post \
       DAWGIE_DB_NAME=ops \
       DAWGIE_DB_PATH=$(cat /proj/sdp/$USER/.pgpass) \
       DAWGIE_DB_PORT=5263 \
       EXCALIBUR_UID=${UID} \
       EXCALIBUR_USER=${USER} \
       EXCALIBUR_SOURCE_PATH=${root}
docker compose -f ${root}/.docker/compose.yaml run tools \
       dawgie.db.tools.post2shelve -O /proj/db -p ${1:-${USER:-undefined}}
       
