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
<<<<<<< HEAD
docker compose -f ${root}/.docker/compose.yaml run --remove-orphans tools \
       dawgie.db.tools.post2shelve -O /proj/db -p ${1:-${USER:-undefined}}
=======
docker compose \
       --file ${root}/.docker/compose.yaml \
       --project-name ${USER,,} \
       run tools \
          dawgie.db.tools.post2shelve -O /proj/db -p ${1:-${USER:-undefined}}
>>>>>>> 91a1722 (19: add project name to allow more than one pipeline per host (#21))
       
