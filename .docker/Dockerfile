FROM ubuntu:24.04 AS base
COPY . /tmp/esp
ENV MPLCONFIGDIR=/tmp
ENV PATH=/usr/local/esp/bin:$PATH
RUN set -ex && \
    export DEBIAN_FRONTEND=noninteractive && \
    apt-get update && \
    apt-get install -y \
            build-essential \
            curl \
            git \
            graphviz \
            haveged \
            libhdf5-serial-dev \
            netcdf-bin \
            libnetcdf-dev \
            postgresql-client \
            python3 \
            python3-dev \
            python3-gdbm \
            python3-pip \
            python3-venv \
            tzdata && \
    apt-get clean && apt-get autoremove && \
    ln -fs /usr/share/zoneinfo/America/Los_Angeles /etc/localtime && \
    dpkg-reconfigure --frontend noninteractive tzdata && \
    python3 -m venv /usr/local/esp && \
    cd /tmp/esp && pip install -r requirements.txt && pip cache purge && \
    mkdir -p /proj/data/ldtk && chmod 777 /proj/data/ldtk && \
    rm -rf /tmp/esp

FROM base AS server
ENV HOME=/tmp
ENTRYPOINT python -m dawgie.pl \
    --context=excalibur.context \
    -l "${EXCALIBUR_USER}.log" -L "${EXCALIBUR_LOG_LEVEL}" -p 8080

FROM base AS tools
ENV DAWGIE_DOCKERIZED_AE_GIT_REVISION=irrelevant
ENV HOME=/tmp
ENTRYPOINT ["python", "-m"]

FROM base AS worker
COPY . /proj/src/ae
ARG AE_VER=undefined
ENV DAWGIE_DOCKERIZED_AE_GIT_REVISION=$AE_VER
ENV HOME=/tmp
ENV NUMEXPR_NUM_THREADS=1
ENV NUMEXPR_MAX_THREADS=1
ENV OMP_NUM_THREADS=1
ENV PYTHONPATH=/proj/src/ae
ENTRYPOINT python -m dawgie.pl.worker \
            -a /proj/src/ae/excalibur \
            -b excalibur \
            -c cluster \
            -g /proj/gnupg \
            -i 0 \
            -n excalibur.jpl.nasa.gov \
            -p ${DAWGIE_WORKER_PORT}
