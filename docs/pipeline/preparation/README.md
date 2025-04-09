# Preparation

If you want to interact with a pipeline, all you need is a broswer. See "How do I [interact with](../interact) a pipeline" for details. For everything else, follow these instructions:

1. [Install](#Install) required software
1. [Clone](#Clone) the repository
1. [Create](#Create) virtual environment
1. [Build](#Build) the docker images

## Install

You will need to have installed:
1. docker
1. git
1. Python 3.12

## Clone

Now that git is installed, clone the repostiory: `git clone <repo>` where `<repo>` is from the big green "Code" button at the [repository site](https://github.com/nasa-jpl-exoplanet/esp).

Now do `cd esp`. This directory - do `pwd` to see what it is - is your repository root directory. Remember this because it will be used often later.

## Create

With the repostory [cloned](#Clone) and Python 3.12 installed, create a virtual environment: `python3 -m venv ${HOME}/.venv/esp`. Anytime you create a new terminal window, you will have to activate this virtual environment: `source ${HOME}/.venv/esp/bin/activate`. The activation is for the bash shell. Use the correct `activate` for your shell.

Activate your virtual environment now. Once the virtual environment is activated, the correct Python command is `python`.

With you virtual environment activated and in your repository root, do `python -m pip install -r requirements.txt`. When pip completes, your environment will contain all of the software that esp requires. For completeness, also do `python -m pip install -e .`.

With the python environment complete, need to create a docker environment as well. Look at the files in the envs directory at the repository root. All of these are setting up virtual docker envirnoments. Pick and choose as desired. If sharing the platform with other users running pipelines, the port number range (DAWGIE_FE_PORT thru and including DAWGIE_SFE_PORT) must always be 6 in length and not overlap with any other user. From here out, this will be referred to as `<environment profile>`.

## Build

If the docker images are not already available, then they have to be built from the repository root:
1. `source envs/<environment profile>`
1. `docker compose -f .docker/compose.yaml build base`
1. `docker compose -f .docker/compose.yaml build server tools workers`

Once complete, you should have four images:
1. esp/worker
1. esp/tools
1. esp/server
1. esp/base

You only need to keep 1-3 inclusive.

