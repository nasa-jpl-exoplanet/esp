#! /usr/bin/env bash

repodir=$(realpath $(dirname $0)/..)
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
pyxbver="$(grep -i pyxb-ctc ${repodir}/requirements.txt | awk -F '==' '{print $2}' | awk '{print $1}')"
if [ -z "$(which pyxbgen)" ]
then
    echo "'pyxbgen' is not defined in your path. It should execute PyXB-CTC ${pysbver}."
    exit -1
fi

if [[ $(pyxbgen --version) == *"${pyxbver}"* ]]
then
    :
else
    echo "pyxbgen is the wrong version"
fi

cd ${repodir}
pyxbgen --schema-location=excalibur/runtime/levers.xsd --module=binding --module-prefix=excalibur.runtime

cd ${repodir}/excalibur/runtime
md5sum levers.xsd binding.py > autogen.md5
