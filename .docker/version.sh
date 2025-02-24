#! /usr/bin/env bash

rdir=$(realpath $(dirname $0)/..)
baseVersion=$(python <<EOF
try:
    import pyblake2 as hashlib
except:
    import hashlib

data = b''
with open ('$rdir/requirements.txt', 'br') as f: data += f.read()
with open ('$rdir/.docker/compose.yaml', 'br') as f: data += f.read()
with open ('$rdir/.docker/Dockerfile.base', 'br') as f: data += f.read()
with open ('$rdir/.docker/Dockerfile.server', 'br') as f: data += f.read()
with open ('$rdir/.docker/Dockerfile.tools', 'br') as f: data += f.read()
with open ('$rdir/.docker/Dockerfile.worker', 'br') as f: data += f.read()
k = hashlib.blake2b (data, digest_size=8)
print (k.hexdigest())
EOF
           )
python <<EOF
import os,sys

with open ('$rdir/.docker/.env', 'rt') as f: config = f.read()
k = 'ESP_VERSION='
i = config.find (k)
if i > -1:
    v = config[i+len(k):config.find('\n',i+len(k))]
    if 'KEEP_CHANGES' in os.environ:
        config = config.replace (v, '$baseVersion')
        v = '$baseVersion'
        with open ('$rdir/.docker/.env', 'tw') as f: f.write(config)
    if v == '$baseVersion': sys.exit(0)
    else:
        print (f'found version $baseVersion but expected version {v}')
        sys.exit(-2)
else:
    print ("failed to find the ESP_VERSION in .docker/.env")
    sys.exit(-1)
EOF
