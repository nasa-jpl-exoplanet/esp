#! /usr/bin/env bash

make_cert () {
    openssl req -new -nodes -newkey rsa:2048 \
            -keyout ${1}.key -out ${1}.csr \
            -subj "/CN=$(id -un)"
    sudo -u sdppiped /proj/sdp/bin/sign.sh ${1}.csr signed.public.pem.$1
    cat ${1}.key /proj/sdp/data/certs/signed.public.pem.$1 > ${1}.ca.signed.pem
    chmod 600 ${1}.ca.signed.pem
}

selfsign() {
    EXT=$(mktemp)
    trap 'rm -f "${EXT}"' EXIT
    cat > "${EXT}" <<EOF
basicConstraints=critical,CA:FALSE
keyUsage=critical,digitalSignature,keyEncipherment
extendedKeyUsage=clientAuth,serverAuth
EOF
    openssl x509 -req -in ${1}.csr -signkey ${1}.key -days 36500 \
            -out ${1}.crt -extfile "${EXT}"
    cat ${1}.key ${1}.crt > ${1}.self.signed.pem
    chmod 600 ${1}.self.signed.pem
}

loc=${1:-${HOME}/.ssh}
user=$(id -un)

mkdir -p ${loc}
cd ${loc}
make_cert ${user}
selfsign ${user}
