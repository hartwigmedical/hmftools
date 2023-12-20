#!/bin/bash

set -x

[[ $# -ne 2 ]] && echo "Provide working directory and version string" && exit 1

REMOTE="gs://pipeline5-jars/${2}.tar"
gsutil stat ${REMOTE} 2>&1 >/dev/null
[[ $? -eq 0 ]] && echo "ERROR: ${REMOTE} already exists" && exit 1

tar --exclude='*-dependencies.jar' --exclude='*hmf-common*' --exclude='*-shaded.jar' -cvf ${1}/${2}.tar $(find ${1} -name '*.jar')
gsutil cp ${1}/${2}.tar ${REMOTE}
