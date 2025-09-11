#!/bin/bash
source /vol/moreratsdata/cactus/cactus_env/bin/activate
TMPDIR=/vol/moreratsdata/tmp/
export TMPDIR
rm -r ./js
cactus-pangenome ./js $@