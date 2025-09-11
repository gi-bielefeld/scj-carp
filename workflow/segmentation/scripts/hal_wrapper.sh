#!/bin/bash
source /vol/moreratsdata/cactus/cactus_env/bin/activate
TMPDIR=/vol/moreratsdata/tmp/
export TMPDIR
rm -r ./hal
cactus-hal2maf ./hal $@ 