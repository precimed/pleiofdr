#!/usr/bin/env bash
# Run the chr21 demo inside a container, offline, with MATLAB's prune indices, and require the
# CSV tables to be byte-identical to the MATLAB reference outputs in tests/fixtures.
#
#   tests/container/demo_parity.sh DEMO_DIR -- CONTAINER_COMMAND...
#
# DEMO_DIR holds the unpacked pleioFDR_demo_data.tar.gz and must be mounted at /data by
# CONTAINER_COMMAND, which must run `pleiofdr` (e.g. `docker run --rm -v DEMO_DIR:/data IMAGE`).
# Needs uv on the host to extract the prune indices from the MATLAB fixture.
set -euo pipefail
demo=$(cd "$1" && pwd)
shift
[ "$1" = "--" ] && shift
repo=$(cd "$(dirname "$0")/../.." && pwd)

uv run --no-project --with scipy python - "$repo/tests/fixtures/golden_conjfdr.mat" "$demo/pruneidx.mat" <<'PY'
import sys
import scipy.io
g = scipy.io.loadmat(sys.argv[1], variable_names=["pruneidx"])
scipy.io.savemat(sys.argv[2], {"pruneidx": g["pruneidx"].astype(bool)})
PY

cat > "$demo/config_parity.txt" <<'CFG'
reffile=/data/ref_1kgPhase3eur_LDr2p1_DEMO.mat
traitfolder=/data
traitfile1=CTG_COG_2018_DEMO.mat
traitname1=COGchr21
traitfiles={'SSGAC_EDU_2016_DEMO.mat'}
traitnames={'EDUchr21'}
stattype=conjfdr
fdrthresh=0.05
randprune=true
randprune_n=20
randprune_file=/data/pruneidx.mat
exclude_chr_pos=[]
reset_pruneidx=true
randprune_repeats=default
pthresh=1
perform_gc=true
use_standard_gc=false
randprune_gc=true
exclude_from_discovery=false
mafthresh=0.005
exclude_ambiguous_snps=true
onscreen=false
outputdir=/data/out_parity
CFG

rm -rf "$demo/out_parity"
"$@" --config /data/config_parity.txt > "$demo/parity.log"

status=0
for ref in "$repo"/tests/fixtures/matlab_conjfdr/*.csv; do
    name=$(basename "$ref")
    if cmp -s "$ref" "$demo/out_parity/$name"; then
        echo "identical to MATLAB: $name"
    else
        echo "DIFFERS from MATLAB: $name"
        status=1
    fi
done
exit $status
