#!/bin/bash

set -e

# ============================================================
# USER CONFIGURATION
# ============================================================

BASE="/nfs/cms/cepeda/FCC/LFV/makeHistos"

#SAMPLE="p8_ee_Ztautau_ecm91_V3"
#CHANNEL="emu"
#TAG="Py8_long_sel"

SAMPLE="Py8_ONLYSIM_2M_LFV_V2"
CHANNEL="emu"
TAG="LFV_long_sel"


INPUT_DIR="${BASE}/../makeTrees/TREE_${SAMPLE}"

OUTPUT_DIR="${BASE}/HIST_${SAMPLE}_${CHANNEL}_${TAG}"

PMAX=50
MMAX=160

SCHEDD="condorsc1.ciemat.es"

# ============================================================

mkdir -p "${OUTPUT_DIR}"
mkdir -p condor_logs

rm -f joblist.txt

N=0

for INPUT in "${INPUT_DIR}"/*.root; do

    [ -e "${INPUT}" ] || continue

    # e.g. tree_157.root -> 157
    BASENAME=$(basename "${INPUT}" .root)
    JOBID="${BASENAME#tree_}"

    echo "${INPUT} ${JOBID} ${OUTPUT_DIR} ${CHANNEL} ${PMAX} ${MMAX}" \
        >> joblist.txt

    N=$((N+1))
done

echo "============================================================"
echo "Sample      : ${SAMPLE}"
echo "Input dir   : ${INPUT_DIR}"
echo "Output dir  : ${OUTPUT_DIR}"
echo "Channel     : ${CHANNEL}"
echo "Jobs        : ${N}"
echo "============================================================"

if [ "${N}" -eq 0 ]; then
    echo "ERROR: no ROOT files found"
    exit 1
fi

condor_submit -name "${SCHEDD}" submit_makeLFVHistos.sub
