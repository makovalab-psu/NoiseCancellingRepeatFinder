#!/usr/bin/env bash
# Evaluate *one* NCRF alignment scoring set, based on its ability (as a
# classifier) to discover repeats known to be embedded in many simulated reads.
#
# Reads:
#   ${seed}_${simulationNum}.noisy.fa.gz
#   ${seed}_${simulationNum}.noisy.truth.dat
#   where ${simulationNum} ranges from 1 to ${numSimulations}
#
# Writes:
#   ${seed}_${simulationNum}.unfiltered.rates.dat
#   ${seed}_${simulationNum}.classifier.dat

ncrfDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

# settings

technology="$1"
seed=$2
numSimulations=$3
motif=$4
minAlignmentLength=$5
maxAlignmentNoise=$6

#=== evaluate the scoring set for its performance on each simulated read

if [[ ${technology} = *"="* ]]; then
    scoring=`${NCRF} ${technology} --report:scoring`
else
    scoring=`${NCRF} --scoring=${technology} --report:scoring`
    fi

simulationNum=0
while [ ${simulationNum} -lt ${numSimulations} ]; do
    simulationNum=$((simulationNum+1))
    echo "=== ${seed}_${simulationNum}.noisy.fa ==="
    #
    gzip -dc ${seed}_${simulationNum}.noisy.fa.gz \
      | ${ncrfDir}/NCRF \
          ${scoring} \
          ${motif} --minlength=${minAlignmentLength} \
          --maxnoise=${maxAlignmentNoise}% \
          --stats=events --positionalevents \
      > ${seed}_${simulationNum}.unfiltered.ncrf
    #
    cat ${seed}_${simulationNum}.unfiltered.ncrf \
      | ${ncrfDir}/ncrf_extract_event_matrix.py --withheader \
      | ${ncrfDir}/event_matrix_to_rates.py \
      > ${seed}_${simulationNum}.unfiltered.rates.dat
    #
    cat ${seed}_${simulationNum}.unfiltered.ncrf \
      | ${ncrfDir}/ncrf_summary.py \
      | ${ncrfDir}/observed_vs_truth.py ${seed}_${simulationNum}.noisy.truth.dat \
          --detail=${seed}_${simulationNum}.classifier.detail.dat \
      > ${seed}_${simulationNum}.classifier.dat
    rm ${seed}_${simulationNum}.unfiltered.ncrf
    #
    cat ${seed}_${simulationNum}.classifier.dat
    cat ${seed}_${simulationNum}.unfiltered.rates.dat
    done
