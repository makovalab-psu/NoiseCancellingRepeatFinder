#!/usr/bin/env bash

ncrfDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

technology=$1     # (or error profile, as for mock_motif_genome)
seed=$2
numSimulations=$3
motif=$4
numRepeats=$5
muLength=$6
sigmaLength=$7
avgFillLength=$8
chromExtra=$9
readLen=${10}
minAlignmentLength=${11}

source ${ncrf_dev}/ncrf_aware.sh

${random_normal} --seed=${seed}.lengths \
    50K --mu=${muLength} --sigma=${sigmaLength} --round \
  | awk '{ if ($0 >= L) print $0 }' L=${minAlignmentLength} \
  | sort -n \
  > ${seed}.len_distrib.dat

simulationNum=0
while [ ${simulationNum} -lt ${numSimulations} ]; do
    simulationNum=$((simulationNum+1))
    echo "=== ${seed}_${simulationNum}.noisy.fa ==="
    #
    ${ncrfDir}/mock_motif_genome.py --seed=${seed}_${simulationNum}.seq \
          ${motif} --errors=${technology} \
          --lengths=${seed}.len_distrib.dat \
          L=${readLen} N=${numRepeats} --name=${seed}_${simulationNum} \
          --catalog=${seed}_${simulationNum}.truth.xxx \
      | gzip \
      > ${seed}_${simulationNum}.noisy.fa.gz
    cat ${seed}_${simulationNum}.truth.xxx \
      > ${seed}_${simulationNum}.noisy.truth.dat
    rm ${seed}_${simulationNum}.truth.xxx
    #
    cat ${seed}_${simulationNum}.noisy.truth.dat \
      | ${ncrfDir}/truth_to_rates}.py \
      > ${seed}_${simulationNum}.noisy.rates.dat
    echo "(embedded error rates: `cat ${seed}_${simulationNum}.noisy.rates.dat`)"
    done
