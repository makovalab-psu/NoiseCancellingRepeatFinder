#!/usr/bin/env bash

ncrfDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

seed=$1
iterationNum=$2
motif=$3
minAlignmentLength=$4
maxAlignmentNoise=$5

source ${ncrf_dev}/ncrf_aware.sh

numScoringSets=`wc -l ${seed}.iter${iterationNum}.scoring_sets | awk '{ print $1 }'`

num=0
time cat ${seed}.iter${iterationNum}.scoring_sets \
  | while read scoring ; do
      num=$((num+1))
      echo "=== #${num} of ${numScoringSets}: ${scoring} ===" | tr "_" " "
      #
      gzip -dc ${seed}.noisy.fa.gz \
        | ${ncrfDir}/NCRF \
            `echo ${scoring} | tr "_" " "` \
            ${motif} --minlength=${minAlignmentLength} \
            --maxnoise=${maxAlignmentNoise}% \
            --stats=events --positionalevents \
        > ${seed}.${scoring}.unfiltered.ncrf
      #
      cat ${seed}.${scoring}.unfiltered.ncrf \
        | ${ncrfDir}/ncrf_summary.py \
        > ${seed}.${scoring}.unfiltered.summary
      #
      cat ${seed}.${scoring}.unfiltered.ncrf \
        | ${ncrfDir}/ncrf_extract_event_matrix.py --withheader \
        | tee ${seed}.${scoring}.unfiltered.events.dat \
        | ${ncrfDir}/event_matrix_to_rates.py \
        > ${seed}.${scoring}.unfiltered.rates.dat
      #
      cat ${seed}.${scoring}.unfiltered.summary \
        | ${ncrfDir}/observed_vs_truth.py ${seed}.noisy.truth.dat \
        > ${seed}.${scoring}.classifier.dat
      done
