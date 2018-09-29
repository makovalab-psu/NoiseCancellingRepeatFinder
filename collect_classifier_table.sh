#!/usr/bin/env bash

seed=$1
iterationNum=$2

scoresHeader=`cat ${seed}.iter${iterationNum}.scoring_sets \
                | head -n 1 \
                | tr "=_" "  " \
                | awk '{ print $1,$3,$5,$7,$9,$11 }'`

cat ${seed}.iter${iterationNum}.scoring_sets \
  | while read scoring ; do
      scoreValues=`echo ${scoring} \
                     | tr "=_" "  " \
                     | awk '{ print $2,$4,$6,$8,$10,$12 }'`
      #
      cat ${seed}.${scoring}.classifier.dat \
        | tr "/" " " \
        | awk '/^TPR/ { TPR=$2/$3; }
               /^FDR/ { FDR=$2/$3; }
               END    {
                      printf("%s\t%.6f\t%.6f\n",scores,100*TPR,100*FDR);
                      }' scores="${scoreValues}"
      done \
  | awk '{ print (++n),$0 }' \
  | sort -nr -k 8 \
  | awk '{
         if (NR==1) print "#num",header,"TPR","FDR";
         print $0;
         }' header="${scoresHeader}" \
  > ${seed}.iter${iterationNum}.classifier.dat
