#!/usr/bin/env bash
# Collect classifier results, for a single NCRF alignment scoring set on many
# reads, into a single table.
#
# The classifier results were computed by evaluate_simulations.
#
# Reads:
#   ${seed}_*.unfiltered.rates.dat
#   ${seed}_*.classifier.dat
#
# Writes:
#   ${seed}.simulations.classifier.dat

ncrfDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

# settings

seed=$1

#=== collect TPR/FDR/etc results

ls ${seed}_*.classifier.dat \
  | sed "s/.*_//" \
  | sed "s/\..*//" \
  | while read simulationNum ; do
      cat ${seed}_${simulationNum}.classifier.dat \
        | tr "/" " " \
        | awk '/^TPR/      { TPR=$2/$3; }
               /^FDR/      { FDR=$2/$3; }
               /^DETECTED/ { DETECTED=$2/$3; }
               END         {
                           printf("%s %.6f %.6f %.1f\n",
                                  num,100*TPR,100*FDR,100*DETECTED);
                           }' num="${simulationNum}"
      done \
  | awk '{
         if (NR==1) print "#num","TPR","FDR","DETECTED";
         print $0;
         }' \
  > farf.${seed}.simulations.classifier.xxx

#=== collect rates

ratesHeader=`cat ${seed}_1.unfiltered.rates.dat | sed "s/=[0-9.]*%//g"`
ls ${seed}_*.unfiltered.rates.dat \
  | grep -v noisy \
  | sed "s/.*_//" \
  | sed "s/\..*//" \
  | while read simulationNum ; do
      cat ${seed}_${simulationNum}.unfiltered.rates.dat \
        | tr -d "%" \
        | tr "\t" " " \
        | sed "s/[a-z]*=//g" \
        | awk '{ print num,$0 }' num="${simulationNum}"
      done \
  | awk '{
         if (NR==1) print "#num",header;
         print $0;
         }' header="${ratesHeader}" \
  > farf.${seed}.simulations.classifier.yyy

#=== collect TPR/FDR/etc with rates

cat farf.${seed}.simulations.classifier.xxx \
    farf.${seed}.simulations.classifier.yyy \
  | awk '/^#/  {
               if (header=="") header = "#num"substr($0,5);
               else            header = header""substr($0,5);
               }
         !/^#/ {
               simNum = $1;  $1=""
               if (!(simNum in sim)) sim[simNum] = simNum" "substr($0,2);
               else                  sim[simNum] = sim[simNum]" "substr($0,2);
               }
         END   {
               print header;
               for (simNum in sim) print sim[simNum]
               }' \
  | awk '/^#/  { print $0; }
         !/^#/ { print $0 | "sort -nr -k 2" }' \
  | tr " " "\t" \
  > ${seed}.simulations.classifier.dat
rm farf.${seed}.simulations.classifier.xxx
rm farf.${seed}.simulations.classifier.yyy
