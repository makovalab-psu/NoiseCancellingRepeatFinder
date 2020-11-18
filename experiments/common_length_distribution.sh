#!/usr/bin/env bash
# Generate samples from two sets of fasta files, so that they have a common
# length distribution
#
# The general idea is to subsample both sets so that each sampled set will have
# identical length distributions. However, the technique used here settles for
# length distributions that are similar rather than identical.
#
# Requiring identical distributions would bias against longer reads. Specific
# long read lengths are rare -- for example, if one set has a single read
# of length 300,000 and the other set has one of length 300,007, neither would
# be chosen in a common *identical* distribution.
#
# To allow for that situation, we bin similar lengths then create distributions
# that are identical with respect to the bin counts. In the long read example
# above, the bin covering lengths 300,000 and 300,007 would have a non-zero
# count in both sets, so the common distribution would contain a non-zero count
# for that bin.
#
# We use a logarithmic binning formula so longer lengths have much wider bins.
# At length 400 the bin size is 5, while at length 300K bin it's about 10K.
# This is accomplished using the formula bin=30.5*log(pos-250) and its inverse
# pos=250+exp(bin/30.5). This binning function was determined by working
# backwards from the stated bin sizes. The formula obviously cannot be applied
# to lengths at or below 250, thus in this experiment we discard any lengths
# below that limit (in fact, we use a higher limit, 375). 
#
# Note that we do not modify sequences. It's possible to achieve distribution
# similarity by truncating some sequences, but we do not do that.
     
# binning formula

minReadLen=375                     # lengths shorter than this are discarded
toBin="30.5*log(pos-250)"          # expression to map a length to its bin
fromBin="250+exp(bin/30.5)"        # expression to map a bin to its length(s)

# define the two sets; each set consists of several fasta files; here
# we have two files in each set

fastaSet1="SRR2036701 SRR2036850"  # for each name in these two lists, we have,
fastaSet2="FAF01127 FAF01169"      # .. name.fa.gz; e.g. SRR2036701.fa.gz

# programs used

fasta_length_distribution="      ../extra/common_length_distribution/fasta_length_distribution.py"
bin_position_counts="            ../extra/common_length_distribution/bin_position_counts.py"
common_length_distribution="     ../extra/common_length_distribution/common_length_distribution.py"
make_length_distribution_spec="  ../extra/common_length_distribution/make_length_distribution_spec.py"
spread_length_distribution="     ../extra/spread_length_distribution.py"
fasta_match_length_distribution="../extra/common_length_distribution/fasta_match_length_distribution.py"

#=== compute the sequence length distributions

echo "${fastaSet1} ${fastaSet2}" | tr " " "\n" \
  | while read f ; do
      echo "computing binned length distribution for ${f}.fa.gz"
      gzip -dc ${f}.fa.gz \
        | ${fasta_length_distribution} \
        | tee ${f}.length_distrib.dat \
        | ${bin_position_counts} "bin=${toBin}" "pos=${fromBin}" \
            --minpos=${minReadLen} \
        > ${f}.length_distrib.binned.dat
      done

echo "computing binned length distribution for fastaSet1"
echo ${fastaSet1} | tr " " "\n" \
  | while read f ; do
      cat ${f}.length_distrib.dat
      done \
  | ${bin_position_counts} "bin=${toBin}" "pos=${fromBin}" \
      --minpos=${minReadLen} \
  > fastaSet1.length_distrib.binned.dat

echo "computing binned length distribution for fastaSet2"
echo ${fastaSet2} | tr " " "\n" \
  | while read f ; do
      cat ${f}.length_distrib.dat
      done \
  | ${bin_position_counts} "bin=${toBin}" "pos=${fromBin}" \
      --minpos=${minReadLen} \
  > fastaSet2.length_distrib.binned.dat

# (cleanup)

echo "${fastaSet1} ${fastaSet2}" | tr " " "\n" \
  | while read f ; do
      rm ${f}.length_distrib.dat
      done

#=== compute the minimum common distribution, reduced by 10%, and apply it to
#=== each set

echo "computing common distribution"
${common_length_distribution} --scale=90% \
    fastaSet1.length_distrib.binned.dat \
    fastaSet2.length_distrib.binned.dat \
  > target.length_distrib.binned.dat

echo "making common distribution spec for fastaSet1"
${make_length_distribution_spec} \
    target.length_distrib.binned.dat \
    fastaSet1.length_distrib.binned.dat \
  > fastaSet1.distribution_spec

echo "making common distribution spec for fastaSet2"
${make_length_distribution_spec} \
    target.length_distrib.binned.dat \
    fastaSet2.length_distrib.binned.dat \
  > fastaSet2.distribution_spec

echo "spreading common distribution over fastaSet1"
echo ${fastaSet1} | tr " " "\n" > temp.components
cat fastaSet1.distribution_spec \
  | ${spread_length_distribution} temp.components \
      --seed=fastaset1.poppy \
      --nowarn \
      --input={component}.length_distrib.binned.dat \
      --output={component}.distribution_spec \
      --progress=written
rm temp.components

echo "spreading common distribution over fastaSet2"
echo ${fastaSet2} | tr " " "\n" > temp.components
cat fastaSet2.distribution_spec \
  | ${spread_length_distribution} temp.components \
      --seed=fastaset2.poppy \
      --nowarn \
      --input={component}.length_distrib.binned.dat \
      --output={component}.distribution_spec \
      --progress=written
rm temp.components

# (cleanup)

rm target.length_distrib.binned.dat
rm fastaSet1.length_distrib.binned.dat
rm fastaSet2.length_distrib.binned.dat
rm fastaSet1.distribution_spec
rm fastaSet2.distribution_spec

echo "${fastaSet1} ${fastaSet2}" | tr " " "\n" \
  | while read f ; do
      rm ${f}.length_distrib.binned.dat
      done

# at this point we have f.distribution_spec for each f.fa.gz

#=== sample reads for each set

echo "sampling reads for fastaSet1"
echo ${fastaSet1} | tr " " "\n" \
  | while read f ; do
      gzip -dc ${f}.fa.gz \
        | ${fasta_match_length_distribution} ${f}.distribution_spec \
            --seed=fastaSet1.watermelon
      done \
  | gzip \
  > fastaSet1.selected.fa.gz

echo "sampling reads for fastaSet2"
echo ${fastaSet2} | tr " " "\n" \
  | while read f ; do
      gzip -dc ${f}.fa.gz \
        | ${fasta_match_length_distribution} ${f}.distribution_spec \
            --seed=fastaSet2.watermelon
      done \
  | gzip \
  > fastaSet2.selected.fa.gz

#== compute sampled distribution of each set, for a sanity check

echo "computing length distribution for sampled fastaSet1"
gzip -dc fastaSet1.selected.fa.gz \
  | ${fasta_length_distribution} \
  | tee fastaSet1.sampled.length_distrib.dat \
  | ${bin_position_counts} "bin=${toBin}" "pos=${fromBin}" \
      --minpos=${minReadLen} \
  > fastaSet1.sampled.length_distrib.binned.dat

echo "computing length distribution for sampled fastaSet2"
gzip -dc fastaSet2.selected.fa.gz \
  | ${fasta_length_distribution} \
  | tee fastaSet2.sampled.length_distrib.dat \
  | ${bin_position_counts} "bin=${toBin}" "pos=${fromBin}" \
      --minpos=${minReadLen} \
  > fastaSet2.sampled.length_distrib.binned.dat

