# Noise Cancelling Repeat Finder

This package finds alignments of short tandem repeats in noisy DNA sequencing
data.  Given a list of known motifs (e.g. "GGCAT", "CCAT", etc.) and DNA
sequences in fasta files (e.g. PacBio or Oxfor Nanopore sequenced reads), it
attempts to align segments of the DNA sequences to repeated copies of the
motifs.

Note that this is *not* intended to be a turnkey solution but more of an
exploratory platform for the user. It will likely require parameter tweaking
and experimentation on the part of the user, as well as some user-written
programs to post-process the output.

### Installation

To install Noise Cancelling Repeat Finder using the source:  
1. Download the latest version of Noise Cancelling Repeat Finder using Github  
```bash  
     git clone https://github.com/makovalab-psu/NoiseCancellingRepeatFinder.git  
```  
2. Compile:  
```bash  
    cd NoiseCancellingRepeatFinder  
    make  
```

### Prerequisites

* gcc or similar C compiler and linker
* python (tested with version 2.7, not likely to work with python 3)

Python is not used by the core program, but is needed for the post-processing
helper programs included in the package.

### Usage Overview

The simplest use, searching reads.fa for the repeated motif GGCAT:

```bash 
cat reads.fa | ./NCRF GGCAT > example.ncrf
```

A more detailed example is shown in example/README.md.  A description of the
output format is included there.  The experiments subdirectory contains
additional examples, provided "as is" without explanation.

The output is usually passed through a series of the ncrf_* post-processing
scripts (e.g. ncrf_summary).


### Usage Details

```bash  
  <fasta>               fasta file containing sequences; read from stdin
  [<name>:]<motif>      dna repeat motif to search for
                        (there can be more than one motif)
  --minmratio=<ratio>   discard alignments with a low frequency of matches;
                        ratio can be between 0 and 1 (e.g. "0.85"), or can be
                        expressed as a percentage (e.g. "85%")
  --maxnoise=<ratio>    (same as --minmratio but with 1-ratio)
  --minlength=<bp>      discard alignments that don't have long enough repeat
                        (default is 500)
  --minscore=<score>    discard alignments that don't score high enough
                        (default is zero)
  --stats=events        show match/mismatch/insert/delete counts
  --positionalevents    show match/mismatch/insert/delete counts by motif
                        position (independent of --stats=events); this may be
                        useful for detecting error non-uniformity, to separate
                        perfect repeats from imperfect
  --help=scoring        show options relating to alignment scoring
  --help=allocation     show options relating to memory allocation
  --help=other          show other, less frequently used options
```

### Additional commands 

#### Additional commands 

¥¥¥
error_nonuniformity_filter
ncrf_sort
ncrf_summary
ncrf_to_bed
ncrf_cat

ncrf_error_clumps
ncrf_error_positions
ncrf_extract_event_matrix
ncrf_extract_mx_matrix
ncrf_msa
ncrf_resolve_overlaps
ncrf_words
event_matrix_to_rates

fasta_length_distribution
sum_length_distributions
bin_position_counts
common_length_distribution
make_length_distribution_spec
spread_length_distribution
fasta_match_length_distribution

sam_to_event_matrix


parameter_ball
truth_to_rates

mock_motif_genome
random_normal
reconstruct_simulated_alignments
map_onto_simulated_reads
observed_vs_truth

minimap2_cs_to_events
harvest_trf_html

ncrf_parse
echydna
interval_dict

¥¥¥

ncrf_words-- Look for prominent "wrong" motifs (words) in the output of Noise
Cancelling Repeat Finder.

```bash  
cat <output_from_NCRF> | ./ncrf_words.py [options]
   minwordratio=<r>      only show words with counts at least r times the
                         motif word's count
                         (default is 1.0)
```

error_nonuniformity_filter-- Filter the output of Noise Cancelling Repeat Finder,
removing alignments in which matches (or errors) aren't uniformly distributed
across the motif positions.

```bash  
cat <output_from_NCRF> | ./error_nonuniformity_filter.py [options]
  --effectsize=<value>  effect size for chi-squared test
                        (default is 0.3)
  --power=<probability> "power of test" for chi-squared test, 1 minus Type II
                        error probability
                        (default is 0.8)
  --discard:good        discard the "good" alignments instead of the "bad" ones
                        (by default we discard the "bad" alignments) 
  --discard:none        don't discard any alignments, just report the test
                        results (see below for how they are reported) 
  --test:matches        perform the test using match counts
                        (this is the default)
  --test:errors         perform the test using error counts; this may increase
                        the likelihood that the test cannot be performed for
                        some alignments
  --head=<number>       limit the number of input alignments
  --batch=<number>      number of input alignments processed by each call to R;
                        our ability to call R fails if the command line we pass
                        it is too long
                        (default is 250)

In a "true" alignment to a given motif unit, we expect the errors to be
distributed randomly and uniformly among the positions in the unit.  (That is
an underlying assumption, but might not itself be true.)  This program
discards alignments that fail a statistical test based on that assumption.

Since error counts may be too small for the statistical test, we use match
counts instead.

The input alignments must include position event information.  This can be
accomplished by using the --positionalevents option of Noise Cancelling Repeat
Finder.

When test results are reported for --discard:none, one of the following lines
is added to each alignment:
# positional chi-squared: match uniformity rejected
# positional chi-squared: match uniformity not rejected
# positional chi-squared: untested
The "untested" case indicates that the chi-squared test could not be
performed, usually because one of the positional match counts is too small.
```

ncrf_sort-- Sort the alignments output by Noise Cancelling Repeat Finder.

```bash  
cat <output_from_NCRF> | ./ncrf_sort.py [options]
   --sortby=mratio[-|+]    sort by decreasing or increasing match ratio
                           (by default we sort by decreasing match ratio)
   --sortby=score[-|+]     sort by decreasing or increasing alignment score
   --sortby=match[-|+]     sort by decreasing or increasing alignment match
                           count
   --sortby=length[-|+]    sort by decreasing or increasing length; length is
                           the number of sequence bases aligned
   --sortby=name           sort by sequence name (and position)
   --sortby=position[-|+]  sort by sequence name (and decreasing or increasing
                           position)
```

ncrf_summary-- Convert the output of Noise Cancelling Repeat Finder to a
summary, a tab-delimited table with one line of stats per alignment.

```bash  
cat <output_from_NCRF> | ./ncrf_summary.py [options]
  --minmratio=<ratio>  discard alignments with a low frequency of matches
                       ratio can be between 0 and 1 (e.g. "0.85"), or can be
                       expressed as a percentage (e.g. "85%")
  --maxnoise=<ratio>   (same as --minmratio but with 1-ratio)

Typical output:
  #line motif seq        start end  strand seqLen querybp mRatio m    mm  i  d
  1     GGAAT FAB41174_6 1568  3021 -      3352   1461    82.6%  1242 169 42 50
  11    GGAAT FAB41174_2 3908  5077 -      7347   1189    82.4%  1009 125 35 55
  21    GGAAT FAB41174_0 2312  3334 -      4223   1060    81.1%  881  115 26 64
   ...
```

ncrf_to_bed-- Convert the output of Noise Cancelling Repeat Finder to bed format.

```bash  
cat <output_from_NCRF> | ./ncrf_to_bed.py [options]
  --minmratio=<ratio>  discard alignments with a low frequency of matches;
                       ratio can be between 0 and 1 (e.g. "0.85"), or can be
                       expressed as a percentage (e.g. "85%")
  --maxnoise=<ratio>   (same as --minmratio but with 1-ratio)

Typical output is shown below.  The 6th column ("score" in the bed spec) is
the match ratio times 1000 (e.g. 826 is 82.6%).
  FAB41174_065680 1568 3021 . - 826
  FAB41174_029197 3908 5077 . - 824
  FAB41174_005950 2312 3334 . - 811
   ...
```

ncrf_error_positions-- Convert the output of Noise Cancelling Repeat Finder to
a position-of-errors summary, a tab-delimited table with one line of error info
per alignment, primarily to show the relative positions of errors along the
alignment.

```bash  
cat <output_from_NCRF> | ./ncrf_error_positions.py [options]
  --minmratio=<ratio>  discard alignments with a low frequency of matches;
                       ratio can be between 0 and 1 (e.g. "0.85"), or can be
                       expressed as a percentage (e.g. "85%")
  --maxnoise=<ratio>   (same as --minmratio but with 1-ratio)

Typical output:
  #line seq             strand start end   querybp mRatio nErrors errors
  1     FAB41174_065680 -      1568  3021  1461    82.6%  261     0.023 0.024 0.025 ... 0.981 0.983
  11    FAB41174_029197 -      3908  5077  1189    82.4%  215     0.021 0.029 0.032 ... 0.989 0.996
  21    FAB41174_005950 -      2312  3334  1060    81.1%  205     0.023 0.027 0.036 ... 0.979 0.995
   ...

Note that lines will usually have different numbers of fields, since they have
different error counts.  "errors" is a vector of the relative positions of
errors along the alignment, with positions ranging from zero to one.
```

ncrf_extract_mx_matrix-- Extract the positional match-and-error-counts matrix
from Noise Cancelling Repeat Finder alignments.

```bash  
cat <output_from_NCRF> | ./ncrf_extract_mx_matrix.py [options]
  --head=<number>         limit the number of input alignments

  The output matrix has R rows and 2M+1 columns, where R is the number of input
  alignments and M is the length of the aligned motif. (It is assumed that all
  alignments are to the same motif, but this is NOT validated).

  The first column is the line number of the alignment in the input file.  The
  next M columns are the positional counts for matches ("m"), and the final M
  columns are the positional counts for errors ("x").  The output is intended
  to be suitable as input to R.
```

ncrf_extract_event_matrix-- Extract the (NON-positional) match-and-error event
counts matrix from Noise Cancelling Repeat Finder alignments.

```bash  
cat <output_from_NCRF> | ./ncrf_extract_event_matrix.py [options]
  --withheader            include a header line in the output
  --sumonly               include only a summation line in the output
                          (by default, we output a separate line for each
                          alignment, and no sum)
  --head=<number>         limit the number of input alignments

  The output matrix has R rows and 9 columns, where R is the number of input
  alignments.

  The first column is the line number of the alignment in the input file. The
  second column is the motif, e.g. GGAAT.  The third column is the match ratio
  ("mRatio"). The remaining columns are, respectively, the counts for matches
  ("m"), mismatches ("mm"), insertion opens ("io"), insertion extensions
  ("ix"), deletion opens ("do"), and deletion extensions ("dx").

  The output is intended to be suitable as input to R, and can be used as input
  to infer_scoring.
```

sam_to_event_matrix-- Extract (NON-positional) match-and-error event counts
from alignments in sam format.

```bash  
cat <sam_file> | ./sam_to_event_matrix.py [options]
  --withheader            include a header line in the output
  --sumonly               include only a summation line in the output
                          (by default, we output a separate line for each
                          alignment, and no sum)
  --head=<number>         limit the number of input alignments
  --progress=<number>     periodically report how many alignments we've read

We expect the sam file to contain alignments of reads to a reference genome,
and it must contain MD tags (if the aligner did not provide these, they can be
added by using samtools calmd).

The output matrix has R rows and 9 columns, where R is the number of input
alignments.

The first column is the line number of the alignment in the input file. The
second column is the read name.  The third column is the match ratio
("mRatio"). The remaining columns are, respectively, the counts for matches
("m"), mismatches ("mm"), insertion opens ("io"), insertion extensions ("ix"),
deletion opens ("do"), and deletion extensions ("dx").

The output is intended to be suitable as input to R, and can be used as input
to infer_scoring.
```

event_matrix_to_rates-- Read match-and-error event counts and report overall
event rates.

```bash  
cat <event_counts_table> | ./event_matrix_to_rates.py [options]
  (currently, there are no options)

The input is usually the output from ncrf_extract_event_matrix.  It has 9
columns, first three of which are ignored here. The remaining columns are,
respectively, counts for matches ("m"), mismatches ("mm"), insertion opens
("io"), insertion extensions ("ix"), deletion opens ("do"), and deletion
extensions ("dx").
```

infer_scoring-- Infer scoring parameters from match-and-error event counts.
_Note that the ematch setting is currently ignored_.
_Note that this program has yet to be shown to work_.

```bash  
cat <event_counts_table> | ./infer_scoring.py [options]
  --ematch=<probability>  expectation that two nucleotides should be the same;
                          this can be computed by fasta_nt_counts
                          (by default, this is 25%)

  The input is usually the output from ncrf_extract_event_matrix.  It has 9
  columns, first three of which are ignored here. The remaining columns are,
  respectively, counts for matches ("m"), mismatches ("mm"), insertion opens
  ("io"), insertion extensions ("ix"), deletion opens ("do"), and deletion
  extensions ("dx").
```

fake_motif_read-- Create a fake "read" with embedded repeat motifs.

```bash  
./fake_motif_read.py [options]
  <motif>                  (cumulative) motif to embed
  --name=<string>          read name
  --length=<bp>       (L=) read length; if this is absent, the repeats will
                           comprise the entire read (no fill DNA)
  --length=<pctg>%    (L=) read length as a percentage of the total repeat
                           length; pctg should be more than 100.
  --length=+<bp>      (L=) read length as delata above the total repeat length
  --repeats=<number>  (N=) number of repeats to embed
                           (default is 1)
  --motif:neighbor=<prob>  probability of embedding motifs that have edit
                           distance 1 from the specified motifs
                           (default is 0)
  --motif:mixture=<prob>   probability of embedding motifs that are 50/50
                           mixtures of a specified motif and one that has edit
                           distance 1
                           (default is 0)
  --lengths=<file>         file containing the repeat length distribution, one
                           length per line; if this is absent, we'll read the
                           repeat length distribution from stdin
  --errors=pacbio          simulate pacbio error profile
                           (by default, errors are not simulated)
  --errors=nanopore        simulate nanopore error profile
  --errors=<pctg>%         simulate simple error profile with the given rate,
                           with mismatch, insertion, and deletion each of equal
                           rates
  --errors=<spec>          simulate error profile with the given spec; <spec>
                           looks like this:
                              mm:1%,i:12%,d:2%
  --catalog=<file>         file to write a catalog of the embedded repeats to
  --wrap=<length>          number of nucleotides per line; 0 means single line
  --seed=<string>          set random seed

Note that if the program is run twice with the same seed, one run with errors,
one without, the same pre-error sequence is generated.  Having the output of
both those runs can be useful.
```

observed_vs_truth-- Compare observed events from alignments to known truth.

```bash  
cat <alignment_summary> | observed_vs_truth.py <truth_catalog> [options]
  <truth_catalog>        File containing aligment "truth" (see below)
  --motif=<motif>        (cumulative) Motifs represented in the summary. Truth
                         intervals for other motifs are discarded. If this is
                         not provided, we use all truth intervals.
  --detection=<portion>  threshold for a truth intervals to be considered as
                         "suffiently" covered; 0<portion<=1, but can be
                         expressed with a % sign
                         (default is 95%)
  --detail=<file>        Report separate true positive rates for each observed
                         interval

We'll consider each base in the genome as an item to be classified. The
alignment summary tells us the classification of each base -- each base in an
alignment interval is classified as a positive, and any other base is
classified as a negative.

The truth catalog tells us the correct classification of each base -- each base
in a truth interval should be classified as a positive (by an ideal
classifier), and any other base should be classified as a negative.

We report
  true postive rate         (TPR) = TP/(TP+FN)
  false negative rate       (FNR) = FN/(TP+FN)
  positive predictive value (PPV) = TP/(TP+FP)   a.k.a. precision
  false discovery rate      (FDR) = FP/(TP+FP)
  detection rate            fraction of truth intervals "suffiently" covered by
                            alignments

Since we don't know the genome size, we can't report anything involving true
negatives.

Note that we AREN'T considering whether the aligner called the correct event
(mm, ins, or del) for a base; only whether it covered that base with an aligned
interval.

The alignment summary is usually the output from ncrf_summary. It has 13
columns but only the aligned intervals and motif are used here ("seq", "start",
"end", and "motif", columns 3, 4, 5, and 2).

The truth catalog is usually the output from the --catalog option of
fake_motif_read. It has 12 columns but only the intervals and motif are used
here ("chrom", "start", "end", and "motif", columns 1, 2, 3, and 4). Intervals
must be distinct; overlaps are not allowed.
```

truth_to_rates-- Read event counts from a truth catalog and report overall
induced event rates.

```bash  
cat <truth_catalog> | ./truth_to_rates.py [options]
  (currently, there are no options)

The truth catalog is usually the output from the --catalog option of
fake_motif_read. It has 12 columns but only the events are used here ("m",
"mm", "i", and "d", columns 9, 10, 11, and 12).
```

parameter_ball-- Vary a set of parameters, inside a "ball".

```bash  
./parameter_ball.py <parameter> [options]
  <parameter>          (cumulative) a parameter and its central value; for
                       example, X=21; the parameter will be varied unless it
                       is in the set of fixed parameters
  --fixed=<parameter>  (cumulative) a parameter NOT to vary
  --radius=<offset>    how much to vary each parameter
                       (default is 1)
  --ball:sparse        the ball is a sparse hypercube
                       (this is the default)
  --ball:hyper         the ball is a complete hypercube; if the radius is more
                       than about 5, --ball:sparse should be used instead
  --ball:spikey        the ball is a "spikey burr", with only one parameter
                       changed relative to the input
  --sample=<number>    number of parameter sets to sample from the ball; note
                       that rejected sets still count toward this limit
                       (by default we report every parameter set in the ball)
  --reject=<formula>   (cumulative) reject a parameter set if the formula is
                       true; e.g. "X>Y" would reject in set in which the
                       paremeter X was greater than the parameter Y
                       (usually the formala has to be enclosed in quotes)
  --nocenter           exclude the central point from the ball
  --seed=<string>      set random seed
```

random_normal-- Generate random numbers, sampled from a normal distribution.

```bash  
random_normal.py <num_values> [options]
  <num_values>             number of values to generate
  --mu=<value>             mean value         (default is 0.0)
  --sigma=<value>          standard deviation (default is 1.0)
  --round                  round to integers
  --floor                  "round down" to integers
  --ceiling                "round up" to integers
  --precision=<digits>     set number of digits after decimal
  --seed=<string>          set random seed
```

fasta_nt_counts-- Report the counts of A, C, G, and T in a fasta file.
_This can be used to derive the ematch setting for scoring inference_.

```bash  
cat <fasta> | ./fasta_nt_counts.py
```

ncrf_error_clumps-- Identify clumps of errors in Noise Cancelling Repeat Finder
alignments.
_This is no longer needed; error clumps are now detected and excised by NCRF_.

```bash  
cat <output_from_NCRF> | ./ncrf_error_clumps.py [options]
  --maxmratio=<ratio>     identify clumps with no more matches than given ratio; 
                          ratio can be between 0 and 1 (e.g. "0.85"), or can be
                          expressed as a percentage (e.g. "85%")
                          (default is 85%)
  --minnoise=<ratio>      (same as --maxmratio but with 1-ratio)
  --mincolumns=<columns>  identified clumps must have at least this many
                          alignment columns
                          (default is 10)
  --head=<number>         limit the number of input alignments
  --report:clumps         report clumps to stderr
                          (by default clumps are just marked on alignments)
```

fasta_length_distribution-- Read a fasta file and report the distribution of
sequence lengths.

```bash  
cat <fasta_file> | ./fasta_length_distribution.py [options]
  --progress=<number>   periodically report how many sequences we've read

The resulting file has two columns -- a sequence length and the number of times
that length was observed. The lines are sorted by increasing length.
```

sum_length_distributions-- Compute the sum of several sequence length
distribution files.

```bash  
cat <length_distribution_files> | ./sum_length_distributions.py [options]
  --progress=<number>   periodically report how many sequences we've read

  --report:totals   report total bp, number of sequences, and averge bp per
                    sequence

Input files are the same as the output of fasta_length_distribution -- two
columns, a sequence length and the number of times that length was observed.
The lines are sorted by increasing length.

Output is the same format.
```

bin_position_counts-- Accumulate per-position counts into binned intervals.

```bash  
cat <counts_file> | ./bin_position_counts.py [options]
  bin=<function>     (required) function that maps a position to its bin;
                     positions that map to the same unit interval are in the
                     same bin; an example is "pos/10"
  pos=<function>     (required) function that maps a bin to its position; this
                     should be the mathematical inverse of the to-bin function;
                     an example is "10*bin"
  --minpos=<number>  positions lower than this are ignored
  --maxpos=<number>  positions higher than this are ignored

_For further detail, see this program's usage report._
```

common_length_distribution-- Given several length distributions, report the
maximum common distibution.

```bash  
./common_length_distribution.py [options]
  <distribution_file>  (cumulative, at least 2 required) length distribution
                       file for one of the input components
  --scale=<number>     factor to multiply the maximum distribution by
                       (0 < number <= 1; default is 1)

The length distribution files consist of two columns -- length, count -- a
length or interval and the number of sequences of that length in that component.

The maximum common distibution is the minimum, over each length, of the
component counts for that length (example below). We call this "maximum"
because it give the largest number that can be sampled from each component to
achieve a common length distribution.

  component 1    component 2    maximum common
  375-378 800    375-378 600    375-378 600
  379-382 900    379-382 750    379-382 750
  383-386 1000   383-386 1040   383-386 1000
  387-391 1050   387-391 1200   387-391 1050
  392-395 1100   392-395 1000   392-395 1000
```

make_length_distribution_spec-- Given a length distribution and a target
distibution, create a spec to be used with fasta_match_length_distribution.

```bash  
./make_length_distribution_spec.py <target> <distribution>
  <target>             (required) length distribution file for the desired
                       fasta output of fasta_match_length_distribution
  <distribution>       (required) length distribution file corresponding to
                       a fasta file (or a set of fasta files) that will be
                       input to fasta_match_length_distribution

The length distribution files consist of two columns -- length, count -- a
length or interval and the number of sequences of that length in that
distribution.

The output is a distribution spec suitable for fasta_match_length_distribution.
```

spread_length_distribution-- Spread a given length distribution over several
component distributions.

```bash  
cat <distribution_spec> | ./spread_length_distribution.py [options]
  <components_file>    (required) file listing the components
  --input=<filespec>   (required) length distribution filename spec, e.g.
                       "{component}.dat"
  --output=<filespec>  (required) distribution spec filename spec, e.g.
                       "{component}.spec"
  --nowarn             don't report a warning for intervals in components that
                       aren't in the distribution spec; such intervals are not
                       processed (whether we report them or not)
  --seed=<string>      random number generator seed

_For further detail, see this program's usage report._
```

fasta_match_length_distribution-- Spread a given length distribution over several
component distributions.

```bash  
cat <fasta_file> | ./fasta_match_length_distribution.py [options]
  <components_file>    (required) file listing the components
  --input=<filespec>   (required) length distribution filename spec, e.g.
                       "{component}.dat"
  --output=<filespec>  (required) distribution spec filename spec, e.g.
                       "{component}.spec"
  --nowarn             don't report a warning for intervals in components that
                       aren't in the distribution spec; such intervals are not
                       processed (whether we report them or not)
  --seed=<string>      random number generator seed

usage: cat <fasta_file> | fasta_match_length_distribution [options]
  <distribution_spec>   (required) file describing the length distribution
  --remainder=<file>    write unfulfilled length distribution to a file
  --wrap=<length>       number of nucleotides per line in output fasta
  --seed=<string>       random number generator seed
  --progress=<number>   periodically report how many sequences we've read

The distribution spec file consists of three columns-- a length (or interval),
the number of sequences of that length we should output, and the number of
sequences of that length we expect to see in input.  An interval is two numbers
connected by a hyphen, e.g. "100-199" indicates sequences at least 100 bp long
but shorter than 200 bp.

The distribution of input lengths allows us to easily select any read with the
correct probability, without having to pre-scan the input to collect that
information. It also allows us to operate over a collection of files, using the
remainder file to communicate from the processing of one file to the next.
```

ncrf_parse.py, echydna.py, interval_dict.py-- _These support the other scripts
and should not be used directly_.


### R plotting functions 

plot_ncrf_event_matrix-- Plot error events vs alignment length, and report the
error event ratios.

in shell:
```bash  
  cat malus.reads.fa \
    | ./NCRF GGAAT \
    | tee malus.unfiltered.ncrf \
    | ./ncrf_extract_event_matrix.py --withheader \
    > malus.unfiltered.events.dat
```

in R:
```bash  
  source("ncrf_plotters.r")
  plot_ncrf_event_matrix("malus reads","malus.unfiltered.events.dat")
```

plot_ncrf_filter_results-- Plot errors vs matches by filtering class.

in shell:
```bash  
  cat malus.reads.fa \
    | ./NCRF GGAAT --positionalevents \
    | tee malus.ncrf \
    | ./error_nonuniformity_filter.py \
    > malus.ncrf
  cat malus.ncrf \
    | ./error_nonuniformity_filter.py --report:matrix \
    > malus.filter_results.dat
```

in R:
```bash  
  source("ncrf_plotters.r")
  plot_ncrf_filter_results("malus reads","malus.filter_results.dat",5)
```

read_ncrf_summary-- Read an ncrf summary.

in shell:
```bash  
  cat malus.reads.fa \
    | ./NCRF GGAAT --positionalevents \
    | ./error_nonuniformity_filter.py \
    | ./ncrf_summary.py \
    > malus.summary
```

in R:
```bash  
  source("ncrf_plotters.r")
  summary = read_ncrf_summary("malus.summary")
  ord = order(summary$querybp)
  plot(summary$querybp[ord],1:nrow(summary),log="x",
       main="GGAAT in malus reads",
       xlab="repeat length (log scale)",
       ylab="repeat instance (sorted by length)")
```

### Contact
For questions regarding usage, please contact Bob Harris <rsharris@bx.psu.edu>. 

### References
As of Oct/2018 this is unpublished work.
