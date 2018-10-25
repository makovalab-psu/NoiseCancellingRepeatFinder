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
cat reads.fa | NCRF GGCAT > example.ncrf
```

A more detailed example is shown in example/README.md.  A description of the
output format is included there.  The experiments subdirectory contains
additional examples, provided "as is" without explanation.

The output is usually passed through a series of the ncrf_* post-processing
scripts (e.g. error_nonuniformity_filter, ncrf_sort, or ncrf_summary).


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

### Additional scripts

#### Primary scripts

ncrf_cat-- Concatenate several output files from Noise Cancelling Repeat Finder.

To feed to output from several runs of NCRF into a script, ncrf_cat is
necessary (as opposed to the usual shell command "cat").

```bash  
ncrf_cat <file1> [<file2> ...]
```

error_nonuniformity_filter-- Filter the output of Noise Cancelling Repeat Finder,
removing alignments in which matches (or errors) aren't uniformly distributed
across the motif positions.

```bash  
ncrf_cat <output_from_NCRF> | error_nonuniformity_filter [options]
  --method=min-max      judge alignments by a "min-max" test  
                        (this is the default)
  --trials=<number>     number of trials for the min-max test
                        (default is 10K)
  --trials=<number>/<number> numbers of successes needed and rials for the
                        min-max test; e.g. "5/10K" will do 10 thousand trials
                        and require at least five success to "pass"
  --discard:good        discard the "good" alignments instead of the "bad" ones
                        (by default we discard the "bad" alignments) 
  --discard:none        don't discard any alignments, just report the test
                        outcomes (see below for how they are reported) 
  --test:matches-insertions perform the test using match counts less insertions
                        (this is the default)
  --test:matches        perform the test using match counts
  --test:errors         perform the test using error counts; this may increase
                        the likelihood that the test cannot be performed for
                        some alignments, since the counts may be too low
  --warn:untested       report each untestable alignment as a warning
  --seed=<string>       set random seed; this is only applicable to the min-max
                        test; two special cases are "default" (use the built-in
                        default seed) and "none" (don't seed the random number
                        generator)
  --subsample=<k>/<n>   Conceptually split the alignments by into <n> groups,
                        and only process the <k>th group; <k> ranges from 1 to
                        <n>
  --head=<number>       limit the number of input alignments
  --progress=<number>   periodically report how many alignments we've tested

In a "true" alignment to a given motif unit, we expect the errors to be
distributed randomly and uniformly among the positions in the unit. (That is
an underlying assumption, possibly not true.)  This program discards alignments
that fail a test based on that assumption.

Error counts are often too small for the statistical test, so a test based on
match counts is used instead. However, an insertion error does not reduce the
match count at any position, so by default we decrease matches by the number of
insertions at that position.

The input alignments must include position event information. This can be
accomplished by using the --positionalevents option of Noise Cancelling Repeat
Finder.

When test outcomes are reported for --discard:none, one of the following lines
is added to each alignment
# positional min-max: match-insert uniformity rejected
# positional min-max: match-insert uniformity not rejected
# positional min-max: untested
The "untested" case indicates that the min-max test could not be performed,
usually because one of the positional match counts is too small.

The user should be aware that the results aren't necessarily determinstic.
When a PRNG is in play (as for min-max), the result for an alignment depends
on the state of the PRNG; and the state of the PRNG depends on all the
alignments that preceded the one being tested.
```

ncrf_sort-- Sort the alignments output by Noise Cancelling Repeat Finder.

```bash  
ncrf_cat <output_from_NCRF> | ncrf_sort [options]
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
ncrf_cat <output_from_NCRF> | ncrf_summary [options]
  --minmratio=<ratio>  discard alignments with a low frequency of matches;
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
ncrf_cat <output_from_NCRF> | ncrf_to_bed [options]
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

#### Less frequently used scripts

ncrf_error_clumps-- Identify clumps of errors in Noise Cancelling Repeat Finder
alignments.
_This is no longer needed; error clumps are now detected and excised by NCRF_.

```bash  
ncrf_cat <output_from_NCRF> | ncrf_error_clumps [options]
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

ncrf_error_positions-- Convert the output of Noise Cancelling Repeat Finder to
a position-of-errors summary, a tab-delimited table with one line of error info
per alignment, primarily to show the relative positions of errors along the
alignment.

```bash  
ncrf_cat <output_from_NCRF> | ncrf_error_positions [options]
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

ncrf_extract_event_matrix-- Extract the (NON-positional) match-and-error event
counts matrix from Noise Cancelling Repeat Finder alignments.

```bash  
ncrf_cat <output_from_NCRF> | ncrf_extract_event_matrix [options]
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

ncrf_extract_mx_matrix-- Extract the positional match-and-error-counts matrix
from Noise Cancelling Repeat Finder alignments.

```bash  
ncrf_cat <output_from_NCRF> | ncrf_extract_mx_matrix [options]
  --head=<number>         limit the number of input alignments

The output matrix has R rows and 2M+1 columns, where R is the number of input
alignments and M is the length of the aligned motif. (It is assumed that all
alignments are to the same motif, but this is NOT validated).

The first column is the line number of the alignment in the input file.  The
next M columns are the positional counts for matches ("m"), and the final M
columns are the positional counts for errors ("x").  The output is intended
to be suitable as input to R.
```

ncrf_msa-- Create an msa (multiple sequence alignment) of each word aligned to
the motif, in the output of Noise Cancelling Repeat Finder.

```bash  
ncrf_cat <output_from_NCRF> | ncrf_msa [options]
  --head=<number>         limit the number of input alignments
```

ncrf_resolve_overlaps-- Resolve overlapping alignments of different motifs.

```bash  
ncrf_resolve_overlaps <alignment_summary..> [options]
  <alignment_summary>    (cumulative) file(s) containing aligment summaries
                         for which overlaps are to be resolved
  --head=<number>        limit the number of input aligment summaries
  --out=<name_template>  file to write overlap groups to; see discussion of
                         name template below; if this option is absent, all
                         output is written to the console

The name template either names a single file or a collection of files.  If
if includes the substring "{motif}", this substring is replaced by a motif name
and any un-overlapped alignments to that motif are written to that file. If
the template name doesn't include "{motif}", all un-overlapped alignments and
overlapping groups are written to one file.

Overlapping groups are either written to the console (if no name template is
given), to the same file with alignments (if the name template doesn't contain
"{motif}"), or to a a file separate from the alignments (with "{motif}"
replaced by "overlaps").

The alignment summaries are usually the output from ncrf_summary. Any file may
contain alignments for more than one motif.

Typical input file:
  #line motif seq        start end  strand seqLen querybp mRatio m    mm  i  d
  1     GGAAT FAB41174_6 1568  3021 -      3352   1461    82.6%  1242 169 42 50
  11    GGAAT FAB41174_2 3908  5077 -      7347   1189    82.4%  1009 125 35 55
  21    GGAAT FAB41174_0 2312  3334 -      4223   1060    81.1%  881  115 26 64
   ...
```

ncrf_words-- Look for prominent "wrong" motifs (words) in the output of Noise
Cancelling Repeat Finder.

```bash  
ncrf_cat <output_from_NCRF> | ncrf_words [options]
  --minwordratio=<r>  only show words that have counts that are at least r
                      times the motif word's count (e.g. r=0.5 would show the
                      words that occur at least half as often as the motif)
                      (default is 1.0)
  --head=<number>     limit the number of input alignments
```

event_matrix_to_rates-- Read match-and-error event counts and report overall
event rates.

```bash  
cat <event_counts_table> | event_matrix_to_rates

The input is usually the output from ncrf_extract_event_matrix.  It has 9
columns, first three of which are ignored here. The remaining columns are,
respectively, counts for matches ("m"), mismatches ("mm"), insertion opens
("io"), insertion extensions ("ix"), deletion opens ("do"), and deletion
extensions ("dx").
```

#### Scripts to sample common read length distributions

fasta_length_distribution-- Read a fasta file and report the distribution of
sequence lengths.

```bash  
cat <fasta_file> | fasta_length_distribution [options]
  --progress=<number>   periodically report how many sequences we've read

The resulting file has two columns -- a sequence length and the number of times
that length was observed. The lines are sorted by increasing length.
```

sum_length_distributions-- Compute the sum of several sequence length
distribution files.

```bash  
cat <length_distribution_files> | sum_length_distributions [options]
  --report:totals   report total bp, number of sequences, and averge bp per
                    sequence

Input files are the same as the output of fasta_length_distribution -- two
columns, a sequence length and the number of times that length was observed.
The lines are sorted by increasing length.

Output is the same format.
```

bin_position_counts-- Accumulate per-position counts into binned intervals.

```bash  
cat <counts_file> | bin_position_counts [options]
    bin=<function>     (required) function that maps a position to its bin;
                       positions that map to the same unit interval are in the
                       same bin; an example is "pos/10"
    pos=<function>     (required) function that maps a bin to its position; this
                       should be the mathematical inverse of the to-bin function;
                       an example is "10*bin"
    --minpos=<number>  positions lower than this are ignored
    --maxpos=<number>  positions higher than this are ignored
  
The counts file consists of two columns-- a position and a count for that
position.  For example, this set of position,count pairs
   375 15    381 26    387 28    393 21    399 23
   376 17    382 27    388 23    394 20    400 27
   377 21    383 26    389 18    395 21    401 24
   378 26    384 23    390 27    396 22    402 24
   379 21    385 14    391 30    397 20    403 28
   380 21    386 19    392 21    398 25    404 30
would be represented in 30 lines in the file:
   375 15
    ...
   404 30

With the binning function bin=pos/10 and its inverse pos=10*bin, these counts
are collected into bins like this. Note that the bins may cover some positions
not present in the input (but no empty bins will be output).
   370-379 100
   380-389 225
   390-399 230
   400-409 133

The reason that 370-379 forms a bin is that the function bin=pos/10 maps each
of them into the unit interval 37 <= bin < 38
   370  37.0
   371  37.1
    ...
   378  37.8
   379  37.9

The bin mapping function should be strictly increasing over the range of
positions presented. It can be more sophisticated than demonstrated above. For
example, bin=30.5*log(pos-250) and inverse pos=250+exp(bin/30.5) will group
small values (pos about 400) into bins of size 5 and large values (pos about
300,000) into much larger bins.
```

common_length_distribution-- Given several length distributions, report the
maximum common distibution.

```bash  
common_length_distribution [options]
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
make_length_distribution_spec <target> <distribution>
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
  <components_file>    (required) file listing the components
  --input=<filespec>   (required) length distribution filename spec, e.g.
                       "{component}.dat"
  --output=<filespec>  (required) distribution spec filename spec, e.g.
                       "{component}.spec"
  --nowarn             don't report a warning for intervals in components that
                       aren't in the distribution spec; such intervals are not
                       processed (whether we report them or not)
  --seed=<string>      random number generator seed

The distribution spec file consists of three columns -- length, outCount,
inCount -- a length (or length interval), the number of sequences of that
length to be output, and the number of sequences of that length we expect would
be seen in input. This is the same format used as input to
fasta_match_length_distribution.

The components file lists component names, one per line. These are mapped to
length distribution files by the input filespec.  Each of the component length
distribution files consists of two columns -- length, inCount -- a length or
interval and the number of sequences of that length in that component.

It is assumed that the lengths in the distribution spec and all the components
are a perfect match, and that a length's inCount in the distribution spec
matches the sum of the length's inCount in all the components.  The goal of
this program is then to spread the outCount of the distribution spec to the
components, and thus create a distribution spec for each component.

For example, suppose we have five lengths/intervals and these three components
and distribution spec:

  component 1    component 2    component 3    distribution spec
  375-378 826    375-378 103    375-378 6804   375-378 4700 7733
  379-382 831    379-382 83     379-382 6893   379-382 2700 7807
  383-386 835    383-386 92     383-386 7043   383-386 4800 7970
  387-391 1026   387-391 131    387-391 8950   387-391 1900 10107
  392-395 820    392-395 90     392-395 7387   392-395 5100 8297

The goal, for the first row, is to spread outCount 4700 among the three
components, essentially randomly sampling from inCount 7733=826+103+6804 with
probabilities corresponding to each component. In the typical result shown
below, that split was 4700=509+66+4125.

  component 1        component 2      component 3
  375-378 509 826    375-378 66 103   375-378 4125 6804
  379-382 283 831    379-382 28 83    379-382 2389 6893
  383-386 507 835    383-386 53 92    383-386 4240 7043
  387-391 199 1026   387-391 25 131   387-391 1676 8950
  392-395 515 820    392-395 55 90    392-395 4530 7387

Note that we make no effort to maintain proportions, but the sampling process
will usually result in something close to proportionality.
```

fasta_match_length_distribution-- Spread a given length distribution over
several component distributions.

```bash  
cat <fasta_file> | fasta_match_length_distribution [options]
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

#### Scripts to derive error models

sam_to_event_matrix-- Extract match, mismatch, insertion, and deletion event
counts from alignments in sam format.

```bash  
cat <sam_file> | sam_to_event_matrix [options]
  --mapq=<num>            ignore alignments with mapping quality (the SAM MAPQ
                          field) below <num>
                          (by default, we accept any mapping quality)
  --withheader            include a header line in the output
  --sumonly               include only a summation line in the output
                          (by default, we output a separate line for each
                          alignment, and no sum)
  --warnandcontinue       warn when alignments violate sanity checks, and
                          discard those alignments
                          (by default, we report it as an error and halt)
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

#### Scripts to derive alignment scores

parameter_ball-- Vary a set of parameters, inside a "ball".

```bash  
parameter_ball <parameter> [options]
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

truth_to_rates-- Read event counts from a truth catalog and report overall
induced event rates.

```bash  
cat <truth_catalog> | truth_to_rates

The truth catalog is usually the output from the --catalog option of
mock_motif_genome. It has 12 columns but only the events are used here ("m",
"mm", "i", and "d", columns 9, 10, 11, and 12).
```

#### Scripts to create simulated genomes and reads

random_normal-- Generate random numbers, sampled from a normal distribution.

```bash  
random_normal <num_values> [options]
  <num_values>             number of values to generate
  --mu=<value>             mean value         (default is 0.0)
  --sigma=<value>          standard deviation (default is 1.0)
  --round                  round to integers
  --floor                  "round down" to integers
  --ceiling                "round up" to integers
  --precision=<digits>     set number of digits after decimal
  --seed=<string>          set random seed
```

mock_motif_genome-- Create a mock genome (or a "read") with embedded repeat motifs.

```bash  
mock_motif_genome <motif> [options]
  <motif>                  (cumulative) motif to embed
  --name=<string>          read name
  --length=<bp>       (L=) read length; if this is absent, the repeats will
                           comprise the entire read (no fill DNA)
  --length=<pctg>%    (L=) read length as a percentage of the total repeat
                           length; pctg should be more than 100.
  --length=+<bp>      (L=) read length as delta above the total repeat length
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
                           repeat length distribution from stdin; if <file>
                           contains "{motif}", we'll use a separate distribution
                           for each motif
  --minfill=<bp>      (F=) minimum fill (random sequence) between repeats
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

reconstruct_simulated_alignments-- Reconstruct alignments between a genome and
simulated reads sampled from the genome.

```bash  
reconstruct_simulated_alignments [options]
  --genome=<filename>      (mandatory) genome file, fasta or gzipped fasta (an
                           input file)
  --reads=<filename>       (mandatory) genome file, fasta or gzipped fasta (an
                           input file)
  --cigars=<filename>      (mandatory) cigar strings file corresponding to the
                           reads (an input file); this may contain cigars for
                           reads not present in the reads file (they are
                           ignored)
  --intervals=<filename>   If this is provided, alignments are truncated to
                           these intervals in the genome (format described
                           below)
  --catalog=<filename>     If this is provided, alignments are truncated to the
                           repeat intervals given by the file, and positional
                           match/mismatch/insert/delete counts are produced. 
                           The format of the input file is the same as produced
                           by mock_motif_genome's --catalog option; positional
                           counts are the same as would be produced by NCRF's
                           --positionalevents option
  --motif=<motif>          (cumulative) motifs of interest; alignments for other
                           motifs are discarded; requires --catalog, and the
                           motif must appear in the same orientation listed in
                           the catalog
                           (if this is not provided, we keep all alignments)
  --chromosome[s]=<names>  (cumulative) only reconstruct alignments on the
                           specified "chromosomes";  <names> is a comma-
                           separated list of sequence names in the genome
                           (default is to report intervals on all chromosomes)
  --minlength=<bp>         discard alignments that aren't long enough on the;
                           genome
                           (but default, we don't filter by length)
  --progress=<number>      periodically report how many reads we've processed

Given a genome and simulated reads sampled by simulate_reads_v4, and the
corresponding cigars file, alignments are reconstructed. Note that this is
*not* an aliger; it is just reconstructing the alignment truth that the
simulate_reads_v4 created.

Note that by default we store the entire genome in memory. If the genome is
large, this could create memory issues. The --chromosomes option can be used to
process alignments on different chromosomes. Also, chromosomes not appearing in
the intervals file (if on is provided) are not stored.

Intervals, if provided, are one per line, <chrom> start> <end>, origin-zero
half-open. Any additional columns are ignored.

Alignment output is in a format compatible with that produced by NCRF.
```

map_onto_simulated_reads-- Map intervals from a "genome" to positions on
simulated reads.

```bash  
cat <intervals_file> | map_onto_simulated_reads [options]
  --cigars=<filename>    (mandatory) cigar strings file (an input file)
  --stranded=<columns>   (cumulative) input columns which are presumed to have
                         strand info (+ or -) as their final character;
                         <columns> is a comma-separated list
  --truncate             truncate mappings at the end of reads; actaully
                         mappings are always truncated, but by default when
                         this happens it is indicated as "<0" or ">1000"
                         (assuming the read length is 1000); this option
                         just removes the "<" and ">" indicators.
  --sortby:reads         sort output by read positions on the genome
                         (by default, output is interval-by-interval in the
                         order intervals are read)
  --separators           print separating lines between different intervals
                         or reads

Given a genome from which simulated reads were sampled by simulate_reads_v4,
and the corresponding cigars file, map intervals (or positions) from the genome
to the corresponding positions on the simulated reads. 

Intervals are one per line, <chrom> start> <end>, origin-zero half-open. Any
additional columns are copied to the output.
```

#### Scripts to evaluate classifier performance

observed_vs_truth-- Compare observed events from alignments to known truth.

```bash  
cat <alignment_summary> | observed_vs_truth <truth_catalog> [options]
  <truth_catalog>        File containing aligment "truth"; the format of this
                         file depends on whether we have alignments to a genome
                         or alignments to reads (see below)
  --genome               alignments are to a genome
                         (this is the default)
  --reads                alignments are to reads
  --motif=<motif>        (cumulative) Motifs represented in the summary. Truth
                         intervals for other motifs are discarded. If this is
                         not provided, we use all truth intervals.
  --detection=<portion>  threshold for a truth interval to be considered as
                         "suffiently" covered; 0<portion<=1, but can be
                         expressed with a % sign
                         (default is 95%)
  --detail=<file>        Report separate true positive rates for each observed
                         interval
  --overcovered          Report over-covered rates too

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

The format of the truth catalog depends on whether we have alignments to a
genome or alignments to reads.

For alignments to a genome, the truth catalog is usually the output from the
--catalog option of mock_motif_genome. It has 12 columns but only the intervals
and motif are used here ("chrom", "start", "end", and "motif", columns 1, 2, 3,
and 4). Intervals must be distinct; overlaps are not allowed.

For alignments to reads, the truth catalog is usually the output from
map_onto_simulated_reads. It has 7 unlabled columns but only the intervals and
motif are used here ("readName", "start", "end", and "motif", columns 4, 5, 6,
and 7). Intervals must be distinct; overlaps are not allowed.
```

harvest_trf_html-- Convert alignments from TRF (Tandem Repeat Finder) html
output to a tabular text format similar to an NcRF summary.

```bash  
cat <trf_html_output> | harvest_trf_html [options]
  --motif=<motif>        (cumulative) motifs of interest; alignments for other
                         motifs are discarded
                         (if this is not provided, we keep all alignments)
  --minlength=<bp>       discard alignments that don't have long enough repeat
                         (but default, we don't filter by length)
  --withheader           include a header line in the output
  --withalignment        include alignment text in the output
```

minimap2_cs_to_events-- Convert the cs tag in minimap2 output to ncrf-style
event counts.

```bash  
cat <output_from_minimap2> | minimap2_cs_to_events [options]
  --minquality=<qual>  discard low quality alignments
  --withheader         include a header line in the output
  --remove:cs          remove the cs tag
  --remove:tags        remove all tags

The minimap2 output should include the cs tag, i.e. minimap2 should have been
run with the "--cs=short" option.
```

#### Other scripts

ncrf_parse, echydna, interval_dict-- _These support the other scripts and
should not be used directly_.

### R plotting functions 

plot_ncrf_event_matrix-- Plot error events vs alignment length, and report the
error event ratios.

in shell:
```bash  
  cat malus.reads.fa \
    | NCRF GGAAT \
    | tee malus.unfiltered.ncrf \
    | ncrf_extract_event_matrix --withheader \
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
    | NCRF GGAAT --positionalevents \
    | tee malus.ncrf \
    | error_nonuniformity_filter \
    > malus.ncrf
  cat malus.ncrf \
    | error_nonuniformity_filter --report:matrix \
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
    | NCRF GGAAT --positionalevents \
    | error_nonuniformity_filter \
    | ncrf_summary \
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
