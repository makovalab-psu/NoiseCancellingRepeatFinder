# Noise Cancelling Repeat Finder

This package finds alignments of short tandem repeats in noisy DNA sequencing
data.  Given a list of known motifs (e.g. "GGCAT", "CCAT", etc.) and DNA
sequences in fasta files (e.g. PacBio or Oxford Nanopore sequenced reads), it
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
3. You may wish to add the path to your copy of the NoiseCancellingRepeatFinder
directory to your PATH variable.

4. You may wish to create symbolic links to the python scripts, so that they
can be used on a command line by name (without the .py part).
```bash  
    cd NoiseCancellingRepeatFinder  
    ./make_symbolic_links.sh  
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
scripts (e.g. ncrf_consensus_filter, ncrf_sort, or ncrf_summary).


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
./ncrf_cat.py <file1> [<file2> ...] [--markend]
  <file1>    an output file from Noise Cancelling Repeat Finder
  <file2>    another output file from Noise Cancelling Repeat Finder
  --markend  assume end-of-file markers are absent in the input, and add an
             end-of-file marker to the output
             (by default we require inputs to have proper end-of-file markers)

Concatenate several output files from Noise Cancelling Repeat Finder.  This
is little more than copying the files and adding a blank line between the
files.

It can also be used to verify that the input files contain end-of-file markers
i.e. that they were not truncated when created.
```

ncrf_consensus_filter-- Filter Noise Cancelling Repeat Finder alignments, discarding alignments that
has a consensus different than the motif unit.

```bash  
./ncrf_cat.py <output_from_NCRF> | ./ncrf_consensus_filter.py [options]
  --consensusonly     just report the consensus motif(s) for each alignment,
                      instead of filtering; these are added to the alignment
                      file with a "# consensus" tag
  --head=<number>     limit the number of input alignments
```

ncrf_sort-- Sort the alignments output by Noise Cancelling Repeat Finder.

```bash  
./ncrf_cat.py <output_from_NCRF> | ./ncrf_sort.py [options]
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
./ncrf_cat.py <output_from_NCRF> | ./ncrf_summary.py [options]
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
./ncrf_cat.py <output_from_NCRF> | ./ncrf_to_bed.py [options]
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

ncrf_resolve_overlaps-- Resolve overlapping alignments of different motifs.

```bash  
./ncrf_resolve_overlaps.py <alignment_summary..> [options]
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

A typical input file is shown below. However, we do not interpret any columns
other than motif, seq, start, and end. This allows, for example, the output
from ncrf_summary_with_consensus.

  #line motif seq        start end  strand seqLen querybp mRatio m    mm  i  d
  1     GGAAT FAB41174_6 1568  3021 -      3352   1461    82.6%  1242 169 42 50
  11    GGAAT FAB41174_2 3908  5077 -      7347   1189    82.4%  1009 125 35 55
  21    GGAAT FAB41174_0 2312  3334 -      4223   1060    81.1%  881  115 26 64
   ...
```

#### Other scripts

ncrf_parse-- _This supports the other scripts and should not be used directly_.

### Contact
For questions regarding usage, please contact Bob Harris <rsharris@bx.psu.edu>. 

### References
As of Oct/2018 this is unpublished work.
