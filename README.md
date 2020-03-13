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
1. Download the latest stable version of Noise Cancelling Repeat Finder using
Github (as of this writing, the latest stable version is v1.01.00):
```bash  
     git clone --branch v1.01.00 https://github.com/makovalab-psu/NoiseCancellingRepeatFinder.git  
```  
Or, you can download the latest release from the releases page:
```bash  
    https://github.com/makovalab-psu/NoiseCancellingRepeatFinder/releases  
```
If you are adventurous and want to just try the latest stuff checked into the
repository, then you probably know how to clone the repo (so I don't need to
show you the command). Beware that it might not be stable.

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

A more detailed example is shown in tutorial/README.md.  A description of the
output format is included there.  The experiments subdirectory contains
additional examples, provided "as is" without explanation.

The output is usually passed through a series of the ncrf_* post-processing
scripts (e.g. ncrf_consensus_filter, ncrf_sort, or ncrf_summary).

### Intended Use Case

NCRF was designed with the idea that it would be used to search noisy reads,
of the type produced by PacBio and Oxford Nanopore. The implementation was
shaped by the expectation that sequences would not be longer than ≈500kbp and
motifs were not longer than ≈200bp. So while it may work for bigger problems
than that, there also may be issues, primarily relating to memory consumption
and speed.

Similarly, we haven't explored it's use in the absence of sequencing noise.
While we expect it would be useful for searching for repeated motifs in an
assembled genome, we have not investigated that use case.

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

And `./NCRF --help=scoring` will produce this list of scoring options

```bash  
  Default scoring is    M=1 MM=5 IO=5 IX=5 DO=5 DX=5
  --scoring=pacbio      use scoring for pacbio reads
                        (M=10 MM=35 IO=33 IX=21 DO=6 DX=28)
  --scoring=nanopore    use scoring for nanopore reads
                        (M=10 MM=63 IO=51 IX=98 DO=27 DX=34)
  --scoring=simple:<M>/<E> simple scoring matrix with match reward <M> and all
                        penalties <E>
  --match=<reward>      (M)  reward when sequence matches motif
  --mismatch=<penalty>  (MM) penalty when sequence mismatches motif
  --iopen=<penalty>     (IO) first penalty when sequence has nt, motif doesn't
  --iextend=<penalty>   (IX) other penalty when sequence has nt, motif doesn't
  --dopen=<penalty>     (DO) first penalty when motif has nt, sequence doesn't
  --dextend=<penalty>   (DX) other penalty when motif has nt, sequence doesn't
  --report:scoring      report scoring parameters and quit

  Mathematically, scaling all scores, including --minscore, by a constant
  factor will produce equivalent alignment results. However, larger M shortens
  the length susceptible to overflow. Overflow may occur if an alignment is
  longer than 2^31/M. For M=100 this is about 20 million. An M larger than 100
  is unlikely to be necessary. The same statements are true for MM, IO, IX, DO,
  and DX
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
usage: ncrf_cat <output_from_NCRF> | ncrf_consensus_filter [options]
  --consensusonly     just report the consensus motif(s) for each alignment,
                      instead of filtering; these are added to the alignment
                      file with a "# consensus" tag; note that the reported
                      consensus will be canonical, the lexigographical minimum
                      of all rotations including reverse complement
  [<name>:]<motif>    dna repeat motif to process; if no motifs are specified,
                      we process all of them (however, see note below)
                      (more than one motif can be specified)
  --head=<number>     limit the number of input alignments
  --progress=<number> periodically report how many alignments we've tested

Any motif that was given a name during the alignment process has to be
specified here, and with the same name. A motif was 'named' if an option of the
form <name>:<motif> was given to NCRF. The nt sequence for named motifs does
not appear in the alignment file produced by NCRF, but this program needs that
sequence.
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

The name template either names a single file or a collection of files. See
below for some examples.

The input alignment summaries are usually the output from ncrf_summary. Any
input file may contain alignments for more than one motif.

A typical input file is shown below. However, we do not interpret any columns
other than motif, seq, start, and end. This allows, for example, the output
from ncrf_summary_with_consensus.

  #line motif seq        start end  strand seqLen querybp mRatio m    mm  i  d
  1     GGAAT FAB41174_6 1568  3021 -      3352   1461    82.6%  1242 169 42 50
  11    GGAAT FAB41174_2 3908  5077 -      7347   1189    82.4%  1009 125 35 55
  21    GGAAT FAB41174_0 2312  3334 -      4223   1060    81.1%  881  115 26 64
   ...

If the output name template includes the substring "{motif}", this substring is
replaced by a motif name and any un-overlapped alignments to that motif are
written to that file. If the template name doesn't include "{motif}", all
un-overlapped alignments and overlapping groups are written to one file.

Overlapping groups are either written to the console (if no name template is
given), to the same file with alignments (if the name template doesn't contain
"{motif}"), or to a a file separate from the alignments (with "{motif}"
replaced by "overlaps").

This is summarized in the table below. We assume for this that the input only
contains two motifs, GGAAT and CATATA.

  name_template    | output
  -----------------+----------------------------------------------------------
  (none)           | un-overlapped and overlap groups written to console
  -----------------+----------------------------------------------------------
  filename         | un-overlapped and overlap groups written to filename
  -----------------+----------------------------------------------------------
  filename.{motif} | un-overlapped GGAAT written to filename.GGAAT
                   | contains un-overlapped CATATA written to filename.CATATA
                   | overlap groups written to filename.overlaps
  -----------------+----------------------------------------------------------

Overlap groups are separated by a single blank line, as shown below (note that
this is a contrived example). When un-overlapped alignments and overlapped ones
are in the same file, the un-overlapped ones are first, with a blank line
separating each alignment.

  #line motif  seq        start end  strand seqLen querybp mRatio m    mm  i  d
  1     GGAAT  FAB41174_6 1568  3021 -      3352   1461    82.6%  1242 169 42 50
  1     CATATA FAB41174_6 1621  2607 -      3352    ...
  (blank line)
  21    GGAAT  FAB41174_0 2312  3334 -      4223   1060    81.1%  881  115 26 64
  41    CATATA FAB41174_0 3276  4098 -      4223    ...
  31    GGAAT  FAB41174_0 3966  4271 +      4223    ...
  (blank line)
   ... (more groups)"""
```

#### Other scripts

ncrf_parse-- _This supports the other scripts and should not be used directly_.

### Contact
For questions regarding usage, please contact Bob Harris <rsharris@bx.psu.edu>. 

### References
As of Oct/2018 this is unpublished work.
