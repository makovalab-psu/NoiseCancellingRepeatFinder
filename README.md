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

Python is not used by the core program, but is needed for a post-processing
helper programs included in the package.

### Usage Overview

The simplest use, searching reads.fa for the repeated motif GGCAT:

```bash 
cat reads.fa | ./NCRF GGCAT > example.ncrf
```

A more detailed example is shown in example/README.md.  A description of the
output format is included there.  The experiments subdirectory contains
additional examples, provided "as is" without explanation.

The output is usually passed through a series of the ncrf* post-processing
scripts.


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

Descriptions to be added.

### Contact
For questions regarding usage, please contact Bob Harris <rsharris@bx.psu.edu>. 

### References
As of Oct/2018 this is unpublished work.
