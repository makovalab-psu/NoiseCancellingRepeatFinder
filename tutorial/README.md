### Brief Tutorial, Finding Repeats in Noisy Sequence Data

This directory contains a toy example in a small fasta file, example.fa, and
expected results files -- expected.unfiltered.ncrf, expected.filtered.ncrf, and
expected.summary.

#### (1) Run Noise Cancelling Repeat Finder to align sequence to repeats.

_Nota Bene: The <--scoring=1.00.XX> option shown below is necessary here to match
the expected tutorial outputs. You would not normally use that option. And, if
you are using a version prior to 1.01.00, this option is not recognized. See
issue #4 for additional details._

```bash 
    cat example.fa \
      | NCRF GGAAT --scoring=1.00.XX \
          --minlength=500 --maxnoise=20% --stats=events \
          --positionalevents \
      | ncrf_sort --sortby=mratio \
      > example.unfiltered.ncrf
```

The output file should match expected.unfiltered.ncrf.

The output consists of a series of alignment blocks, each showing an alignment
between a subsequence of the input and a repeated version of a motif. A typical
alignment block is shown here:

```bash 
# score=1193 querybp=771 mRatio=91.7% m=744 mm=9 i=40 d=18 ==========x=======...
NCRF_EXAMPLE                       50000 793bp 31854-32647 GAATGGAATGCAATGGAA...
GGAAT+                                   771bp score=1193  GAATGGAATGgAATGGAA...
# position 0 [G] mRatio=88.8% m=143 mm=2 i=7 d=9 mmA=1 mmC=1 mmG=0 mmT=0 x=18
# position 1 [G] mRatio=91.5% m=151 mm=4 i=10 d=0 mmA=2 mmC=1 mmG=0 mmT=1 x=14
# position 2 [A] mRatio=93.8% m=150 mm=1 i=6 d=3 mmA=0 mmC=0 mmG=1 mmT=0 x=10
# position 3 [A] mRatio=93.3% m=154 mm=0 i=11 d=0 mmA=0 mmC=0 mmG=0 mmT=0 x=11
# position 4 [T] mRatio=91.3% m=146 mm=2 i=6 d=6 mmA=0 mmC=1 mmG=1 mmT=0 x=14
```

The first line gives the alignment score, the number of bases from the repeated
motif sequence (the "query") that were aligned, the "match ratio" (defined
below), the number of matched bases ("m"), mismatched bases ("mm"), bases
inserted into the repeated motif ("i"), and bases deleted from the repeated
motif ("d"). The last part of the line shows where there are errors in the
alignment (mismatches or indels), an "=" indicating a match and an "x"
indicating an indel. In the example, we have an alignment of a 771 bp query
(the repeated motif) with 9 mismatches, 40 insertions and 18 deletions. Note
that the actual line is much longer than shown here (the same is true for the
second and third lines).

Note: we regard indels as having been introduced as sequencing errors, and the
repeated motif as having existed, perfectly repeated, in the molecule that was
sequenced. Thus insertions are bases in the sequence but not in the repeated
motif, and deletions are bases in the repeated motif but not in the sequence.

The second line shows the name of the aligned sequence, its full length, the
length of the aligned segment, the coordinates of the aligned segment, and the
segment's nucleotides. Any deletions are shown as dashes. In this example we
aligned a 793 bp segment of the 50,000 bp sequence, from position 31,854 up to,
but not including, position 32,647. Note that the coordinates are zero-based
and exclude the end position.

The third line shows the repeat motif, the orientation in which it was aligned
(+ or -), the length of the aligned repeat sequence (same as "querybp" as from
line 1), the alignment score (same as from line 1), and the aligned repeat
sequence nucleotides. Any insertions are shown as dashes, and any mismatches
are shown in lower case. In the example, the motif GGAAT aligned in forward (+)
orientation. The motif was repeated enough times to create a 771 bp sequence.
Positions 0 through 9 match the aligned segment; position 10 has a G where the
segment has a C.

Match ratio is the percentage of alignment columns that are matches. 100%
indicates a perfect alignment and lower numbers indicate more errors. The
formula is mRatio = m/(m+mm+i+d).

After the first three lines, alignment statistics are broken down by position
within the repeat motif. mRatio, m, mm, i, and d have the same meaning as for
the first line. mmA, mmC, mmG, and mmT are mismatch counts for each substituted
nucleotide. x is the sum of mm, i, and d. In the example above we see the
errors are a little less prevalent in position 0 (lower mRatio, higher x),
which might suggest that another motif is be a better match for this query
segment (or for parts of it). Reporting of these statistics can be disabled by
removing "--positionalevents" from the NCRF command line.

#### (2) False alignments to similar motifs.

One of the difficulties that arises, because we allow so much noise in the
alignments, is a motif will align to segments that are better matches for
other, similar, motifs. This is demonstrated by the following example, which
begins at line 145.

```bash 
# score=789 querybp=801 mRatio=84.7% m=698 mm=19 i=23 d=84 ========x====x==x=...
NCRF_EXAMPLE                       50000 740bp 37074-37814 TGGAATGG-ATGGAAAAG...
GGAAT+                                   801bp score=789   TGGAATGGAATGG-AAtG...
# position 0 [G] mRatio=89.7% m=148 mm=3 i=5 d=9 mmA=0 mmC=2 mmG=0 mmT=1 x=17
# position 1 [G] mRatio=94.0% m=157 mm=3 i=7 d=0 mmA=0 mmC=1 mmG=0 mmT=2 x=10
# position 2 [A] mRatio=51.5% m=85 mm=4 i=5 d=71 mmA=0 mmC=1 mmG=3 mmT=0 x=80
# position 3 [A] mRatio=92.1% m=152 mm=8 i=5 d=0 mmA=0 mmC=3 mmG=1 mmT=4 x=13
# position 4 [T] mRatio=96.3% m=156 mm=1 i=1 d=4 mmA=1 mmC=0 mmG=0 mmT=0 x=6
```

We have an alignment of GGAAT to a segment that might actually a better match
for GGAT. This is evident from the fact that the alignment errors are highly
skewed toward motif position 2 (mRatio=51.5%), where a large number (71) of
deletions were observed. In fact, this segment probably contains a mix of those
two motifs (GGAAT and GGAT).

#### (3) Filtering for consensus.

ncrf_consensus_filter.py is a post processor that automatically discards
alignments like the one in the previous example. It derives a consensus
from the segments of the sequence that aligned to the motif. If the consensus
doesn't match the motif, the alignment is discarded.

Run Noise Cancelling Repeat Finder again, passing the output through the
consensus filter.

```bash 
    cat example.fa \
      | NCRF GGAAT \
          --minlength=500 --maxnoise=20% \
          --stats=events --positionalevents \
      | ncrf_consensus_filter \
      | ncrf_sort --sortby=mratio \
      > example.filtered.ncrf
```

The output file should match expected.filtered.ncrf. Comparing to the earlier
unfiltered output, the suspect alignment is not included in the filtered output.

#### (4) Alignment summary.

Pass the output through ncrf_summary to get a tabular listing of the aligned
segments.

```bash 
    ncrf_cat example.filtered.ncrf \
      | ncrf_summary \
      > example.summary
```

The tab-delimited output should look like this:

```bash 
#line motif seq          start end   strand seqLen querybp mRatio m    mm i  d
1     GGAAT NCRF_EXAMPLE 31854 32647 +      50000  771     91.7%  744  9  40 18
10    GGAAT NCRF_EXAMPLE 14029 15525 +      50000  1453    91.4%  1400 17 79 36
19    GGAAT NCRF_EXAMPLE 24105 24972 +      50000  852     90.9%  814  10 43 28
28    GGAAT NCRF_EXAMPLE 33044 33866 +      50000  796     90.8%  763  15 44 18
37    GGAAT NCRF_EXAMPLE 45378 46122 -      50000  723     90.5%  692  10 42 21
46    GGAAT NCRF_EXAMPLE 40889 41729 -      50000  825     90.4%  785  12 43 28
55    GGAAT NCRF_EXAMPLE 47778 48345 -      50000  563     90.2%  532  8  27 23
64    GGAAT NCRF_EXAMPLE 661   1703  +      50000  1011    90.1%  965  17 60 29
73    GGAAT NCRF_EXAMPLE 46204 46731 +      50000  518     90.1%  493  5  29 20
82    GGAAT NCRF_EXAMPLE 20556 21422 -      50000  846     90.0%  805  13 48 28
91    GGAAT NCRF_EXAMPLE 48657 49997 -      50000  1310    89.9%  1243 24 73 43
100   GGAAT NCRF_EXAMPLE 41829 42471 +      50000  627     89.5%  594  11 37 22
109   GGAAT NCRF_EXAMPLE 42762 43707 +      50000  918     89.5%  874  12 59 32
118   GGAAT NCRF_EXAMPLE 22894 24093 -      50000  1148    89.2%  1098 18 83 32
127   GGAAT NCRF_EXAMPLE 16059 17456 -      50000  1335    89.1%  1277 21 99 37
```

The first column is the line number of the alignment in the alignment file. The
second column is the motif, e.g. GGAAT. Columns 3 through 5, "seq", "start",
and "end", report the sequence and position of the aligned segment (coordinates
are zero-based and exclude the end position). Column 6, "strand", shows the
orientation in which the motif was aligned. Column 7, "seqLen", is the full
length of the sequence containing the alignment. Column 8, "querybp", is the
number of motif bases in the alignment. The remaining five columns, "mRatio",
"m", "mm", "i", and "d", are as defined earlier.
