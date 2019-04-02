# Scripts supporting tests described in the manuscript

#### Scripts to create simulated genomes and reads

random_normal-- Generate random numbers, sampled from a normal distribution.

```bash  
./random_normal.py <num_values> [options]
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
./mock_motif_genome.py <motif> [options]
  --arrays=<filename>      specific arrays to embed; each line is of the form
                           <length> <motif>[<strand>]; <length> is in bp;
                           <strand> is ignored;
                           this cannot be used with any command line <motif>s,
                           nor with --repeats, --lengths, --motif:neighbor, or
                           --motif:mixture
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
./reconstruct_simulated_alignments.py [options]
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
                           (by default, we don't filter by length)
  --noisygenome            tolerate *some* noise in the genome's repeat arrays
                           (by default, we expect pristine repeat arrays in the
                           genome)
  --progress=<number>      periodically report how many reads we've processed

Given a genome and simulated reads sampled by ncrf_read_simulator, and the
corresponding cigars file, alignments are reconstructed. Note that this is
*not* an aliger; it is just reconstructing the alignment truth that the
ncrf_read_simulator created.

Note that by default we store the entire genome in memory. If the genome is
large, this could create memory issues. The --chromosomes option can be used to
process alignments on different chromosomes. Also, chromosomes not appearing in
the intervals file (if on is provided) are not stored.

Intervals, if provided, are one per line, <chrom> start> <end>. Coordinates are
zero-based and exclude the end position. Any additional columns are ignored.

Alignment output is in a format compatible with that produced by NCRF.
```

map_onto_simulated_reads-- Map intervals from a "genome" to positions on
simulated reads.

```bash  
cat <intervals_file> | ./map_onto_simulated_reads.py [options]
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

Given a genome from which simulated reads were sampled by ncrf_read_simulator,
and the corresponding cigars file, map intervals (or positions) from the genome
to the corresponding positions on the simulated reads. 

Intervals are one per line, <chrom> start> <end>. Coordinates are zero-based
and exclude the end position. Any additional columns are copied to the
output.
```

ncrf_read_simulator-- Convert a fasta file to sampled "reads" in fasta or fastq
format.

```bash  
[cat <fasta_file> |] ./ncrf_read_simulator.py [options]
  <num>x<length>[,<length>]  number and length of reads to generate; the second
                             length is used for paired reads when the mates
                             have different lengths; if <length> is "stdin",
                             lengths are read from stdin
  [<weight>:]<fasta_file>    (cumulative) sample reads from the specified
                             genome file; if no files are given, the genome is
                             read from stdin; <weights> can be used to
                             control the mixture of reads from different
                             genomes, and are internally scaled by the length
                             of the corresponding genome; if the weight is
                             absent, it is 1.0 by default
  --insert=<avg>,<stdev>     length of inserts for paired reads
                             (by default, reads are not paired)
  --orientation=<T2T|H2H|..> orientation for paired reads; the orientation
                             can be one of the following:
                               T2T, RF, PE (these are equivalent)
                               H2H, FR, MP (these are equivalent)
                               H2T, FF     (these are equivalent)
                               T2H, RR     (these are equivalent)
                             (default is tail-to-tail)
  --noise=<probability>      inject random sequencing errors (substitutions);
                             each base suffers a substitution error with the
                             given probability 
  --indel=<open>,<extend>    inject random sequencing errors (indels); <open>
                             is the probability of starting an indel at a
                             particular base; <extend> is the probability of
                             extended the gap after each base in the gap;
                             insertions and deletions are equally probable
  --errors=pacbio            simulate pacbio error profile
  --errors=nanopore          simulate nanopore error profile
  --errors=<spec>            simulate error profile with the given spec; <spec>
                             looks like this:
                                mm:1%,i:12%,d:2%
  --lower                    show substitutions in lowercase
  --prohibitn[=<length>]     don't sample reads containing a run of Ns; if
                             <length> is not given, no Ns are allowed
                             (by default, we do not check for Ns)
  --name=<template>          name for reads (see description below)
                             (default is FAKE_READ_[6], but with 6 replaced by
                             the smallest number sufficient)
  --width=<characters>       number of characters for each line of output dna;
                             this is only relevant for fasta output
                             (default is 100)
  --fastq                    output in fastq format
                             (default is fasta format)
  --intervals                output intervals, without sequence
  --output[+]=<filename>     write output to the specified file; if <filename>
                             contains {mate} or {zmate}, paired reads are
                             written to two files; if the plus sign is used,
                             we append to the file(s)
                             (default is stdout)
  --quality=<character>      set fastq qualities to a constant string of the
                             specified character (nothing better offered yet)
                             (default is J)
  --sam=<filename>           write a sam file equivalent of the generated reads
                             (currently not available for paired reads)
  --step=<number>            rather than sampling randomly, start at the
                             beginning and output a read starting at every Nth
                             position
  --strand=<+|->             only create reads from the specified strand
                             (by default we random choose strands)
  --cigars=<filename>        write cigar strings to the specified file
  --seed=<string>            set random seed
  --progress=<number>        periodically report how many reads we've generated

<template> is like this: BASE{4}_[6] where {4} is replaced by four random
letters/numbers and [6] is replaced by the read number (in six digits). [*]
is the read number in 'just enough' digits. Other recognized fields are
   {chrom}   the name of the sequence the read was drawn from
   {uchrom}  the name of the sequence the read was drawn from, in uppercase
   {start}   the starting position on that sequence (origin 1)
   {zstart}  the starting position on that sequence (origin 0)
   {end}     the ending position on that sequence
   {pstart}  the pair's starting position on that sequence (origin 1)
   {zpstart} the pair's starting position on that sequence (origin 0)
   {pend}    the pair's ending position on that sequence
   {strand}  the orientation on that sequence
   {mate}    for paired reads (1 or 2)
   {zmate}   for paired reads (0 or 1)
```

#### Scripts to evaluate classifier performance

observed_vs_truth-- Compare observed events from alignments to known truth.

```bash  
cat <alignment_summary> | ./observed_vs_truth.py <truth_catalog> [options]
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
cat <trf_html_output> | ./harvest_trf_html.py [options]
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
cat <output_from_minimap2> | ./minimap2_cs_to_events.py [options]
  --minquality=<qual>  discard low quality alignments
  --withheader         include a header line in the output
  --remove:cs          remove the cs tag
  --remove:tags        remove all tags

The minimap2 output should include the cs tag, i.e. minimap2 should have been
run with the "--cs=short" option.
```

#### Scripts to derive error models

sam_to_event_matrix-- Extract match, mismatch, insertion, and deletion event
counts from alignments in sam format.

```bash  
cat <sam_file> | ./sam_to_event_matrix.py [options]
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

truth_to_rates-- Read event counts from a truth catalog and report overall
induced event rates.

```bash  
cat <truth_catalog> | ./truth_to_rates

The truth catalog is usually the output from the --catalog option of
mock_motif_genome. It has 12 columns but only the events are used here ("m",
"mm", "i", and "d", columns 9, 10, 11, and 12).
```

#### Less frequently used scripts

event_matrix_to_rates-- Read match-and-error event counts and report overall
event rates.

```bash  
cat <event_counts_table> | ./event_matrix_to_rates

The input is usually the output from ncrf_extract_event_matrix.  It has 9
columns, first three of which are ignored here. The remaining columns are,
respectively, counts for matches ("m"), mismatches ("mm"), insertion opens
("io"), insertion extensions ("ix"), deletion opens ("do"), and deletion
extensions ("dx").
```

#### Other scripts

echydna and prob_table-- _These support the other scripts and should not be
used directly_.

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
    | ./ncrf_consensus_filter.py \
    > malus.ncrf
  cat malus.ncrf \
    | ./ncrf_consensus_filter.py --report:matrix \
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
    | ./ncrf_consensus_filter.py \
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
