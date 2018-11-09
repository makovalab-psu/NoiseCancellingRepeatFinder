# Stuff for sampling common read length distributions

_add general description here_

#### Scripts to sample common read length distributions

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
cat <fasta_file> | ./fasta_match_length_distribution.py [options]
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

#### Other scripts

interval_dict-- _This support the other scripts and should not be used
directly_.
