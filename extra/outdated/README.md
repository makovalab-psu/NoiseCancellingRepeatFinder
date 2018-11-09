# Outdated stuff

error_nonuniformity_filter-- Filter the output of Noise Cancelling Repeat Finder,
removing alignments in which matches (or errors) aren't uniformly distributed
across the motif positions.

```bash  
./ncrf_cat.py <output_from_NCRF> | ./error_nonuniformity_filter.py [options]
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

ncrf_words-- Look for prominent "wrong" motifs (words) in the output of Noise
Cancelling Repeat Finder.

```bash  
./ncrf_cat.py <output_from_NCRF> | ./ncrf_words.py [options]
  --minwordratio=<r>  only show words that have counts that are at least r
                      times the motif word's count (e.g. r=0.5 would show the
                      words that occur at least half as often as the motif)
                      (default is 1.0)
  --head=<number>     limit the number of input alignments
```

ncrf_extract_mx_matrix-- Extract the positional match-and-error-counts matrix
from Noise Cancelling Repeat Finder alignments.

```bash  
./ncrf_cat.py <output_from_NCRF> | ./ncrf_extract_mx_matrix.py [options]
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
./ncrf_cat.py <output_from_NCRF> | ./ncrf_extract_event_matrix.py [options]
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

ncrf_error_positions-- Convert the output of Noise Cancelling Repeat Finder to
a position-of-errors summary, a tab-delimited table with one line of error info
per alignment, primarily to show the relative positions of errors along the
alignment.

```bash  
./ncrf_cat.py <output_from_NCRF> | ./ncrf_error_positions.py [options]
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

ncrf_error_clumps-- Identify clumps of errors in Noise Cancelling Repeat Finder
alignments.
_This is no longer needed; error clumps are now detected and excised by NCRF_.

```bash  
./ncrf_cat.py <output_from_NCRF> | ./ncrf_error_clumps.py [options]
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

