# Noise Cancelling Repeat Finder -- Mock Genomes

Simulated pacbio and nanopore reads from fake genomes embedded with some tandem
repeats.

# Mock Genome 1

(to be described later)

# Mock Genome 2

(to be described later)

# Mock Genome 3

The genome is a mix of random sequence and period-5 and period-32 tandem
repeats. The primary repeat elements are have periods of 5, 10, 20, 40, 80,
and 171 bp, shown below. GGAAT is a highly studied 5-mer, and the 171-mer is
alpha satellite. The other four are randomly generated sequences.

Motifs of various sizes were embedded into a ≈2Mbp genome.

```bash  
motif5 GGAAT
motif10 CACTGCTGGT
motif20 GTTAGCGGTCGTGCTGATGG
motif40 GGCTCCTATCTCGCCTGTTCCCGGGTTCCTCTTATTCTCA
motif80 GAGATTGGAGTCCAAGAAATTCAGTCACCTTTCAGCGGTTCCAGTCACGGCGCTAAGTGCCTATTGACCCGCTACTGTTT
motif171 TCTGTCTAGTTTTTATATGAAGATATTCCCTTTTCCACCATAATCCTCAAAGCGCTCCAAATATCCACTTGCAGATTCTACAAAAAGAGTGTTTCCAAACTGCTCTATCAAAAGAAATGTTCAACTCTGTGAGTTGAATACACACATCACAAAGAAGTTTCTGAGAATGCT
```bash  


The genome is a mix of random sequence and period-5 and period-32 tandem
repeats. The primary repeat elements are GGAAT and
TAGTCGAGGAATTATGGATGGCGAAGATAAAT. Some elements are minor variants of those
elements, e.g. CGAAT or TAGTCGAGGAATTATGGATGGCGCAGATAAAT. All of the repeats
have a total length, in the genome, of at least 500 bp. FakeTRGenome.truth.dat
shows the placement of these repeats.

Fake reads were sampled from this genome, simulating a pacbio sequencing error
profile derived from alignments published in the genome-in-a-bottle project --
1.7% mismatch rate, 8.9% single-base insertion, 4.3% single-base deletion.

Overall coverage (of the fake genome) is intended to be about 20X. The
sampling positions within the genome are encoded in the read names. For
example, FAKE_READ_R_43796_47012 means the read was sampled from the genomic
interval 43796..47012 (origin-zero, half open); the "R" means it was reverse
complemented.

