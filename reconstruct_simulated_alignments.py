#!/usr/bin/env python
"""
Reconstruct alignments between a genome and simulated reads sampled from the
genome.
"""

from sys        import argv,stdin,stdout,stderr,exit
from os         import path as os_path
from gzip       import open as gzip_open
from ncrf_parse import reverse_complement


def usage(s=None):
	message = """
usage: reconstruct_simulated_alignments [options]
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
                           below).
  --chromosome[s]=<names>  (cumulative) only reconstruct alignments on the
                           specified "chromosomes";  <names> is a comma-
                           separated list of sequence names in the genome
                           (default is to report intervals on all chromosomes)

Given a genome and simulated reads sampled by simulate_reads_v4, and the
corresponding cigars file, alignments are reconstructed.  Note that this is
*not* an aliger; it is just reconstructing the alignment truth that the
simulate_reads_v4 created.

Note that by default we store the entire genome in memory. If the genome is
large, this could create memory issues. The --chromosomes option can be used
to process alignments on different chromosomes.  Also, chromosomes not
appearing in the intervals file (if on is provided) are not stored.

Intervals, if provided, are one per line, <chrom> start> <end>, origin-zero
half-open. Any additional columns are ignored."""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global debug

	# parse the command line

	genomeFilename    = None
	readsFilename     = None
	cigarFilename     = None
	intervalsFilename = None
	chromsOfInterest  = None
	debug             = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--genome=")):
			genomeFilename = argVal
		elif (arg.startswith("--reads=")) or (arg.startswith("--read=")):
			readsFilename = argVal
		elif (arg.startswith("--cigars=")) or (arg.startswith("--cigar=")):
			cigarFilename = argVal
		elif (arg.startswith("--intervals=")) or (arg.startswith("--interval=")):
			intervalsFilename = argVal
		elif (arg.startswith("--chromosome=")) or (arg.startswith("--chromosomes=")) \
		  or (arg.startswith("--chrom="))      or (arg.startswith("--chroms=")):
			if (chromsOfInterest == None): chromsOfInterest = set()
			for chrom in argVal.split(","):
				chromsOfInterest.add(chrom)
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	if (genomeFilename == None):
		usage("you need to give me a genome file")
	if (readsFilename == None):
		usage("you need to give me a reads file")
	if (cigarFilename == None):
		usage("you need to give me a cigar strings file")

	# read the intervals
	#
	# nota bene: this can modify chromsOfInterest, restricting it to the
	# chromosomes in the intervals list

	chromToIntervals = None

	if (intervalsFilename != None):
		chromToIntervals = {}

		intervalsF = file(intervalsFilename,"rt")

		for (chrom,gStart,gEnd,tags) in read_intervals(intervalsF):
			if (chromsOfInterest != None) and (chrom not in chromsOfInterest): continue
			if (chrom not in chromToIntervals): chromToIntervals[chrom] = []
			chromToIntervals[chrom] += [(gStart,gEnd)]

		intervalsF.close()

		if (chromsOfInterest == None):
			chromsOfInterest = set(chromToIntervals)
		else:
			for chrom in chromsOfInterest:
				if (chrom not in chromToIntervals):
					chromsOfInterest.remove(chrom)

	# read the genome

	chromToSequence = {}

	if (genomeFilename.endswith(".gz")) or (genomeFilename.endswith(".gzip")):
		genomeF = gzip_open(genomeFilename,"rt")
	else:
		genomeF = file(genomeFilename,"rt")

	for (chrom,seq) in read_fasta_sequences(genomeF,chromsOfInterest):
		if (chrom in chromToSequence):
			exit("%s: \"%s\" appears more than once in \"%s\""
			   % (os_path.basename(argv[0]),chrom,genomeFilename))
		chromToSequence[chrom] = seq

	genomeF.close()

	if (chromsOfInterest != None):
		for chrom in chromsOfInterest:
			if (chrom not in chromToSequence):
				exit("%s: \"%s\" doesn't appear in \"%s\""
				   % (os_path.basename(argv[0]),chrom,genomeFilename))

	# read the cigar strings

	cigarF = file(cigarFilename,"rt")

	readNameToCigar = {}

	for (readName,chrom,strand,gStart,gEnd,cigar) in read_cigars(cigarF):
		if (chromsOfInterest != None) and (chrom not in chromsOfInterest): continue
		(gLength,rLength) = cigar_lengths(cigar)
		readNameToCigar[readName] = (chrom,gStart,gEnd,gLength,strand,rLength,cigar)
		assert (gLength == gEnd-gStart)  # ... yank this

	cigarF.close()

	# process the reads

	chromW    = max([len(chrom)    for chrom    in chromsOfInterest])
	readNameW = max([len(readName) for readName in readNameToCigar])
	nameW     = max(chromW+1,readNameW)

	if (readsFilename.endswith(".gz")) or (readsFilename.endswith(".gzip")):
		readsF = gzip_open(readsFilename,"rt")
	else:
		readsF = file(readsFilename,"rt")

	isFirst = True
	for (readName,rNucs) in read_fasta_sequences(readsF):
		if (readName not in readNameToCigar):
			exit("%s: \"%s\" doesn't appear in \"%s\""
			   % (os_path.basename(argv[0]),readNameToCigar,cigarFilename))

		(chrom,gStart,gEnd,gLength,strand,rLength,cigar) = readNameToCigar[readName]
		assert (rLength == len(rNucs))  # ... yank this
		gNucs = chromToSequence[chrom][gStart:gEnd]

		if (strand == "-"):
			gNucs = reverse_complement(gNucs)

		(gText,rText) = reconstruct_alignment(gNucs,rNucs,cigar)
		# $$$ if we have intervals, truncate alignment to the intersecting intervals

		if (isFirst): isFirst = False
		else:         print

		print "%-*s %s" % (nameW,chrom+strand,rText)
		print "%-*s %s" % (nameW,readName,gText)

	readsF.close()


def read_intervals(f):
	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		if (line == ""): continue
		if (line.startswith("#")): continue

		fields = line.split()
		if (len(fields) < 3):
			exit("%s: not enough fields at line %d (%d, expected at least %d)"
			   % (os_path.basename(argv[0]),lineNumber,len(fields),3))

		try:
			chrom  = fields[0]
			gStart = int(fields[1])
			gEnd   = int(fields[2])
			if (gEnd < gStart): raise ValueError
			tags   = None if (len(fields) == 3) else fields[3:]
		except ValueError:
			exit("%s: bad interval (at line %d): %s"
			   % (os_path.basename(argv[0]),lineNumber,line))

		yield (chrom,gStart,gEnd,tags)


def read_fasta_sequences(f,namesOfInterest=None):
	name     = None
	ignoring = False
	lines    = None

	lineNum = 0
	for line in f:
		lineNum += 1
		line = line.rstrip()

		if (line.startswith(">")):
			if (lines != None):
				seq = "".join(lines)
				yield (name,seq)
			name = line[1:].strip()
			if (namesOfInterest == None) or (name in namesOfInterest):
				lines = []
			else:
				ignoring = True
				lines    = None
		elif (ignoring):
			pass
		elif (name == None):
			exit("%s: first sequence has no header"
			   % os_path.basename(argv[0]))
		else:
			lines += [line.upper()]

	if (lines != None):
		seq = "".join(lines)
		yield (name,seq)


def read_cigars(f):
	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		if (line == ""): continue
		if (line.startswith("#")): continue

		fields = line.split()
		if (len(fields) < 5):
			exit("%s: not enough fields at line %d (%d, expected at least %d)"
			   % (os_path.basename(argv[0]),lineNumber,len(fields),5))

		try:
			name   = fields[0]
			chrom  = fields[1]
			gStart = int(fields[2])
			gEnd   = int(fields[3])
			cigar  = " ".join(fields[4:])
			if (gEnd < gStart): raise ValueError
		except ValueError:
			exit("%s: bad cigar line (at line %d): %s"
			   % (os_path.basename(argv[0]),lineNumber,line))

		if   (chrom.endswith("+")): (chrom,strand) = (chrom[:-1],"+")
		elif (chrom.endswith("-")): (chrom,strand) = (chrom[:-1],"-")
		else:                        strand = "+"

		try:
			cigar = list(cigar_ops(cigar))
		except ValueError:
			exit("%s: unparsable cigar string (at line %d): %s"
			   % (os_path.basename(argv[0]),lineNumber,cigar))

		yield (name,chrom,strand,gStart,gEnd,cigar)


# cigar_ops--
#	Convert cigar string into a series of (count,op) pairs

def cigar_ops(cigar):
	count = ""
	for ch in cigar:
		if (ch in "0123456789"):
			count += ch
			if (count == "0"): raise ValueError
		elif (ch in "MID"):
			if (count == ""): raise ValueError
			yield (int(count),ch)
			count = ""
		elif (count == "") and (ch in [" ","\t"]):
			pass  # allow whitespace before count
		else:
			raise ValueError

	if (count != ""): raise ValueError


# cigar_lengths--

def cigar_lengths(cigar):
	gLength = rLength = 0
	for (count,op) in cigar:
		if (op == "I"):
			rLength += count
		elif (op == "D"):
			gLength += count
		else: # if (op == "M"):
			gLength += count
			rLength += count

	return (gLength,rLength)

# reconstruct_alignment--

def reconstruct_alignment(gNucs,rNucs,cigar):
	gText = []
	rText = []

	gPos = rPos = 0
	for (count,op) in cigar:
		if (op == "I"):
			gText += ["-"*count]
			rText += [rNucs[rPos:rPos+count]]
			rPos += count
		elif (op == "D"):
			gText += [gNucs[gPos:gPos+count]]
			rText += ["-"*count]
			gPos += count
		else: # if (op == "M"):
			gText += [gNucs[gPos:gPos+count]]
			rText += [rNucs[rPos:rPos+count]]
			gPos += count
			rPos += count

	return ("".join(gText),"".join(rText))


if __name__ == "__main__": main()
