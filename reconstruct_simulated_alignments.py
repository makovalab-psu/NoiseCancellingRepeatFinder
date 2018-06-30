#!/usr/bin/env python
"""
Reconstruct alignments between a genome and simulated reads sampled from the
genome.
"""

from sys        import argv,stdin,stdout,stderr,exit
from os         import path as os_path
from gzip       import open as gzip_open
from ncrf_parse import reverse_complement

class Alignment: pass


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
	global nameFieldW,lengthFieldW,countFieldW,rangeFieldW
	global debug

	# parse the command line

	genomeFilename    = None
	readsFilename     = None
	cigarFilename     = None
	intervalsFilename = None
	chromsOfInterest  = None
	nameFieldW        = 1
	lengthFieldW      = 1
	countFieldW       = 1
	rangeFieldW       = 1
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
		elif (arg.startswith("--fields=")) or (arg.startswith("F=")):
			(nameFieldW,lengthFieldW,countFieldW,rangeFieldW) = argVal.split(",",4)
			nameFieldW   = max(int(nameFieldW),1)
			lengthFieldW = max(int(lengthFieldW),1)
			countFieldW  = max(int(countFieldW),1)
			rangeFieldW  = max(int(rangeFieldW),1)
		elif (arg.startswith("--namefield=")) or (arg.startswith("F1=")):
			nameFieldW = max(int(argVal),1)
		elif (arg.startswith("--lengthfield=")) or (arg.startswith("F2=")):
			lengthFieldW = max(int(argVal),1)
		elif (arg.startswith("--countfield=")) or (arg.startswith("F3=")):
			countFieldW = max(int(argVal),1)
		elif (arg.startswith("--intervalfield=")) or (arg.startswith("F4=")):
			rangeFieldW = max(int(argVal),1)
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

		for chrom in chromToIntervals:
			chromToIntervals[chrom].sort()

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

	for (lineNumber,line,readName,chrom,strand,gStart,gEnd,cigar) in read_cigars(cigarF):
		if (chromsOfInterest != None) and (chrom not in chromsOfInterest): continue
		(rLength,gLength) = cigar_lengths(cigar)
		readNameToCigar[readName] = (chrom,gStart,gEnd,gLength,strand,rLength,cigar)
		if (gLength != gEnd-gStart):
			exit("%s: bad cigar line (at line %d); cigar doesn't match interval length (%d vs %d)\n%s"
			   % (os_path.basename(argv[0]),lineNumber,gLength,gEnd-gStart,line))

	cigarF.close()

	# process the reads

	if (readsFilename.endswith(".gz")) or (readsFilename.endswith(".gzip")):
		readsF = gzip_open(readsFilename,"rt")
	else:
		readsF = file(readsFilename,"rt")

	for (readName,rNucs) in read_fasta_sequences(readsF):
		if (readName not in readNameToCigar):
			exit("%s: \"%s\" doesn't appear in \"%s\""
			   % (os_path.basename(argv[0]),readNameToCigar,cigarFilename))

		(chrom,gStart,gEnd,gLength,strand,rLength,cigar) = readNameToCigar[readName]
		gNucs = chromToSequence[chrom][gStart:gEnd]

		if (strand == "-"):
			gNucs = reverse_complement(gNucs)

		a = Alignment()
		a.readName = readName
		a.rStart   = 0
		a.rEnd     = rLength
		a.rLength  = rLength
		a.rNucs    = rNucs
		a.chrom    = chrom
		a.strand   = strand
		a.gStart   = gStart
		a.gEnd     = gEnd
		a.gNucs    = gNucs
		a.score    = 0
		a.motif    = "%s:%d-%d%s" % (chrom,a.gStart,a.gEnd,strand)

		(a.rText,a.gText) = reconstruct_alignment(rNucs,gNucs,cigar)

		if (chromToIntervals == None):
			print_alignment(a)
		else:
			intervals = chromToIntervals[chrom]
			for (s,e) in intersecting_intervals(intervals,gStart,gEnd):
				aSliced = slice_alignment(a,s,e)
				print_alignment(aSliced)

				if ("intervalsanity" in debug):
					rText    = remove_gaps(a.rText)
					realText = rNucs[a.rStart:a.rEnd]
					if (realText != rText):
						exit("%s: sanity check failed for read:\n\"%s\"\n\"%s\""
						   % (os_path.basename(argv[0]),rText,realText))

					gText    = remove_gaps(a.gText).upper()
					realText = chromToSequence[chrom][a.gStart:a.gEnd]
					if (strand == "-"): realText = reverse_complement(realText)
					if (realText != gText):
						exit("%s: sanity check failed for genome:\n\"%s\"\n\"%s\""
						   % (os_path.basename(argv[0]),gText,realText))

	readsF.close()


# print_alignment--

havePrintedAlignments = False

def print_alignment(a):
	global havePrintedAlignments

	(nMatch,nMismatch,nInsO,nInsX,nDelO,nDelX) = extract_events(a.rText,a.gText)
	(eText,rText,gText) = alignment_to_error_text(a.rText,a.gText)

	mCount     = nMatch
	mmCount    = nMismatch
	iCount     = nInsO+nInsX
	dCount     = nDelO+nDelX
	matchRatio = float(mCount) / (mCount + mmCount + iCount + dCount)

	rangeStr        = "%d-%d" % (a.rStart,a.rEnd)
	scoreStr        = "score=%d" % a.score
	motifWithStrand = a.motif
	rLengthStr      = "%d" %   a.rLength
	rBaseCountStr   = "%dbp" % len(a.rNucs)
	gBaseCountStr   = "%dbp" % len(a.gNucs)

	nameW   = max(nameFieldW,len(a.readName),len(motifWithStrand))
	lengthW = max(lengthFieldW,len(rLengthStr))
	countW  = max(countFieldW,len(rBaseCountStr),len(gBaseCountStr))
	rangeW  = max(rangeFieldW,len(rangeStr),len(scoreStr))

	statsStr = "# score=%d querybp=%d mRatio=%.1f%% m=%d mm=%d i=%d d=%d" \
	         % (a.score,len(a.gNucs),100*matchRatio,mCount,mmCount,iCount,dCount)
	if (len(statsStr) > nameW+lengthW+countW+rangeW+3):
		nameW = len(statsStr) - (lengthW+countW+rangeW+3)

	if (not havePrintedAlignments):
		havePrintedAlignments = True
	else:
		print

	print "%-*s %s" % (nameW+lengthW+countW+rangeW+3,statsStr,eText)

	print "%-*s %-*s %-*s %-*s %s" \
	    % (nameW,   a.readName,
	       lengthW, rLengthStr,
	       countW,  rBaseCountStr,
	       rangeW,  rangeStr,
	       a.rText)
	print "%-*s %-*s %-*s %-*s %s" \
	    % (nameW,   motifWithStrand,
	       lengthW, "",
	       countW,  gBaseCountStr,
	       rangeW,  scoreStr,
	       a.gText)


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

		yield (lineNumber,line,name,chrom,strand,gStart,gEnd,cigar)


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
			rLength += count
			gLength += count

	return (rLength,gLength)


# reconstruct_alignment--

def reconstruct_alignment(rNucs,gNucs,cigar):
	rText = []
	gText = []

	gPos = rPos = 0
	for (count,op) in cigar:
		if (op == "I"):
			rText += [rNucs[rPos:rPos+count]]
			gText += ["-"*count]
			rPos += count
		elif (op == "D"):
			rText += ["-"*count]
			gText += [gNucs[gPos:gPos+count]]
			gPos += count
		else: # if (op == "M"):
			rText += [rNucs[rPos:rPos+count]]
			gText += [gNucs[gPos:gPos+count]]
			rPos += count
			gPos += count

	assert (rPos == len(rNucs))
	assert (gPos == len(gNucs))

	return ("".join(rText),"".join(gText))


# alignment_to_error_text--

def alignment_to_error_text(rText,gText):
	if (len(rText) != len(gText)):
		exit(("%s: internal error, alignment text lengths aren't the same"
		       + " for alignment at line %d")
		   % (os_path.basename(argv[0]),a.lineNumber))

	eText = []
	gNew  = []

	for (ix,(rCh,gCh)) in enumerate(zip(rText,gText)):
		if (rCh == gCh):
			eText += ["="]
			gNew  += [gCh]
		elif (rCh == "-") or (gCh == "-"):
			eText += ["x"]
			gNew  += [gCh]
		else:
			eText += ["x"]
			gNew  += [gCh.lower()]

	return ("".join(eText),rText,"".join(gNew))


# extract_events--

def extract_events(rText,gText):
	if (len(rText) != len(gText)):
		exit(("%s: internal error, alignment text lengths aren't the same"
		       + " for alignment at line %d")
		   % (os_path.basename(argv[0]),a.lineNumber))

	nMatch = nMismatch = nInsO = nInsX = nDelO = nDelX = 0

	prevEvent = None
	for (ix,(rCh,gCh)) in enumerate(zip(rText,gText)):
		if (rCh == "-") and (gCh == "-"):
			exit(("%s: internal error, gap aligned to gap in column %d" \
			       + " of alignment at line %d")
			   % (os_path.basename(argv[0]),ix+1,a.lineNumber))

		if (rCh == "-"):
			if (prevEvent != "d"):
				nDelO += 1
				prevEvent = "d"
			else:
				nDelX += 1
		elif (gCh == "-"):
			if (prevEvent != "i"):
				nInsO += 1
				prevEvent = "i"
			else:
				nInsX += 1
		elif (rCh == gCh):
			nMatch += 1
			prevEvent = "m"
		else: # if (rCh != gCh):
			nMismatch += 1
			prevEvent = "mm"

	return (nMatch,nMismatch,nInsO,nInsX,nDelO,nDelX)


# intersecting_intervals--
# 
# nota bene: we assume intervals have been sorted by increases start, but they
#            may be overlapping

def intersecting_intervals(intervals,start,end):
	for (s,e) in intervals:
		if (s >= end):   continue  # tempting to break, but that would be incorrect
		if (e <= start): continue
		yield (max(s,start),min(e,end))


# slice_alignment--

def slice_alignment(a,gStart,gEnd):
	aSliced = Alignment()
	aSliced.readName = a.readName
	aSliced.rLength  = a.rLength
	aSliced.chrom    = a.chrom
	aSliced.strand   = a.strand
	aSliced.score    = 1

	if (a.strand == "+"):
		(rText,gText) = (a.rText,a.gText)
	else:
		(rText,gText) = (reverse_complement(a.rText),reverse_complement(a.gText))

	startIx = endIx = None
	rSliceStart = rSliceEnd = None
	gSliceStart = gSliceEnd = None
	(rPos,gPos) = (0,a.gStart)
	for (ix,(rCh,gCh)) in enumerate(zip(rText,gText)):
		if (gPos >= gStart) and (startIx == None):
			startIx = ix
			(rSliceStart,gSliceStart) = (rPos,gPos)
		if (gPos >= gEnd):
			endIx = ix
			(rSliceEnd,gSliceEnd) = (rPos,gPos)
			break

		if (gCh == "-"):     # insertion
			rPos += 1
		elif (rCh == "-"):   # deletion
			gPos += 1
		else:                # match or mismatch
			rPos += 1
			gPos += 1

	if (endIx == None):
		endIx = len(a.gText)
		(rSliceEnd,gSliceEnd) = (rPos,gPos)

	if (a.strand == "+"):
		aSliced.rText  = rText[startIx:endIx]
		aSliced.gText  = gText[startIx:endIx]
	else:
		aSliced.rText  = reverse_complement(rText[startIx:endIx])
		aSliced.gText  = reverse_complement(gText[startIx:endIx])
		(rSliceStart,rSliceEnd) = (a.rEnd-rSliceEnd,a.rEnd-rSliceStart)

	aSliced.rNucs  = remove_gaps(aSliced.rText)
	aSliced.rStart = rSliceStart
	aSliced.rEnd   = rSliceEnd

	aSliced.gNucs  = remove_gaps(aSliced.gText).upper()
	aSliced.gStart = gSliceStart
	aSliced.gEnd   = gSliceEnd
	aSliced.motif  = "%s:%d-%d%s" \
	               % (aSliced.chrom,aSliced.gStart,aSliced.gEnd,aSliced.strand)

	return aSliced


# remove_gaps--

def remove_gaps(text):
	return "".join([nuc for nuc in text if (nuc != "-")])


if __name__ == "__main__": main()
