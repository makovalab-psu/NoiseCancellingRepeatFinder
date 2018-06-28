#!/usr/bin/env python
"""
Map intervals from a "genome" to positions on simulated reads.
"""

from sys  import argv,stdin,stdout,stderr,exit


def usage(s=None):
	message = """
usage: [cat <intervals_file> |] map_onto_simulated_reads [options]
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

Given a genome from which simulated reads were sampled by simulate_reads_v4,
and the corresponding cigars file, map intervals (or positions) from the genome
to the corresponding positions on the simulated reads.

Intervals are one per line, <chrom> start> <end>, origin-zero half-open. Any
additional columns are copied to the output."""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global debug

	# parse the command line

	cigarFilename      = None
	strandedTags       = None
	indicateTruncation = True
	sortByReads        = False
	separateIntervals  = False
	debug              = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--cigars=")) or (arg.startswith("--cigar=")):
			cigarFilename = argVal
		elif (arg.startswith("--stranded=")):
			if (strandedTags == None): strandedTags = set()
			for col in map(int,argVal.split(",")):
				strandedTags.add(col)
		elif (arg == "--truncate"):
			indicateTruncation = False
		elif (arg == "--sortby:reads"):
			sortByReads = True
		elif (arg == "--separators"):
			separateIntervals = True
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	if (cigarFilename == None):
		usage("you need to give me a cigar strings file")

	if (strandedTags != None):
		strandedTags = [col-4 for col in strandedTags]
		strandedTags.sort()

	# read the cigar strings

	cigarF = file(cigarFilename,"rt")

	chroms = []
	chromToCigars = {}
	nameToGenome = {}

	for (name,chrom,strand,gStart,gEnd,cigar) in read_cigars(cigarF):
		if (strand == "-"): cigar = cigar[::-1]  # reverse order of cigar ops
		(gLength,rLength) = cigar_lengths(cigar)
		if (chrom not in chromToCigars):
			chroms += [chrom]
			chromToCigars[chrom] = []
		chromToCigars[chrom] += [(gStart,gEnd,gLength,name,strand,rLength,cigar)]
		nameToGenome[name] = (gStart,gEnd)

	cigarF.close()

	for chrom in chromToCigars:
		chromToCigars[chrom].sort()

	# process the intervals

	oppositeStrand = {"+":"-", "-":"+"}

	chromToMappings = {}
	for chrom in chroms: chromToMappings[chrom] = []

	haveOutput = False
	for (chrom,gStart,gEnd,tags) in read_intervals(stdin):
		if (chrom not in chromToCigars): continue
		cigarInfo = chromToCigars[chrom]

		needSeparator = separateIntervals
		for (name,strand,rStart,rEnd) in map_interval(cigarInfo,gStart,gEnd):
			if (indicateTruncation):
				if (type(rStart) == tuple): rStart = "%s%d" % rStart
				if (type(rEnd)   == tuple): rEnd   = "%s%d" % rEnd
			else:
				if (type(rStart) == tuple): rStart = rStart[1]
				if (type(rEnd)   == tuple): rEnd   = rEnd[1]

			if (tags == None):
				oTags = ""
			else:
				oTags = list(tags)
				if (strand == "-") and (strandedTags != None):
					for col in strandedTags:
						if (col >= len(oTags)): continue
						tailCh = oTags[col][-1]
						if (tailCh in "+-"):
							oTags[col] = oTags[col][:-1] + oppositeStrand[tailCh]
				oTags = "\t" + "\t".join(oTags)

			mappedStr = "%s\t%d\t%d\t%s\t%s\t%s\t%s" \
			          % (chrom,gStart,gEnd,name,rStart,rEnd,oTags)

			if (sortByReads):
				(s,e) = nameToGenome[name]
				chromToMappings[chrom] += [(s,e,rStart,rEnd,mappedStr)]
			else:
				if (haveOutput) and (needSeparator): print ""
				print mappedStr
				haveOutput    = True
				needSeparator = False

	if (sortByReads):
		haveOutput = False
		for chrom in chroms:
			chromToMappings[chrom].sort()
			needSeparator = separateIntervals
			for (_,_,_,_,mappedStr) in chromToMappings[chrom]:
				if (haveOutput) and (needSeparator): print ""
				print mappedStr
				haveOutput    = True
				needSeparator = False


def map_interval(cigarInfo,gStart,gEnd):
	# Note that insertions are nucleotides that are in the read but not in the
	# genome; deletions are nucleotides that are in the genome but not in the
	# read

	# Also note that the cigar operations list has already been reversed if
	# the read was pulled from revcomp of genome

	if ("mapping" in debug):
		print >>stderr, "mapping %d..%d" % (gStart,gEnd)

	for (s,e,gLength,name,strand,rLength,cigar) in cigarInfo:
		if (e <= gStart): continue
		if (s >= gEnd): break

		if ("mapping" in debug):
			print >>stderr, "  intersects with %d..%d" % (s,e)

		(gPos,rPos) = (s,0)
		rStart = rEnd = None

		if (gStart < s): rStart = ("<",0)

		for (count,op) in cigar:
			if ("mapping" in debug):
				print >>stderr, "  g=%d r=%d" % (gPos,rPos)

			if (rStart == None) and (gPos == gStart):
				rStart = rPos
			if (gPos == gEnd):
				rEnd = rPos
				break

			if (op == "I"):
				rPos += count
			elif (op == "D"):
				gPos += count
			else: # if (op == "M"):
				if (rStart == None) and (gPos < gStart < gPos+count):
					rStart = rPos + gStart-gPos
				if (gPos < gEnd < gPos+count):
					rEnd = rPos + gEnd-gPos
					break
				gPos += count
				rPos += count

		if ("mapping" in debug):
			print >>stderr, "  g=%d r=%d" % (gPos,rPos)

		if (rEnd == None):
			assert (rPos == rLength)
			if (gPos == gEnd): rEnd = rPos
			else:              rEnd = (">",rPos)

		assert (rStart != None)

		# if read was pulled from revcomp of genome, we need to reverse
		# the positions here

		if (strand == "-"):
			if (type(rEnd)   == tuple): reverseStart = ("<",0)
			else:                       reverseStart = rLength-rEnd
			if (type(rStart) == tuple): reverseEnd   = (">",rLength)
			else:                       reverseEnd   = rLength-rStart
			(rStart,rEnd) = (reverseStart,reverseEnd)

		yield (name,strand,rStart,rEnd)


def read_intervals(f):
	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		if (line == ""): continue
		if (line.startswith("#")): continue

		fields = line.split()
		assert (len(fields) >= 3), \
		       "not enough fields at line %d (%d, expected at least %d)" \
		     % (lineNumber,len(fields),3)

		try:
			chrom  = fields[0]
			gStart = int(fields[1])
			gEnd   = int(fields[2])
			if (gEnd < gStart): raise ValueError
			tags   = None if (len(fields) == 3) else fields[3:]
		except ValueError:
			assert (False), "bad line (%d): %s" % (lineNumber,line)

		yield (chrom,gStart,gEnd,tags)


def read_cigars(f):
	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()
		if (line == ""): continue
		if (line.startswith("#")): continue

		fields = line.split()
		assert (len(fields) >= 5), \
		       "not enough fields at line %d (%d, expected at lest %d)" \
		     % (lineNumber,len(fields),5)

		try:
			name   = fields[0]
			chrom  = fields[1]
			gStart = int(fields[2])
			gEnd   = int(fields[3])
			cigar  = " ".join(fields[4:])
			if (gEnd < gStart): raise ValueError
		except ValueError:
			assert (False), "bad cigar line (%d): %s" % (lineNumber,line)

		if   (chrom.endswith("+")): (chrom,strand) = (chrom[:-1],"+")
		elif (chrom.endswith("-")): (chrom,strand) = (chrom[:-1],"-")
		else:                        strand = "+"

		try:
			cigar = list(cigar_ops(cigar))
		except ValueError:
			assert (False), "unparsable cigar string (line %d): %s" % (lineNumber,cigar)

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


if __name__ == "__main__": main()
