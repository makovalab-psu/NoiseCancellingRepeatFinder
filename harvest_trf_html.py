#!/usr/bin/env python
"""
Convert alignments from TRF (Tandem Repeat Finder) html output to a tabular
text format similar to an NcRF summary.
"""

from sys    import argv,stdin,stdout,stderr,exit
from os     import path as os_path
from re     import compile as re_compile
from string import maketrans

class Alignment: pass


def usage(s=None):
	message = """
usage: cat <trf_html_output> | harvest_trf_html [options]
  --motif=<motif>        (cumulative) motifs of interest; alignments for other
                         motifs are discarded
                         (if this is not provided, we keep all alignments)
  --minlength=<bp>       discard alignments that don't have long enough repeat
                         (but default, we don't filter by length)
  --withheader           include a header line in the output
  --withalignment        include alignment text in the output"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global debug

	# parse the command line

	motifs             = None
	minLength          = None
	writeHeader        = False
	writeAlignmentText = False
	debug              = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--motif=")):
			if (motifs == None): motifs = set()
			motifs.add(argVal)
		elif (arg.startswith("--minlength=")) or (arg.startswith("--minlen=")):
			try:
				minLength = int(argVal)
				if (minLength < 0): raise ValueError
				if (minLength == 0): minLength = None
			except ValueError:
				usage("bad length in \"%s\"" % arg)
		elif (arg in ["--withheader","--with=header","--with:header"]):
			writeHeader = True
		elif (arg in ["--withalignment","--with=alignment","--with:alignment"]):
			writeAlignmentText = True
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	if (motifs != None):
		motifMap = {}
		for motif in motifs:
			(canonicalMotif,strand) = canonical_motif(motif)
			if (canonicalMotif in motifMap):
				print >>stderr, "WARNING: motif %s is the same as %s (we'll report it as %s)" \
				              % (motif,motifMap[canonicalMotif][0],motifMap[canonicalMotif][0])
			else:
				motifMap[canonicalMotif] = (motif,strand)

	# process the alignments

	for a in read_trf_html_alignments(stdin):
		if (writeHeader):
			header = ["line","trfScore","motif",
			          "seq","start","end","strand","seqLen",
			          "querybp",
			          "mRatio","m","mm","io","ix","do","dx"]
			if (writeAlignmentText):
				header += ["seqAligned","motifAligned"]

			print "#" + "\t".join(header)
			writeHeader = False

		if (motifs != None):
			if (a.motif not in motifMap): continue
			(a.motif,strand) = motifMap[a.motif]
			if (strand == "-"):
				a.strand = "-" if (a.strand == "+") else "+"

		(nMatch,nMismatch,nInsO,nInsX,nDelO,nDelX) = extract_events(a)
		nEvents = nMatch + nMismatch + nInsO + nInsX + nDelO + nDelX
		queryBp = nMatch + nMismatch + nDelO + nDelX

		if (minLength != None):
			if (queryBp < minLength): continue

		vec = [a.lineNumber,a.score,a.motif,
		       a.seqName,a.start,a.end,a.strand,a.seqLength,
		       queryBp,
		       "%.1f" % (100.0*nMatch / nEvents),
		       nMatch,nMismatch,nInsO,nInsX,nDelO,nDelX]
		if (writeAlignmentText):
			vec += [a.seqText,a.motifText]
		print "\t".join(map(str,vec))


# read_trf_html_alignments--
#   yields the next alignment from a TRF html file

def read_trf_html_alignments(f):
	(lineNumber,line) = skip_til_line(f,"Sequence: ")
	if (line == None): return

	# $$$ need to finish converting asserts to error messages

	try:
		(seqName,seqLength) = parse_sequence_line(line)
	except ValueError:
		exit(("%s: expected \"Sequence\" line in TRF html"
		       + " (in alignment at line %d)")
		   % (os_path.basename(argv[0]),lineNumber))

	while (True):
		# consume e.g. "Found at i:3332 original size:5 final size:5"

		(lineCount,line) = skip_til_line(f,"Found at ")
		lineNumber += lineCount
		if (line == None): return

		a = Alignment()
		a.lineNumber = lineNumber
		a.seqName    = seqName
		a.seqLength  = seqLength

		# consume e.g. "Indices: 3322--3565  Score: 104"

		(lineCount,line) = skip_til_line(f,"Indices: ",stripPrefix=True)
		if (line == None):
			exit(("%s: expected \"Indices\" line in TRF html"
			       + " (in alignment at line %d)")
			   % (os_path.basename(argv[0]),a.lineNumber))
		lineNumber += lineCount

		(start,end,score) = parse_indices_line(line)
		a.start = start
		a.end   = end
		a.score = score

		# consume e.g. "Period size: 5  Copynumber: 74.4  Consensus size: 5"

		(lineCount,line) = skip_til_line(f,"Period size: ",stripPrefix=True)
		if (line == None):
			exit(("%s: expected \"Period size\" line in TRF html"
			       + " (at or after line %d, in alignment at line %d)")
			   % (os_path.basename(argv[0]),lineNumber+1,a.lineNumber))
		if (lineCount != 1):
			exit(("%s: unexpected lines before \"Period size\" line in TRF html"
			       + " (at or after line %d, in alignment at line %d)")
			   % (os_path.basename(argv[0]),lineNumber+1,a.lineNumber))
		lineNumber += lineCount

		(period,copyNumber,consensusSize) = parse_period_size_line(line)
		a.period        = period
		a.copyNumber    = copyNumber
		a.consensusSize = consensusSize

		# consume the alignment text

		if ("triplets" in debug):
			print >>stderr, "parse_alignment_text at %d" % (lineNumber+1)

		(lineCount,a.errorText,a.seqText,a.motifText) = parse_alignment_text(f,dbgLineNumber=lineNumber)
		lineNumber += lineCount
		assert (a.seqText != None)

		# consume e.g.
		#   "Statistics"
		#   "Matches: 191,  Mismatches: 26, Indels: 44"
		#
		# nota bene:
		#   I'd like to sanity check these TRF stats vs what I count in the
		#   alignment, but I can't make any sense out of them. This line from
		#     https://tandem.bu.edu/trf/trf.definitions.html#alignment
		#   is a clue:
		#     "Statistics refers to the matches, mismatches and indels overall
		#     between adjacent copies in the sequence, not between the sequence
		#     and the consensus pattern."

		(lineCount,line) = skip_til_line(f,"Statistics")
		lineNumber += lineCount
		if (line == None):
			exit(("%s: alignment appears to be truncated at line %d"
			       + " (in alignment at line %d);"
			       + "\n  .. expected \"Statistics\"")
			   % (os_path.basename(argv[0]),lineNumber,a.lineNumber+1))
		if (lineCount != 1):
			exit(("%s: unexpected lines before \"Statistics\"at line %d"
			       + " (in alignment at line %d)")
			   % (os_path.basename(argv[0]),lineNumber,a.lineNumber+1))

		#lineNumber += 1
		#line = f.readline()
		#assert (line != "")
		#(matches,mismatches,indels) = parse_match_stats(line)
		#(eMatches,eMismatches,eIndels) = count_match_stats(a.errorText)
		#assert (matches == eMatches)  this won't pass sanity check
		#assert (mismatches == eIndels)
		#assert (indels == eMismatches+eIndels)

		# consume e.g.
		#   "Consensus pattern (5 bp):   "
		#   "TTCCA"
		# note that the consensus can extend for more than one line; it seems
		# to wrap at 65 bp.

		(lineCount,line) = skip_til_line(f,"Consensus pattern ")
		lineNumber += lineCount
		assert (line != None)

		nucs = ""
		while (len(nucs) < a.consensusSize):
			lineNumber += 1
			line = f.readline()
			if (line == ""): break
			nucs += line.rstrip()

		(a.motif,a.strand) = canonical_motif(nucs)
		assert (len(a.motif) == a.consensusSize)

		yield a


# parse_alignment_text--

def parse_alignment_text(f,dbgLineNumber=0):
	if (dbgLineNumber != None):
		inAlignmentAt = " (in alignment at line %d)" % (dbgLineNumber+1)
	else:
		inAlignmentAt = ""

	# consume the flanking prefix; note that we're not guaranteed to always
	# have a flanking prefix -- if the alignment starts at position 1, then the
	# first line is the mmText line of the first triplet

	(lineCount,line) = skip_blank_lines(f)
	if (line == None):
		exit("%s: expected alignment text at line %d, but line was blank%s"
		   % (os_path.basename(argv[0]),dbgLineNumber+lineCount,inAlignmentAt))

	flankingPrefix = line

	# consume one blank line after the flanking prefix; no that this *might*
	# actually be the seqText line of the first triplet

	lineCount += 1
	line = f.readline()

	accidentalMmText = accidentalSeqText = None
	if (line.strip() != ""):
		accidentalMmText  = flankingPrefix.rstrip("\n")
		accidentalSeqText = line.rstrip()

	# consume alignment triplets, e.g.
    #   "          *      *                            *             ** "
    #   "3824 TTCGCG TTCGCG TTCCA TTCCA TTCCA TTCCA TTCTA TT-CA GTGTTAAA"
    #   "   1 TTC-CA TTC-CA TTCCA TTCCA TTCCA TTCCA TTCCA TTCCA ---TTCCA"
    # each followed by one blank line.

	allErrorText = []
	allSeqText   = []
	allMotifText = []

	while (True):
		if (accidentalMmText != None):
			mmTextLine  = accidentalMmText
			seqTextLine = accidentalSeqText
			accidentalMmText = accidentalSeqText = None
		else:
			# consume mismatch/indel line

			lineCount += 1
			line = f.readline()
			if (line == ""):
				exit(("%s: alignment appears to be truncated at line %d%s;"
				       + "\n  .. expected mismatch/indel string")
				   % (os_path.basename(argv[0]),dbgLineNumber+lineCount,inAlignmentAt))
			mmTextLine = line.rstrip("\n")
			if (mmTextLine == ""): break  # alignment has no flanking suffix
			if ("triplets" in debug):
				print >>stderr, "  mmText at %d" % lineCount

			# consume aligned sequence line

			lineCount += 1
			line = f.readline()
			if (line == ""):
				exit(("%s: alignment appears to be truncated at line %d%s;"
				       + "\n  .. expected aligned sequence string")
				   % (os_path.basename(argv[0]),dbgLineNumber+lineCount,inAlignmentAt))
			seqTextLine = line.rstrip()

		# consume aligned repeat element line

		lineCount += 1
		line = f.readline()
		if (line == ""):
			exit(("%s: alignment appears to be truncated at line %d%s;"
			       + "\n  .. expected aligned repeat element string")
			   % (os_path.basename(argv[0]),dbgLineNumber+lineCount,inAlignmentAt))
		motifTextLine = line.rstrip()
		if (motifTextLine == ""): break  # mismatch/indel line was the flanking suffix

		# consume blank line

		lineCount += 1
		line = f.readline()
		if (line == ""):
			exit("%s: alignment appears to be truncated at line %d%s"
			   % (os_path.basename(argv[0]),dbgLineNumber+lineCount,inAlignmentAt))
		line = line.strip()
		if (line != ""):
			exit("%s: expected blank line at line %d%s, but got:\n%s"
			   % (os_path.basename(argv[0]),dbgLineNumber+lineCount,inAlignmentAt,line))

		# parse the triplet

		(indent,mmText) = (mmTextLine[:12],mmTextLine[12:])
		if (indent.rstrip() != ""):
			exit("%s: expected 12 character indentation in mismatch/indel string%s, but got:\n%s"
			   % (os_path.basename(argv[0]),inAlignmentAt,mmTextLine))

		# $$$ sanity check the seqStart values (and the motifStart values too)
		(seqStart,spacer,seqText) = (seqTextLine[:11],seqTextLine[11],seqTextLine[12:])
		try:
			if (spacer != " "): raise ValueError
			seqStart = int(seqStart)
		except ValueError:
			exit("%s: problem parsing position in aligned sequence string%s:\n%s"
			   % (os_path.basename(argv[0]),inAlignmentAt,seqTextLine))

		(motifStart,spacer,motifText) = (motifTextLine[:11],motifTextLine[11],motifTextLine[12:])
		try:
			if (spacer != " "): raise ValueError
			motifStart = int(motifStart)
			#if (motifStart != 1): raise ValueError  test fails when alignments wrap (after 65 characters)
		except ValueError:
			exit("%s: problem parsing position in aligned repeat element string%s:\n%s"
			   % (os_path.basename(argv[0]),inAlignmentAt,motifTextLine))

		if (len(motifText) != len(seqText)):
			exit("%s: aligned sequence and aligned repeat element strings have different lengths%s"
			   % (os_path.basename(argv[0]),inAlignmentAt))
		if (len(mmText) < len(seqText)): mmText += " " * (len(seqText)-len(mmText))

		errorText    = []
		newSeqText   = []
		newMotifText = []
		for (mmCh,seqCh,motifCh) in zip(mmText,seqText,motifText):
			if (seqCh == " "):
				assert (motifCh == " ")
				assert (mmCh    == " ")
				continue

			assert (motifCh != " ")
			newSeqText   += [seqCh]
			newMotifText += [motifCh]

			if (seqCh == "-") or (motifCh == "-"):
				errorText += ["x"]
			elif (seqCh == motifCh):
				assert (mmCh == " ")
				errorText += ["="]
			else:
				assert (mmCh == "*")
				errorText += ["*"]

		allErrorText += ["".join(errorText)]
		allSeqText   += ["".join(newSeqText)]
		allMotifText += ["".join(newMotifText)]

	return (lineCount,"".join(allErrorText),"".join(allSeqText),"".join(allMotifText))


# parse_sequence_line--
#	e.g. "Sequence: SRR2036394.5614"
#	OR   "Sequence: SRR2036394.5614 length=12345"

def parse_sequence_line(s):
	sequenceRe1 = re_compile("^Sequence: +(?P<name>[^ ]+)$")
	sequenceRe2 = re_compile("^Sequence: +(?P<name>[^ ]+)"
	                       + " +length=(?P<length>[0-9]+)$")
	m = sequenceRe1.match(s)
	if (m != None):
		name = m.group("name")
		return (name,None)
	m = sequenceRe2.match(s)
	if (m != None):
		name   =     m.group("name")
		length = int(m.group("length"))
		return (name,length)
	raise ValueError


# parse_indices_line--
#	e.g. "Indices: 3322--3565  Score: 104"

def parse_indices_line(s):
	indicesRe = re_compile("^Indices: +(?P<start>[0-9]+)--(?P<end>[0-9]+)"
	                     + " +Score: +(?P<score>[0-9]+)$")
	m = indicesRe.match(s)
	if (m == None): raise ValueError
	start = int(m.group("start")) - 1  # subtract 1 because TRF start is origin
	end   = int(m.group("end"))        # .. one and we want origin zero
	score = int(m.group("score"))
	return (start,end,score)


# parse_period_size_line--
#	e.g. "Period size: 5  Copynumber: 74.4  Consensus size: 5"

def parse_period_size_line(s):
	periodSizeRe = re_compile("^Period size: +(?P<period>[0-9]+)"
	                        + " +Copynumber: +(?P<copyNumber>[0-9]+\.[0-9]+)"
	                        + " +Consensus size: +(?P<consensusSize>[0-9]+)$")
	m = periodSizeRe.match(s)
	if (m == None): raise ValueError
	period        = int  (m.group("period"))
	copyNumber    = float(m.group("copyNumber"))
	consensusSize = int  (m.group("consensusSize"))
	return (period,copyNumber,consensusSize)


# parse_match_stats--
#	e.g. "Matches: 191,  Mismatches: 26, Indels: 44"

def parse_match_stats(s):
	matchStatsRe = re_compile("^Matches: +(?P<matches>[0-9]+),"
	                        + " +Mismatches: +(?P<mismatches>[0-9]+),"
	                        + " +Indels: +(?P<indels>[0-9]+)$")
	m = matchStatsRe.match(s)
	if (m == None): raise ValueError
	matches    = int(m.group("matches"))
	mismatches = int(m.group("mismatches"))
	indels     = int(m.group("indels"))
	return (matches,mismatches,indels)


# skip_til_line--
#	Read lines from a file until we find one that starts with a specified
#	trigger.

def skip_til_line(f,trigger,stripPrefix=False):
	lineCount = 0
	while (True):
		lineCount += 1
		line = f.readline()
		if (line == ""): return (lineCount,None)
		if (stripPrefix): line = line.strip()
		else:             line = line.rstrip()

		if (line.startswith(trigger)):
			return (lineCount,line)


# skip_blank_lines--
#	Read lines from a file until we find one that is not blank.

def skip_blank_lines(f,stripPrefix=False):
	lineCount = 0
	while (True):
		lineCount += 1
		line = f.readline()
		if (line == ""): return (lineCount,None)
		if (stripPrefix): line = line.strip()
		else:             line = line.rstrip("\n")

		if (line != ""): return (lineCount,line)


# count_match_stats--

def count_match_stats(errorText):
	matches = mismatches = indels = 0
	for errorCh in errorText:
		if   (errorCh == "="): matches    += 1
		elif (errorCh == "*"): mismatches += 1
		elif (errorCh == "x"): indels     += 1
		else:
			exit("%s: unexpected \"%c\" in error text"
			   % (os_path.basename(argv[0]),errorCh))
	return (matches,mismatches,indels)


# extract_events--

def extract_events(a):
	if (len(a.seqText) != len(a.motifText)):
		exit(("%s: internal error, alignment text lengths aren't the same"
		       + " for alignment at line %d")
		   % (os_path.basename(argv[0]),a.lineNumber))

	nMatch = nMismatch = nInsO = nInsX = nDelO = nDelX = 0

	prevEvent = None
	for (ix,(seqCh,motifCh)) in enumerate(zip(a.seqText,a.motifText)):
		if (seqCh == "-") and (motifCh == "-"):
			exit(("%s: internal error, gap aligned to gap in column %d" \
			       + " of alignment at line %d")
			   % (os_path.basename(argv[0]),ix+1,a.lineNumber))

		if (seqCh == "-"):
			if (prevEvent != "d"):
				nDelO += 1
				prevEvent = "d"
			else:
				nDelX += 1
		elif (motifCh == "-"):
			if (prevEvent != "i"):
				nInsO += 1
				prevEvent = "i"
			else:
				nInsX += 1
		elif (seqCh == motifCh):
			nMatch += 1
			prevEvent = "m"
		else: # if (seqCh != motifCh):
			nMismatch += 1
			prevEvent = "mm"

	return (nMatch,nMismatch,nInsO,nInsX,nDelO,nDelX)


# canonical_motif--

def canonical_motif(s):
	(minNucs,strand) = (s,"+")
	rev = reverse_complement(s)
	if (rev < minNucs): (minNucs,strand) = (rev,"-")
	for rotIx in xrange(1,len(s)):
		rotNucs = s[rotIx:] + s[:rotIx]
		if (rotNucs < minNucs): (minNucs,strand) = (rotNucs,"+")
		rotNucs = rev[rotIx:] + rev[:rotIx]
		if (rotNucs < minNucs): (minNucs,strand) = (rotNucs,"-")
	return (minNucs,strand)


# reverse_complement--

complementMap = maketrans("ACGTSWRYMKBDHVNacgtswrymkbdhvn",
                          "TGCASWYRKMVHDBNtgcaswyrkmvhdbn")

def reverse_complement(nukes):
	return nukes[::-1].translate(complementMap)


if __name__ == "__main__": main()
