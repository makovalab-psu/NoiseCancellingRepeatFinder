#!/usr/bin/env python
"""
Extract (NON-positional) match-and-error event counts from alignments in sam 
format.
"""

from sys        import argv,stdin,stdout,stderr,exit
from ncrf_parse import int_with_unit,commatize


def usage(s=None):
	message = """

usage: cat <sam_file> | sam_to_event_matrix [options]
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
to infer_scoring."""

	if (s == None): exit (message)
	else:           exit ("%s%s" % (s,message))


def main():
	global warnOnError

	# parse the command line

	writeHeader    = False
	writeWhat      = "per alignment"
	warnOnError    = False
	headLimit      = None
	reportProgress = None

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg in ["--withheader","--with=header","--with:header"]):
			writeHeader = True
		elif (arg in ["--sumonly","--sum=only","--sum:only"]):
			writeWhat = "sum only"
		elif (arg == "--warnandcontinue"):
			warnOnError = True
		elif (arg.startswith("--head=")):
			headLimit = int_with_unit(argVal)
		elif (arg.startswith("--progress=")):
			reportProgress = int_with_unit(argVal)
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	# process the sam records

	sum = {"m":0, "mm":0, "io":0, "ix":0, "do":0, "dx":0}

	recordNum = alignmentNum = 0
	for a in read_sam_plain(stdin):
		recordNum += 1
		if (reportProgress != None) and (recordNum % reportProgress == 0):
			sum["events"] = (sum["m"] + sum["mm"] + sum["io"] + sum["ix"] + sum["do"] + sum["dx"])
			mRatio = float(sum["m"]) / sum["events"]
			vec = [mRatio,sum["m"],sum["mm"],sum["io"],sum["ix"],sum["do"],sum["dx"]]
			print >>stderr, "progress: processing sam record %s (mRatio=%.3f m=%d mm=%d io=%d ix=%d do=%d dx=%d)" \
			              % (commatize(recordNum),
			                 mRatio,sum["m"],sum["mm"],sum["io"],sum["ix"],sum["do"],sum["dx"])

		if (headLimit != None) and (recordNum > headLimit):
			print >>stderr, "limit of %s sam records reached" % commatize(headLimit)
			break

		if (a.rName == "*"): continue  # read did not align

		alignmentNum += 1
		events = sam_to_events(a)
		if (type(events) == str):
			print >>stderr, events
			continue
		(nMatch,nMismatch,nInsO,nInsX,nDelO,nDelX) = events

		if (writeHeader):
			print "\t".join(["line","read","mRatio","m","mm","io","ix","do","dx"])
			writeHeader = False

		if (writeWhat == "per alignment"):
			mRatio = float(nMatch) / (nMatch+nMismatch+nInsO+nInsX+nDelO+nDelX)
			mRatio = "%.3f" % mRatio
			vec = [a.lineNumber,a.qName,mRatio,nMatch,nMismatch,nInsO,nInsX,nDelO,nDelX]
			print "\t".join(map(str,vec))

		sum["m"]  += nMatch
		sum["mm"] += nMismatch
		sum["io"] += nInsO
		sum["ix"] += nInsX
		sum["do"] += nDelO
		sum["dx"] += nDelX

	sum["events"] = (sum["m"] + sum["mm"] + sum["io"] + sum["ix"] + sum["do"] + sum["dx"])

	if (alignmentNum == 0):
		print >>stderr, "WARNING: input contained no alignments"
	elif (writeWhat == "sum only"):
		alignmentNumStr = "(%d)" % alignmentNum
		mRatio = float(sum["m"]) / sum["events"]
		mRatio = "%.3f" % mRatio
		vec = ["all",alignmentNumStr,mRatio,sum["m"],sum["mm"],sum["io"],sum["ix"],sum["do"],sum["dx"]]
		print "\t".join(map(str,vec))


# sam_to_events --

def sam_to_events(a):

	# count indels in the cigar string
	#
	# note that in SAM cigar,
	#   I indicates bases in the read that are not in the reference
	#   D indicates bases in the reference that are not in the read
	# and in NCRF
	#   I indicates bases in the sequence that are not in the motif
	#   D indicates bases in the motif that are not in the sequence
	# If we assume NCRF was used to align motifs to reads (i.e. "sequence" is
	# "read") then I and D have the same meaning in both contexts

	nInsO = nInsX = nDelO = nDelX = 0
	cigarMatch = 0

	try:
		for (count,op) in cigar_ops(a.cigar):
			if (op == "I"):
				nInsO += 1
				nInsX += count-1
			elif (op == "D"):
				nDelO += 1
				nDelX += count-1
			elif (op == "M"):
				cigarMatch += count
	except ValueError:
		warning = "problem with cigar string at line %d\n%s" % (a.lineNumber,a.cigar)
		if (warnOnError): return warning
		raise ValueError, warning

	# count matches and mismatches in the md tag

	nMatch = nMismatch = 0
	mdTagDelO = mdTagDelX = 0

	try:
		for (count,op) in md_tag_ops(a.mdTag):
			if (op in "ACGT"):
				nMatch    += count
				nMismatch += 1
			elif (op == "."):
				nMatch += count
			else: # if (op.startswith("^"))
				nMatch    += count
				mdTagDelO += 1
				mdTagDelX += len(op)-2
	except ValueError:
		warning = "problem with MD tag at line %d\n%s" % (a.lineNumber,a.mdTag)
		if (warnOnError): return warning
		raise ValueError, warning

	# sanity check

	if (nMatch + nMismatch != cigarMatch):
		warning = "problem at line %d nMatch+nMismatch!=cigarMatch (%d+%d!=%d)\n%s\n%s" \
		        % (a.lineNumber,nMatch,nMismatch,cigarMatch,a.cigar,a.mdTag)
		if (warnOnError): return warning
		raise ValueError, warning

	if (mdTagDelO != nDelO):
		warning = "problem at line %d mdTagDelO!=nDelO (%d!=%d)\n%s\n%s" \
		        % (a.lineNumber,mdTagDelO,nDelO,a.cigar,a.mdTag)
		if (warnOnError): return warning
		raise ValueError, warning

	if (mdTagDelX != nDelX):
		warning = "problem at line %d mdTagDelX!=nDelX (%d!=%d)\n%s\n%s" \
		        % (a.lineNumber,mdTagDelX,nDelX,a.cigar,a.mdTag)
		if (warnOnError): return warning
		raise ValueError, warning

	return (nMatch,nMismatch,nInsO,nInsX,nDelO,nDelX)


# cigar_ops--

def cigar_ops(cigar):
	count = ""
	for ch in cigar:
		if (ch in "0123456789"):
			count += ch
			if (count == "0"): raise ValueError
		elif (ch in "HSMID"):
			if (count == ""): raise ValueError
			yield (int(count),ch)
			count = ""
		elif (ch in "=X"):
			# for historical reasons, we count these as M
			if (count == ""): raise ValueError
			yield (int(count),"M")
			count = ""
		else:
			raise ValueError

	if (count != ""): raise ValueError


# md_tag_ops--

def md_tag_ops(mdTag):
	prefix = "MD:Z:"
	if (not mdTag.startswith(prefix)): raise ValueError
	count    = ""
	deletion = ""
	for ch in mdTag[len(prefix):]:
		if (ch in "0123456789"):
			if (deletion != ""):
				if (deletion == "^"): raise ValueError
				yield (int(count),deletion)
				count = ""
				deletion = ""
			if (count == "0"): raise ValueError
			count += ch
		elif (ch in "ACGTN"): # $$$ it's not clear that N should be treated like ACGT
			if (deletion != ""):
				deletion += ch
			else:
				if (count == ""): raise ValueError
				yield (int(count),ch)
				count = ""
		elif (ch == "^"):
			if (deletion != ""): raise ValueError
			deletion = "^"
		else:
			raise ValueError

	if (deletion != ""):
		if (deletion == "^"): raise ValueError
		yield (int(count),deletion)
	elif (count != ""):
		yield (int(count),".")


# read_sam_plain--

class Alignment: pass

SAM_QNAME_COLUMN = 0
SAM_RNAME_COLUMN = 2
SAM_CIGAR_COLUMN = 5
SAM_MIN_COLUMNS  = 11

def read_sam_plain(f):
	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip()

		if (line.startswith("@")):
			continue

		fields = line.split()
		numFields = len(fields)
		assert (numFields >= SAM_MIN_COLUMNS), \
		      "not enough columns at line %d (%d, expected %d)" \
		    % (lineNumber,numFields,SAM_MIN_COLUMNS)

		a = Alignment()
		a.lineNumber = lineNumber
		a.rName      = fields[SAM_RNAME_COLUMN]
		a.qName      = fields[SAM_QNAME_COLUMN]
		a.cigar      = fields[SAM_CIGAR_COLUMN]

		if (a.rName == "*"):
			yield a
			continue

		mdFound = False
		for field in fields[SAM_MIN_COLUMNS:]:
			if (not field.startswith("MD:")): continue
			assert (not mdFound), \
			       "SAM record at line %d contains more than one MD tag\n%s" \
			     % (lineNumber,line)
			mdFound = True
			a.mdTag = field

		assert (mdFound), \
		       "SAM record at line %d lacks MD tag\n%s" \
		     % (lineNumber,line)

		yield a


if __name__ == "__main__": main()
