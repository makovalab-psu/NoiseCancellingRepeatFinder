#!/usr/bin/env python
"""
Create an msa (multiple sequence alignmnet) of each word aligned to the motif,
in the output of Noise Cancelling Repeat Finder.
"""

from sys        import argv,stdin,stdout,stderr,exit
from ncrf_parse import alignments,reverse_complement,int_with_unit


def usage(s=None):
	message = """
usage: cat <output_from_NCRF> | ncrf_msa [options]
  --head=<number>     limit the number of input alignments"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global debug

	# parse the command line

	headLimit  = None
	requireEof = True
	debug      = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--head=")):
			headLimit = int_with_unit(argVal)
		elif (arg in ["--noendmark","--noeof","--nomark"]):   # (unadvertised)
			requireEof = False
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		else:
			usage("unrecognized option: %s" % arg)

	# process the alignments

	alignmentNum = 0
	for a in alignments(stdin,requireEof):
		alignmentNum += 1 

		if (headLimit != None) and (alignmentNum > headLimit):
			print >>stderr, "limit of %d alignments reached" % headLimit
			break

		if (alignmentNum > 1): print
		print "\n".join(a.lines)

		motifText = a.motifText
		seqText   = a.seqText
		if ("noflip" in debug):
			pass
		elif (a.strand == "-") and (a.start < a.end):
			# alignment was reported in reverse complement of motif, so flip it
			motifText = reverse_complement(motifText)
			seqText   = reverse_complement(seqText)

		seqChunks = chunkify(a.motif,motifText,seqText)

		motifLen = len(a.motif)
		positionLength = [1] * motifLen
		for chunk in seqChunks:
			for (motifIx,seqNucs) in enumerate(chunk):
				if (seqNucs == None): continue
				positionLength[motifIx] = max(positionLength[motifIx],len(seqNucs))

		line = []
		for (motifIx,motifNuc) in enumerate(a.motif):
			line += [motifNuc.ljust(positionLength[motifIx],".")]
		print "# msa.query %s" % "".join(line)

		for chunk in seqChunks:
			line = []
			for (motifIx,seqNucs) in enumerate(chunk):
				if (seqNucs == None):
					line += ["."*positionLength[motifIx]]
				elif (seqNucs == a.motif[motifIx]):
					line += ["="*positionLength[motifIx]]
				else:
					line += [seqNucs.ljust(positionLength[motifIx],".")]
			print "# msa.seq   %s" % "".join(line)

	if (requireEof):
		print "# ncrf end-of-file"


# chunkify--
#	Returns a list of sublists. Each sublist is as long as the motif, consisting
#	of the string (the nt or nts) matched to each position in the motif.

def chunkify(motif,motifText,seqText):
	motifText = motifText.upper()
	motifLen  = len(motif)
	(motifPos,direction) = position_in_motif(motif,motifText)
	if ("motifpos" in debug):
		print >>stderr, "motifPos=%d%s" % (motifPos,direction)

	chunks = []
	motifIx = motifPos
	chunk = [None] * motifIx
	for (ix,textNuc) in enumerate(motifText):
		seqNuc = seqText[ix]
#...
#		print >>stderr, seqNuc,textNuc
#...
		if (textNuc == "-"):     # insertion, extra character in seqText
			if (chunk != []):    # add inserted character to latest token
				if (chunk[-1] != None): chunk[-1] += seqNuc
			elif (chunks != []): # add inserted character to final token in latest chunk
				latestChunk = chunks[-1]
				latestChunk[-1] += seqNuc
				chunks[-1] = latestChunk
			continue

		chunk += [seqNuc]        # match or deletion
#...
#		print >>stderr, chunk
#...
		motifIx += 1
		if (motifIx == motifLen):
			chunks += [chunk]
			chunk  =  []
			motifIx = 0

	if (chunk != []):
		if (len(chunk) < motifLen): chunk += [None] * (motifLen - len(chunk))
		chunks += [chunk]

	if (direction == "-"):
		for chunk in chunks:
			chunk.reverse()

	if ("chunks" in debug):
		for (chunkIx,chunk) in enumerate(chunks):
			print >>stderr, "chunk[%d] = %s" % (chunkIx,str(chunk))

	return chunks


# position_in_motif--
#	Returns the index in 'motif' at which the first position of 'text" aligns

def position_in_motif(motif,text):
	motifLen = len(motif)

	for rotIx in xrange(motifLen):
		rotMotif = motif[rotIx:] + motif[:rotIx]
		pos = 0
		motifIx = 0
		isMatch = True
		while (pos < len(text)):
			nuc = text[pos]
			pos += 1
			if (nuc == '-'): continue
			if (nuc != rotMotif[motifIx]):
				isMatch = False
				break
			motifIx += 1
			if (motifIx == motifLen):
				break
		if (isMatch):
			return (rotIx,"+")

	motif = reverse_complement(motifLen)
	for rotIx in xrange(motif):
		rotMotif = motif[rotIx:] + motif[:rotIx]
		pos = 0
		motifIx = 0
		isMatch = True
		while (pos < len(text)):
			nuc = text[pos]
			pos += 1
			if (nuc == '-'): continue
			if (nuc != rotMotif[motifIx]):
				isMatch = False
				break
			motifIx += 1
			if (motifIx == motifLen):
				break
		if (isMatch):
			return (rotIx,"-")

	# no match found, just report position 0,forward

	return (0,"+")


if __name__ == "__main__": main()
