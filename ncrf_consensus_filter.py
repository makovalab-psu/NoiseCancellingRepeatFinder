#!/usr/bin/env python
"""
Filter Noise Cancelling Repeat Finder alignments, discarding alignments that
has a consensus different than the motif unit.
"""

from sys        import argv,stdin,stdout,stderr,exit
from math       import log,exp
from ncrf_parse import alignments,reverse_complement,parse_probability,int_with_unit


def usage(s=None):
	message = """
usage: ncrf_cat <output_from_NCRF> | ncrf_consensus_filter [options]
  --consensusonly     just report the consensus motif(s) for each alignment,
                      instead of filtering; these are added to the alignment
                      file with a "# consensus" tag
  --head=<number>     limit the number of input alignments"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global debug

	# parse the command line

	filterToKeep    = "consensus"
	reportConsensus = False
	reportMsa       = False
	pThreshold      = 0.75  # (see derive_consensuses)
	winnerThreshold = 0.65  # (see derive_consensuses)
	headLimit       = None
	requireEof      = True
	debug           = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg == "--nonconsensus"):   # (unadvertised)
			filterToKeep    = "non consensus"
			reportMsa       = False
			reportConsensus = True
		elif (arg == "--nonconsensus,msa"):   # (unadvertised)
			filterToKeep    = "non consensus"
			reportMsa       = True
			reportConsensus = True
		elif (arg == "--consensusonly"):
			filterToKeep    = "no filter"
			reportMsa       = False
			reportConsensus = True
		elif (arg == "--filter,consensus"):   # (unadvertised)
			filterToKeep    = "consensus"
			reportMsa       = False
			reportConsensus = True
		elif (arg == "--msa"):   # (unadvertised)
			filterToKeep    = "no filter"
			reportMsa       = True
			reportConsensus = True
		elif (arg.startswith("--winner=")):   # (unadvertised)
			winnerThreshold = parse_probability(argVal)
		elif (arg.startswith("--head=")):
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
	alignmentsWritten = 0
	for a in alignments(stdin,requireEof):
		alignmentNum += 1 

		if (headLimit != None) and (alignmentNum > headLimit):
			print >>stderr, "limit of %d alignments reached" % headLimit
			break

		motifText = a.motifText
		seqText   = a.seqText
		if ("noflip" in debug):
			pass
		elif (a.strand == "-") and (a.start < a.end):
			# alignment was reported in reverse complement of motif, so flip it
			motifText = reverse_complement(motifText)
			seqText   = reverse_complement(seqText)

		# derive consensus; this can actually be a set of more than one
		# candidate

		seqChunks   = chunkify(a.motif,motifText,seqText)


		if ("consensus" in debug):
			print >>stderr
			print >>stderr, "%d score=%d" % (a.lineNumber,a.score)

		consensuses = derive_consensuses(seqChunks,
		                                 expectedMotif=a.motif,pThreshold=pThreshold,
		                                 winnerThreshold=winnerThreshold)
		consensuses = list(consensuses)

		# discard the alignment if it meets the filtering criterion (if there
		# is any such criterion)

		if (filterToKeep == "consensus"):
			if (a.motif not in consensuses): continue  # (discard it)
		elif (filterToKeep == "non consensus"):
			if (a.motif in consensuses): continue  # (discard it)
		else: # if (filterToKeep == "no filter"):
			pass

		# copy the (unfiltered) alignment to the output

		if (alignmentsWritten > 0): print
		alignmentsWritten += 1

		print "\n".join(a.lines)

		# report the consensus and msa, if we're supposed to

		if (reportConsensus):
			if (consensuses == []):
				print "# consensus (none)"
			else:
				for word in consensuses:
					print "# consensus %s" % word

		if (reportMsa):
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
		if (textNuc == "-"):     # insertion, extra character in seqText
			if (chunk != []):    # add inserted character to latest token
				if (chunk[-1] != None): chunk[-1] += seqNuc
			elif (chunks != []): # add inserted character to final token in latest chunk
				latestChunk = chunks[-1]
				latestChunk[-1] += seqNuc
				chunks[-1] = latestChunk
			continue

		chunk += [seqNuc]        # match or deletion
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


# derive_consensuses--
#	Yields a series of potential consensus motifs for a given alignment. The
#	input is a list of sublists (as created by chunkify).  Each sublist is as
#	long as the motif, consisting of the string (the nt or nts) matched to each
#	position in the motif.

def derive_consensuses(seqChunks,
                       expectedMotif=None,pThreshold=0.75,
                       winnerThreshold=0):
	if (seqChunks == []): return
	motifLen = len(seqChunks[0])

	# count the number of times each motif position is observed to match a
	# particular 'token'

	ixToTokens = {}
	for motifIx in xrange(motifLen):
		ixToTokens[motifIx] = {}

	for chunk in seqChunks:
		for (motifIx,seqNucs) in enumerate(chunk):
			if (seqNucs == None): seqNucs = ""
			if (seqNucs not in ixToTokens[motifIx]):
				ixToTokens[motifIx][seqNucs] = 1
			else:
				ixToTokens[motifIx][seqNucs] += 1

	# at each position, sort the tokens from most-observed to least-observed,
	# then reduce them to best (or almost tied for best)

	ixToWinners = {}
	for motifIx in xrange(motifLen):
		tokensSeen = ixToTokens[motifIx]
		ixToTokens[motifIx] = [(ixToTokens[motifIx][seqNucs],seqNucs) for seqNucs in tokensSeen]
		ixToTokens[motifIx].sort()
		ixToTokens[motifIx].reverse()
		tokensCount = sum([count for (count,_) in ixToTokens[motifIx]])

		if ("consensus" in debug):
			s = []
			for (count,seqNucs) in ixToTokens[motifIx]:
				s += ["%d:\"%s\"" % (count,seqNucs)]
			print >>stderr, "# [%d] tokensCount=%d %s" % \
			                (motifIx,tokensCount," ".join(s))

		ixToWinners[motifIx] = [ixToTokens[motifIx][ix]
		                          for (ix,(count,_)) in enumerate(ixToTokens[motifIx])
		                          if (count >= winnerThreshold*tokensCount)]

	# determine whether the 'expected' consensus is sufficiently represented
	# in the observations

	motifsReported = set()

	if (expectedMotif != None):
		logOfP = 0.0
		for (motifIx,expectedNuc) in enumerate(expectedMotif):
			tokensSeen = {seqNucs:count for (count,seqNucs) in ixToTokens[motifIx]}
			if (expectedNuc not in tokensSeen):
				logOfP = None
				if ("consensus" in debug):
					print >>stderr, "# %s p%d=%d/%d=%.3f" % \
					                (expectedNuc,motifIx,0,float(tokensCount),0)
				break

			tokensCount = sum([tokensSeen[count] for count in tokensSeen])
			logOfP += log(tokensSeen[expectedNuc] / float(tokensCount))

			if ("consensus" in debug):
				print >>stderr, "# %s p%d=%d/%d=%.3f log=%.3f" % \
				                (expectedNuc,
				                 motifIx,tokensSeen[expectedNuc],float(tokensCount),
				                 tokensSeen[expectedNuc] / float(tokensCount),
				                 log(tokensSeen[expectedNuc] / float(tokensCount)))

	if (logOfP != None):
		if ("consensus" in debug):
			print >>stderr, "# logOfP=%.3f avg=%.3f avgP=%.3f" % \
			                (logOfP,logOfP/motifLen,exp(logOfP/motifLen))
		if (logOfP >= motifLen*log(pThreshold)):
			# exp(logOfP/motifLen) >= pThreshold
			yield expectedMotif
			motifsReported.add(expectedMotif)

	# generate the possible consensus motifs
	# $$$ change this to include the almost-tied ones

	if (winnerThreshold != None):
		motif = []
		for motifIx in xrange(motifLen):
			tokensSeen = ixToWinners[motifIx]
			if (tokensSeen == []):   # (no consensus can be formed, because nothing
				return               #  .. in this column is a clear winner)
			(_,bestSeqNucs) = tokensSeen[0]
			motif += [bestSeqNucs]
		motif = "".join(motif)

		if (motif not in motifsReported):
			yield motif
			motifsReported.add(motif)


if __name__ == "__main__": main()
