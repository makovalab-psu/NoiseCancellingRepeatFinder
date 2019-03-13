#!/usr/bin/env python
"""
Create a mock genome (or a "read") with embedded repeat motifs.
"""

from sys        import argv,stdin,stderr,exit
from random     import seed as random_seed,randint,choice,random as unit_random,shuffle
from copy       import deepcopy
from echydna    import EchyDna,reverse_complement
from ncrf_parse import int_with_unit,parse_probability

class CatalogEntry: pass


# error profiles
#
# pacbio   (v3) mm=1.30% i=6.37% d=3.62% is from Guiblet et al (submitted)
# pacbio   (v1) mm=1.7%  i=8.9%  d=4.3%  is derived from GIAB HG002 blasr alignments
# pacbio        mm=1%    i=12%   d=2%    is from readsim-1.6 (ref [1] below)
# nanopore (v3) mm=4.60% i=3.78% d=7.73% is derived from GIAB HG002 minimap2 alignments
# nanopore (v1) mm=7.4%  i=2%    d=10.2% is derived from data from Jain et al (ref [2] below)
# nanopore      mm=3%    i=3%    d=3%    is from readsim-1.6 (ref [1] below)
#
# References
#   [1] readsim-1.6
#         https://sourceforge.net/p/readsim
#   [2] Miten Jain et al. "Nanopore sequencing and assembly of a human genome
#       with ultra-long reads."
#         https://www.biorxiv.org/content/early/2017/04/20/128835
#       Data at
#         https://github.com/nanopore-wgs-consortium/NA12878/blob/master/Genome.md
#       Later version, peer-reviewed, at
#         https://www.nature.com/articles/nbt.4060

errorProfilePacbioV3        = {"mm":0.0130, "i":0.0637, "d":0.0362}
errorProfilePacbioV1        = {"mm":0.0170, "i":0.0890, "d":0.0430}
errorProfilePacbioReadsim   = {"mm":0.01,   "i":0.12,   "d":0.02  }
errorProfileNanoporeV3      = {"mm":0.0460, "i":0.0378, "d":0.0773}
errorProfileNanoporeV1      = {"mm":0.074,  "i":0.02,   "d":0.102 }
errorProfileNanoporeReadSim = {"mm":0.03,   "i":0.03,   "d":0.03  }


def usage(s=None):
	message = """
usage: mock_motif_genome <motif> [options]
  --arrays=<filename>      specific arrays to embed; each line is of the form
                           <length> <motif>[<strand>]; <length> is in bp;
                           <strand> is ignored;
                           this cannot be used with any command line <motif>s,
                           nor with --repeats, --lengths, --motif:neighbor, or
                           --motif:mixture
  <motif>                  (cumulative) motif to embed
  --name=<string>          read name
  --length=<bp>       (L=) read length; if this is absent, the repeats will
                           comprise the entire read (no fill DNA)
  --length=<pctg>%    (L=) read length as a percentage of the total repeat
                           length; pctg should be more than 100.
  --length=+<bp>      (L=) read length as delta above the total repeat length
  --repeats=<number>  (N=) number of repeats to embed
                           (default is 1)
  --motif:neighbor=<prob>  probability of embedding motifs that have edit
                           distance 1 from the specified motifs
                           (default is 0)
  --motif:mixture=<prob>   probability of embedding motifs that are 50/50
                           mixtures of a specified motif and one that has edit
                           distance 1
                           (default is 0)
  --lengths=<file>         file containing the repeat length distribution, one
                           length per line; if this is absent, we'll read the
                           repeat length distribution from stdin; if <file>
                           contains "{motif}", we'll use a separate distribution
                           for each motif
  --minfill=<bp>      (F=) minimum fill (random sequence) between repeats
  --errors=pacbio          simulate pacbio error profile
                           (by default, errors are not simulated)
  --errors=nanopore        simulate nanopore error profile
  --errors=<pctg>%         simulate simple error profile with the given rate,
                           with mismatch, insertion, and deletion each of equal
                           rates
  --errors=<spec>          simulate error profile with the given spec; <spec>
                           looks like this:
                              mm:1%,i:12%,d:2%
  --catalog=<file>         file to write a catalog of the embedded repeats to
  --wrap=<length>          number of nucleotides per line; 0 means single line
  --seed=<string>          set random seed

Note that if the program is run twice with the same seed, one run with errors,
one without, the same pre-error sequence is generated.  Having the output of
both those runs can be useful."""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():

	# parse the command line

	arraysFilename   = None
	motifs          = []
	sequenceName    = None
	sequenceLen     = 0
	numRepeats      = None
	genNeighbors    = 0.0
	genMixture      = 0.0
	lengthsFilename = None
	minFill         = None
	errorProfile    = None
	catalogFilename = None
	wrapLength      = 100

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--arrays=")):
			arraysFilename = argVal
		elif (arg.startswith("--name=")):
			sequenceName = argVal
		elif (arg.startswith("--length=")) or (arg.startswith("--len=")) or (arg.startswith("L=")):
			if (argVal.endswith("%")):
				sequenceLen = float(argVal[:-1]) / 100.0
				assert (sequenceLen >= 1.0)
				sequenceLen = ("%",sequenceLen)
			elif (argVal.startswith("+")):
				sequenceLen = int_with_unit(argVal[1:])
				assert (sequenceLen >= 0)
				sequenceLen = ("+",sequenceLen)
			else:
				sequenceLen = int_with_unit(argVal)
				assert (sequenceLen >= 0)
		elif (arg.startswith("--repeats=")) or (arg.startswith("N=")):
			numRepeats = int_with_unit(argVal)
			assert (numRepeats > 0)
		elif (arg.startswith("--motif:neighbor=")):
			genNeighbors = parse_probability(argVal)
		elif (arg.startswith("--motif:mixture=")):
			genMixture = parse_probability(argVal)
		elif (arg.startswith("--lengths=")):
			lengthsFilename = argVal
		elif (arg.startswith("--minfill=")) or (arg.startswith("F=")):
			minFill = int(argVal)
			if (minFill < 0):
				print >>stderr, "WARNING: \"%s\" interpreted as no minimum fill" % argVal
				minFill = None
			if (minFill == 0):
				minFill = None
		elif (arg.startswith("--errors=")):
			errorProfile = None
			if (argVal in ["pacbio","pacbio.v3","pacbio.GIAB","pacbio.giab"]):
				errorProfile = errorProfilePacbioV3
			elif (argVal == "pacbio.v2"):  # for historical reasons, v2 is an alias for v3
				errorProfile = errorProfilePacbioV3
			elif (argVal in ["pacbio.v1","pacbio.Guiblet","pacbio.guiblet"]):
				errorProfile = errorProfilePacbioV1
			elif (argVal in ["pacbio.readsim"]):
				errorProfile = errorProfilePacbioReadsim
			elif (argVal in ["nanopore","nanopore.v3","nanopore.GIAB","nanopore.giab"]):
				errorProfile = errorProfileNanoporeV3
			elif (argVal == "nanopore.v2"):  # for historical reasons, v2 is an alias for v3
				errorProfile = errorProfileNanoporeV3
			elif (argVal in ["nanopore.v1","nanopore.Jain","nanopore.jain"]):
				errorProfile = errorProfileNanoporeV1
			elif (argVal in ["nanopore.readsim"]):
				errorProfile = errorProfileNanoporeReadSim
			elif (":" in argVal):
				try:
					errorProfile = parse_error_spec(argVal)
				except ValueError:
					pass
			else:
				p = parse_probability(argVal)
				errorProfile = {"mm":p, "i":p, "d":p }
			if (errorProfile == None):
				usage("\"%s\" is not a valid error spec" % argVal)
			subProb       = errorProfile["mm"]
			insOpenProb   = errorProfile["i"]
			delOpenProb   = errorProfile["d"]
			insExtendProb = delExtendProb = 0.0
		elif (arg.startswith("--catalog=")):
			catalogFilename = argVal
		elif (arg.startswith("--wrap=")):
			wrapLength = int(argVal)
			if (wrapLength <= 0): wrapLength = None
		elif (arg.startswith("--seed=")):
			# nota bene: if the seed is a number, use it as a number, since
			#            string seeds can produce different sequences on
			#            different versions/builds of python
			seed = argVal
			try:
				seed = int(seed)
			except ValueError:
				try:               seed = float(seed)
				except ValueError: pass
			random_seed(seed)
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		elif (is_nucleotide_string(arg)):
			motifs += [arg.upper()]
		else:
			usage("unrecognized option: %s" % arg)

	if (arraysFilename != None):
		if (motifs != []):
			usage("command line <motif>s cannot be used with --arrays")
		if (numRepeats != None):
			usage("--repeats cannot be used with --arrays")
		if (lengthsFilename != None):
			usage("--lengths cannot be used with --arrays")
		if (genNeighbors != 0.0):
			usage("--motif:neighbor cannot be used with --arrays")
		if (genMixture != 0.0):
			usage("--motif:mixture cannot be used with --arrays")
	elif (motifs == []):
		usage("you have to give me at least one motif")

	if (numRepeats == None) and (arraysFilename != None):
		numRepeats = 1
	
	# read the arrays file, if we have one

	repeatLengths = {}
	haveSpecificArrays = False

	if (arraysFilename != None):
		haveSpecificArrays = True
		f = file(arraysFilename,"rt")
		numRepeats = 0
		for (length,motif,_) in read_arrays(f,arraysFilename):
			numRepeats += 1
			motifs += [(motif)]
			if (motif not in repeatLengths): repeatLengths[motif] =  [length]
			else:                            repeatLengths[motif] += [length]
		f.close()

		if (motifs == []):
			usage("array file \"%s\" contains no arrays" % arraysFilename)

	# read the lengths file

	if (repeatLengths == {}):
		if (lengthsFilename == None):
			lengths = read_integers(stdin)
			for motif in motifs:
				repeatLengths[motif] = lengths
		elif ("{motif}" not in lengthsFilename):
			f = file(lengthsFilename,"rt")
			lengths = read_integers(f,lengthsFilename)
			f.close()
			for motif in motifs:
				repeatLengths[motif] = lengths
		else:
			for motif in motifs:
				motifLengthsFilename = lengthsFilename.replace("{motif}",motif)
				f = file(motifLengthsFilename,"rt")
				lengths = read_integers(f,motifLengthsFilename)
				f.close()
				repeatLengths[motif] = lengths

	# generate the number and type of motifs we'll embed
	#
	# note: to satisfy the requirement that the same seed generates the same
	#       pre-error sequence, we should have no variance in the use of the
	#       PRNG until after we've generated that sequence; see "point A" below

	embeddings = []

	if (haveSpecificArrays):
		for motif in motifs:
			for length in repeatLengths[motif]:
				strand = choice(["+","-"])
				offset = choice(xrange(len(motif)))
				embeddings += [(1.0,motif,motif,strand,offset,length)]
		shuffle(embeddings)
	else:
		for _ in xrange(numRepeats):
			motif = choice(motifs)
			length = choice(repeatLengths[motif])
			u = unit_random()
			if (genNeighbors > 0) and (u < genNeighbors):
				motif = motif_neighbor(motif)
				(mix,motif2) = (1.0,motif)
			elif (genMixture > 0) and (u < genNeighbors+genMixture):
				(mix,motif2) = (0.5,motif_neighbor(motif))
			else:
				(mix,motif2) = (1.0,motif)
			strand = choice(["+","-"])
			offset = choice(xrange(len(motif)))
			embeddings += [(mix,motif,motif2,strand,offset,length)]

	totalRepeatBp = sum([length for (_,_,_,_,_,length) in embeddings])

	# assign each repeat a position within the "fill" sequence;  note that we
	# might have more than one repeat assigned to the same position, in which
	# case they will be back-to-back with no fill between them

	if (type(sequenceLen) == tuple):
		(op,sequenceLen) = sequenceLen
		if (op == "%"):
			sequenceLen = int(round(totalRepeatBp*sequenceLen))
		else: # if (op == "+"):
			sequenceLen = totalRepeatBp + sequenceLen

	if (totalRepeatBp > sequenceLen):
		fillBp = 0
		if (sequenceLen > 0):
			print >>stderr, "WARNING: length of embedded repeats (%d) exceeds specified" % totalRepeatBp
			print >>stderr, "         sequence length (%d); there will be no fill DNA"   % sequenceLen
	elif (minFill != None):
		fillBp = sequenceLen - totalRepeatBp
		totalMinFill = (numRepeats+1) * minFill
		if (totalMinFill > fillBp):
			print >>stderr, "WARNING: minimum fill of %d cannot be achieved"           % minFill
			print >>stderr, "         total minimum fill (%d) exceeds total fill (%d)" % (totalMinFill,fillBp)
			minFill = fillBp / (numRepeats+1)
		fillBp -= minFill * (numRepeats+1)
	else:
		fillBp = sequenceLen - totalRepeatBp

	positions = [randint(0,fillBp) for _ in xrange(numRepeats)]
	positions.sort()

	if (minFill != None):
		fillBp += minFill * (numRepeats+1)
		for rptNum in xrange(numRepeats):
			positions[rptNum] += (rptNum+1) * minFill

	# generate the sequence

	catalog = None
	if (catalogFilename != None):
		catalog = []

	fillSeq = str(EchyDna(fillBp))
	seq = []
	seqPos  = 0
	prevEnd = 0
	fillPos = 0
	for (ix,pos) in enumerate(positions):
		if (fillPos < pos):
			seq += [fillSeq[fillPos:pos]]
			seqPos  += pos - fillPos
			fillPos =  pos

		(mix,motif,motif2,strand,offset,length) = embeddings[ix]
		if (catalog != None):
			c = CatalogEntry()
			c.start        = seqPos
			c.end          = seqPos+length
			c.mix          = mix
			c.motif        = motif
			c.motif2       = motif2
			c.strand       = strand
			c.repeatLength = length
			c.offset       = offset
			catalog += [c]

		enoughCopies = (length+offset+len(motif)-1) / len(motif)
		if (strand == "-"): motif = reverse_complement(motif)

		if (mix >= 1.0):
			repeat = motif * enoughCopies
		else:
			repeat = []
			for _ in xrange(enoughCopies):
				if (unit_random() < mix): repeat += [motif]
				else:                     repeat += [motif2]
			repeat = "".join(repeat)

		seq += repeat[offset:offset+length]
		seqPos += length
		prevEnd = seqPos

	if (fillPos < fillBp):
		seq += [fillSeq[fillPos:fillBp]]

	seq = "".join(seq)

	#=== point A: it's now safe to make additional use of the PRNG ===

	# apply error profile

	events = profile = None
	if (argVal in ["pacbio","pacbio.v3","pacbio.GIAB","pacbio.giab"]):
		errorProfile = errorProfilePacbioV3
	elif (argVal == "pacbio.v2"):  # for historical reasons, v2 is an alias for v3
		errorProfile = errorProfilePacbioV3
	elif (argVal in ["pacbio.v1","pacbio.Guiblet","pacbio.guiblet"]):
		errorProfile = errorProfilePacbioV1
	elif (argVal in ["pacbio.readsim"]):
		errorProfile = errorProfilePacbioReadsim
	elif (argVal in ["nanopore","nanopore.v3","nanopore.GIAB","nanopore.giab"]):
		errorProfile = errorProfileNanoporeV3
	elif (argVal == "nanopore.v2"):  # for historical reasons, v2 is an alias for v3
		errorProfile = errorProfileNanoporeV3
	elif (argVal in ["nanopore.v1","nanopore.Jain","nanopore.jain"]):
		errorProfile = errorProfileNanoporeV1
	elif (argVal in ["nanopore.readsim"]):
		errorProfile = errorProfileNanoporeReadSim
	elif (type(errorProfile) == float):
		eRate = errorProfile / 3.0;
		profile = {"mm":eRate, "i":eRate, "d":eRate }
	elif (type(errorProfile) == dict):
		profile = dict(errorProfile)

	if (profile != None):
		print >>stderr, "(applying error profile mm=%.2f%% i=%.2f%% d=%.2f%%)" \
		             % (100*profile["mm"],100*profile["i"],100*profile["d"])
		(seq,catalog,events) = apply_errors(profile,seq,catalog)

	# write the sequence

	if (sequenceName != None):
		print ">%s" % sequenceName

	if (wrapLength == None):
		print seq
	else:
		for i in range(0,len(seq),wrapLength):
			print seq[i:i+wrapLength]

	# write the catalog

	if (catalogFilename != None):
		catalogF = file(catalogFilename,"wt")
		if (sequenceName in [None,""]): seqNameForCatalog = "seq"
		else:                           seqNameForCatalog = sequenceName

		if (events == None):
			print >>catalogF, "#%s\t%s\t%s\t%s\t%s\t%s\t%s" \
			                % ("chrom","start","end","motif","rptLen","len","fill")
		else:
			print >>catalogF, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" \
			                % ("chrom","start","end","motif","rptLen","len","fill",
			                   "mRatio","m","mm","i","d")

		prevEnd = 0
		for (catIx,c) in enumerate(catalog):
			motifStr = c.motif
			if (c.mix < 1.0): motifStr += "," + c.motif2
			motifStr += ".%s%s" % (c.offset,c.strand)
			if (events == None):
				print >>catalogF, "%s\t%s\t%s\t%s\t%s\t%s\t%s" \
				                % (seqNameForCatalog,c.start,c.end,motifStr,
				                   c.repeatLength,c.end-c.start,c.start-prevEnd)
			else:
				if (catIx in events):
					(m,mm,i,d) = events[catIx]
					mRatio = "%.1f%%" % (100.0*m/(m+mm+i+d))
				else:
					mRatio = m = mm = i = d = "NA"
				print >>catalogF, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" \
				                % (seqNameForCatalog,c.start,c.end,motifStr,
				                   c.repeatLength,c.end-c.start,c.start-prevEnd,
				                   mRatio,m,mm,i,d)
			prevEnd = c.end

		catalogF.close()


# apply_errors--

def apply_errors(profile,seq,catalog):
	pMm = 0.01
	pI  = 0.12
	pD  = 0.02
	pMm = profile["mm"]
	pI  = profile["i"]
	pD  = profile["d"]

	if (catalog == None):
		newCatalog = None
		events     = None
	else:
		newCatalog = deepcopy(catalog)
		startToIx = {}
		endToIx   = {}
		for (catIx,c) in enumerate(catalog):
			startToIx[c.start] = catIx
			endToIx  [c.end  ] = catIx
			c.start = c.end = None  # (so we'll know if we failed to change them)
		events = {}

	newSeq = []
	newPos = 0
	m = None
	for pos in xrange(len(seq)+1):
		if (newCatalog != None):
			# nota bene: we assume catalog intervals don't overlap, but they
			#            may abut
			if (pos in endToIx):
				catIx = endToIx[pos]
				newCatalog[catIx].end = newPos
				events[catIx] = (m,mm,i,d)
				m = None
			if (pos in startToIx):
				catIx = startToIx[pos]
				newCatalog[catIx].start = newPos
				m = mm = i = d = 0

		if (pos == len(seq)):
			break

		nuc = seq[pos]
		r = unit_random()
		if (r < pMm):
			newSeq += [choice(mismatchLookup[nuc])]
			newPos += 1
			if (m != None): mm += 1
		elif (r < pMm+pI):
			newSeq += [choice("ACGT")+nuc]
			newPos += 2
			if (m != None): i += 1
		elif (r < pMm+pI+pD):
			if (m != None): d += 1
		else:
			newSeq += [nuc]
			newPos += 1
			if (m != None): m += 1

	return ("".join(newSeq),newCatalog,events)


# motif_neighbor--
#	Choose a motif with edit distance 1 from a specified motif.

def motif_neighbor(motif):
	newMotif = [nuc for nuc in motif]
	motifLen = len(motif)
	possibleEdits = 3*motifLen + (motifLen+1) + motifLen
	r = randint(0,possibleEdits-1)
	if (r < 3*motifLen):
		# substitution
		ix    = r / 3
		subIx = r % 3
		nuc = motif[ix]
		newMotif[ix] = mismatchLookup[nuc][subIx]
	elif (r-(3*motifLen) < motifLen+1):
		# insertion
		r -= 3*motifLen
		newMotif.insert(r,choice("ACGT"))
	else:
	# deletion
		r -= 3*motifLen + (motifLen+1)
		del newMotif[r]

	return "".join(newMotif)


# read_arrays--
#	Yield the next (length,motif) pair from a file

def read_arrays(f,filename=None):
	if (filename == None): filename = "input"

	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip("\n")
		if (line == ""): continue

		try:
			fields = line.split(" ")
			if (len(fields) != 2): raise ValueError
			length = int(fields[0])
			motif  = fields[1]
			strand = None
			if (motif.endswith("+")) or (motif.endswith("-")):
				(motif,strand) = (motif[:-1],motif[-1])
			if (motif == ""): raise ValueError
			if (length < len(motif)): raise ValueError
		except ValueError:
			raise ValueError, \
			      "bad array (line %d in %s):\n%s" \
			    % (lineNumber,filename,line)

		yield (length,motif,strand)


# read_integers--
#	Read integers from a file, one per line

def read_integers(f,filename=None):
	if (filename == None): filename = "input"

	numbers = []

	lineNumber = 0
	for line in f:
		lineNumber += 1
		line = line.strip("\n")
		if (line == ""): continue

		try:
			numbers += [int(line)]
		except ValueError:
			raise ValueError, \
			      "bad number (line %d in %s):\n%s" \
			    % (lineNumber,filename,line)

	return numbers


# is_nucleotide_string--

def is_nucleotide_string(s):
	for nuc in s:
		if (nuc not in "ACGTacgt"): return False
	return True

# mismatchLookup

mismatchLookup = {"A":[    "C","G","T"],
                  "C":["A",    "G","T"],
	              "G":["A","C",    "T"],
	              "T":["A","C","G"    ]}


# parse_error_spec--
#	A typical spec looks like this: mm:1%,i:12%,d:2%

def parse_error_spec(spec):
	fields = spec.split(",")
	if (len(fields) != 3): raise ValueError

	errorProfile = {}
	for field in fields:
		field = field.strip()
		if (":" not in field): raise ValueError
		(event,p) = field.split(":",1)
		if (event not in ["mm","i","d"]): raise ValueError
		if (event in errorProfile): raise ValueError  # (event specified twice)
		errorProfile[event] = parse_probability(p)

	return errorProfile


if __name__ == "__main__": main()
