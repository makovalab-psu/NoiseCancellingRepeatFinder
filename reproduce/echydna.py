#!/usr/bin/env python

from math      import log,floor
from string    import maketrans
from itertools import izip
from random    import Random,random as unit_random,choice,shuffle


# EchyDna--
#	Class to deal with strings of DNA.

class EchyDna(object):
	background="ACGT"    # default background distribution for generated sequence

	def __init__(self,lengthOrString,name=None,background=None,seed=None):
		"""
		if 1st arg is a length:   generate random sequence of length N
		if 1st arg is a string:   use the string as the sequence
		if 1st arg is a EchyDna:  use the EchyDna's seq as the sequence
		if 1st arg is a list:     concatenate as the sequence; list elements
		                          must be strings or EchyDnas
		name:                     string to use for name in a few cases (such as
		                          .. output to a fasta file)
		background:               letters to choose from; multiplicities are OK,
		                          .. eg "AAACCGGTTT"
		if seed is None:          use the standard prng, unseeded
		if seed is a Random:      use it as the prng
		if seed is string or int: use a fresh prng with that seed
		"""

		if (name == None): self.name = ""
		else:              self.name = name

		# if the caller has given us a list of sequences, concatenate it

		if (isinstance(lengthOrString,list)):
			seq = []
			for item in lengthOrString:
				seq += [str(item)]
			self.seq = "".join(seq)
			return

		# if the caller has given us a sequence, use it

		if (isinstance(lengthOrString,EchyDna)) or (isinstance(lengthOrString,str)):
			self.seq = str(lengthOrString)
			return

		# otherwise, generate a random sequence

		if (isinstance(lengthOrString,int)):
			length = lengthOrString
		elif (isinstance(lengthOrString,float)):
			length = int(floor(lengthOrString))
		else:
			raise ValueError("first argument is neither a string nor valid length")

		(_,chooser,_) = resolve_prng(seed)

		if (length <= 0):
			self.seq = ""
		else:
			if (background == None): background = EchyDna.background
			self.seq = "".join([chooser(background) for _ in xrange(length)])

	@classmethod
	def from_fasta(EchyDna,filename,requireFastaHeader=False):
		"""generator; yields each sequence read from the file"""
		f = file(filename,"rt")
		seqName = None
		seqNucs = None
		for line in f:
			line = line.strip()
			if (line.startswith(">")):
				if (seqName != None):
					yield EchyDna("".join(seqNucs),name=seqName)
				seqName = line[1:].strip()
				seqNucs = []
			elif (seqName == None):
				if (requireFastaHeader):
					raise ValueError("The firsrt sequence in \"" + filename + "\" has no fasta header")
				seqName = ""
				seqNucs = [line]
			else:
				seqNucs += [line]
		if (seqName != None):
			yield EchyDna("".join(seqNucs),name=seqName)
		f.close()

	def fasta(self,filename=None,name=None,wrap=None,append=False):
		"""
		if filename is a string, we'll write a file with that name and close it
		if it is a file object, we'll write to it and leave it open
		append is only relevant if filename is a string
		"""
		if (name == None): name = self.name
		fasta = [">" + name]
		if (wrap == None) or (wrap == 0):
			fasta += [self.seq]
		else:
			if (not isinstance(wrap,int)) or (wrap < 0):
				raise ValueError("wrap argument is not a valid count")
			for ix in range(0,len(self.seq),wrap):
				fasta += [self.seq[ix:ix+wrap]]
		fasta = "\n".join(fasta)
		if (filename == ""): filename = None
		if (filename == None):
			return fasta
		closeFile = False
		if (isinstance(filename,str)):
			if (append): f = file(filename,"at")
			else:        f = file(filename,"wt")
			closeFile = True
		elif (isinstance(filename,file)):
			f = filename
		else:
			raise ValueError("filename argument is not a valid string")
		print >>f, fasta
		if (closeFile):
			f.close()
		return fasta

	def __str__(self):
		return self.seq

	def __repr__(self):
		return "EchyDna(" + self.seq + ")"

	def __len__(self):
		return len(self.seq)

	def __getitem__(self,slice):
		return EchyDna(self.seq.__getitem__(slice))

	def __add__(self,other):
		"""a+b concatentates a and b"""
		if (not isinstance(other,EchyDna)) and (not isinstance(other,str)):
			raise ValueError("operand is not a valid string type")
		return EchyDna(self.seq+str(other))

	def __iadd__(self,other):
		"""a+=b concatentates b to a"""
		if (not isinstance(other,EchyDna)) and (not isinstance(other,str)):
			raise ValueError("operand is not a valid string type")
		self.seq += str(other)
		return self

	def __sub__(self,other):
		"""a-b concatentates a and reverse-complement of b"""
		if (not isinstance(other,EchyDna)) and (not isinstance(other,str)):
			raise ValueError("operand is not a valid string type")
		return EchyDna(self.seq+reverse_complement(str(other)))

	def __isub__(self,other):
		"""a-=b concatentates reverse-complement of b to a"""
		if (not isinstance(other,EchyDna)) and (not isinstance(other,str)):
			raise ValueError("operand is not a valid string type")
		self.seq += reverse_complement(str(other))
		return self

	def __mul__(self,other):
		"""a*n concatentates n copies of a; if n is negative the reverse complement of a is used"""
		if (isinstance(other,int)):
			pass
		elif (isinstance(other,float)):
			other = int(floor(other))
		else:
			raise ValueError("operand is not a valid count")
		if (other >= 0):
			return EchyDna(self.seq * other)
		else:
			return EchyDna(reverse_complement(self.seq) * -other)

	def __imul__(self,other):
		"""a*=n replaces a with n-1 copies of itself; if n is negative the reverse complement of a is used"""
		if (isinstance(other,int)):
			pass
		elif (isinstance(other,float)):
			other = int(floor(other))
		else:
			raise ValueError("operand is not a valid count")
		if (other >= 0):
			self.seq *= other
		else:
			self.seq = reverse_complement(self.seq) * -other
		return self

	def __lshift__(self,other):
		"""a<<n discards n nts on the left (abcde<<1 is bcde)"""
		if (isinstance(other,int)):
			n = other
		elif (isinstance(other,float)):
			n = int(floor(other))
		else:
			raise ValueError("operand is not a valid count")
		if (n <= 0): return EchyDna(self.seq)
		else:        return EchyDna(self.seq[n:])

	def __rshift__(self,other):
		"""a>>n discards n nts on the right (abcde>>1 is abcd)"""
		if (isinstance(other,int)):
			n = other
		elif (isinstance(other,float)):
			n = int(floor(other))
		else:
			raise ValueError("operand is not a valid count")
		if (n <= 0): return EchyDna(self.seq)
		else:        return EchyDna(self.seq[:-n])

	def __neg__(self):
		"""-a gives reverse-complement of a"""
		return EchyDna(reverse_complement(self.seq))

	def __invert__(self):
		"""~a gives un-reversed complement of a"""
		return EchyDna(forward_complement(self.seq))

	def __iter__(self):
		self.iterIx = 0
		return self

	def next(self):  # for python2;  for python3, this should be __next__
		if (self.iterIx == len(self.seq)):
			raise StopIteration
		self.iterIx += 1
		return self.seq[self.iterIx-1]

	def reverse(self):
		"""a.reverse gives un-complemented reverse of a"""
		return EchyDna(self.seq[::-1])

	def upper(self):
		return EchyDna(self.seq.upper())

	def lower(self):
		return EchyDna(self.seq.lower())

	def __eq__(self,other):
		if (not isinstance(other,EchyDna)) and (not isinstance(other,str)):
			raise ValueError("operand is not a valid string type")
		return (self.seq == str(other))

	def __ne__(self,other):
		if (not isinstance(other,EchyDna)) and (not isinstance(other,str)):
			raise ValueError("operand is not a valid string type")
		return (self.seq != str(other))

	def __lt__(self,other):
		if (not isinstance(other,EchyDna)) and (not isinstance(other,str)):
			raise ValueError("operand is not a valid string type")
		return (self.seq < str(other))

	def __le__(self,other):
		if (not isinstance(other,EchyDna)) and (not isinstance(other,str)):
			raise ValueError("operand is not a valid string type")
		return (self.seq <= str(other))

	def __gt__(self,other):
		if (not isinstance(other,EchyDna)) and (not isinstance(other,str)):
			raise ValueError("operand is not a valid string type")
		return (self.seq > str(other))

	def __ge__(self,other):
		if (not isinstance(other,EchyDna)) and (not isinstance(other,str)):
			raise ValueError("operand is not a valid string type")
		return (self.seq >= str(other))

	def __contains__(self,other):
		"""b in a is true iff b is a subtring of a; *is* sensitive to upper/lowercase"""
		if (not isinstance(other,EchyDna)) and (not isinstance(other,str)):
			raise ValueError("operand is not a valid string type")
		return (str(other) in self.seq)

	def replace(self,old,new,maxReplace=None):
		if (maxReplace != None):
			return EchyDna(self.seq.replace(old,new,maxReplace))
		else:
			return EchyDna(self.seq.replace(old,new))

	def shuffle(self,seed=None):
		(prng,_,_) = resolve_prng(seed)
		seq = list(self.seq)
		if (prng == None): shuffle(seq)
		else:              prng.shuffle(seq)
		return EchyDna("".join(seq))

	def kmers(self,k):
		for ix in xrange(k,len(self.seq)+1):
			yield self.seq[ix-k:ix]

	def mutate(self,subRate=None,indelRate=None,indels=None,background=None,seed=None):
		"""
		subRate:    the per-base probability of a substitution
		indelRate:  the per-base probability of an insertion or deletion event
		indels:     generally, a function that maps a random number (in 0..1) to
		            .. an insertion or deletion length; insertions are positive,
		            ..   deletions negative
		            .. indels can also be a list of lengths, to be chosen
		            ..   with uniform probability
		            .. indels can also be an integer N, representing a uniform
		            ..   distribution of -N..+N without 0
		            .. if this is not provided, indels of length 1 are used
		background: letters to choose from; multiplicities are OK;
		seed:       same as for __init__
		"""

		# allow the caller to pass indel information as a pair in either
		# indel argument

		if (indelRate != None) \
		  and (indels == None) \
		  and (type(indelRate) == tuple) \
		  and (len(indelRate) == 2):
			(indelRate,indels) = indelRate
		elif (indels != None) \
		  and (indelRate == None) \
		  and (type(indels) == tuple) \
		  and (len(indels) == 2):
			(indelRate,indels) = indels

		# if the mutation rates are both zero, just return a copy

		if (subRate   == None) or (subRate   <= 0): subRate   = 0.0
		if (indelRate == None) or (indelRate <= 0): indelRate = 0.0

		if (subRate == 0.0) and (indelRate == 0.0):
			return EchyDna(self.seq)

		# sanity check to make sure the sum of the rates is not more than .99;
		# because we will allow back-to-back indels, an indel rate of 1 would
		# create an infinite loop

		maxRate   = .99
		eventRate = subRate + indelRate

		if (subRate > maxRate):
			(subRate,indelRate) = (maxRate,0.0)
			eventRate = maxRate
		elif (eventRate > maxRate):
			indelRate = maxRate - subRate
			eventRate = maxRate

		# make sure we have an indel length function, if we need one
		#  none provided: give 1 bp insert or delete with equal probability;
		#  int provided:  choose 1..N and -1..-N equal probabilities
		#  list provided: choose from the list with equal probabilities

		if (indelRate > 0):
			if (indels == None):
				indels = lambda u: -1 if (u<0.5) else 1
			elif (type(indels) == int) and (indels != 0):
				maxLength = abs(indels)
				scale = 2*maxLength - 1e-5
				indels = lambda u: -(1+int(floor(u*scale))) if (u<0.5) \
				                else 1+int(floor((u-0.5)*scale))
			elif (type(indels) in [list,tuple]) and (len(indels) > 0):
				indelList = indels
				indels = lambda u: indelList[int(floor(u*len(indelList)))]

		# build nucleotide substitution tables

		if (background == None): background = EchyDna.background

		if (subRate != None):
			ntToSubs = {}
			for nt in "ACGT":
				ntToSubs[nt] = [subNt for subNt in background if (subNt != nt)]
				if (len(ntToSubs[nt]) == 0):
					raise ValueError("background contains no substitutions for \"%s\"" % nt)

		# perform mutation

		(_,chooser,spinner) = resolve_prng(seed)

		seq = []
		ix = 0
		while (ix < len(self.seq)):
			nt = self.seq[ix]

			spin = spinner()
			if (spin >= eventRate):      # no event
				seq += [nt]
				ix  += 1
				continue

			if (spin >= indelRate):      # substitution event
				if (nt in ntToSubs): seq += [chooser(ntToSubs[nt])]
				else:                seq += [chooser(background)]
				ix += 1
				continue

			spin /= indelRate            # remap spin into 0..1
			indelLen = indels(spin)
			if (indelLen == 0):          # no event
				seq += [nt]
				ix  += 1
				continue

			if (indelLen < 0):           # deletion event
				#print "delete %d at %d" % (-indelLen,ix)
				ix += (-indelLen)
			else:                        # insertion event
				#print "insert %d at %d" % (indelLen,ix)
				seq += [chooser(background) for _ in xrange(indelLen)]
				# we don't increase ix -- back-to-back indels are allowed

		# check for insertion(s) after the end of the sequence

		if (indelRate > 0):
			while (True):
				spin = spinner()
				if (spin >= indelRate):  # no event
					break
				spin /= indelRate        # remap spin into 0..1
				indelLen = indels(spin)
				if (indelLen <= 0):      # no event
					break
				#print "insert %d at %d" % (indelLen,len(self.seq))
				seq += [chooser(background) for _ in xrange(indelLen)]

		return EchyDna("".join(seq))

	def differences(self,other,symbols="-x"):
		"""a.differences(b) shows the positions where a and b differ"""
		if (not isinstance(other,EchyDna)) and (not isinstance(other,str)):
			raise ValueError("operand is not a valid string type")
		other = str(other)
		if (len(other) != len(self.seq)):
			raise ValueError("operand does not match string length")
		if (not isinstance(symbols,str)) or (len(symbols) != 2):
			raise ValueError("symbols argument is not a valid string")
		diff = lambda (a,b): symbols[0] if (a==b) else symbols[1]
		return "".join(map(diff,izip(self.seq,other)))

# geometric_indel--
#	geometric_indel(p) can be used the indels argument for EchyDna.mutate()

def geometric_indel(pExtend):
	return lambda u: -geometric_distribution(pExtend,2*u) if (u<0.5) \
	             else geometric_distribution(pExtend,2*u-1)

def geometric_distribution(pExtend,u):
	return int(floor(1+log(1-u)/log(pExtend)))

# resolve_prng--
#	Figure out which prng to use, as per the needs of EchyDna.__init__

def resolve_prng(seed=None):
	""""""
	prng = None
	if (isinstance(seed,Random)):
		prng = seed
		chooser = prng.choice
		spinner = prng.random
	elif (seed != None):
		prng = Random()
		prng.seed(seed)
		chooser = prng.choice
		spinner = prng.random
	else:
		chooser = choice
		spinner = unit_random
	return (prng,chooser,spinner)

# reverse_complement--

complementMap = maketrans("ACGTSWRYMKBDHVNacgtswrymkbdhvn",
                          "TGCASWYRKMVHDBNtgcaswyrkmvhdbn")

def reverse_complement(nts):
	return nts[::-1].translate(complementMap)

def forward_complement(nts):
	return nts.translate(complementMap)
