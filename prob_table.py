#!/usr/bin/env python
"""
The ProbabilityTable class supports the fast generation of a finite random
variable from a table of probabilties or weights.  Precomputation builds a
table that allows us to convert a random value in the unit interval to an
instance of the random variable with a single table lookup and one comparison.

I believe the earliest reference to this method is A. J. Walker "An efficient
method for generating discrete random variables with general distributions,"
ACM Transactions on Mathematical Software, vol. 3, pp. 253-256, 1977.  It has
also appeared in papers by Marsaglia where it is referred to as "squaring the
histogram".

The method divides the probability space into N unit bars (where N is the
number of items to choose from).  Each bar is divided into two pieces, with the
left piece assigned to one item and the right piece assigned to another (some
bars may have no right piece).

The division boundaries and item assignments can be computed by a reasonably
simple algorithm.  The implementation here repeatedly chooses the min and max
remaining probabilities aided by a heap.  I don't recall where I found this
algorithm.  A potentially faster algorithm is described in R.A. Kronmal et al,
"On the alias method for generating random variables from a discrete
distribution," The American Statistician 33.4 (1979): 214-218.

Below we show an example probability vector and corresponding table.  To look
up the selection corresponding to the unit random number r, we scale r by 10,
lookup the corresponding entry in the table (rounding down), then choose the
left or right symbol depending on whether the scaled random number is less or
more than the indicated value.  For example, r=.503 scales to 5.03;  entry [5]
is "D 5.91 H", and since 5.03 is less than 5.91, we select D.

	                 +---------------------------------------------------+
	     probability | .077 .093 .044 .091 .126 .147 .016 .168 .169 .069 |
	                 +---------------------------------------------------+
	table            |    A    B    C    D    E    F    G    H    I    J |
	-----------------+---------------------------------------------------+
	[0]  G  0.16  I  |    -    -    -    -    -    - 0.16    - 0.84    - |
	[1]  C  1.44  H  |    -    - 0.44    -    -    -    - 0.56    -    - |
	[2]  J  2.69  F  |    -    -    -    -    - 0.31    -    -    - 0.69 |
	[3]  A  3.77  E  | 0.77    -    -    - 0.23    -    -    -    -    - |
	[4]  I  4.85  F  |    -    -    -    -    - 0.15    -    - 0.85    - |
	[5]  D  5.91  H  |    -    -    - 0.91    -    -    - 0.09    -    - |
	[6]  B  6.93  H  |    - 0.93    -    -    -    -    - 0.07    -    - |
	[7]  H  7.96  E  |    -    -    -    - 0.04    -    - 0.96    -    - |
	[8]  E  8.99  F  |    -    -    -    - 0.99 0.01    -    -    -    - |
	[9]  F 10.00     |    -    -    -    -    - 1.00    -    -    -    - |
	-----------------+---------------------------------------------------+
	selection weight | 0.77 0.93 0.44 0.91 1.26 1.47 0.16 1.68 1.69 0.69 |
	                 +---------------------------------------------------+

==========

Using standard python stuff, the same thing could be done like this (except I
don't think x can be characters; it has to be integers):

	from scipy.stats import rv_discrete
	x  = ["A" ,"B" ,"C" ,"D" ,"E" ,"F" ,"G" ,"H" ,"I" ,"J"]
	px = [.077,.093,.044,.091,.126,.147,.016,.168,.169,.069]
	sample = rv_discrete(values=(x,px)).rvs(size=1)

==========

Another solution can be found here:
	http://code.google.com/p/lea

"""

__author__ = "Bob Harris (rsharris@bx.psu.edu)"


from sys    import argv,stdin,stderr,exit
from heapq  import heapify,heappop
from random import seed as random_seed,random as unit_random

floatSlop = 1e-7


def usage(s=None):
	message = """
usage: prob_table [options]
  <symbol>:<count>   relative likelihood of the specified symbol
  --length=<number>  number of values to generate
                     (default is 10 thousand)
  --wrap=<number>    output line length
                     (by default each symbol is written on a separate line)
  --seed=<string>    set random seed"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():

	# parse the command line

	symToCount  = None
	sequenceLen = None
	wrap        = None
	debug       = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--length=")) or (arg.startswith("--len=")):
			sequenceLen = int(argVal)
			assert (sequenceLen > 0)
		elif (arg.startswith("--wrap=")):
			wrap = int(argVal)
			assert (wrap > 0)
		elif (arg.startswith("--seed=")):
			random_seed(argVal)
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		elif (len(arg) > 2) and (arg[1] == ":"):
			(sym,count) = (arg[0],arg[2:])
			if (symToCount == None): symToCount = {}
			assert (sym not in symToCount), \
			       "%s has more than one count (%d and %s)" \
			     % (sym,symToCount[sym],count)
			try:
				count = int(count)
				if (count < 0): raise ValueError
			except ValueError:
				assert (False), \
				       "%s as not a valid count (%s:*d)" \
				     % (count,sym,count)
			symToCount[sym] = count
		else:
			usage("unrecognized option: %s" % arg)

	if (sequenceLen == None):
		sequenceLen = 10000

	# convert counts to probabilities, or use defaults if no counts given

	if (symToCount != None):
		totalCount = sum([symToCount[sym] for sym in symToCount])
		assert (sum > 0), \
		       "you have to give me SOME positive counts!"

		pMap = []
		for sym in symToCount:
			count = symToCount[sym]
			if (count == 0): continue
			pMap += [(count/float(totalCount),sym)]
	else:
		pMap = [(.077,"A"), (.093,"B"), (.044,"C"), (.091,"D"), (.126,"E"),
		        (.147,"F"), (.016,"G"), (.168,"H"), (.169,"I"), (.069,"J")]

	# run the test

	t = ProbabilityTable(pMap)
	if ("table" in debug):
		print >>stderr, t

	if (wrap != None): line = []
	for i in xrange(sequenceLen):
		sym = t.choice()
		if (wrap == None):
			print sym
		else:
			line += [sym]
			if (len(line) >= wrap):
				print "".join(line)
				line = []

	if (wrap != None) and (len(line) > 0):
		print "".join(line)


class ProbabilityTable(object):
	"""
	selector = ProbabilityTable(distribution)
		distribution is a list of (p,key)
		sum of p is usually 1
		keys are usually ints or strings, usually distinct

	selector.choice() selects an item from the distribution

	selector[item] is overloaded
		selector[x], with 0<=x<1 returns the item associated with x;  if x is
		.. uniformly random, the result is a selection from the distribution
		selector[item] is the probability that a particular item (key) will be
		.. selected
	"""

	def __init__(self,pMap,scaleToOne=False):
		"""
		pMap is a list of (p,sym) with p adding up to 1.0
		internally, the sum of p is N (the size of the list)
		"""
		self.set_pmap(pMap,scaleToOne=scaleToOne)
		self.build_lookup()


	def __len__(self):
		return len(self.pMap)


	def __getitem__(self,key):
		if (type(key) == float) and (0 <= key < 1):
			return self.choice(randVal=key)
		else:
			return self.symToP[key]


	def __str__(self):
		wP    = max([len("%.6f"%p)  for (p,sym1,sym2) in self.table])
		wSym1 = max([len("%s"%sym1) for (p,sym1,sym2) in self.table])
		wSym2 = max([len("%s"%sym2) for (p,sym1,sym2) in self.table if (sym2 != None)])
		wIx   = len("%d" % (len(self.table)-1))
		s = []
		for (ix,(p,sym1,sym2)) in enumerate(self.table):
			if (sym2 == None): sym2 = ""
			s += ["[%*d] %*.6f %-*s %-*s" % (wIx,ix,wP,p,wSym1,sym1,wSym2,sym2)]
		return "\n".join(s)


	def set_pmap(self,pMap,scaleToOne=False):
		if (len(pMap) < 1): raise ValueError

		symToP = {}
		syms   = []
		for (p,sym) in pMap:
			if (p == 0): continue
			if (sym not in symToP):
				symToP[sym] = p
				syms += [sym]
			else:
				symToP[sym] = p
		pMap = [(symToP[sym],sym) for sym in syms]
		if (len(pMap) < 1): raise ValueError

		pSum = sum([p for (p,sym) in pMap])
		if (scaleToOne):
			if (min([p for (p,sym) in pMap]) < 0.0): raise ValueError
			scale = len(pMap) / float(pSum)
		else:
			for (p,sym) in pMap:
				if (p <= floatSlop): raise ValueError
			if (not (1-floatSlop <= pSum <= 1+floatSlop)): raise ValueError
			scale = len(pMap)
		self.pMap = [(p*scale,sym) for (p,sym) in pMap]

		symToP = {}
		for (p,sym) in self.pMap:
			symToP[sym] = p / len(self.pMap)
		self.symToP = symToP


	def build_lookup(self):
		self.table = None

		q = list(self.pMap)
		heapify(q)

		table = []
		while (q != []):
			(p1,sym1) = heappop(q)
			#assert (p1 <= 1.0)
			if (p1 >= 1-floatSlop):		# (were it not for roundoff error,
				(p1,sym2) = (1.0,None)	#  .. this test should be (p1 == 1.0)
			else:
				(p,sym2) = max(q)
				ix = q.index((p,sym2))
				assert (p >= 1.0)
				#assert (p >= 1-floatSlop)
				#if (p < 1.0): p = 1.0
				q[ix] = (p-(1-p1),sym2)
				reheapify(q,ix)
			table += [(len(table)+p1,sym1,sym2)]

		self.table = table


	def choice(self,count=None,randVal=None):
		if (count == None):
			if (randVal == None): randVal = unit_random()
			randVal *= len(self.table)
			(p,sym1,sym2) = self.table[int(randVal)]
			if (randVal < p): return sym1
			else:             return sym2
		else:
			if (randVal != None): raise ValueError
			if (count < 0):       raise ValueError
			choices = []
			for _ in xrange(count):
				randVal = unit_random() * len(self.table)
				(p,sym1,sym2) = self.table[int(randVal)]
				if (randVal < p): choices += [sym1]
				else:             choices += [sym2]
			return choices


def reheapify(q,ix):
	while (ix > 0):
		parent = (ix-1) / 2
		if (q[ix] > q[parent]): break
		(q[ix],q[parent]) = (q[parent],q[ix])
		ix = parent


if __name__ == "__main__": main()
