#!/usr/bin/env python
"""
Vary a set of parameters, inside a "ball".
"""

from sys        import argv,stdin,stderr,exit
from math       import *
from sets       import Set
from random     import seed as random_seed,randint
from ncrf_parse import int_with_unit


def usage(s=None):
	message = """
usage: parameter_ball <parameter> [options]
  <parameter>          (cumulative) a parameter and its central value; for
                       example, X=21; the parameter will be varied unless it
                       is in the set of fixed parameters
  --fixed=<parameter>  (cumulative) a parameter NOT to vary
  --radius=<offset>    how much to vary each parameter
                       (default is 1)
  --ball:sparse        the ball is a sparse hypercube
                       (this is the default)
  --ball:hyper         the ball is a complete hypercube; if the radius is more
                       than about 5, --ball:sparse should be used instead
  --ball:spikey        the ball is a "spikey burr", with only one parameter
                       changed relative to the input
  --sample=<number>    number of parameter sets to sample from the ball; note
                       that rejected sets still count toward this limit
                       (by default we report every parameter set in the ball)
  --reject=<formula>   (cumulative) reject a parameter set if the formula is
                       true; e.g. "X>Y" would reject in set in which the
                       paremeter X was greater than the parameter Y
                       (usually the formala has to be enclosed in quotes)
  --nocenter           exclude the central point from the ball
  --seed=<string>      set random seed"""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():

	# parse the command line

	names          = []
	nameToVal      = {}
	fixedNames     = Set()
	radius         = 1
	ballKind       = "sparse hypercube"
	sampleSize     = None
	rejectCriteria = []
	excludeCenter  = False

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg.startswith("--fixed=")):
			for name in argVal.split(","):
				fixedNames.add(name)
		elif (arg.startswith("--radius=")):
			if ("by" not in argVal):
				radius = abs(int(argVal))
			else:
				(radius,step) = argVal.split("by",1)
				radius = abs(int(radius))
				step   = abs(int(step))
				assert (radius % step == 0)
				if (step != 1):
					radius = (radius/step,step)
		elif (arg in ["--ball:sparse","--ball=sparse"]):
			ballKind = "sparse hypercube"
		elif (arg in ["--ball:hyper","--ball=hyper"]):
			ballKind = "hypercube"
		elif (arg in ["--ball:spikey","--ball=spikey","--spikey"]):
			ballKind = "spikey burr"
		elif (arg.startswith("--sample=")):
			sampleSize = int_with_unit(argVal)
		elif (arg.startswith("--reject=")):
			rejectCriteria += [argVal]
		elif (arg == "--nocenter"):
			excludeCenter = True
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
		elif ("=" in arg):
			name = arg.split("=",1)[0]
			val  = int(argVal)
			if (name in nameToVal) and (nameToVal[name] != val):
				usage("you have given me more than one value for %s, %d and %d"
				    % (name,nameToVal[name],val))
			if (name not in nameToVal):
				names += [name]
				nameToVal[name] = val
		else:
			usage("unrecognized option: %s" % arg)
	
	if (names == []):
		usage("you have to give me at least one parameter to vary")

	for name in fixedNames:
		if (name not in nameToVal):
			print >>stderr, "WARNING: no value was provided for \"%s\"" % name

	# separate fixed and varying names

	variables     = []
	variableToVal = {}

	for name in names:
		if (name not in fixedNames):
			variables += [name]
			variableToVal[name] = nameToVal[name]

	# generate the parameter sets

	if (ballKind == "spikey burr"):
		ball = SpikeyBurr(variables,variableToVal,radius)
	elif (ballKind == "sparse hypercube"):
		if (sampleSize == None): sampleSize = 1
		ball = SparseHyperCube(variables,variableToVal,radius,sampleSize,excludeCenter)
	else: # if (ballKind == "hypercube"):
		ball = HyperCube(variables,variableToVal,radius)

	if (sampleSize != None):
		ballSize = ball.size(excludeCenter)
		if (sampleSize >= ballSize):
			sampleSize = None

	if (sampleSize == None):
		for params in ball.ball():
			if (excludeCenter):
				if (params_are_same(params,variableToVal)): continue

			reject = False
			for formula in rejectCriteria:
				if (evaluate(formula,params) == True):
					reject = True
					break
			if (reject): continue

			for name in names:
				if (name not in params): params[name] = nameToVal[name]
			print " ".join(["%s=%d" % (name,params[name]) for name in names])
	else:
		leftToSample = sampleSize
		leftInBall   = ballSize
		for params in ball.ball():
			if (excludeCenter):
				if (params_are_same(params,variableToVal)): continue

			reject = False
			for formula in rejectCriteria:
				if (evaluate(formula,params) == True):
					reject = True
					break
			if (reject): continue

			if (randint(0,leftInBall-1) < leftToSample):
				for name in names:
					if (name not in params): params[name] = nameToVal[name]
				print " ".join(["%s=%d" % (name,params[name]) for name in names])
				leftToSample -= 1
			leftInBall -= 1


def params_are_same(params1,params2):
	for name in params1:
		if (name not in params2): return False
	for name in params2:
		if (name not in params1): return False
		if (params1[name] != params2[name]): return False
	return True


class HyperCube(object):

	def __init__(self,names,variableToVal,radius):
		self.names     = list(names)
		self.variableToVal = dict(variableToVal)
		if (type(radius) == tuple):
			(self.radius,self.step) = radius
		else:
			self.radius = radius
			self.step   = 1

		if (excludeCenter): ballSize -= 1

	def size(self,excludeCenter):
		size = 1
		for _ in self.names:
			size *= 2*self.radius+1
		if (excludeCenter): size -= 1
		return size

	def ball(self):
		nameToMin = {}
		nameToMax = {}
		params    = {}
		for name in self.names:
			nameToMin[name] = self.variableToVal[name] - self.radius*self.step
			nameToMax[name] = self.variableToVal[name] + self.radius*self.step
			params   [name] = nameToMin[name]

		while (True):
			yield params
			stop = True
			for name in self.names:
				if (params[name] == nameToMax[name]):
					params[name] = nameToMin[name]
				else:
					stop = False
					break
			if (stop): break
			params[name] += self.step


class SparseHyperCube(HyperCube):

	def __init__(self,names,nameToVal,radius,sampleSize,excludeCenter):
		assert (sampleSize > 0)
		self.names     = list(names)
		self.nameToVal = dict(nameToVal)
		if (type(radius) == tuple):
			(self.radius,self.step) = radius
		else:
			self.radius = radius
			self.step   = 1
		self.sampleSize    = sampleSize
		self.excludeCenter = excludeCenter

	def size(self,excludeCenter):
		size = 1
		for _ in self.names:
			size *= 2*self.radius+1
		if (excludeCenter): size -= 1
		return min(size,self.sampleSize)

	def ball(self):
		history = Set()

		if (self.excludeCenter):
			params = dict(self.nameToVal)
			key = self.to_key(params)
			history.add(key)

		for _ in xrange(self.sampleSize):
			while (True):
				params = dict(self.nameToVal)
				for name in self.names:
					offset = randint(-self.radius,self.radius)
					params[name] = self.nameToVal[name] + offset*self.step

				key = self.to_key(params)
				if (key not in history): break
			history.add(key)
			yield params

	def to_key(self,params):
		return (params[name] for name in self.names)


class SpikeyBurr(HyperCube):

	def __init__(self,names,nameToVal,radius):
		self.names     = list(names)
		self.nameToVal = dict(nameToVal)
		if (type(radius) == tuple):
			(self.radius,self.step) = radius
		else:
			self.radius = radius
			self.step   = 1

	def size(self,excludeCenter):
		size = 0
		for _ in self.names:
			size += 2*self.radius+1
		if (excludeCenter): size -= 1
		return size

	def ball(self):
		params = dict(self.nameToVal)
		yield params
		for name in self.names:
			params = dict(self.nameToVal)
			for offset in xrange(-self.radius,self.radius+1):
				if (offset == 0): continue
				params[name] = self.nameToVal[name] + offset*self.step
				yield params

# evaluate--
#	Let the python interpreter evaluate a formula.
#
#	See "Using eval() safely in python" (lybniz2.sourceforge.net/safeeval.html)

def evaluate(formula,params):
	context = dict(safeDict)
	for param in params:
		context[param] = params[param]
	return eval(formula,{"__builtins__":None},context)


def round_int(x): return int(round(x))

def floor_int(x): return int(floor(x))

def ceil_int(x):  return int(ceil(x))

def log2(x):      return log(x) / log(2.0)

safeList = ["e", "exp", "log", "log10", "pi", "sqrt"]
safeDict = dict([(k,locals().get(k,None)) for k in safeList])
safeDict["abs"]   = abs
safeDict["int"]   = int
safeDict["float"] = float
safeDict["ceil"]  = ceil_int
safeDict["floor"] = floor_int
safeDict["round"] = round_int
safeDict["log2"]  = log2
safeDict["max"]   = max
safeDict["min"]   = min


if __name__ == "__main__": main()
