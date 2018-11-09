#!/usr/bin/env python
"""
Dictionary that maps non-overlapping intervals to an Interval instance.  By
"intervals" we mean segments of the number line, *not* genomic intervals.

This isn't perfectly anologous to a python dict object.  Instead of the caller
being able to use __setitem_ to establish a key,value pair, the caller uses
add(minX,maxX) to add a key to the object.  The value for that key is an object
i, created internally, with i.minX and i.maxX.  The caller is expected to add
other instance variables to this object.

Intervals are added and used like this:

	intervals = IntervalDict()
	 ...
	interval = intervals.add(minX,maxX)  # minX and maxX are a "closed" interval
	if (interval == None):               # new interval overlaps previous interval
		report error

	 ...
	try: 
	    interval = intervals[(x)]        # we're asking what interval x is in
	except KeyError:                     # x is not contained in any interval    
		handle error
	do whatever                          # interval.minX <= x <= interval.maxX

	 ... or ...
	interval = intervals.get(x)          # we're asking what interval x is in
	if (interval != None):               # x is contained in "interval"
		do whatever                      # interval.minX <= x <= interval.maxX
	else:                                # x is not contained in any interval
		do whatever
	 ...

The current implementation simply maintains the intervals in a sorted list
using insertion sort.  This works well in the typical case where the incoming
intervals are nearly sorted.  Lookup is a simple binary search.

The implementation could clearly be inproved in the general case -- for example
using a self-balancing binary tree (like AVL) on blocks of intervals.
"""

from sys import stderr

class IntervalDict(object):
	def __init__(self):
		self.intervals = []

	def add(self,minX,maxX):
		if (self.intervals == []):
			node = Interval(minX,maxX)
			self.intervals += [node]
			return node

		ix = self.find(minX)
		if (ix < 0):
			if (maxX >= self.intervals[0].minX): return None
			node = Interval(minX,maxX)
			self.intervals = [node] + self.intervals
			return node

		if (ix == len(self.intervals)-1):
			if (minX <= self.intervals[ix].maxX): return None
			node = Interval(minX,maxX)
			self.intervals += [node]
			return node

		if (minX <= self.intervals[ix].maxX): return None
		if (maxX >= self.intervals[ix+1].minX): return None
		node = Interval(minX,maxX)
		self.intervals = self.intervals[:ix+1] + [node] + self.intervals[ix+1:]
		return node

	def __getitem__(self,x):
		node = self.get(x)
		if (node == None): raise KeyError
		return node

	def get(self,x):
		ix = self.find(x)
		if (ix < 0): return None
		node = self.intervals[ix]
		return node if (x <= node.maxX) else None

	def overlapper(self,minX,maxX):
		if (self.intervals == []):
			return None

		ix = self.find(minX)
		if (ix < 0):
			node = self.intervals[0]
			return node if (maxX >= node.minX) else None

		node = self.intervals[ix]
		if (ix == len(self.intervals)-1):
			return node if (minX <= node.maxX) else None

		if (minX <= node[ix].maxX): return node
		node = self.intervals[ix+1]
		return node if (maxX >= node.minX) else None

	def find(self,x):
		# generally returns ix such that
		#   intervals[ix].minX <= x < intervals[ix+1].minX
		# special cases:
		#   returns -1 if intervals is empty
		#   returns -1 if x < intervals[0].minX
		#   returns L-1 if x >= intervals[L-1].minX, where L=len(intervals)

		if (self.intervals == []): return -1

		ixLo = 0
		xLo  = self.intervals[ixLo].minX
		if (x < xLo): return -1

		ixHi = len(self.intervals)-1
		xHi  = self.intervals[ixHi].minX
		if (x >= xHi): return ixHi

		while (ixHi > ixLo+1):
			# invariants:
			#   ixHi-ixLo >= 2
			#   xLo = intervals[ixLo].minX
			#   xHi = intervals[ixHi].minX
			#   xLo <= x < xHi

			ixMid = (ixLo + ixHi) / 2
			xMid  = self.intervals[ixMid].minX
			if (x < xMid): (ixHi,xHi) = (ixMid,xMid)
			else:          (ixLo,xLo) = (ixMid,xMid)

		return ixLo

	def __iter__(self):
		for node in self.intervals:
			yield node

	def __str__(self):
		return " ".join([str(node) for node in self.intervals])


class Interval(object):
	def __init__(self,minX,maxX):
		self.minX   = minX
		self.maxX   = maxX

	def __str__(self):
		if (self.minX == self.maxX): return "%s"    %  self.minX
		else:                        return "%s-%s" % (self.minX,self.maxX)
