#!/usr/bin/env python

__author__ = 'H.L.'

import sys
import os
import pysam
from time import time

class peak(object):
	"""
	the related input peak calling algorithms. input filename should be sam or bam file.
	"""
	def __init__(self, filename, type='rb'):
		self.filename = filename
		self.path = '/'.join(filename.split('/')[:-1])
		file = pysam.AlignmentFile(filename, type)
		chrs_forward = {}
		chrs_reverse = {}
		for ref in file.references:
			chrs_forward[ref] = []
			chrs_reverse[ref] = []
		itr = file.fetch(until_eof=True)
		t0 = time()
		print 'start reading file...'
		for i in itr:
			if i.is_unmapped:
				print 'unexpected unmapped read'
				break
			else:
				if i.is_reverse:
					chrs_reverse[file.references[i.reference_id]].append(i.pos)
				else:
					chrs_forward[file.references[i.reference_id]].append(i.pos)
		print 'start sorting file...'
		for ref in file.references:
			chrs_forward[ref].sort()
			chrs_reverse[ref].sort()
		print 'time spent to read and sort input file:', round((time()-t0)/60, 3), 'mins'
		self.reads_forward = chrs_forward
		self.reads_reverse = chrs_reverse

	def dbscan(self, eps, min_samples, output_name):
		pass







