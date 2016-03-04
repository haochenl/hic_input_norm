#!/usr/bin/env python

__author__ = 'H.L.'

import sys
import os
import pysam
from time import time
import genome

class callPeaks(object):
	"""
	the related input peak calling algorithms. input filename should be sam or bam file.
	"""
	def __init__(self, filename, type='rb'):
		self.filename = filename
		self.path = '/'.join(filename.split('/')[:-1])
		file = pysam.AlignmentFile(filename, type)
		chrs_forward = {}
		chrs_reverse = {}
		self.references = genome.reSites._chrSort(file.references)
		for ref in self.references:
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

	def writeBed(self, eps, min_samples):
		out_forward = open(os.path.splitext(self.filename)[0]+'.forward.bed', 'w')
		out_reverse = open(os.path.splitext(self.filename)[0]+'.reverse.bed', 'w')
		for ref in self.references:
			forwards = self.reads_forward[ref]
			reverses = self.reads_reverse[ref]
			clusters_forward = self.dbscan(sorted_list=forwards, eps, min_samples)['clusters']
			clusters_reverse = self.dbscan(sorted_list=reverses, eps, min_samples)['clusters']
			for i in clusters_forward:
				coords = [forwards[elem] for elem in i] 
				out_forward.write('\t'.join([ref] + [str(min(coords))] + [str(max(coords))]) + '\n')
			for i in clusters_reverse:
				coords = [reverses[elem] for elem in i] 
				out_reverse.write('\t'.join([ref] + [str(min(coords))] + [str(max(coords))]) + '\n')
		out_forward.close()
		out_reverse.close()

	def dbscan(self, sorted_list, eps, min_samples):
		clusters = []
		visited = [0] * len(sorted_list)
		noise = [0] * len(sorted_list)
		for i in range(len(sorted_list)):
			if visited[i] == 1:
				pass
			else:
				visited[i] = 1
				neighbors = self._regionQuery(sorted_list, index=i, eps)
				if len(neighbors) < min_samples:
					noise[i] = 1
				else:
					C = self._expandCluster(sorted_list, index=i, visited, neighbors, eps, min_samples, clusters)
					clusters.append(C)
		return {'clusters': clusters, 'noises': noise}

	def _expandCluster(self, sorted_list, index, visited, neighbors, eps, min_samples, clusters):
		C = [index]
		for i in neighbors:
			if visited[i] == 0:
				visited[i] = 1
				neighbors_prime = self._regionQuery(sorted_list, index=i, eps)
				if len(neighbors_prime) >= min_samples:
					new_neighbors = [elem for elem in neighbors_prime if elem not in neighbors]
					neighbors += new_neighbors
			if len(clusters) == 0:
				C.append(i)
			else:
				if i not in clusters[-1]:
					C.append(i)
		return C

	def _regionQuery(self, sorted_list, index, eps):
		points = [index]
		plus = index + 1
		minus = index - 1
		while plus <= len(sorted_list):
			if sorted_list[plus] - sorted_list[index] <= eps:
				points.append(plus)
				plus += 1
			else:
				break
		while minus >= 0:
			if sorted_list[index] - sorted_list[minus] <= eps:
				points.append(minus)
				minus -= 1
			else:
				break
		return points








