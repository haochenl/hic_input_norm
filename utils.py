#!/usr/bin/env python

__author__ = 'H.L.'

import sys
import os
import pysam
from time import time

#-----------------------------------------------------------------
class callClusters(object):
    """
    the related input cluster calling algorithms. input filename should be sam or bam file.
    """
    def __init__(self, filename, type='rb'):
        """
        Initialization by reading the alignment file and generating an array where the coordinates which are sorted

        Parameters
        ----------
        filename : string
            sequencing alignment (bam) file
        type : string
            alignment file format (by default is bam file)

        """
        self.filename = filename
        self.path = '/'.join(filename.split('/')[:-1])
        self.file = pysam.AlignmentFile(filename, type)
        chrs_forward = {}
        chrs_reverse = {}
        self.references = chrSort(self.file.references)
        for ref in self.references:
            chrs_forward[ref] = []
            chrs_reverse[ref] = []
        itr = self.file.fetch(until_eof=True)
        t0 = time()
        print 'start reading input file...'
        for i in itr:
            if i.is_unmapped:
                print 'unexpected unmapped read'
                break
            else:
                if i.is_reverse:
                    chrs_reverse[self.references[i.reference_id]].append(i.pos)
                else:
                    chrs_forward[self.references[i.reference_id]].append(i.pos)
        print 'start sorting input file...'
        for ref in self.references:
            chrs_forward[ref].sort()
            chrs_reverse[ref].sort()
        print 'time spent to read and sort input file:', round((time()-t0)/60, 3), 'mins'
        self.reads_forward = chrs_forward
        self.reads_reverse = chrs_reverse

    def DBSCANwriteBed(self, eps, min_samples):
        """
        Write the dbscan input(bias) peaks/clusters into bed files for different (+/-) strands

        Parameters
        ----------
        eps : int
            the density based scanning radius
        min_samples : int
            the minimum amount of samples in each peak/cluster

        Returns
        -------
        out : bed files (coordinates are string sorted)

        """
        out_forward = open(os.path.splitext(self.filename)[0]+'.dbscan_peaks_eps%d_min%d.forward.bed' % (eps, min_samples), 'w')
        out_reverse = open(os.path.splitext(self.filename)[0]+'.dbscan_peaks_eps%d_min%d.reverse.bed' % (eps, min_samples), 'w')
        for ref in self.references:
            forwards = self.reads_forward[ref]
            reverses = self.reads_reverse[ref]
            print 'building clusters for', ref, 'forward strand...'
            clusters_forward = dbscan(forwards).clusters(eps, min_samples)
            print 'building clusters for', ref, 'reverse strand...'
            clusters_reverse = dbscan(reverses).clusters(eps, min_samples)
            print 'writing clusters for', ref, 'into bed files...'
            for i in clusters_forward:
                coords = [forwards[elem] for elem in i] 
                out_forward.write('\t'.join([ref] + [str(min(coords))] + [str(max(coords))] + [str(len(coords))]) + '\n')
            for i in clusters_reverse:
                coords = [reverses[elem] for elem in i] 
                out_reverse.write('\t'.join([ref] + [str(min(coords))] + [str(max(coords))] + [str(len(coords))]) + '\n')
        out_forward.close()
        out_reverse.close()

    def BINwriteBed(self, bin_size):
        """
        Write the binned input(bias) into bed files for different (+/-) strands

        Parameters
        ----------
        bin_size : int
            the binning for tiling the input(bias) reads

        Returns
        -------
        out : bed files (coordinates are string sorted)
        
        """
        out_forward = open(os.path.splitext(self.filename)[0]+'.bin%d.forward.bed' % bin_size, 'w')
        out_reverse = open(os.path.splitext(self.filename)[0]+'.bin%d.reverse.bed' % bin_size, 'w')
        reference_dict = self.buildRefDict()
        print 'building input hash table for forward strand...'
        forward_dict = self.buildBiasDict(self.reads_forward, bin_size)
        print 'building input hash table for reverse strand...'
        reverse_dict = self.buildBiasDict(self.reads_reverse, bin_size)
        for chr in self.references:
            for key in forward_dict[reference_dict[chr]]:
                line = '\t'.join([chr] + [str(int(key)*bin_size)] + [str((int(key)+1)*bin_size)] + [str(forward_dict[reference_dict[chr]][key])]) + '\n'
                out_forward.write(line)
        for chr in self.references:
            for key in reverse_dict[reference_dict[chr]]:
                line = '\t'.join([chr] + [str(int(key)*bin_size)] + [str((int(key)+1)*bin_size)] + [str(reverse_dict[reference_dict[chr]][key])]) + '\n'
                out_reverse.write(line)
        out_forward.close()
        out_reverse.close()

    def buildBiasDict(self, reads, bin_size):
        """
        Generate input(bias) reads hash table for Hi-C contact assignment

        Parameters
        ----------
        reads : list
            the input(bias) reads array generated during initialization
        bin_size : int
            the binning for tiling the input(bias) reads

        Returns
        -------
        out : dict

        """
        bias_dict = []
        for chr in self.references:
            chr_reads = reads[chr]
            chr_dict = {}
            pointer = 0
            for i in range(chr_reads[0]/bin_size, chr_reads[-1]/bin_size):
                count = 0
                while chr_reads[pointer] < bin_size * (i + 1):
                    count += 1
                    pointer += 1
                chr_dict[str(i)] = count
            chr_dict[str(chr_reads[-1]/bin_size)] = len(chr_reads[pointer:])
            bias_dict.append(chr_dict)
        return bias_dict

    def buildRefDict(self):
        """
        Build hash talbe for chromosome assignment
        """
        ref_dict = {}
        for i in range(len(self.references)):
            ref_dict[self.references[i]] = i
        return ref_dict





#-----------------------------------------------------------------------
class dbscan(object):
    """
    DBSCAN algorithm for one dimensional sorted list
    """
    def __init__(self, sorted_list):
        """
        Parameters
        ----------
        sorted_list : list
            sorted input(bias) reads list generated during class callClusters initialization

        """
        self.visited = [0] * len(sorted_list) 
        self.noises = [0] * len(sorted_list) 
        self.sorted_list = sorted_list
        
    def clusters(self, eps, min_samples):
        """
        Generated the clusters list of 1D DBSCAN

        Parameters
        ----------
        eps : int
            the density based scanning radius
        min_samples : int
            the minimum amount of samples in each peak/cluster

        Returns
        -------
        out : list (indeces of the reads in the sorted list)

        """
        clusters = []
        t0 = time()
        for i in range(len(self.sorted_list)):
            if self.visited[i] == 1:
                pass
            else:
                self.visited[i] = 1
                neighbors = self._regionQuery(self.sorted_list, i, eps)
                if len(neighbors) < min_samples:
                    self.noises[i] = 1
                else:
                    C = self._expandCluster(self.sorted_list, i, neighbors, eps, min_samples, clusters)
                    clusters.append(C)
        print 'time spent for scanning clusters', round((time() - t0)/60, 2), 'mins'
        return clusters

    def noises(self):
        """
        Returns indeces of noise reads in the sorted list
        """
        return self.noises

    def _expandCluster(self, sorted_list, index, neighbors, eps, min_samples, clusters):
        """
        Expand from the current candidate sample to generate a DBSCAN cluster

        Parameters
        ----------
        sorted_list : list
            sorted input(bias) reads list generated during class callClusters initialization
        index : int
            index of the current candidate sample in the sorted list
        neighbors : list
            indeces of neighbors of the current candidate sample returned from the _regionQuery method
        eps : int
            the density based scanning radius
        min_samples : int
            the minimum amount of samples in each peak/cluster
        clusters : list
            list of the clusters found so far

        Returns
        -------
        out : list (the next cluster to add)

        """
        C = []
        neighbors_dict = {}
        for i in neighbors:
            neighbors_dict[str(i)] = 1
        clusters_dict = {}
        if len(clusters) == 0:
            pass
        else:
            for i in clusters[-1]:
                clusters_dict[str(i)] = 1
        for i in neighbors:
            if self.visited[i] == 0:
                self.visited[i] = 1
                neighbors_prime = self._regionQuery(sorted_list, i, eps)
                if len(neighbors_prime) >= min_samples:
                    for elem in neighbors_prime:
                        try: 
                            neighbors_dict[str(elem)]
                            pass
                        except KeyError:
                            neighbors.append(elem)
                            neighbors_dict[str(elem)] = 1
            try:
                clusters_dict[str(i)]
                pass
            except KeyError:
                C.append(i)
        return C

    def _regionQuery(self, sorted_list, index, eps):
        """
        Find the indeces of the neighbors of a given sample

        Parameters
        ----------
        sorted_list : list
            sorted input(bias) reads list generated during class callClusters initialization
        index : int
            index of the givne sample in the sorted list
        eps : int
            the density based scanning radius

        Returns
        -------
        out : list (indeces of neighbors)

        """
        points = [index]
        plus = index + 1
        minus = index - 1
        while plus < len(sorted_list):
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


##############################################################################
## other utility functions
##############################################################################

def chrSort(chr_list):
    """
    sort a list of chromosomes in the regular order
    """
    num2chr = {}
    others = []
    chrs = []
    for i in chr_list:
        if i.startswith('chr'):
            digits = i[3:]
            try:
                chrs.append(int(digits))
                num2chr[digits] = i
            except ValueError:
                if digits == 'X':
                    num2chr['23'] = i
                    chrs.append(23)
                elif digits == 'Y':
                    num2chr['24'] = i
                    chrs.append(24)
                elif digits == 'M':
                    num2chr['25'] = i
                    chrs.append(25)
                else:
                    print 'unexpected chromosome name.'
                    others.append(i)
        else:
            others.append(i)
    others.sort()
    chrs.sort()
    return [ num2chr[str(i)] for i in chrs ] + others



