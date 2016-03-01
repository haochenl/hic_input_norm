#!/usr/bin/env python

__author__ = 'H.L.'

import sys
import re
import os

class reSites(object):
    """
    Generate restriction enzyme cutting site locations across the reference. With initial input of fasta files directory.
    """
    def __init__(self, ref_fa_path):
        f = os.walk(ref_fa_path).next()
        self.reference_path = f[0]
        fa_files = [ elem for elem in f[2] if os.path.splitext(elem)[1] == '.fa' ]
        chrs = []
        for _ in fa_files:
            f = open(self.reference_path + '/' + _, 'r')
            header = f.readline()
            chrs.append(header[1:].strip())
            f.close()
        file_key = {}
        for i in range(len(chrs)):
            file_key[chrs[i]] = fa_files[i]
        self.file_key = file_key
        self.reference_names = self._chrSort(chrs)
        self.fa_files = [ file_key[i] for i in self.reference_names ]

    def _chrSort(self, files):
        num2chr = {}
        others = []
        chrs = []
        for i in files:
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

    def getREsites(self, site, output_filename):
        """
        read through the fasta files and generate a table that shows all the fragments. With inputs of the RE site and output filename.
        """
        header = ['index', 'chrom', 'start', 'end']
        out = open(output_filename, 'w')
        out.write('\t'.join(header) + '\n')
        index = 0
        for i in self.reference_names:
            f = open(self.reference_path + '/' + self.file_key[i], 'r')
            line = f.readline()
            line = f.readline()
            seq = ''
            while line:
                seq += line.strip().upper()
                line = f.readline()
            f.close()
            start = 0
            while start >= 0:
                end = seq.find(site.upper(), start+1)
                index += 1
                if end >= 0:
                    out.write(str(index) + '\t' + i + '\t' + str(start+1) + '\t' + str(end) + '\n')
                else:
                    out.write(str(index) + '\t' + i + '\t' + str(start+1) + '\t' + str(len(seq)) + '\n')
                start = end
        out.close()

                    


                
            




