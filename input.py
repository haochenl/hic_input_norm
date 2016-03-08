#!/usr/bin/env python

__author__ = 'H.L.'

import sys
import os
import re
import pysam


class pairBam(object):
    """
    The object will manipulate Hi-C raw sequencing alignment files and generate disired outputs.
    """
    _gatc = {'G': 'C', 'A': 'T', 'T': 'A', 'C': 'G', 'g': 'c', 'a': 't', 't': 'a', 'c': 'g'}

    def __init__(self, read1_filename, read2_filename, type='rb'):
        self.read1          = pysam.AlignmentFile(read1_filename, type)
        self.read2          = pysam.AlignmentFile(read2_filename, type)
        self.read1_filename = read1_filename
        self.read2_filename = read2_filename
        self.path           = '/'.join(self.read1_filename.split('/')[:-1])

    def separatePair(self, site, prefix, cutoff=1000, mapq=30):
        self.read1.reset()
        self.read2.reset()
        itr1 = self.read1.fetch(until_eof=True)
        itr2 = self.read2.fetch(until_eof=True)
        sig1_output  = pysam.AlignmentFile(os.path.join(self.path, prefix + '.sig.type1.bam'), 'wb', template=self.read1)
        dist1_output = pysam.AlignmentFile(os.path.join(self.path, prefix + '.dist.type1.bam'), 'wb', template=self.read2)
        sig2_output  = pysam.AlignmentFile(os.path.join(self.path, prefix + '.sig.type2.bam'), 'wb', template=self.read1)
        dist2_output = pysam.AlignmentFile(os.path.join(self.path, prefix + '.dist.type2.bam'), 'wb', template=self.read2)
        hic1_output  = pysam.AlignmentFile(os.path.join(self.path, prefix + '.hic1.bam'), 'wb', template=self.read1)
        hic2_output  = pysam.AlignmentFile(os.path.join(self.path, prefix + '.hic2.bam'), 'wb', template=self.read2)
        junk1_output = pysam.AlignmentFile(os.path.join(self.path, prefix + '.jk1.bam'), 'wb', template=self.read1)
        junk2_output = pysam.AlignmentFile(os.path.join(self.path, prefix + '.jk2.bam'), 'wb', template=self.read2)
        for r1 in itr1:
            r2 = itr2.next()
            if r1.qname == r2.qname:
                if self._getContact(r1, r2, cutoff, mapq):
                    hic1_output.write(r1)
                    hic2_output.write(r2)
                elif self._typeIcandidate(r1, r2, cutoff, mapq):
                    if self._sigEnd(r1, site):
                        sig1_output.write(r1)
                        dist1_output.write(r2)
                    if self._sigEnd(r2, site):
                        sig1_output.write(r2)
                        dist1_output.write(r1)
                elif self._typeIIcandidate(r1, r2, mapq):
                    if r1.is_unmapped and self._sigEnd(r1, site) and not self._findJunction(r1, site):
                        sig2_output.write(r1)
                        dist2_output.write(r2)
                    if r2.is_unmapped and self._sigEnd(r2, site) and not self._findJunction(r2, site):
                        sig2_output.write(r2)
                        dist2_output.write(r1)
                else:
                    junk1_output.write(r1)
                    junk2_output.write(r2)
            else:
                print 'unexpected header unmatch. truncated files.'
                break
        sig1_output.close()
        dist1_output.close()
        sig2_output.close()
        dist2_output.close()
        hic1_output.close()
        hic2_output.close()

    def _sigEnd(self, read, site):
        """
        This method will judge if a read has the cutting site signature
        """
        l = len(site)
        if read.is_unmapped:
            if read.seq[:l] == site.upper():
                return True
            else:
                return False
        else:
            if not read.is_reverse and read.seq[:l] == site.upper():
                return True
            elif read.is_reverse and read.seq[-l:] == site.upper():
                return True
            else:
                return False

    def _findJunction(self, read, site):
        """
        This method only identify ligation junctions generated by regular six cutters and four cutters.
        """
        l = len(site)
        if l == 4:
            junction = str(site.upper())*2
            if re.search(junction, read.seq):
                return True
            else:
                return False
        elif l == 5:
            junction = _gatc[str(site.upper())[4]] + str(site.upper())[:4] + str(site.upper())
            if re.search(junction, read.seq):
                return True
            else:
                return False
        else:
            print 'unusual cutting site signature.'
            sys.exit()

    def _typeIcandidate(self, read1, read2, cutoff, mapq):
        """
        This method will read through the pair and judge if the pair is a type one input candidate or not, which means they will satisfy all constrains except for the cutting site signature.
        """
        if (not read1.is_unmapped) and (not read2.is_unmapped):
            if read1.reference_id != read2.reference_id:
                return False
            else:
                if not read1.is_reverse and read2.is_reverse:
                    if (read1.mapq >= mapq) and (read2.mapq >= mapq) and (0 <= read2.pos-read1.pos < cutoff):
                        return True
                    else:
                        return False
                elif read1.is_reverse and not read2.is_reverse:
                    if (read1.mapq >= mapq) and (read2.mapq >= mapq) and (0 <= read1.pos-read2.pos < cutoff):
                        return True
                    else:
                        return False
                else:
                    return False
        else:
            return False

    def _typeIIcandidate(self, read1, read2, mapq):
        if read1.is_unmapped and (not read2.is_unmapped) and read2.mapq >= mapq:
            return True
        elif read2.is_unmapped and (not read1.is_unmapped) and read1.mapq >= mapq:
            return True
        else:
            return False

    def _getContact(self, read1, read2, cutoff, mapq):
        if (not read1.is_unmapped) and (not read2.is_unmapped) and (read1.mapq >= mapq) and (read2.mapq >= mapq):
            if read1.reference_id != read2.reference_id:
                return True
            else:
                if abs(read1.pos-read2.pos) >= cutoff:
                    return True
                else:
                    return False
        else:
            return False

    def setCutoff(self, fragments_filename, n=10000):
        """
        Get n number of alignment pairs and plot the distance sum to RE sites distribution and determine the cutoff to identify input pairs.
        """
        self.read1.reset()
        self.read2.reset()
        itr1 = self.read1.fetch(until_eof=True)
        itr2 = self.read2.fetch(until_eof=True)

    def combineTypes(self, ):
        pass


if __name__ == '__main__':
    read1_filename = str(sys.argv[1])
    read2_filename = str(sys.argv[2])
    site = str(sys.argv[3])
    prefix = read1_filename.split('/')[-1][:10]
    rao = pairBam(read1_filename, read2_filename)
    rao.separatePair(site, prefix)

