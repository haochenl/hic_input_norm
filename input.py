#!/usr/bin/env python

__author__ = 'H.L.'

import sys
import os
import re
import pysam
import utils


class pairBam(object):
    """
    The object will manipulate Hi-C raw sequencing alignment files and generate disired outputs.
    """
    _gatc = {'G': 'C', 'A': 'T', 'T': 'A', 'C': 'G', 'g': 'c', 'a': 't', 't': 'a', 'c': 'g'}

    def __init__(self, read1_filename, read2_filename, type='rb'):
        """
        Parameters
        ----------
        read1_filename : string
            read 1 alignment file name
        read2_filename : string
            read 2 alignment file name
        type : string
            file format of alignmet files

        """
        self.read1          = pysam.AlignmentFile(read1_filename, type)
        self.read2          = pysam.AlignmentFile(read2_filename, type)
        self.read1_filename = read1_filename
        self.read2_filename = read2_filename
        self.path           = '/'.join(self.read1_filename.split('/')[:-1])

    def separatePair(self, site, prefix, cutoff=1000, mapq=30):
        """
        separate the paired alignment files into Hi-C contacts, typeI and typeII inputs, and junk files

        Parameters
        ----------
        site : string
            restriction end sequence (4/6 cutters)
        prefix : string
            file name prefix of the output
        cutoff : int
            distance cutoff separate input and Hi-C contacts
        mapq : int
            mapping quality score cutoff

        Returns
        -------
        out : bam files

        """
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
        total = 0
        type1 = 0
        type2 = 0
        hic   = 0
        junk  = 0
        print 'ligation junction formed by the RE is:', self._returnJunction(site)
        for r1 in itr1:
            r2 = itr2.next()
            total += 1
            if r1.qname == r2.qname:
                if self._getContact(r1, r2, cutoff, mapq):
                    hic += 1
                    hic1_output.write(r1)
                    hic2_output.write(r2)
                elif self._typeIcandidate(r1, r2, cutoff):
                    is_r1TypeI = self._sigEnd(r1, site) and r2.mapq >= mapq
                    is_r2TypeI = self._sigEnd(r2, site) and r1.mapq >= mapq
                    if is_r1TypeI:
                        type1 += 1
                        sig1_output.write(r1)
                        dist1_output.write(r2)
                    if is_r2TypeI:
                        type1 += 1
                        sig1_output.write(r2)
                        dist1_output.write(r1)
                    if not is_r1TypeI and not is_r2TypeI:
                        junk += 1
                        junk1_output.write(r1)
                        junk2_output.write(r2)
                elif self._typeIIcandidate(r1, r2, mapq):
                    is_r1TypeII = r1.is_unmapped and self._sigEnd(r1, site) and not self._findJunction(r1, site)
                    is_r2TypeII = r2.is_unmapped and self._sigEnd(r2, site) and not self._findJunction(r2, site)
                    if is_r1TypeII:
                        type2 += 1
                        sig2_output.write(r1)
                        dist2_output.write(r2)
                    if is_r2TypeII:
                        type2 += 1
                        sig2_output.write(r2)
                        dist2_output.write(r1)
                    if not is_r1TypeII and not is_r2TypeII:
                        junk += 1
                        junk1_output.write(r1)
                        junk2_output.write(r2)
                else:
                    junk += 1
                    junk1_output.write(r1)
                    junk2_output.write(r2)
            else:
                print 'unexpected header unmatch. truncated files.'
                break
        input = type1 + type2
        print 'total number of reads processed:', total
        print 'number of hic contacts:', hic
        print 'number of input reads:', input
        print 'number of typeI input reads:', type1
        print 'number of typeII input reads:', type2
        print 'number of junk reads:', junk
        sig1_output.close()
        dist1_output.close()
        sig2_output.close()
        dist2_output.close()
        hic1_output.close()
        hic2_output.close()
        junk1_output.close()
        junk2_output.close()

    def _sigEnd(self, read, site):
        """
        This method will judge if a read has the cutting site signature. This method only support regular 4 cutter such as MboI and 6 cutter such as HindIII.

        Parameters
        ----------
        read : <pysam.calignmentfile.AlignedSegment>
            each input read
        site : string
            restriction end sequence (4/6 cutters)

        Returns
        -------
        out : bool

        """
        l = len(site)
        if l == 4:
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
        elif l == 5:
            site_reverse = self._gatc[site[-1].upper()] + site[:-1].upper()
            if read.is_unmapped:
                if read.seq[:l] == site.upper():
                    return True
                else:
                    return False
            else:
                if not read.is_reverse and read.seq[:l] == site.upper():
                    return True
                elif read.is_reverse and read.seq[-l:] == site_reverse:
                    return True
                else:
                    return False
        else:
            print 'unusual cutting site signature.'
            sys.exit()


    def _findJunction(self, read, site):
        """
        This method only identify ligation junctions generated by regular six cutters and four cutters.

        Parameters:
        -----------
        read : <pysam.calignmentfile.AlignedSegment>
            each input read
        site : string
            restriction end sequence (4/6 cutters)

        Returns
        -------
        out : bool

        """
        l = len(site)
        if l == 4:
            junction = self._returnJunction(site)
            if re.search(junction, read.seq):
                return True
            else:
                return False
        elif l == 5:
            junction = self._returnJunction(site)
            if re.search(junction, read.seq):
                return True
            else:
                return False
        else:
            print 'unusual cutting site signature.'
            sys.exit()

    def _returnJunction(self, site):
        """
        This method only identify ligation junctions generated by regular six cutters and four cutters.

        Parameters
        ----------
        site : string
            restriction end sequence (4/6 cutters)

        Returns
        -------
        out : string

        """
        l = len(site)
        if l == 4:
            junction = str(site.upper())*2
            return junction
        elif l == 5:
            junction = self._gatc[str(site.upper())[4]] + str(site.upper())[:4] + str(site.upper())
            return junction
        else:
            print 'unusual cutting site signature.'
            sys.exit()

    def _typeIcandidate(self, read1, read2, cutoff):
        """
        This method will read through the pair and judge if the pair is a type one input candidate or not, which means they will satisfy all constrains except for the cutting site signature.

        Parameters
        ----------
        read1 : <pysam.calignmentfile.AlignedSegment>
            input read 1 
        read2 : <pysam.calignmentfile.AlignedSegment>
            input read 2
        cutoff : int
            distance cutoff to separate input and Hi-C contacts

        Returns
        -------
        out : bool

        """
        if (not read1.is_unmapped) and (not read2.is_unmapped):
            if read1.reference_id != read2.reference_id:
                return False
            else:
                if not read1.is_reverse and read2.is_reverse:
                    if 0 <= read2.pos-read1.pos < cutoff:
                        return True
                    else:
                        return False
                elif read1.is_reverse and not read2.is_reverse:
                    if 0 <= read1.pos-read2.pos < cutoff:
                        return True
                    else:
                        return False
                else:
                    return False
        else:
            return False

    def _typeIIcandidate(self, read1, read2, mapq):
        """
        This method will read through the pair and judge if the pair is a type one input candidate or not. 5' cutting site signature and ligation junction filtering is to be determined.

        Parameters
        ----------
        read1 : <pysam.calignmentfile.AlignedSegment>
            input read 1
        read2 : <pysam.calignmentfile.AlignedSegment>
            input read 2
        mapq : int
            mapping quality score cutoff

        Returns
        -------
        out : bool

        """
        if read1.is_unmapped and (not read2.is_unmapped) and read2.mapq >= mapq:
            return True
        elif read2.is_unmapped and (not read1.is_unmapped) and read1.mapq >= mapq:
            return True
        else:
            return False

    def _getContact(self, read1, read2, cutoff, mapq):
        """
        This method will read through the pair and judge if the pair is a Hi-C contact

        Parameters
        ----------
        read1 : <pysam.calignmentfile.AlignedSegment>
            input read 1
        read2 : <pysam.calignmentfile.AlignedSegment>
            input read 2
        cutoff : int 
            distance cutoff to separate input and Hi-C contacts: int
        mapq : int
            mapping quality score cutoff: int

        Returns
        -------
        out : bool

        """
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

    def buildTagPeaks(self, peaks_forward_file, peaks_reverse_file, eps, output_filename):
        """
        This method gets reference from the input peak bed files on +/- strands and generate contact tag files with bias (input) info

        Parameters
        ----------
        peaks_forward_file : string
            input peak bed file name for forward (+) strand
        peaks_reverse_file : string
            input peak bed file name for reverse (-) strand
        eps : int
            total flexibility (eps) of peak assignment (same eps if use the dbscan method)
        output_filename : string
            output file name of the binary contacts

        Returns
        -------
        out : tag file

        """
        references = self.read1.references
        chr_dict = {}
        step = eps/2
        for i in range(len(references)):
            chr_dict[references[i]] = i
        print 'building bias reference of forward strand...'
        peaks_forward_dict = self._peaksDict(peaks_forward_file, eps)
        print 'building bias reference of reverse strand...'
        peaks_reverse_dict = self._peaksDict(peaks_reverse_file, eps)
        self.read1.reset()
        self.read2.reset()
        itr1 = self.read1.fetch(until_eof=True)
        itr2 = self.read2.fetch(until_eof=True)
        out = open(output_filename, 'w')
        print 'writing Hi-C contact pairs with bias info...'
        for r1 in itr1:
            r2 = itr2.next()
            if r1.qname == r2.qname:
                chr_r1 = references[r1.reference_id]
                chr_r2 = references[r2.reference_id]
                pos_r1 = r1.pos
                pos_r2 = r2.pos
                if r1.is_reverse:
                    r1_strand = '-'
                    r1_bias = self._findBias(chr_dict, peaks_reverse_dict, chr_r1, pos_r1, step)
                else:
                    r1_strand = '+'
                    r1_bias = self._findBias(chr_dict, peaks_forward_dict, chr_r1, pos_r1, step)
                if r2.is_reverse:
                    r2_strand = '-'
                    r2_bias = self._findBias(chr_dict, peaks_reverse_dict, chr_r2, pos_r2, step)
                else:
                    r2_strand = '+'
                    r2_bias = self._findBias(chr_dict, peaks_forward_dict, chr_r2, pos_r2, step)
                out.write('\t'.join([chr_r1]+[r1_strand]+[str(pos_r1)]+[str(r1_bias)]+[chr_r2]+[r2_strand]+[str(pos_r2)]+[str(r2_bias)])+'\n')
            else:
                print 'unexpected header unmatch. truncated files.'
                break
        out.close()

    def _findBias(self, chr_dict, peaks_dict, chr, pos, step):
        """
        Get the bias (input) info for each Hi-C contact
        """
        try:
            bias = peaks_dict[chr_dict[chr]][str(pos/step)]
        except KeyError:
            bias = '*'
        return bias

    def _peaksDict(self, peaks_file, eps):
        """
        Generate peak info hash table from the input(bias) bed file

        Parameters
        ----------
        peaks_file : string
            input(bias) bed file (+/- strand)
        eps: int
            total flexibility (eps) of peak assignment (same eps if use the dbscan method)

        Returns
        -------
        out : dict

        """
        references = self.read1.references
        step = eps/2
        chr_dict = {}
        peaks_dict = []
        for i in range(len(references)):
            chr_dict[references[i]] = i
            peaks_dict.append({})
        file = open(peaks_file, 'r')
        line = file.readline().strip().split()
        while line:
            chr = line[0]
            start = int(line[1])
            end = int(line[2])
            count = int(line[3])
            for i in range(start/step - 1, end/step + 2):
                peaks_dict[chr_dict[chr]][str(i)] = count*500.0/(end-start)
            line = file.readline().strip().split()
        file.close()
        return peaks_dict

    def buildTagBins(self, bias_filename, bin_size, output_filename):
        """
        This method gets reference from the binned input bed files on +/- strands and generate contact tag files with bias (input) info

        Parameters
        ----------
        bias_filename : string
            input(bias) alignment bam file (combined with typeI and typeII)
        bin_size : int
            size of the bin for tiling the input(bias) alignment files: int
        output_filename : string
            output file name of the binary contacts: string

        Returns
        -------
        out : tag file

        """
        bias_file = utils.callClusters(bias_filename)
        references = self.read1.references
        chr_dict = bias_file.buildRefDict()
        print 'building input hash table for forward strand...'
        forward_dict = bias_file.buildBiasDict(bias_file.reads_forward, bin_size)
        print 'building input hash table for reverse strand...'
        reverse_dict = bias_file.buildBiasDict(bias_file.reads_reverse, bin_size)
        self.read1.reset()
        self.read2.reset()
        itr1 = self.read1.fetch(until_eof=True)
        itr2 = self.read2.fetch(until_eof=True)
        out = open(output_filename, 'w')
        print 'writing Hi-C contact pairs with bias info...'
        for r1 in itr1:
            r2 = itr2.next()
            if r1.qname == r2.qname:
                chr_r1 = references[r1.reference_id]
                chr_r2 = references[r2.reference_id]
                pos_r1 = r1.pos
                pos_r2 = r2.pos
                if r1.is_reverse:
                    r1_strand = '-'
                    try:
                        r1_bias = reverse_dict[chr_dict[chr_r1]][str(pos_r1/bin_size)] 
                    except KeyError:
                        r1_bias = 0
                else:
                    r1_strand = '+'
                    try:
                        r1_bias = forward_dict[chr_dict[chr_r1]][str(pos_r1/bin_size)] 
                    except KeyError:
                        r1_bias = 0
                if r2.is_reverse:
                    r2_strand = '-'
                    try:
                        r2_bias = reverse_dict[chr_dict[chr_r2]][str(pos_r2/bin_size)] 
                    except KeyError:
                        r1_bias = 0
                else:
                    r2_strand = '+'
                    try:
                        r2_bias = forward_dict[chr_dict[chr_r2]][str(pos_r2/bin_size)] 
                    except KeyError:
                        r1_bias = 0
                out.write('\t'.join([chr_r1]+[r1_strand]+[str(pos_r1)]+[str(r1_bias)]+[chr_r2]+[r2_strand]+[str(pos_r2)]+[str(r2_bias)])+'\n')
            else:
                print 'unexpected header unmatch. truncated files.'
                break
        out.close()


if __name__ == '__main__':
    read1_filename = str(sys.argv[1])
    read2_filename = str(sys.argv[2])
    site = str(sys.argv[3])
    prefix = read1_filename.split('/')[-1][:10]
    rao = pairBam(read1_filename, read2_filename)
    rao.separatePair(site, prefix)

