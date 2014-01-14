"""hapcount.py -- Infer an unknown number of haplotypes from SNP data
and read fragments.

62,147,723
"""
import pysam
import numpy as np
from collections import deque, namedtuple
from math import log10
import sys
from itertools import groupby
import pdb

# pysam's CIGAR codes
MATCH, INS, DEL, REF_SKIP, SOFT_CLIP, HARD_CLIP, PAD, EQUAL, DIFF = range(9)

SNP = namedtuple("SNP", ("tid", "pos", "ref", "alt", "prob"))
Insertion = namedtuple("Insertion", ("tid", "pos", "seq"))
Deletion = namedtuple("Deletion", ("tid", "pos", "width"))
Region = namedtuple("Region", ("chrom", "start", "end"))

def tagfilter(tags, tag):
    for key, value in tags:
        if key == tag:
            return value

def qual2prob(qual):
    """
    Assumes Sanger-encoded quality.
    """
    val = ord(qual) - 33
    return 10**(val/-10.0)

def contains_indels(aln):
    return any([op == INS or op == DEL for op, length in aln.cigar])

def trim_softclipped(aln):
    """Trim off any soft-clipped regions from the alignment ends; assert
    there are no soft clips in middle.
    """
    cigar = deque(aln.cigar) # copy; don't tinker with internal attrs
    if aln.cigar[0][0] == SOFT_CLIP:
        cigar.popleft()
    if aln.cigar[-1][0] == SOFT_CLIP:
        cigar.pop()
    assert(not any([op == SOFT_CLIP for op, _ in cigar]))
    return cigar

def MDTag2Variants(aln):
    """
    Parse an alignment's MD tag, extracting all variants.
    """
    pdb.set_trace()

def find_variants(aln):
    """Find variants in MD tag; check they match CIGAR string and NM tag.
    """
    
    variants = list()
    read_seq = aln.seq[aln.qstart:aln.qend]

    # remove soft clips from CIGAR, since our alignment excludes these
    trim_cigar = trim_softclipped(aln)
    
    # TODO
    return variants
            
class AlignmentProcessor(object):
    def __init__(self, bam_file):
        self.bam = pysam.Samfile(bam_file, "r")
        pdb.set_trace()
        self.bam_file = bam_file

        # storing current state information
        self._curr_tid = None
        self._curr_seq = None
        
    def block_reset(self, pos):
        """Reset a current block, starting with position `pos`.
        """
        self._block_variants = set()
        self._block_read_variants = deque()
        self._block_start = pos
        self._block_end = None
        self._block_variant_last_pos = None
        
    def process_alignment(self, aln):
        """Process an alignment: 
        (1) find variants in a single read
        (2) register a read's variants in the current block
        """
        if aln.is_unmapped:
            return None
        variants = find_variants(aln, self._curr_seq)
        if len(variants) > 0:
            self._block_read_variants.append(variants)
            
            # load in set for quick membership testing
            for var in variants:
                self._block_variants.add(var)
                if var.pos > self._block_variant_last_pos:
                    self._block_variant_last_pos = var.pos
            #pdb.set_trace()
    
    def is_block_end(self, aln):
        """Check if we can terminate the current block. We terminate if we hit
        the end of the chromosome, or if we hit a position-sorted read
        that does not overlap the last (in terms of position) of
        current alleles. Note too that we can only terminate if we've
        hit variants.
        """
        if len(self._block_variants) == 0:
            return False
        same_chrom = aln.tid == self._curr_tid
        if not same_chrom:
            return True
        past_last_variant = aln.pos > self._block_variant_last_pos
        return past_last_variant

    def run(self, region=None):
        start_pos = 0 if region is None else region.start
        self.block_reset(start_pos)
        bamstream = self.bam
        last_pos = 0
        if region is not None:
            bamstream = bamstream.fetch(region.chrom, region.start, region.end)
        
        stop = False
        for alignedread in bamstream:
            if alignedread.is_unmapped:
                continue

            assert(alignedread.pos >= last_pos)
            if self._curr_tid != alignedread.tid:
                # done with the last chromosome, or start of first
                if self._curr_tid is not None:
                    # process end of last chromosome TODO
                    pass
                self._curr_tid = alignedread.tid
            
            # only check if stop critera has been met once we process
            # all alignments at a position (either 0, or start of
            # region)
            if alignedread.pos > start_pos:
                stop = self.is_block_end(alignedread)
            if stop:
                # now, take the candidate variants and process each
                pdb.set_trace()
                self.process_variants()
                self.block_reset(alignedread.pos)
                stop = False
            else:
                # alignment is within block bounds, so process
                end = self.process_alignment(alignedread)
            last_pos = alignedread.pos

    def process_variants(self):
        """After alignments have been processed and stopping point as been
        reached, process all alignments, inferring haplotypes from these.

        Our primary inference at this step is whether our site is
        polymorphic, but this depends on how many underlying
        haplotypes there are.

        """
        pass

if __name__ == "__main__":
    bam_file = sys.argv[1]

    ap = AlignmentProcessor(bam_file)
    ap.run(region=Region("1", 265744795, 265747800))
