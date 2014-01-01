"""hapcount.py -- Infer an unknown number of haplotypes from SNP data
and read fragments.


"""
import pysam
from collections import deque, namedtuple
from math import log10
import sys

import pdb

# pysam's CIGAR codes
MATCH, INS, DEL, REF_SKIP, SOFT_CLIP, HARD_CLIP, PAD, EQUAL, DIFF = range(9)

AlleleCandidate = namedtuple("AlleleCandidate", ("tid", "pos", "ref", "alt", "prob"))
Allele = namedtuple("Allele", ("chrom", "pos", "ref", "alt"))
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

def mismatches(aln, ref_seq):
    """
    TODO: handle indels.
    """
    candidates = list()
    nmismatches = 0
    assert(len(ref_seq) == aln.alen)
    read_seq = aln.seq[aln.qstart:aln.qend]
    # our aligned sequence excludes soft-clipped regions, so trim
    # these off CIGAR string.
    trim_cigar = trim_softclipped(aln)
    if contains_indels(aln):
        pass
    else:
        if len(trim_cigar) > 1:
            pdb.set_trace()
    current_op, op_len = trim_cigar.popleft()
    cigar_pos = op_len
    readiter = enumerate(read_seq)
    for read_pos, read_base in readiter:
        if read_pos > cigar_pos:
            current_op, op_len = trim_cigar.popleft()
            cigar_pos += op_len
        if current_op == MATCH:
            if read_base != ref_seq[read_pos]:
                qual = aln.qual[read_pos]
                candidates.append(AlleleCandidate(aln.tid,
                                                  aln.pos+read_pos,
                                                  ref_seq[read_pos],
                                                  read_base, qual2prob(qual)))
                nmismatches += 1
        elif current_op == INS:
            ins_seq = list()
            for i in range(op_len):
                read_pos, read_base = next(readiter)
                ins_seq.append(read_base)
            ins_seq = "".join(ins_seq)
            pdb.set_trace()
        elif current_op == DEL:
            # deletion with respect to reference, remove bases from
            # reference seq so later parts match up.
            
            for i in range(op_len):
                read_pos, read_base = next(readiter)
                ins_seq.append(read_base)
        else:
            pdb.set_trace()
    assert(nmismatches == tagfilter(aln.tags, "NM"))
    return candidates

class Reference(object):
    def __init__(self, ref_file):
        self.ref = pysam.Fastafile(sys.argv[2])
        self.ref_file = ref_file
    def get_region(self, chrom, start=None, end=None):
        return self.ref.fetch(chrom, start, end)

class AlignmentProcessor(object):
    def __init__(self, bam_file, ref):
        self.bam = pysam.Samfile(bam_file, "rb")
        self.ref = ref
        self.bam_file = bam_file
        self._alignments = deque()
        self._alleles = deque()

        # storing current state information
        self._curr_tid = None
        self._curr_seq = None
        
    def block_reset(self, pos):
        """Reset a current block, starting with position `pos`.
        """
        self._block_alleles = set()
        self._block_start = pos
        self._block_end = None
        self._block_allele_last_pos = None
        
    def process_alignment(self, aln):
        """Process an alignment: 
        (1) compare against reference sequeunce, considering CIGAR string
        (2) identify mismatches as allele candidates
        """
        if aln.is_unmapped:
            return None
        ref_seq = self._curr_seq[aln.pos:aln.aend] # TODO check handling of clipping
        candidates = mismatches(aln, ref_seq)
        if len(candidates) > 0:
            for candidate in candidates:
                self._block_alleles.add(candidate)
                if candidate.pos > self._block_allele_last_pos:
                    self._block_allele_last_pos = candidate.pos
            #pdb.set_trace()
    
    def is_block_end(self, aln):
        """Check if we can terminate the current block. We terminate if we hit
        the end of the chromosome, or if we hit a position-sorted read
        that does not overlap any of our current alleles. Note too
        that we can only terminate if we've hit variants.
        """
        if len(self._block_alleles) == 0:
            return False
        same_chrom = aln.tid == self._curr_tid
        if not same_chrom:
            return True
        past_last_allele = aln.pos > self._block_allele_last_pos
        return past_last_allele

    def run(self, region=None):
        start_pos = 0 if region is None else region.start
        self.block_reset(start_pos)
        bamstream = self.bam
        last_pos = 0
        if region is not None:
            bamstream = bamstream.fetch(region.chrom, region.start, region.end)
        
        stop = False
        for alignedread in bamstream:
            assert(alignedread.pos >= last_pos)
            if self._curr_tid != alignedread.tid:
                # done with the last chromosome, or start of first
                if self._curr_tid is not None:
                    # process end of last chromosome TODO
                    pass
                self._curr_tid = alignedread.tid
                # fetch chromosome sequence
                self._curr_seq = self.ref.get_region(self.bam.getrname(self._curr_tid))
            end = self.process_alignment(alignedread)
            stop = self.is_block_end(alignedread)
            if stop:
                # now, take the candidate alleles and process each
                self.process_alleles()
                self.block_reset(alignedread.pos)
            last_pos = alignedread.pos

    def process_alleles(self):
        """After alignments have been processed and stopping point as been
        reached, process all alignments, inferring haplotypes from these.

        Our primary inference at this step is whether our site is
        polymorphic, but this depends on how many underlying
        haplotypes there are.

        """
        pass

if __name__ == "__main__":
    bam_file = sys.argv[1]
    ref = Reference(sys.argv[2])
    
    ap = AlignmentProcessor(bam_file, ref)
    ap.run(region=Region("1", 265745000, 265747800))
