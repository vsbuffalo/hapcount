"""hapcount.py -- Infer an unknown number of haplotypes from SNP data
and read fragments.


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

def find_mismatches(read_seq, ref_seq, tid, pos):
    """Find mismatches between sequences. Further versions could
    categorize these to SNPs and MNPs.
    """
    assert(len(read_seq) == len(ref_seq))
    mismatches = list()
    for i, read_base in enumerate(read_seq):
        if read_base != ref_seq[i]:
            snp = SNP(tid, pos+i, ref=ref_seq[i], alt=read_base, prob=None)
            mismatches.append(snp)
    return mismatches

def find_alleles(aln, ref_seq):
    """TODO: rework so that we keep indices to active position in both
    reference sequence and read sequence. This will allow for better
    handling of indels.

    There's three bits of state we need to maintain: (1) position in
    reference sequence, (2) position in read sequence, and (3)
    corresponding CIGAR operation (as this affects (1) and (2), via
    indels).
    """
    candidates = list()
    read_seq = aln.seq[aln.qstart:aln.qend]

    # remove soft clips from CIGAR, since our alignment excludes these
    trim_cigar = trim_softclipped(aln)

    # break up read into CIGAR parts and corresponding ref parts
    start, ref_start = 0, aln.pos
    cigar_parts = list()
    for op, length in trim_cigar:
        if op == MATCH:
            refseq_frag = ref_seq[ref_start:ref_start+length]
            readseq_frag = read_seq[start:start+length]
            mm = find_mismatches(readseq_frag, refseq_frag, aln.tid, aln.pos+start)
            candidates.extend(mm)
            # full match; both incremented equally
            start += length
            ref_start += length
        elif op == INS:
            readseq_frag = read_seq[start:start+length]
            ins = Insertion(aln.tid, ref_start, readseq_frag)
            candidates.append(ins)
            # insertion wrt ref means increment read, not ref
            start += length
        elif op == DEL:
            dl = Deletion(aln.tid, ref_start, length)
            candidates.append(dl)
            # deletion wrt ref means increment ref, not read
            ref_start += length
            pdb.set_trace()
        else:
            pdb.set_trace()
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
        candidates = find_alleles(aln, self._curr_seq)
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
