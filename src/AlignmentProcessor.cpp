#include <vector>
#include <deque>
#include "AlignmentProcessor.h"
#include "api/BamReader.h"

#define DEBUG(msg) std::cerr << msg << endl;

using namespace BamTools;

void AlignmentProcessor::blockReset() {
  // Block initialization
  assert(block_start == 0);
  blockReset(0);
}

void AlignmentProcessor::blockReset(pos_t pos) {
  // Destroy all of the last block's data. Note that we can use the
  // destructor since we're storing values, not pointers in the deque.
  if (block_variants != nullptr) 
    delete block_variants;
  if (block_reads_variants != nullptr)
    delete block_reads_variants;
  if (block_reads != nullptr)
    delete block_reads;
  
  block_variants = new set<Variant>;
  block_reads_variants = new deque< vector<Variant> >;
  block_reads = new deque<BamAlignment>;

  block_start = pos;
  block_end = -1;
  block_last_variant_pos = -1;
}

pos_t AlignmentProcessor::processAlignment(const BamAlignment& alignment) {
  /* 
     For each alignment, extract the MD and NM tags, validate against
     CIGAR string, and store alleles. All reads for a block are stored
     in a deque, and processed again to create candidate haplotypes.

     Returns the start position of this alignment (TODO correct?)
  */  
  
  if (!alignment.HasTag("NM") || !alignment.HasTag("MD")) {
    std::cerr << "error: BamAlignment '" << alignment.Name << 
      "' does not have either NM or MD tags" << std::endl;
  }
  
  string nm_td; 
  string md_td;
  unsigned int aln_len = alignment.GetEndPosition() - alignment.Position;

  alignment.GetTag("MD", md_td);
  
  std::cout << "MD: " << md_td << " - len: " << aln_len << endl;

  // (1) Create a ReadHaplotype from a BamAlignment.

  // (2) Store which variants are in this ReadHaplotype in a
  // block-level set, which allows us to see all variant position.
  
}

bool AlignmentProcessor::isBlockEnd(const BamAlignment& alignment) {
  
}

int AlignmentProcessor::run() {

  int nmapped = 0;
  pos_t last_aln_pos = 0;
  bool stop;
  BamAlignment al; // TODO: mem copying issues?
  
  if (!reader.IsOpen()) {
    std::cerr << "error: BAM file '" << filename << "' is not open." << std::endl;
  }

  while (reader.GetNextAlignment(al)) {
    if (!al.IsMapped()) continue;

    // TODO add chromosome checking code here

    assert(al.Position >= last_aln_pos); // ensure is sorted

    if (al.pos > block_start) {
      // only check if we can stop if we've moved in the block
      stop = isBlockEnd(al);
    }

    if (stop) {
      // end of block; process all variants in this block and output
      // haplotype count statistics
      // TODO
      
      // reset block
      blockReset((pos_t) al.Position);
      stop = false;
    } else {
      // process read
      last_aln_pos = processAlignment(al);
    }
    
    nmapped++;
  }
  return nmapped;
}
