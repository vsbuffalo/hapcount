#include <vector>
#include <deque>
#include "AlignmentProcessor.h"
#include "api/BamReader.h"

#define DEBUG(msg) std::cerr << msg << endl;

using namespace BamTools;


void AlignmentProcessor::blockReset(pos_t pos) {
  // Destroy all of the last block's data. Note that we can use the
  // destructor since we're storing values, not pointers in the deque.
  if (block_variants != nullptr) 
    delete block_variants;
  if (block_reads_variants != nullptr)
    delete block_reads_variants;
  if (block_reads != nullptr)
    delete block_reads;
  
  block_start = block_end;
  block_end = 0;
  block_last_variant_pos = 0;
}

void AlignmentProcessor::processAlignment(const BamAlignment& alignment) {
  /* 
     For each alignment, extract the MD and NM tags, validate against
     CIGAR string, and store alleles. All reads for a block are stored
     in a deque, and processed again to create candidate haplotypes.
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
  
}

void AlignmentProcessor::isBlockEnd(const BamAlignment& alignment) {
  
}

int AlignmentProcessor::run() {

  int nmapped = 0;
  pos_t last_aln_pos = 0;

  BamAlignment al; // TODO: mem copying issues?
  
  if (!reader.IsOpen()) {
    std::cerr << "error: BAM file '" << filename << "' is not open." << std::endl;
  }
  
  while (reader.GetNextAlignment(al)) {

    if (!al.IsMapped()) continue;

    assert(al.Position >= last_aln_pos); // ensure is sorted
    
    this->processAlignment(al);

    last_aln_pos = al.Position;
    nmapped++;
  }
  return nmapped;
}
