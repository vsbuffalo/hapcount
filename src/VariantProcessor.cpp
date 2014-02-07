#include <algorithm>
#include <iostream>
#include <cctype>
#include <vector>
#include <string>
#include <deque>
#include <string>
#include <memory>

#include "VariantProcessor.h"
#include "Variant.h"
#include "api/BamReader.h"

using namespace BamTools;
using namespace std;

#define DEBUG(msg) std::cerr << "[debug] " << msg << endl;

// Block processing:
// (1) First process all alignments until stop condition is met, gathering all
// variant positions
// (2) Iterate through each alignment (stored in memory) in a block,
// identifying which alleles it carries. At this point, any necessary allele
// filtering can be carried out.
// (3) Build haplotypes.

void VariantProcessor::blockReset() {
  // Block initialization
  assert(block_start == 0);
  blockReset(0);
}

void VariantProcessor::blockReset(pos_t pos) {
  // Clear all temporary data structures; these use smart pointers so
  // members will be deleted.
  block_variants.clear();
  block_read_variants.clear();
  block_alignments.clear();

  block_start = pos;
  block_end = -1;
  block_last_variant_pos = -1;
}

void VariantProcessor::processBlockAlignments() {
  // Reprocess all alignments in a block, identifying the alleles
  // cared by different reads
  pos_t last_aln_pos = -1;
  while (block_read_variants.size()) {
    BamAlignment alignment = block_alignments.front();
    last_aln_pos = alignment.Position;


    block_alignments.pop_front();
  }
}

bool VariantProcessor::isBlockEnd(const BamAlignment& alignment) {
  if (alignment.RefID != curr_refid) return true;
  if (block_variants.size() == 0) return false;
  if (alignment.GetEndPosition() > block_last_variant_pos) return true;
  return false;
}

int VariantProcessor::run() {

  int nmapped = 0;
  last_aln_pos = 0;
  bool stop;
  BamAlignment al; // TODO: mem copying issues?

  if (!reader.IsOpen()) {
    std::cerr << "error: BAM file '" << filename << "' is not open." << std::endl;
  }

  while (reader.GetNextAlignment(al)) {
    if (!al.IsMapped()) continue;

    // TODO add chromosome checking code here

    assert(al.Position >= last_aln_pos); // ensure is sorted

    if (al.Position > block_start) {
      // only check if we can stop if we've moved in the block
      stop = isBlockEnd(al);
    }

    if (stop) {
      // end of block; process all variants in this block and output
      // haplotype count statistics
      processBlockAlignments();

      // reset block
      blockReset((pos_t) al.Position);
      stop = false;
    } else {
      // process read
      last_aln_pos = processAlignment(al);
    }

    // for debug TODO
    for (set<VariantPtr>::const_iterator it = block_variants.begin(); it != block_variants.end(); ++it) {
      cout << **it << endl;
    }

    nmapped++;
  }

  return nmapped;
}
