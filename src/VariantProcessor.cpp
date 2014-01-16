#include <algorithm>
#include <iostream>
#include <cctype>
#include <vector>
#include <string>
#include <deque>
#include "VariantProcessor.h"
#include "Variant.h"
#include "api/BamReader.h"

using namespace BamTools;
using namespace std;

#define DEBUG(msg) std::cerr << "[debug] " << msg << endl;

struct MDToken {
  enum type {isSNP, isMatch, isDel};
  string seq;
  int length;
  MDToken(enum MDToken::type type, string seq, unsigned int length): seq(seq), length(length) {
    type = type;
  };
};

int TokenizeMD(const string& md, vector<MDToken>& tokens) {
  // Tokenize the MD tag, through the vector vector<MDTokens>. Returns
  // width of MD tag in bases.
  int length;
  int total=0;
  string::const_iterator it = md.begin(), end = md.end(), pos;
  assert(!std::isalpha(*it)); // should not start with character

  while (it != end) {
    // http://stackoverflow.com/a/15039964/147427
    pos = std::find_if_not(it, md.end(), static_cast<int(*)(int)>(std::isdigit));
    length = std::stoi(string(it, pos));
    
    if (length >= 0) {
      // handle end matches
      if (length > 0) {
	tokens.push_back(MDToken(MDToken::isMatch, "", length));
	total += length;
      }
      //DEBUG("adding " << length << " matches at from MD tag: " << md);
      it = pos;
    }

    if (std::isalpha(*pos)) {
      // Base mismatch, add to tokens
      length = std::distance(it, pos);
      assert(length == 0); // no back to back SNPs; MD always uses 0
      tokens.push_back(MDToken(MDToken::isSNP, string(1, *pos), 1));
      total++;
      //DEBUG("adding " << length << " SNP '" << string(1, *pos) << "' from MD tag: " << md);
      it = ++pos;
    }
    
    if (*pos == '^') {
      // Deletion; gather following alpha characters
      it = ++pos;
      pos = std::find_if_not(pos, md.end(), static_cast<int(*)(int)>(std::isalpha));
      length = std::distance(it, pos);
      tokens.push_back(MDToken(MDToken::isDel, string(it, pos), length));
      total += length;
      //DEBUG("adding " << length << " insertion '" << string(it, pos) << "' from MD tag: " << md);
      it = pos;
    }
  }
  return total;
}


void VariantProcessor::blockReset() {
  // Block initialization
  assert(block_start == 0);
  blockReset(0);
}

void VariantProcessor::blockReset(pos_t pos) {
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

pos_t VariantProcessor::processAlignment(const BamAlignment& alignment) {
  /* 
     For each alignment, extract the MD and NM tags, validate against
     CIGAR string, and create Variants and ReadHaplotypes. All reads
     for a block are stored in a deque, and processed again to create
     candidate haplotypes.

     Returns the start position of this alignment (TODO correct?)
  */  
  
  if (!alignment.HasTag("NM") || !alignment.HasTag("MD")) {
    std::cerr << "error: BamAlignment '" << alignment.Name << 
      "' does not have either NM or MD tags" << std::endl;
  }
  
  int nm_td; 
  string md_td;
  unsigned int aln_len = alignment.GetEndPosition() - alignment.Position;

  alignment.GetTag("MD", md_td);
  alignment.GetTag("NM", nm_td);
  
  // Find all variants using MD tag and CIGAR string. Note that MD
  // only indicates deletions, since this is information that would
  // require information from the reference.
  int md_pos = 0; // where we are in MD tag
  int cigar_pos = 0; // where we are in CIGAR string
  vector<MDToken> tokens;
  int md_length = TokenizeMD(md_td, tokens);

  for (std::vector<CigarOp>::const_iterator op = alignment.CigarData.begin(); op != alignment.CigarData.end(); ++op) {
    if (op->Type == 'S') continue; // soft clipped regions skipped
    
    if (op->Type == 'M') {
      // can be match or mismatch; use MD tag to search for mismatches
      
    }

    

  }

  // Load any and all variants in reads in to block_variants
  

  // Old note:
  // (1) Create a ReadHaplotype from a BamAlignment.

  // (2) Store which variants are in this ReadHaplotype in a
  // block-level set, which allows us to see all variant position.
  
  return 0; // TODO
}

bool VariantProcessor::isBlockEnd(const BamAlignment& alignment) {
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
