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

// Todo(vsbuffalo)
// - stop rule
// - block alignment processing

enum class MDType {
  isSNP, 
  isMatch, 
  isDel
};

struct MDToken {
  MDType type;
  string seq;
  int length;
  MDToken(MDType type, string seq, unsigned int length): type(type), seq(seq), length(length) {};
};

int TokenizeMD(const string& md, vector<MDToken>& tokens) {
  // Tokenize the MD tag, through the vector vector<MDTokens>. Returns
  // the number of bases matching or mismatching in reference and read
  // (e.g. no insertions or deletions) -- the number of aligned bases.
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
	tokens.push_back(MDToken(MDType::isMatch, "", length));
	total += length;
      }
      //DEBUG("adding " << length << " matches at from MD tag: " << md);
      it = pos;
    }

    if (std::isalpha(*pos)) {
      // Base mismatch, add to tokens
      length = std::distance(it, pos);
      assert(length == 0); // no back to back SNPs; MD always uses 0
      tokens.push_back(MDToken(MDType::isSNP, string(1, *pos), 1));
      total++;
      //DEBUG("adding " << length << " SNP '" << string(1, *pos) << "' from MD tag: " << md);
      it = ++pos;
    }
    
    if (*pos == '^') {
      // Deletion; gather following alpha characters
      it = ++pos;
      pos = std::find_if_not(pos, md.end(), static_cast<int(*)(int)>(std::isalpha));
      length = std::distance(it, pos);
      tokens.push_back(MDToken(MDType::isDel, string(it, pos), length));
      //DEBUG("adding " << length << " insertion '" << string(it, pos) << "' from MD tag: " << md);
      it = pos;
    }
  }
  return total;
}

string createReferenceSequence(const BamAlignment& alignment) {
  // Recreate a reference sequence for a particular alignment. This is
  // the reference sequence that is identical to the reference at this
  // spot. This means skipping insertions or soft clipped regions in
  // reads, adding deletions back in, and keeping read matches.
  const vector<CigarOp> cigar = alignment.CigarData;
  const string querybases = alignment.QueryBases;
  string md_tag;
  alignment.GetTag("MD", md_tag);
  
  vector<MDToken> tokens;
  string refseq, alignedseq; // final ref bases; aligned portion of ref bases
  int md_len = TokenizeMD(md_tag, tokens);

  // Create reference-aligned sequence of read; doesn't contain soft
  // clips or insertions. Then, deletions and reference alleles are
  // added onto this.
  int pos=0;
  for (vector<CigarOp>::const_iterator op = cigar.begin(); op != cigar.end(); ++op) {
    if (!(op->Type == 'S' || op->Type == 'I')) {
      alignedseq.append(querybases.substr(pos, op->Length));
      pos += op->Length;
    } else {
      pos += op->Length; // increment read position past skipped bases
    }
  }

  // the size of the aligned sequence MUST equal what is returned from
  // TokenizeMD: the number of aligned bases. Not the real reference
  // sequence is this length + deletions, which we add in below.
  assert(alignedseq.size() == md_len);

  pos = 0;
  for (vector<MDToken>::const_iterator it = tokens.begin(); it != tokens.end(); ++it) {
    if (it->type == MDType::isMatch) {
      refseq.append(alignedseq.substr(pos, it->length));
      pos += it->length;
    } else if (it->type == MDType::isSNP) {
      assert(it->length == it->seq.size());
      refseq.append(it->seq);
      pos += it->length;
    } else if (it->type == MDType::isDel) {
      // does not increment position in alignedseq
      assert(it->length == it->seq.size());
      refseq.append(it->seq);
    } else {
      assert(false);
    }
  }
  return refseq;
}

bool hasInDel(const BamAlignment& alignment) {
  auto& ops = alignment.CigarData;
  return ops.end() != find_if(ops.begin(), ops.end(), [](const CigarOp& op) { 
      return op.Type == 'I' || op.Type == 'D'; });
}

bool hasIns(const BamAlignment& alignment) {
  auto& ops = alignment.CigarData;
  return ops.end() != find_if(ops.begin(), ops.end(), [](const CigarOp& op) { 
      return op.Type == 'I'; });
}

bool hasDel(const BamAlignment& alignment) {
  auto& ops = alignment.CigarData;
  return ops.end() != find_if(ops.begin(), ops.end(), [](const CigarOp& op) { 
      return op.Type == 'D'; });
}


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

pos_t VariantProcessor::processMatchOrMismatch(const BamAlignment& alignment, 
					       vector<VariantPtr>& read_variants, 
					       const uint32_t& op_length, const string& refseq, 
					       const pos_t& refpos, const pos_t& readpos) {
  // Process a matching or mismatching sequence in the CIGAR string,
  // adding any SNP variants present.
  int endpos = alignment.GetEndPosition();
  for (int i = 0; i < op_length; i++) {
    assert(alignment.Position + i < endpos);
    char query_base = alignment.QueryBases[readpos + i];
    assert(refpos + i < refseq.size());
    char ref_base = refseq[refpos + i];
    if (ref_base != query_base) {
      // SNP
      string ref(1, ref_base), alt(1, query_base);
      char qual_base = alignment.Qualities[refpos + i]; // TODO check

      VariantPtr snp(new Variant(VariantType::SNP, alignment.RefID, 
				 alignment.Position+i, 1, 0, ref, alt));
      block_variants.insert(snp);
      read_variants.push_back(snp);
      //cout << "mismatch at " << alignment.Position + i <<" refbase: " << ref_base << " querybase: " << query_base << endl;      
    }
  }
}

pos_t VariantProcessor::processInsertion(const BamAlignment& alignment, 
					 vector<VariantPtr>& read_variants, 
					 const uint32_t& op_length, const string& refseq, 
					 const pos_t& refpos, const pos_t& readpos) {
  string ref;
  string alt = alignment.QueryBases.substr(readpos, op_length);
  VariantPtr ins(new Variant(VariantType::Insertion, alignment.RefID, 
			     alignment.Position + refpos, 1, 0, ref, alt));
  block_variants.insert(ins);
  read_variants.push_back(ins);
}

pos_t VariantProcessor::processDeletion(const BamAlignment& alignment, 
					vector<VariantPtr>& read_variants, 
					const uint32_t& op_length, const string& refseq, 
					const pos_t& refpos, const pos_t& readpos) {
  string alt;
  string ref = refseq.substr(refpos, op_length);
  VariantPtr del(new Variant(VariantType::Deletion, alignment.RefID, 
			     alignment.Position + refpos, 1, 0, ref, alt));
  block_variants.insert(del);
  read_variants.push_back(del);
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
  
  int nm_tag; 
  string md_tag;
  unsigned int aln_len = alignment.GetEndPosition() - alignment.Position;

  alignment.GetTag("MD", md_tag);
  alignment.GetTag("NM", nm_tag);
  
  // Reconstruct reference sequence using MD tags
  string refseq = createReferenceSequence(alignment);

  // With reconstructed reference sequence and query sequence, look
  // for variants. It's a bit roundabout to reconstruct reference from
  // MD, then use it to find variants (already in MD) but keeping
  // state between CIGAR and MD is tricky. This also is a good
  // validation; variants found must much the number of variants in
  // CIGAR/MD.
  vector<VariantPtr> variants;
  vector<VariantPtr> read_variants;
  const vector<CigarOp>& cigar = alignment.CigarData;
  int refpos = 0, readpos = 0;
  
  for (vector<CigarOp>::const_iterator op = cigar.begin(); op != cigar.end(); ++op) {
    if (op->Type == 'S') {
      readpos += op->Length;
    } else if (op->Type == 'M') {
      // match or SNP
      processMatchOrMismatch(alignment, read_variants, op->Length, refseq, refpos, readpos);
      readpos += op->Length;
      refpos += op->Length;
    } else if (op->Type == 'I') {
      processInsertion(alignment, read_variants, op->Length, refseq, refpos, readpos);
      readpos += op->Length;
    } else if (op->Type == 'D') {
      processDeletion(alignment, read_variants, op->Length, refseq, refpos, readpos);
      refpos += op->Length; // deletion w.r.t reference; skip ref length
    } else {
      cerr << "error: unidentified CIGAR type: " << op->Type << endl;
      exit(1);
    }
  }

  // Add to alignments list
  block_alignments.push_back(alignment);
  return 0; // TODO(vsbuffalo)
}

void VariantProcessor::processBlockAlignments() {
  // Re-process all alignments in a block, identifying the alleles
  // cared by different reads
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
      (*it)->print();
    }

    nmapped++;
  }

  return nmapped;
}
